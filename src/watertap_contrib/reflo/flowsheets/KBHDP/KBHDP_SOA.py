from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    check_optimal_termination,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed

from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.pressure_changer import Pump

from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.costing import TreatmentCosting

from watertap_contrib.reflo.flowsheets.KBHDP.components import *
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_system",
    "add_connections",
    "add_constraints",
    "add_costing",
    "set_inlet_conditions",
    "set_operating_conditions",
    "report_MCAS_stream_conc",
    "display_system_stream_table",
    "display_system_build",
    "init_system",
]


def build_system():
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.MCAS_properties = MCASParameterBlock(
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Cl_-",
            "Mg_2+",
            "K_+",
            "SiO2",
            "Na_+",
            "SO4_2-",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )

    m.fs.RO_properties = NaClParameterBlock()
    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])

    build_treatment(m)

    return m


def build_treatment(m):

    m.fs.feed = Feed(property_package=m.fs.MCAS_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.sludge = Product(property_package=m.fs.MCAS_properties)
    m.fs.UF_waste = Product(property_package=m.fs.UF_properties)

    # Define the Unit Models
    m.fs.softener = FlowsheetBlock(dynamic=False)
    m.fs.UF = FlowsheetBlock(dynamic=False)
    m.fs.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.RO = FlowsheetBlock(dynamic=False)
    m.fs.DWI = FlowsheetBlock(dynamic=False)

    # Add translator blocks
    m.fs.MCAS_to_ZO_TDS_translator = TranslatorMCAStoZO(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.ZO_TDS_to_NaCl_translator = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_softener(m.fs.softener, prop_package=m.fs.MCAS_properties)
    build_UF(m.fs.UF, prop_package=m.fs.UF_properties)
    build_ro(m.fs.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(m.fs.DWI, prop_package=m.fs.RO_properties)

    # Set default scaling
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e-2, index=("H2O"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1, index=("tds"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e5, index=("tss"))

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def add_connections(m):

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.feed.inlet,
    )

    m.fs.softener_to_translator = Arc(
        source=m.fs.softener.product.outlet,
        destination=m.fs.MCAS_to_ZO_TDS_translator.inlet,
    )

    m.fs.softener_to_sludge = Arc(
        source=m.fs.softener.disposal.outlet,
        destination=m.fs.sludge.inlet,
    )

    m.fs.MCAS_ZO_translator_to_UF = Arc(
        source=m.fs.MCAS_to_ZO_TDS_translator.outlet,
        destination=m.fs.UF.feed.inlet,
    )

    m.fs.UF_to_translatorZONaCl = Arc(
        source=m.fs.UF.product.outlet,
        destination=m.fs.ZO_TDS_to_NaCl_translator.inlet,
    )

    m.fs.UF_to_waste = Arc(
        source=m.fs.UF.disposal.outlet,
        destination=m.fs.UF_waste.inlet,
    )

    m.fs.translatorZONaCl_to_pump = Arc(
        source=m.fs.ZO_TDS_to_NaCl_translator.outlet,
        destination=m.fs.pump.inlet,
    )

    m.fs.pump_to_ro = Arc(
        source=m.fs.pump.outlet,
        destination=m.fs.RO.feed.inlet,
    )

    m.fs.RO_to_product = Arc(
        source=m.fs.RO.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.RO_to_disposal = Arc(
        source=m.fs.RO.disposal.outlet,
        destination=m.fs.DWI.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_constraints(m):

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )


def apply_scaling(m):
    add_UF_scaling(m.fs.UF)
    add_ro_scaling(m.fs.RO)
    iscale.calculate_scaling_factors(m)


def add_costing(m):
    m.fs.costing = TreatmentCosting()
    elec_cost = value(
        pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)
    )
    m.fs.costing.electricity_cost.fix(elec_cost)

    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    add_softener_costing(m.fs.softener, costing_block=m.fs.costing)
    add_UF_costing(m.fs.UF, costing_block=m.fs.costing)
    add_ro_costing(m.fs.RO, costing_block=m.fs.costing)
    add_DWI_costing(m.fs.DWI, costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )

    m.fs.costing.initialize()


def set_inlet_conditions(m, Qin=4):

    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    Qin = Qin * pyunits.Mgallon / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3
    print('\n=======> SETTING FEED CONDITIONS <======="\n')

    inlet_dict = {
        "Ca_2+": 0.61 * pyunits.kg / pyunits.m**3,
        "Mg_2+": 0.161 * pyunits.kg / pyunits.m**3,
        "Alkalinity_2-": 0.0821 * pyunits.kg / pyunits.m**3,
        "SiO2": 0.13 * pyunits.kg / pyunits.m**3,
        "Cl_-": 5.5 * pyunits.kg / pyunits.m**3,
        "Na_+": 5.5 * pyunits.kg / pyunits.m**3,
        "K_+": 0.016 * pyunits.kg / pyunits.m**3,
        "SO4_2-": 0.23 * pyunits.kg / pyunits.m**3,
    }

    # initialize feed
    m.fs.feed.pressure[0].fix(101325)
    m.fs.feed.temperature[0].fix(293)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(Qin * rho)

    for solute, solute_conc in inlet_dict.items():
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            Qin * solute_conc
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute]),
            index=("Liq", solute),
        )

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )


def set_operating_conditions(m, RO_pressure=20e5):

    set_inlet_conditions(m, Qin=4)
    set_softener_op_conditions(m.fs.softener)
    set_UF_op_conditions(m.fs.UF)
    m.fs.pump.efficiency_pump.fix(0.8)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m.fs.RO, mem_area=10000)


def init_system(m):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_softener)

    init_softener(m.fs.softener)
    propagate_state(m.fs.softener_to_translator)
    propagate_state(m.fs.softener_to_sludge)
    m.fs.sludge.initialize()

    m.fs.MCAS_to_ZO_TDS_translator.initialize()
    propagate_state(m.fs.MCAS_ZO_translator_to_UF)
    init_UF(m.fs.UF)
    propagate_state(m.fs.UF_to_translatorZONaCl)
    propagate_state(m.fs.UF_to_waste)
    m.fs.UF_waste.initialize()

    m.fs.ZO_TDS_to_NaCl_translator.initialize()

    propagate_state(m.fs.translatorZONaCl_to_pump)
    m.fs.pump.initialize()

    propagate_state(m.fs.pump_to_ro)

    init_ro_system(m.fs.RO)
    propagate_state(m.fs.RO_to_product)
    propagate_state(m.fs.RO_to_disposal)

    m.fs.product.initialize()
    init_DWI(m.fs.DWI)


def optimize(
    m,
    water_recovery=0.8,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))
    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        m.fs.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        m.fs.pump.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        m.fs.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for _, stage in m.fs.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for _, stage in m.fs.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)


def display_costing_breakdown(m, w=30):
    title = "KBHDP Unit Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(f'{"LCOW":<{w}s}{f"${m.fs.costing.LCOW():<{w}.3f}"}{"$/m3":<{w}s}')
    print(
        f'{"Feed Flow":<{w}s}{f"{value(pyunits.convert(m.fs.feed.properties[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w},.1f}"}{"MGD":<{w}s}'
    )
    print(
        f'{"Product Flow":<{w}s}{f"{value(pyunits.convert(m.fs.product.properties[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w},.1f}"}{"MGD":<{w}s}'
    )
    print(
        f'{"Pump Capital Cost":<{w}s}{f"${value(m.fs.pump.costing.capital_cost):<{w},.0f}"}{pyunits.get_units(m.fs.pump.costing.capital_cost)}'
    )
    print_RO_costing_breakdown(m.fs.RO, w=w)
    print_softening_costing_breakdown(m.fs.softener, w=w)
    print_UF_costing_breakdown(m.fs.UF, w=w)
    print_DWI_costing_breakdown(m.fs.DWI, w=w)
    title = "System Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n\n{header}\n")
    print(f'{"LCOW":<{w}s}{f"${m.fs.costing.LCOW():<{w}.3f}"}{"$/m3":<{w}s}\n')
    print(
        f'{"Total Annualized Cost":<{w}s}{f"${m.fs.costing.total_annualized_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.total_annualized_cost)}'
    )
    print(
        f'{"Total CAPEX":<{w}s}{f"${m.fs.costing.total_capital_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.total_capital_cost)}'
    )
    print(
        f'{"Total OPEX":<{w}s}{f"${m.fs.costing.total_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.total_operating_cost)}\n'
    )
    print(
        f'{f"Total Fixed OPEX":<{w}s}{f"${m.fs.costing.total_fixed_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.total_fixed_operating_cost)}'
    )
    print(
        f'{f"Agg Fixed OPEX":<{w}s}{f"${m.fs.costing.aggregate_fixed_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.aggregate_fixed_operating_cost)}'
    )
    print(
        f'{f"MLC OPEX":<{w}s}{f"${m.fs.costing.maintenance_labor_chemical_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.maintenance_labor_chemical_operating_cost)}'
    )
    print(
        f'{f"Total Variable OPEX":<{w}s}{f"${m.fs.costing.total_variable_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.total_variable_operating_cost)}'
    )
    print(
        f'{f"Agg Variable OPEX":<{w}s}{f"${m.fs.costing.aggregate_variable_operating_cost():<{w},.3f}"}{pyunits.get_units(m.fs.costing.aggregate_variable_operating_cost)}'
    )
    print(f"{'.' * (3 * w)}")
    for flow in m.fs.costing.aggregate_flow_costs:
        print(
            f'{f"{flow.upper()} Cost":<{w}s}{f"${m.fs.costing.aggregate_flow_costs[flow]():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_flow_costs[flow])}'
        )
    print("\n\n")


# def build_sweep():

#     m = build_system()
#     add_connections(m)
#     add_constraints(m)
#     set_operating_conditions(m)
#     apply_scaling(m)
#     init_system(m)
#     add_costing(m)
#     optimize(m, ro_mem_area=None, water_recovery=0.8)

#     return m


def report_all_results(m):
    display_RO_flow_table(m.fs.RO)
    report_softener(m.fs.softener)
    report_UF(m.fs.UF)
    report_RO(m.fs.RO)
    display_costing_breakdown(m)


def main():

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    optimize(m)
    results = solve(m)
    assert_optimal_termination(results)
    report_all_results(m)

    return m


if __name__ == "__main__":
    m = main()
