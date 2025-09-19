from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed

from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)

from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

from watertap_contrib.reflo.flowsheets.KBHDP.components import *
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve


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

    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])
    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = SteamParameterBlock()

    build_treatment(m)
    build_energy(m)

    return m


def build_treatment(m):
    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.liquid_prop)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.LTMED = FlowsheetBlock(dynamic=False)
    treatment.DWI = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_ZOTDS_translator = TranslatorMCAStoZO(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    treatment.ZOTDS_to_TDS_translator = TranslatorZOtoSW(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.liquid_prop,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_EC(treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(treatment.UF, prop_package=m.fs.UF_properties)
    build_LTMED(
        treatment.LTMED, liquid_prop=m.fs.liquid_prop, vapor_prop=m.fs.vapor_prop
    )
    build_DWI(treatment.DWI, prop_package=m.fs.liquid_prop)

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e-2, index=("H2O"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1, index=("tds"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e5, index=("tss"))

    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def build_energy(m):
    m.fs.energy = Block()
    m.fs.energy.FPC = FlowsheetBlock(dynamic=False)
    build_FPC(m.fs.energy.FPC)


def add_connections(m):
    treatment = m.fs.treatment

    treatment.feed_to_translator = Arc(
        source=treatment.feed.outlet,
        destination=treatment.MCAS_to_ZOTDS_translator.inlet,
    )

    treatment.translator_to_EC = Arc(
        source=treatment.MCAS_to_ZOTDS_translator.outlet,
        destination=treatment.EC.feed.inlet,
    )

    treatment.EC_to_UF = Arc(
        source=treatment.EC.product.outlet,
        destination=treatment.UF.feed.inlet,
    )

    treatment.EC_to_sludge = Arc(
        source=treatment.EC.disposal.outlet,
        destination=treatment.sludge.inlet,
    )

    treatment.UF_to_translator3 = Arc(
        source=treatment.UF.product.outlet,
        destination=treatment.ZOTDS_to_TDS_translator.inlet,
    )

    treatment.UF_to_waste = Arc(
        source=treatment.UF.disposal.outlet,
        destination=treatment.UF_waste.inlet,
    )

    treatment.translator_to_LTMED = Arc(
        source=treatment.ZOTDS_to_TDS_translator.outlet,
        destination=treatment.LTMED.feed.inlet,
    )

    treatment.LTMED_to_product = Arc(
        source=treatment.LTMED.product.outlet,
        destination=treatment.product.inlet,
    )

    treatment.LTMED_to_dwi = Arc(
        source=treatment.LTMED.disposal.outlet,
        destination=treatment.DWI.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_constraints(m):
    treatment = m.fs.treatment

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=treatment.feed.properties[0].flow_vol * m.fs.water_recovery
        == treatment.product.properties[0].flow_vol
    )


def add_treatment_costing(blk):
    blk.costing = TreatmentCosting()

    add_EC_costing(blk.EC, costing_block=blk.costing)
    add_UF_costing(blk.UF, costing_block=blk.costing)
    add_LTMED_costing(blk.LTMED, costing_block=blk.costing)
    add_DWI_costing(blk.DWI, costing_block=blk.costing)

    blk.costing.ultra_filtration.capital_a_parameter.fix(500000)
    blk.costing.total_investment_factor.fix(1)

    blk.costing.cost_process()
    blk.costing.initialize()


def add_energy_costing(blk):
    blk.costing = EnergyCosting()

    add_FPC_costing(blk.FPC, costing_block=blk.costing)

    blk.costing.cost_process()
    blk.costing.add_LCOH()
    blk.costing.initialize()


def add_costing(m):

    add_treatment_costing(m.fs.treatment)
    add_energy_costing(m.fs.energy)

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()

    m.fs.costing.add_annual_water_production(
        m.fs.treatment.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    # Adding SEC
    feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
    )

    m.fs.treatment.costing._add_flow_component_breakdowns(
        "heat", "SEC_th", feed_m3h, period=pyunits.hr
    )

    m.fs.treatment.costing._add_flow_component_breakdowns(
        "electricity", "SEC_elec", feed_m3h, period=pyunits.hr
    )
    m.fs.costing.initialize()


def apply_system_scaling(m):
    iscale.set_scaling_factor(
        m.fs.treatment.sludge.properties[0.0].flow_mass_comp["tss"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )

    iscale.set_scaling_factor(
        m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3
    )
    iscale.set_scaling_factor(
        m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3
    )


def apply_scaling(m):

    add_ec_scaling(m.fs.treatment.EC)
    add_UF_scaling(m.fs.treatment.UF)
    add_LTMED_scaling(m.fs.treatment.LTMED)
    add_FPC_scaling(m.fs.energy.FPC)
    apply_system_scaling(m)
    iscale.calculate_scaling_factors(m)


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
    m.fs.treatment.feed.pressure[0].fix(101325)
    m.fs.treatment.feed.temperature[0].fix(293)
    m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(Qin * rho)

    for solute, solute_conc in inlet_dict.items():
        m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            Qin * solute_conc
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1
            / value(
                m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute]
            ),
            index=("Liq", solute),
        )

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )


def set_operating_conditions(m, system_capacity=50):
    treatment = m.fs.treatment
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4)
    set_EC_operating_conditions(treatment.EC)
    set_UF_op_conditions(treatment.UF)
    set_LTMED_operating_conditions(treatment.LTMED)
    # Set system capacity for FPC
    set_FPC_op_conditions(m.fs.energy.FPC, system_capacity=system_capacity)


def init_treatment(blk):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(blk)}")
    blk.feed.initialize()
    propagate_state(blk.feed_to_translator)

    blk.MCAS_to_ZOTDS_translator.initialize()
    propagate_state(blk.translator_to_EC)

    init_EC(blk.EC)
    propagate_state(blk.EC_to_UF)

    init_UF(blk.UF)
    propagate_state(blk.UF_to_translator3)
    propagate_state(blk.UF_to_waste)

    blk.ZOTDS_to_TDS_translator.initialize()
    propagate_state(blk.translator_to_LTMED)

    init_LTMED(blk.LTMED)
    propagate_state(blk.LTMED_to_product)
    propagate_state(blk.LTMED_to_dwi)

    blk.product.initialize()
    init_DWI(blk.DWI)


def init_system(m):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    init_treatment(m.fs.treatment)


def optimize_system(
    m,
    water_recovery=None,
    grid_frac_heat=None,
    heat_price=None,
):
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    m.fs.obj = Objective(expr=m.fs.costing.LCOT)

    if water_recovery is not None:
        # Fix water recovery
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.treatment.LTMED.unit.recovery_vol_phase[0.0, "Liq"].unfix()
        m.fs.water_recovery.fix(water_recovery)

    if grid_frac_heat is not None:
        m.fs.energy.FPC.unit.system_capacity.unfix()
        m.fs.costing.frac_heat_from_grid.fix(grid_frac_heat)

    if heat_price is not None:
        m.fs.costing.heat_cost_buy.fix(heat_price)
        m.fs.energy.FPC.unit.system_capacity.unfix()
        m.fs.costing.frac_heat_from_grid.unfix()


def display_costing_breakdown(m, w=30):
    title = "KBHDP Unit Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print_EC_costing_breakdown(m.fs.treatment.EC, w=w)
    print_UF_costing_breakdown(m.fs.treatment.UF, w=w)
    print_LTMED_costing_breakdown(m.fs.treatment.LTMED, w=w)
    print_DWI_costing_breakdown(m.fs.treatment.DWI, w=w)
    print(f"{'*' * (3 * w)}")
    title = "System Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n\n{header}\n")
    print(f'{"LCOW":<{w}s}{f"${m.fs.costing.LCOW():<{w}.3f}"}{"$/m3":<{w}s}')
    print(f'{"LCOT":<{w}s}{f"${m.fs.costing.LCOT():<{w}.3f}"}{"$/m3":<{w}s}\n')
    sec_th = value(m.fs.treatment.costing.SEC_th_component["fs.treatment.LTMED.unit"])
    print(f'{"SEC (thermal)":<{w}s}{f"{sec_th :<{w}.3f}"}{"kWh/m3":<{w}s}\n')
    print(
        f'{"Total CAPEX":<{w}s}{f"${m.fs.costing.total_capital_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.total_capital_cost)}'
    )
    print(
        f'{"Total OPEX":<{w}s}{f"${m.fs.costing.total_operating_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.total_operating_cost)}\n'
    )
    print(
        f'{f"Total Fixed OPEX":<{w}s}{f"${m.fs.costing.total_fixed_operating_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.total_fixed_operating_cost)}'
    )
    print(
        f'{f"Agg Fixed OPEX":<{w}s}{f"${m.fs.costing.aggregate_fixed_operating_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_fixed_operating_cost)}'
    )
    print(
        f'{f"Total Variable OPEX":<{w}s}{f"${m.fs.costing.total_variable_operating_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.total_variable_operating_cost)}'
    )
    print(
        f'{f"Agg Variable OPEX":<{w}s}{f"${m.fs.costing.aggregate_variable_operating_cost():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_variable_operating_cost)}'
    )
    print(f"{'.' * (3 * w)}")
    for flow in m.fs.treatment.costing.aggregate_flow_costs:
        if flow == "electricity":
            print(
                f'{f"{flow.upper()} Cost":<{w}s}{f"${m.fs.costing.aggregate_flow_electricity_purchased():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_flow_electricity_purchased)}'
            )
        else:
            print(
                f'{f"{flow.upper()} Cost":<{w}s}{f"${m.fs.treatment.costing.aggregate_flow_costs[flow]():<{w},.0f}"}{pyunits.get_units(m.fs.treatment.costing.aggregate_flow_costs[flow])}'
            )
    print("\n\n")


def report_all_results(m, w=25):
    print(f"{'*' * (3 * w)}")
    print(f"{'*' * (3 * w)}")
    title = "KBHDP RPT 2 Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Feed Flow":<{w}s}{f"{value(pyunits.convert(m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w},.1f}"}{"MGD":<{w}s}'
    )
    print(
        f'{"Product Flow":<{w}s}{f"{value(pyunits.convert(m.fs.treatment.product.properties[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w},.1f}"}{"MGD":<{w}s}'
    )
    salinity_out = (
        value(m.fs.treatment.product.properties[0].conc_mass_phase_comp["Liq", "TDS"])
        * 1000
    )
    print(f'{"Product Salinity":<{w}s}{f"{salinity_out:<{w},.2f}"}{"mg/L":<{w}s}')
    print(
        f'{"System Recovery":<{w}s}{f"{value(m.fs.water_recovery)*100:<{w},.1f}"}{"%":<{w}s}'
    )
    report_EC(m.fs.treatment.EC, w=w)
    report_UF(m.fs.treatment.UF, w=w)
    report_LTMED(m.fs.treatment.LTMED, w=w)
    report_FPC(m.fs.energy.FPC, w=w)
    display_costing_breakdown(m, w=w)


# def build_sweep(
#     grid_frac_heat=None,
#     heat_price=None,
#     water_recovery=None,
# ):
#     m = build_system()
#     add_connections(m)
#     add_constraints(m)
#     set_operating_conditions(m)
#     apply_scaling(m)
#     init_system(m)
#     add_costing(m)
#     optimize_system(
#         m,
#         water_recovery=water_recovery,
#         grid_frac_heat=grid_frac_heat,
#         heat_price=heat_price,
#     )
#     return m


def main():

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)

    optimize_system(m, water_recovery=0.35, heat_price=0.005)
    results = solve(m)
    assert_optimal_termination(results)
    report_all_results(m)

    return m


if __name__ == "__main__":

    m = main()
