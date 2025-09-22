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
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.pressure_changer import Pump
from watertap.core.zero_order_properties import WaterParameterBlock

from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap_contrib.reflo.flowsheets.KBHDP.components import *
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve


def build_system(RE=None):
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

    build_energy(m)

    return m


def build_treatment(m):
    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.RO_properties)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.RO = FlowsheetBlock(dynamic=False)
    treatment.DWI = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_TDS_translator = TranslatorMCAStoZO(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    treatment.TDS_to_NaCl_translator = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_EC(treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(treatment.DWI, prop_package=m.fs.RO_properties)

    m.fs.treatment.product.properties[0].conc_mass_phase_comp

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


def build_energy(m):
    m.fs.energy = Block()
    m.fs.energy.pv = FlowsheetBlock()
    build_pv(m.fs.energy.pv)


def add_connections(m):
    treatment = m.fs.treatment

    treatment.feed_to_translator = Arc(
        source=treatment.feed.outlet,
        destination=treatment.MCAS_to_TDS_translator.inlet,
    )

    treatment.translator_to_EC = Arc(
        source=treatment.MCAS_to_TDS_translator.outlet,
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
        destination=treatment.TDS_to_NaCl_translator.inlet,
    )

    treatment.UF_to_waste = Arc(
        source=treatment.UF.disposal.outlet,
        destination=treatment.UF_waste.inlet,
    )

    treatment.translator_to_pump = Arc(
        source=treatment.TDS_to_NaCl_translator.outlet,
        destination=treatment.pump.inlet,
    )

    treatment.pump_to_ro = Arc(
        source=treatment.pump.outlet,
        destination=treatment.RO.feed.inlet,
    )

    treatment.ro_to_product = Arc(
        source=treatment.RO.product.outlet,
        destination=treatment.product.inlet,
    )

    treatment.ro_to_dwi = Arc(
        source=treatment.RO.disposal.outlet,
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

    blk.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=blk.costing,
    )
    add_EC_costing(blk.EC, costing_block=blk.costing)
    add_UF_costing(blk.UF, costing_block=blk.costing)
    add_ro_costing(blk.RO, costing_block=blk.costing)
    add_DWI_costing(blk.DWI, costing_block=blk.costing)

    blk.costing.ultra_filtration.capital_a_parameter.fix(500000)
    blk.costing.total_investment_factor.fix(1)

    blk.costing.cost_process()
    blk.costing.initialize()


def add_energy_costing(blk):
    blk.costing = EnergyCosting()

    add_pv_costing(blk.pv, costing_block=blk.costing)

    blk.costing.cost_process()
    blk.costing.add_LCOE()
    blk.costing.initialize()


def add_costing(m):

    add_treatment_costing(m.fs.treatment)
    add_energy_costing(m.fs.energy)

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()
    elec_cost = value(
        pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)
    )
    m.fs.costing.electricity_cost.set_value(elec_cost)

    m.fs.costing.add_annual_water_production(
        m.fs.treatment.product.properties[0].flow_vol
    )
    m.fs.treatment.costing.add_specific_energy_consumption(
        m.fs.treatment.product.properties[0].flow_vol, name="SEC"
    )
    m.fs.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOE()

    m.fs.costing.initialize()


def apply_scaling(m):

    add_ec_scaling(m.fs.treatment.EC)
    add_UF_scaling(m.fs.treatment.UF)
    add_ro_scaling(m.fs.treatment.RO)
    add_pv_scaling(m.fs.energy.pv)
    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )
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


def set_operating_conditions(m, RO_pressure=20e5, system_capacity=50):
    treatment = m.fs.treatment
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4)
    set_EC_operating_conditions(treatment.EC)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(0.8)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(treatment.RO, mem_area=10000)
    # Set system capacity for PV
    set_pv_op_conditions(m.fs.energy.pv, system_capacity=system_capacity)


def init_treatment(blk):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    blk.feed.initialize()
    propagate_state(blk.feed_to_translator)
    blk.MCAS_to_TDS_translator.initialize()
    propagate_state(blk.translator_to_EC)

    init_EC(blk.EC)
    propagate_state(blk.EC_to_UF)

    init_UF(blk.UF)
    propagate_state(blk.UF_to_translator3)
    propagate_state(blk.UF_to_waste)

    blk.TDS_to_NaCl_translator.initialize()
    propagate_state(blk.translator_to_pump)

    blk.pump.initialize()

    propagate_state(blk.pump_to_ro)

    init_ro_system(blk.RO)
    propagate_state(blk.ro_to_product)

    blk.product.initialize()
    init_DWI(blk.DWI)


def init_energy(blk):
    print("\n\n-------------------- INITIALIZING ENERGY --------------------\n\n")
    init_pv(blk.pv)


def init_system(m):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    init_treatment(m.fs.treatment)
    # init_energy(m.fs.energy)


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Pump Pressure":<{w}s}{value(pyunits.convert(blk.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<{w}.1f}{"bar"}'
    )
    print(
        f'{"Pump Work":<{w}s}{value(pyunits.convert(blk.control_volume.work[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )


def optimize_system(
    m,
    water_recovery=0.8,
    fixed_pressure=None,
    ro_mem_area=None,
    grid_frac=None,
    elec_price=None,
):
    treatment = m.fs.treatment
    print(f"\n\nDOF before optimization: {degrees_of_freedom(m)}")

    m.fs.obj = Objective(expr=m.fs.costing.LCOT)

    if water_recovery is not None:
        # Fix water recovery
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        # Optimize water recovery
        lower_bound = 0.01
        upper_bound = 0.99
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        # Membrane area is unfixed, pressure is fixed
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        treatment.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        # Membrane area is fixed, optimize pressure
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        treatment.pump.control_volume.properties_out[0].pressure.unfix()
        treatment.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        treatment.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        # Membrane area is fixed, optimize pressure
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        treatment.RO.total_membrane_area.fix(ro_mem_area)
        for _, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)
    else:
        # Membrane area is unfixed, pressure is fixed
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for _, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)

    if grid_frac is not None:
        # Grid fraction is fixed, optimize system capacity
        print(f"\n------- Fixed Grid Fraction at {100*grid_frac}% -------\n")
        m.fs.costing.frac_elec_from_grid.fix(grid_frac)
        m.fs.energy.pv.unit.system_capacity.unfix()
        m.fs.energy.pv.unit.electricity_annual.unfix()

    if elec_price is not None:
        # Electricity price is fixed, optimize grid fraction and system capacity
        print(f"\n------- Fixed Electricity Price at ${elec_price}/kWh -------\n")
        m.fs.costing.electricity_cost_buy.set_value(elec_price)
        m.fs.costing.frac_elec_from_grid.unfix()
        m.fs.energy.pv.unit.system_capacity.unfix()
        m.fs.energy.pv.unit.electricity_annual.unfix()

    for _, stage in treatment.RO.stage.items():
        stage.module.width.setub(5000)
        stage.module.feed_side.velocity[0, 0].unfix()
        stage.module.feed_side.velocity[0, 1].setlb(0.0)
        stage.module.feed_side.K.setlb(1e-6)
        stage.module.feed_side.friction_factor_darcy.setub(50)
        stage.module.flux_mass_phase_comp.setub(1)
        stage.module.feed_side.cp_modulus.setub(10)
        stage.module.rejection_phase_comp.setlb(1e-4)
        stage.module.feed_side.N_Re.setlb(1)
        stage.module.recovery_mass_phase_comp.setlb(1e-7)


def display_costing_breakdown(m, w=30):
    title = "KBHDP Unit Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print_EC_costing_breakdown(m.fs.treatment.EC, w=w)
    print_UF_costing_breakdown(m.fs.treatment.UF, w=w)
    title = "Pump Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{"Pump Capital Cost":<{w}s}{f"${value(m.fs.treatment.pump.costing.capital_cost):<{w},.0f}"}{pyunits.get_units(m.fs.treatment.pump.costing.capital_cost)}'
    )
    print_RO_costing_breakdown(m.fs.treatment.RO, w=w)
    print_DWI_costing_breakdown(m.fs.treatment.DWI, w=w)
    print(f"{'*' * (3 * w)}")
    title = "System Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n\n{header}\n")
    print(f'{"LCOW":<{w}s}{f"${m.fs.costing.LCOW():<{w}.3f}"}{"$/m3":<{w}s}\n')
    print(
        f'{"SEC":<{w}s}{f"${m.fs.treatment.costing.SEC():<{w}.3f}"}{"kWh/m3":<{w}s}\n'
    )
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
    title = "KBHDP RPT 1 Report"
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
        value(m.fs.treatment.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"])
        * 1000
    )
    print(f'{"Product Salinity":<{w}s}{f"{salinity_out:<{w},.2f}"}{"mg/L":<{w}s}')
    print(
        f'{"System Recovery":<{w}s}{f"{value(m.fs.water_recovery)*100:<{w},.1f}"}{"%":<{w}s}'
    )
    report_EC(m.fs.treatment.EC, w=w)
    report_UF(m.fs.treatment.UF, w=w)
    report_RO(m.fs.treatment.RO, w=w)
    report_pump(m.fs.treatment.pump, w=w)
    display_RO_flow_table(m.fs.treatment.RO, w=w)
    report_PV(m.fs.energy.pv, w=w)
    display_costing_breakdown(m, w=w)


def main():

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)

    optimize_system(
        m,
        ro_mem_area=20000,
        fixed_pressure=None,
        water_recovery=0.8,
        grid_frac=None,
        elec_price=0.15,
    )
    # optimize_system(
    #     m,
    #     ro_mem_area=None,
    #     fixed_pressure=62e5,
    #     water_recovery=None,
    #     grid_frac=None,
    #     elec_price=0.15,
    # )
    results = solve(m)
    assert_optimal_termination(results)
    report_all_results(m)

    return m


if __name__ == "__main__":
    m = main()
