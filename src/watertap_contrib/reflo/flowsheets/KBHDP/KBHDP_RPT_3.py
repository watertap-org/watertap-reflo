import pathlib
from pyomo.environ import (
    ConcreteModel,
    Param,
    TransformationFactory,
    Objective,
    Block,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state

import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

from watertap_contrib.reflo.flowsheets.KBHDP.components import *
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/kbhdp_case_study.yaml"


def build_system(Qin=4, Cin=12, water_recovery=0.8):

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.treatment = Block()

    m.inlet_flow_rate = pyunits.convert(
        Qin * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = water_recovery
    m.fs.water_recovery = Param(initialize=water_recovery, mutable=True)

    m.fs.treatment.costing = TreatmentCosting()

    # Property package
    m.fs.properties = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.treatment.feed = Feed(property_package=m.fs.properties)
    m.fs.treatment.product = Product(property_package=m.fs.properties)
    # m.fs.disposal = Product(property_package=m.fs.properties)

    # Create MD unit model at flowsheet level
    m.fs.treatment.md = FlowsheetBlock()

    build_md(m.fs.treatment.md, prop_package=m.fs.properties)
    m.fs.treatment.dwi = FlowsheetBlock()
    build_DWI(m.fs.treatment.dwi, m.fs.properties)

    m.fs.energy = Block()
    m.fs.energy.FPC = FlowsheetBlock()
    m.fs.energy.costing = EnergyCosting()
    build_FPC(m.fs.energy.FPC)

    return m


def add_connections(m):

    treatment = m.fs.treatment

    treatment.feed_to_md = Arc(
        source=treatment.feed.outlet, destination=treatment.md.feed.inlet
    )

    treatment.md_to_product = Arc(
        source=treatment.md.permeate.outlet, destination=treatment.product.inlet
    )

    treatment.md_to_dwi = Arc(
        source=treatment.md.concentrate.outlet,
        destination=treatment.dwi.unit.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_costing(m):

    m.fs.treatment.costing = TreatmentCosting()
    m.fs.energy.costing = EnergyCosting()

    add_FPC_costing(m.fs.energy.FPC, costing_block=m.fs.energy.costing)

    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)

    add_DWI_costing(m.fs.treatment.dwi, costing_block=m.fs.treatment.costing)

    # Treatment costing
    m.fs.treatment.costing.cost_process()

    m.fs.treatment.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)

    # Energy costing
    m.fs.energy.costing.cost_process()
    m.fs.energy.costing.add_LCOH()
    # System costing
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(0.01)
    elec_cost = value(
        pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)
    )
    m.fs.costing.electricity_cost_buy.set_value(elec_cost)
    m.fs.costing.cost_process()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")

    m.fs.treatment.costing.initialize()
    m.fs.energy.costing.initialize()
    m.fs.costing.initialize()

    m.fs.costing.add_annual_water_production(
        m.fs.treatment.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOH()
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

    # For reporting purposes
    m.fs.treatment.md.unit.capital_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.capital_cost), mutable=True
    )
    m.fs.treatment.md.unit.fixed_operating_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.fixed_operating_cost),
        mutable=True,
    )
    m.fs.treatment.md.unit.module_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.module_cost), mutable=True
    )
    m.fs.treatment.md.unit.other_capital_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.other_capital_cost),
        mutable=True,
    )


def apply_scaling(m):

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "TDS"))

    iscale.set_scaling_factor(m.fs.energy.FPC.unit.system_capacity, 1e-5)
    iscale.set_scaling_factor(m.fs.energy.FPC.unit.electricity_annual, 1e-4)

    iscale.calculate_scaling_factors(m)


def set_inlet_conditions(m):

    print(f'\n{"=======> SETTING FEED CONDITIONS <=======":^60}\n')

    m.fs.treatment.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.inlet_flow_rate,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.inlet_salinity,
            ("temperature", None): 298.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def set_operating_conditions(m, system_capacity=50):
    set_inlet_conditions(m)
    set_FPC_op_conditions(m.fs.energy.FPC, system_capacity=system_capacity)


def init_system(m):

    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    treatment.feed.initialize()

    init_md(treatment.md)

    propagate_state(treatment.md_to_product)
    treatment.product.initialize()

    propagate_state(treatment.md_to_dwi)

    init_DWI(treatment.dwi)

    init_FPC(m.fs.energy.FPC)


def report_all_results(m, w=30):

    title = "KBHDP RPT 3 Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.h,
    )

    product_m3h = pyunits.convert(
        m.fs.treatment.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.h,
    )
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(f'{"System recovery":<{w}s}{product_m3h() / feed_m3h() * 100:<{w}.2f}{"%"}')
    print(
        f'{"Electricity demand":<{w}s}{pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.MW * pyunits.h / pyunits.year)():<{w}.2f}{"MWh/year"}'
    )
    print(
        f'{"Heat demand":<{w}s}{pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat, to_units=pyunits.MW * pyunits.h / pyunits.year)():<{w}.2f}{"MWh/year"}'
    )
    print(
        f'{"Heat demand":<{w}s}{pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat, to_units=pyunits.MW * pyunits.h / pyunits.year)():<{w}.2f}{"MWh/year"}'
    )
    report_MD(m.fs.treatment.md)
    report_FPC(m.fs.energy.FPC)
    display_costing_breakdown(m, w=w)


def display_costing_breakdown(m, w=30):

    title = "System Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print_MD_costing_breakdown(m.fs.treatment.md, w=w)
    print_DWI_costing_breakdown(m.fs.treatment.dwi, w=w)
    print_FPC_costing_breakdown(m.fs.energy.FPC, w=w)

    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"LCOT":<{w}s}{value(m.fs.costing.LCOT):<{w},.2f}{pyunits.get_units(m.fs.costing.LCOT)}'
    )
    print(
        f'{"LCOW":<{w}s}{value(m.fs.treatment.costing.LCOW):<{w}.2f}{pyunits.get_units(m.fs.treatment.costing.LCOW)}\n'
    )
    print(
        f'{"Energy LCOH":<{w}s}{value(m.fs.energy.costing.LCOH):<{w}.2f}{pyunits.get_units(m.fs.energy.costing.LCOH)}'
    )
    print(
        f'{"Fraction of heat from grid":<{w}s}{value(m.fs.costing.frac_heat_from_grid)*100:<{w}.2f}{pyunits.get_units(m.fs.costing.frac_heat_from_grid)}\n'
    )
    print(
        f'{"Total CAPEX":<{w}s}{value(m.fs.costing.total_capital_cost):<{w},.2f}{pyunits.get_units(m.fs.costing.total_capital_cost)}'
    )
    print(
        f'{"Total OPEX":<{w}s}{value(m.fs.costing.total_operating_cost):<{w},.2f}{pyunits.get_units(m.fs.costing.total_operating_cost)}'
    )


def main():

    m = build_system()

    add_connections(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)

    results = solve(m.fs.treatment.md)
    assert_optimal_termination(results)
    results = solve(m, raise_on_failure=False)
    assert_optimal_termination(results)

    add_costing(m)

    results = solve(m)
    assert_optimal_termination(results)

    m.fs.lcot_objective = Objective(expr=m.fs.costing.LCOT)

    results = solve(m)
    assert_optimal_termination(results)

    report_all_results(m)

    return m


if __name__ == "__main__":

    m = main()
