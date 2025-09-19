from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed, StateJunction

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.surrogate import LTMEDSurrogate
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_LTMED",
    "set_LTMED_operating_conditions",
    "init_LTMED",
    "add_LTMED_costing",
    "add_LTMED_scaling",
    "report_LTMED",
    "print_LTMED_costing_breakdown",
]


def build_LTMED(blk, liquid_prop=None, vapor_prop=None, number_effects=12):
    print(f'\n{"=======> BUILDING LTMED <=======":^60}\n')

    if liquid_prop is None:
        m = blk.model()
        liquid_prop = m.fs.liquid_prop
    if vapor_prop is None:
        m = blk.model()
        vapor_prop = m.fs.vapor_prop

    blk.feed = StateJunction(property_package=liquid_prop)
    blk.product = StateJunction(property_package=liquid_prop)
    blk.disposal = StateJunction(property_package=liquid_prop)

    blk.unit = LTMEDSurrogate(
        property_package_liquid=liquid_prop,
        property_package_vapor=vapor_prop,
        number_effects=number_effects,
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.feed,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.distillate,
        destination=blk.product.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.brine,
        destination=blk.disposal.inlet,
    )


def set_LTMED_operating_conditions(blk, recovery_ratio=0.45, steam_temperature=80):

    blk.unit.steam_props[0].temperature.fix(steam_temperature + 273.15)
    blk.unit.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    blk.unit.recovery_vol_phase[0, "Liq"].setub(None)


def add_LTMED_scaling(blk):
    iscale.set_scaling_factor(blk.unit.specific_energy_consumption_thermal, 1e-2)
    iscale.set_scaling_factor(blk.unit.feed_props[0].flow_vol_phase["Liq"], 10)
    iscale.set_scaling_factor(blk.unit.cooling_out_props[0].flow_vol_phase["Liq"], 10)
    iscale.set_scaling_factor(blk.unit.distillate_props[0].flow_vol_phase["Liq"], 1e2)
    iscale.set_scaling_factor(
        blk.unit.distillate_props[0].flow_mass_phase_comp["Liq", "H2O"], 1e-3
    )
    iscale.set_scaling_factor(blk.unit.brine_props[0].flow_vol_phase["Liq"], 100)
    iscale.constraint_scaling_transform(blk.unit.eq_specific_area_per_m3_day, 1e-2)
    iscale.constraint_scaling_transform(blk.unit.eq_specific_area_kg_s, 1e-2)
    iscale.constraint_scaling_transform(blk.unit.eq_steam_mass_flow, 1e-2)
    iscale.constraint_scaling_transform(
        blk.unit.eq_specific_thermal_energy_consumption, 1e-2
    )
    iscale.constraint_scaling_transform(blk.unit.eq_thermal_power_requirement, 1e-4)
    iscale.constraint_scaling_transform(blk.unit.eq_feed_cool_vol_flow, 1e-4)
    iscale.constraint_scaling_transform(blk.unit.eq_feed_cool_mass_flow, 1e-2)

    for stream in [
        blk.unit.eq_distillate_temp,
        blk.unit.eq_brine_temp,
        blk.unit.eq_cooling_temp,
    ]:
        for e in stream:
            # print(stream[e])
            iscale.constraint_scaling_transform(stream[e], 1e-2)
        for pressure_stream in [
            blk.unit.eq_feed_to_distillate_isobaric,
            blk.unit.eq_feed_to_brine_isobaric,
            blk.unit.eq_feed_to_cooling_isobaric,
        ]:
            for e in pressure_stream:
                iscale.constraint_scaling_transform(pressure_stream[e], 1e-5)

    iscale.set_scaling_factor(blk.unit.feed_props[0.0].temperature, 1e-2)


def add_LTMED_costing(blk, costing_block=None):

    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def init_LTMED(blk):
    print("\n\n-------------------- INITIALIZING LT-MED --------------------\n\n")

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize()

    propagate_state(blk.unit_to_product)
    propagate_state(blk.unit_to_disposal)

    blk.product.initialize()
    blk.disposal.initialize()


def set_system_operating_conditions(m):
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(171.76)
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.6423)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(25 + 273.15)


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.liquid_prop)
    m.fs.disposal = Product(property_package=m.fs.liquid_prop)
    m.fs.product = Product(property_package=m.fs.liquid_prop)

    m.fs.LTMED = FlowsheetBlock(dynamic=False)
    build_LTMED(m.fs.LTMED, liquid_prop=m.fs.liquid_prop, vapor_prop=m.fs.vapor_prop)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.LTMED.feed.inlet,
    )
    m.fs.unit_to_disposal = Arc(
        source=m.fs.LTMED.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )
    m.fs.unit_to_product = Arc(
        source=m.fs.LTMED.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def init_system(m):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)

    init_LTMED(m.fs.LTMED)

    propagate_state(m.fs.unit_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.unit_to_disposal)
    m.fs.disposal.initialize()


def report_LTMED(blk, w=30):
    title = "LT-MED Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{f"Stream":<{w}}{f"Mass Flow Water (kg/s)":<{w}}{f"Mass Flow TDS (kg/s)":<{w}}'
    )
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Feed":<{w}}{value(blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<{w}.2f}{value(blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]):.2f}'
    )
    print(
        f'{"Product":<{w}}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<{w}.2f}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "TDS"]):.2f}'
    )
    print(
        f'{"Disposal":<{w}}{value(blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<{w}.2f}{value(blk.disposal.properties[0].flow_mass_phase_comp["Liq", "TDS"]):.2f}'
    )
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"STEC":<{w}}{value(blk.unit.specific_energy_consumption_thermal):<{w}.2f}kWh/m3'
    )
    print(
        f'{"GOR":<{w}}{value(blk.unit.gain_output_ratio):<{w}.2f}{pyunits.get_units(blk.unit.gain_output_ratio)}'
    )
    print(
        f'{"Thermal Power Req.":<{w}}{value(blk.unit.thermal_power_requirement):<{w},.2f}{pyunits.get_units(blk.unit.thermal_power_requirement)}'
    )
    print(
        f'{"Water Recovery":<{w}}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]/blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]) * 100:<{w}.2f}%'
    )
    print("\n\n")


def print_LTMED_costing_breakdown(blk, w=25):
    unit = blk.unit
    title = "LT-MED Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Capital Cost":<{w}s}{f"${unit.costing.capital_cost.value:<{w},.0f}"}{pyunits.get_units(unit.costing.capital_cost)}'
    )
    print(
        f'{"     Membrane System":<{w}s}{f"${unit.costing.membrane_system_cost.value:<{w},.0f}"}{pyunits.get_units(unit.costing.membrane_system_cost)}'
    )
    print(
        f'{"     Evaporator System":<{w}s}{f"${unit.costing.evaporator_system_cost.value:<{w},.0f}"}{pyunits.get_units(unit.costing.evaporator_system_cost)}'
    )
    print(
        f'{"Operating Cost":<{w}s}{f"${unit.costing.fixed_operating_cost.value:<{w},.0f}"}{pyunits.get_units(unit.costing.fixed_operating_cost)}'
    )


def main():

    m = build_system()
    set_system_operating_conditions(m)
    set_LTMED_operating_conditions(m.fs.LTMED)
    add_LTMED_costing(m.fs.LTMED)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.LTMED.product.properties[0].flow_vol)

    add_LTMED_scaling(m.fs.LTMED)
    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "TDS")
    )

    iscale.calculate_scaling_factors(m)

    assert degrees_of_freedom(m) == 0
    init_system(m)
    results = solve(m)
    assert_optimal_termination(results)

    report_LTMED(m.fs.LTMED)
    print_LTMED_costing_breakdown(m.fs.LTMED)

    return m


if __name__ == "__main__":
    m = main()
