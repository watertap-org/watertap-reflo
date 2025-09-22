from pyomo.environ import (
    ConcreteModel,
    Expression,
    Objective,
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

from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock as CrystParameterBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)

from watertap_contrib.reflo.unit_models import MultiEffectCrystallizer
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_MEC",
    "set_MEC_op_conditions",
    "add_MEC_costing",
    "init_MEC",
    "scale_MEC",
    "rescale_MEC",
    "report_MEC",
    "print_MEC_costing_breakdown",
]


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.costing = TreatmentCosting()
    m.fs.cryst_properties = CrystParameterBlock()
    m.fs.vapor_properties = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.cryst_properties)
    m.fs.product = Product(property_package=m.fs.cryst_properties)

    m.fs.mec = FlowsheetBlock()

    build_MEC(m.fs.mec)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.mec.feed.inlet,
    )
    m.fs.unit_to_product = Arc(
        source=m.fs.mec.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_MEC(blk, prop_package=None, vapor_prop_package=None):
    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.cryst_properties
    if vapor_prop_package is None:
        m = blk.model()
        vapor_prop_package = m.fs.vapor_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)

    blk.unit = MultiEffectCrystallizer(
        property_package=prop_package, property_package_vapor=vapor_prop_package
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )
    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_MEC_op_conditions(
    blk,
    operating_pressures=[0.45, 0.25, 0.208, 0.095],
    feed_H2O=111.14146,
    feed_NaCl=28.4782,
    nacl_yield=0.9,
    heat_transfer_coeff=1300,
):
    """
    For the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """

    # Guessed values for initialization
    feed_pressure = 101325
    feed_temperature = 273.15 + 20

    # Total going in to MEC
    flow_mass_phase_water_total = feed_H2O
    flow_mass_phase_salt_total = feed_NaCl

    # Total into each effect initial
    flow_mass_phase_water_per = (
        flow_mass_phase_water_total
        / (flow_mass_phase_water_total + flow_mass_phase_salt_total)
        * pyunits.kg
        / pyunits.s
    )
    flow_mass_phase_salt_per = (
        flow_mass_phase_salt_total
        / (flow_mass_phase_water_total + flow_mass_phase_salt_total)
        * pyunits.kg
        / pyunits.s
    )

    saturated_steam_pressure = 101325 * pyunits.Pa + pyunits.convert(
        3 * pyunits.bar, to_units=pyunits.Pa
    )

    # Fix unit model parameters
    for (_, eff), op_pressure in zip(blk.unit.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            flow_mass_phase_water_per
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            flow_mass_phase_salt_per
        )

        eff.effect.properties_in[0].pressure.fix(feed_pressure)
        eff.effect.properties_in[0].temperature.fix(feed_temperature)

        eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        # We will need this later...
        eff.effect.properties_in[0].conc_mass_phase_comp[...]

        eff.effect.crystallization_yield["NaCl"].fix(nacl_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(heat_transfer_coeff)

    first_effect = blk.unit.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(heat_transfer_coeff)
    first_effect.heating_steam[0].pressure_sat
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # Fix control volume properties except the liquid mass flow rates
    # The liquid phase flow rates are the sum of the flow rates in all effects
    blk.unit.control_volume.properties_in[0].pressure.fix(feed_pressure)
    blk.unit.control_volume.properties_in[0].temperature.fix(feed_temperature)
    blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    for _, eff in blk.unit.effects.items():
        assert degrees_of_freedom(eff.effect) == 0

    blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].unfix()


def scale_MEC(blk, cryst_prop_pack=None, vapor_prop_pack=None):

    if cryst_prop_pack is None:
        m = blk.model()
        cryst_prop_pack = m.fs.cryst_properties
    if vapor_prop_pack is None:
        m = blk.model()
        vapor_prop_pack = m.fs.vapor_properties

    m = blk.model()

    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )
    vapor_prop_pack.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    vapor_prop_pack.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    for _, eff in blk.unit.effects.items():
        iscale.set_scaling_factor(
            eff.effect.properties_solids[0].flow_vol_phase["Vap"], 1e5
        )
        iscale.set_scaling_factor(
            eff.effect.properties_out[0].flow_vol_phase["Vap"], 1e5
        )


def rescale_MEC(
    blk,
    flow_mass_phase_water_total=111.1414,
    flow_mass_phase_salt_total=28.478,
    cryst_prop_pack=None,
    vapor_prop_pack=None,
):
    """
    Rescaling may be needed for extremely large feed flow rates
    """

    if cryst_prop_pack is None:
        m = blk.model()
        cryst_prop_pack = m.fs.cryst_properties
    if vapor_prop_pack is None:
        m = blk.model()
        vapor_prop_pack = m.fs.vapor_properties

    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
        index=("Liq", "NaCl"),
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Vap", "H2O"),
    )
    cryst_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
        index=("Sol", "NaCl"),
    )
    vapor_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Vap", "H2O"),
    )
    vapor_prop_pack.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )


def add_MEC_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
        costing_method_arguments={"cost_work_as": "heat"},
    )


def set_system_op_conditions(m, feed_H2O=111.14146, feed_NaCl=28.4782):

    m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].set_value(feed_H2O)
    m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].set_value(feed_NaCl)
    m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].unfix()
    m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].set_value(0)
    # m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].unfix()
    # m.fs.mec.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].set_value(0)
    m.fs.mec.unit.inlet.temperature[0].set_value(273.15 + 30.51)
    m.fs.mec.unit.inlet.pressure[0].set_value(101325)

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_H2O)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_NaCl)
    m.fs.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    # feed steam flow is dictated by unit.heating_steam[0] so we don't fix the inlet steam flow here
    # m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
    m.fs.feed.properties[0].temperature.fix(273.15 + 30.51)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].conc_mass_phase_comp[...]

    set_MEC_op_conditions(m.fs.mec, feed_H2O=feed_H2O, feed_NaCl=feed_NaCl)


def init_MEC(blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    # Initialize each effect sequentially
    # NOTE: DO NOT use blk.unit.initialize() here
    for n, eff in blk.unit.effects.items():
        eff.effect.initialize()
        if n > 1:
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    results = solve(blk)
    assert_optimal_termination(results)

    # Release 1st effect flow rate
    first_effect = blk.unit.effects[1].effect
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    blk.unit.inlet.temperature[0].unfix()
    blk.unit.inlet.pressure[0].unfix()

    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def init_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)
    init_MEC(m.fs.mec)

    propagate_state(m.fs.unit_to_product)
    m.fs.product.initialize()


def report_MEC(blk, w=25):
    title = "MEC Report"
    side = int(((4 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    conv = pyunits.convert(
        1 * pyunits.m**3 / pyunits.s, to_units=pyunits.Mgallons / pyunits.day
    )()
    print(f"\n{header}\n")
    print(f'{f"Phase, Component":<{w}}{f"Inlet (kg/s)":<{w}}{f"Outlet (kg/s)":<{w}}')
    print(f"{'-' * (3 * w)}")
    pcs = [(0, "Liq", "H2O"), (0, "Liq", "NaCl"), (0, "Sol", "NaCl"), (0, "Vap", "H2O")]
    for pc in pcs:
        print(
            f'{f"{pc[1]}, {pc[2]}":<{w}}{value(blk.unit.inlet.flow_mass_phase_comp[pc]):<{w}.2f}{value(blk.unit.outlet.flow_mass_phase_comp[pc]):<{w}.2f}'
        )

    print(
        f'\n{f"Node":<{w}}{f"Water @ Inlet (MGD)":<{w}}{f"Water @ Outlet (MGD)":<{w}}{f"Solids (kg/s)":<{w}}{f"Steam (kg/s)":<{w}}'
    )
    print(f"{'-' * (5 * w)}")
    cv = blk.unit.control_volume
    flow_vol_out = value(
        blk.unit.recovery_vol_phase["Liq"] * blk.unit.total_flow_vol_in
    )
    print(
        f'{f"System Total":<{w}}{value(blk.unit.total_flow_vol_in)*conv:<{w}.2f}{flow_vol_out*conv:<{w}.2f}{value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp["Sol", "NaCl"]):<{w}.2f}{value(cv.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]):<{w}.2f}'
    )
    for n, eff in blk.unit.effects.items():
        print(
            f'{f"Effect {n}":<{w}}{value(eff.effect.properties_in[0].flow_vol_phase["Liq"])*conv:<{w}.2f}{value(eff.effect.properties_pure_water[0].flow_vol_phase["Liq"])*conv:<{w}.2f}{value(eff.effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"]):<{w}.2f}{value(eff.effect.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]):<{w}.2f}'
        )
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{f"Inlet Salinity":<{w}}{value(blk.unit.effects[1].effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]):<{w}.2f}{pyunits.get_units(blk.unit.effects[1].effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"])}'
    )

    for n, eff in blk.unit.effects.items():
        print(
            f'{f"Effect {n} HX Area":<{w}}{value(eff.effect.heat_exchanger_area):<{w}.2f}{pyunits.get_units(eff.effect.heat_exchanger_area)}'
        )
        print(
            f'{f"Effect {n} Work":<{w}}{value(eff.effect.work_mechanical[0]):<{w}.2f}{pyunits.get_units(eff.effect.work_mechanical[0])}'
        )


def print_MEC_costing_breakdown(blk, w=30):
    title = "MEC Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"MEC Capital Cost":<{w}s}{f"${value(blk.unit.costing.capital_cost):<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost)}'
    )
    print("\n\n")


def main():

    # feed_H2O = 111.14146762116363 # kg/s
    # feed_NaCl = 28.47821365277779 # kg/s

    m = build_system()

    set_system_op_conditions(m)
    scale_MEC(m.fs.mec)
    iscale.calculate_scaling_factors(m)
    init_system(m)

    # rescale_MEC(m.fs.mec)
    # iscale.calculate_scaling_factors(m)
    add_MEC_costing(m.fs.mec)

    m.fs.costing.cost_process()

    flow_out = pyunits.convert(
        m.fs.mec.unit.recovery_vol_phase["Liq"] * m.fs.mec.unit.total_flow_vol_in,
        to_units=pyunits.m**3 / pyunits.hr,
    )
    m.fs.costing.add_LCOW(flow_out)
    m.fs.costing.add_specific_energy_consumption(flow_out, name="SEC")
    m.fs.costing.SEC_th = Expression(
        expr=pyunits.convert(
            m.fs.costing.aggregate_flow_heat / flow_out,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    assert degrees_of_freedom(m) == 0
    results = solve(m)
    assert_optimal_termination(results)
    report_MEC(m.fs.mec)
    print_MEC_costing_breakdown(m.fs.mec)

    return m


if __name__ == "__main__":
    m = main()
