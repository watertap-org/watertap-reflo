import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

# from analysisWaterTAP.utils.flowsheet_utils import *
# from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
#     calculate_operating_pressure,
# )

# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
    NaClStateBlock,
)
from watertap_contrib.reflo.unit_models.multi_effect_crystallizer import (
    MultiEffectCrystallizer,
)

from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect

rho = 1000 * pyunits.kg / pyunits.m**3
feed_pressure = 101325 * pyunits.Pa
atm_pressure = 101325 * pyunits.Pa
feed_temperature = 273.15 + 20

__all__ = [
    "build_mec",
    "init_mec",
    "display_mec_streams",
    "set_mec_scaling",
    "set_mec_initial_scaling",
    "display_mec_dof",
]


def build_system():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.solids = Product(property_package=m.fs.properties)

    m.fs.MEC = mec = FlowsheetBlock(dynamic=False)

    build_mec(m, m.fs.MEC)

    m.fs.feed_to_unit = Arc(source=m.fs.feed.outlet, destination=mec.unit.inlet)

    m.fs.mec_to_product = Arc(source=mec.product.outlet, destination=m.fs.product.inlet)

    m.fs.mec_to_solids = Arc(source=mec.solids.outlet, destination=m.fs.solids.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mec(
    m,
    blk,
    prop_package=None,
    prop_package_vapor=None,
    number_effects=4,
):

    if prop_package is None:
        prop_package = m.fs.properties
    if prop_package_vapor is None:
        prop_package_vapor = m.fs.vapor_properties

    blk.product = StateJunction(property_package=prop_package)
    blk.solids = StateJunction(property_package=prop_package)

    blk.unit = MultiEffectCrystallizer(
        property_package=prop_package,
        property_package_vapor=prop_package_vapor,
        number_effects=number_effects,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    blk.unit_to_solids = Arc(
        source=blk.unit.solids,
        destination=blk.solids.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_system_operating_conditions(
    m,
    blk,
    Qin=2.8,  # MGD
    tds=12,  # g/L
    operating_pressures=[0.45, 0.25, 0.208, 0.095],
    upstream_recovery=0.9,  # assumed
    crystallizer_yield=0.5,
    saturated_steam_pressure_gage=3,
    heat_transfer_coefficient=0.1,
    eps=1e-8,
    **kwargs,
):
    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.gram / pyunits.liter
    tds = tds / (1 - upstream_recovery)
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
    )

    m.flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    m.flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)
    m.upstream_recovery = upstream_recovery
    m.operating_pressures = operating_pressures
    m.crystallizer_yield = crystallizer_yield
    m.heat_transfer_coefficient = heat_transfer_coefficient
    m.saturated_steam_pressure = saturated_steam_pressure
    m.saturated_steam_pressure_gage = saturated_steam_pressure_gage

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(m.flow_mass_water)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(m.flow_mass_tds)
    m.fs.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)


def init_system(m, blk):

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)

    init_mec(m, blk)

    propagate_state(m.fs.MEC.unit_to_product)
    m.fs.MEC.product.initialize()

    propagate_state(m.fs.mec_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.MEC.unit_to_solids)
    m.fs.MEC.solids.initialize()

    propagate_state(m.fs.mec_to_solids)
    m.fs.solids.initialize()

    add_mec_costing(m, m.fs.MEC)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.initialize()


def init_mec(m, blk, feed_props=None, verbose=True, solver=None):

    if feed_props is None:
        feed_props = m.fs.feed.properties[0]

    feed_props.conc_mass_phase_comp
    feed_props.flow_vol_phase
    feed_props.parent_block().initialize()

    tds = feed_props.conc_mass_phase_comp["Liq", "NaCl"]
    # tds = 120 * pyunits.kg / pyunits.m**3

    mec = blk.unit
    assert len(m.operating_pressures) == mec.config.number_effects

    set_mec_initial_scaling(m, blk)

    """
    Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """

    flow_mass_phase_water_per = rho / (rho + tds) * 1 * pyunits.kg / pyunits.s
    flow_mass_phase_salt_per = tds / (rho + tds) * 1 * pyunits.kg / pyunits.s
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        m.saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
    )

    ### FIX UNIT MODEL PARAMETERS
    for (_, eff), op_pressure in zip(mec.effects.items(), m.operating_pressures):

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
        eff.effect.properties_in[0].conc_mass_phase_comp[...]

        eff.effect.crystallization_yield["NaCl"].fix(m.crystallizer_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(m.heat_transfer_coefficient)

    first_effect = mec.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(m.heat_transfer_coefficient)
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

    ### FIX CV PROPERTIES EXCEPT FOR THE LIQUID FLOW RATES
    mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
    mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    for n, eff in mec.effects.items():
        eff.effect.initialize()

    ### UNFIX THE INLET FLOW RATES OF EACH EFFECT
    for n, eff in mec.effects.items():
        if n > 1:
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    solver = get_solver()
    results = solver.solve(mec)
    assert_optimal_termination(results)
    results = solver.solve(blk)
    assert_optimal_termination(results)

    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        m.flow_mass_water
    )
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        m.flow_mass_tds
    )

    set_mec_scaling(m, blk)

    results = solver.solve(mec)
    assert_optimal_termination(results)
    results = solver.solve(blk)
    assert_optimal_termination(results)

    for k, v in blk.unit.control_volume.properties_in[0].define_port_members().items():
        if k == "flow_mass_phase_comp":
            for i, vv in v.items():
                vv.unfix()
        else:
            v.unfix()


def set_mec_initial_scaling(m, blk):
    """
    Note:
    Keep in mind that we assumes feed flow rate to the 1st effect to be 1 kg/s,
    and the total feed flow rate is to be determined
    """
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )
    m.fs.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    m.fs.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )


def set_mec_scaling(m, blk):

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.flow_mass_water),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.flow_mass_tds),
        index=("Liq", "NaCl"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 10, index=("Vap", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Sol", "NaCl")
    )
    m.fs.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Vap", "H2O")
    )
    m.fs.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    calculate_scaling_factors(blk)


def add_mec_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )


def display_mec_dof(m, blk, where=None):
    if where is not None:
        print(f"\n{where}")
    else:
        print()
    for n, eff in blk.unit.effects.items():
        print(f"DOF effect {n}: {degrees_of_freedom(eff.effect)}")
    print(f"Degrees of Freedom model: {degrees_of_freedom(m)}")
    print(f"Degrees of Freedom MEC blk: {degrees_of_freedom(blk)}")
    print(f"Degrees of Freedom MEC unit: {degrees_of_freedom(blk.unit)}")


def display_mec_streams(m, blk):

    fe = blk.unit.effects[1].effect  # first effect

    print("\nm.fs.feed")
    print(
        f"\t{'Feed Mass Flow In Liq, H2O:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s"
    )
    print(
        f"\t{'Feed Mass Flow In Liq, TDS:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s"
    )
    print(
        f"\t{'Feed Mass Flow In Sol, TDS:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s"
    )
    print(
        f"\t{'Feed Mass Flow In Vap, H2O:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s"
    )
    # print("\nMEC Feed blk")
    # print(f"\t{'MEC Feed Mass Flow In Liq, H2O:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    # print(f"\t{'MEC Feed Mass Flow In Liq, TDS:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
    # print(f"\t{'MEC Feed Mass Flow In Sol, TDS:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
    # print(f"\t{'MEC Feed Mass Flow In Vap, H2O:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    for n, eff in blk.unit.effects.items():
        print(f"\nMEC Effect {n}")
        print(
            f"\t{f'Effect {n} TDS Conc:':<50} {value(blk.unit.effects[n].effect.properties_in[0].conc_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/m3"
        )
        print(
            f"\t{f'Effect {n} Mass Flow In Liq, H2O:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s"
        )
        print(
            f"\t{f'Effect {n} Mass Flow In Liq, TDS:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s"
        )
        print(
            f"\t{f'Effect {n} Mass Flow In Sol, TDS:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s"
        )
        print(
            f"\t{f'Effect {n} Mass Flow In Vap, H2O:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s"
        )
        if n == 1:
            print(
                f"\t{f'Effect {n} Mass Flow In Vap, H2O:':<50} {value(blk.unit.effects[n].effect.heating_steam[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s"
            )
    # print("\nm.fs.steam")
    # print(f"\t{'Steam Mass Flow In Liq, H2O:':<50} {value(m.fs.steam.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    # print(f"\t{'Steam Mass Flow In Vap, H2O:':<50} {value(m.fs.steam.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    # print("\nMEC Steam blk")
    # print(f"\t{'MEC Steam Mass Flow In Liq, H2O:':<50} {value(blk.steam.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    # print(f"\t{'MEC Steam Mass Flow In Vap, H2O:':<50} {value(blk.steam.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nMEC CV IN")
    print(
        f"\t{'MEC Control Volume Mass Flow In Liq, H2O:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow In Liq, TDS:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow In Sol, TDS:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow In Vap, H2O:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s"
    )
    print("\nMEC CV OUT")
    print(
        f"\t{'MEC Control Volume Mass Flow Out Liq, H2O:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow Out Liq, TDS:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow Out Sol, TDS:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s"
    )
    print(
        f"\t{'MEC Control Volume Mass Flow Out Vap, H2O:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s"
    )

    total_mass_flow_water = sum(
        value(
            blk.unit.effects[n]
            .effect.properties_in[0]
            .flow_mass_phase_comp["Liq", "H2O"]
        )
        for n, _ in blk.unit.effects.items()
    )
    total_mass_flow_tds = sum(
        value(
            blk.unit.effects[n]
            .effect.properties_in[0]
            .flow_mass_phase_comp["Liq", "NaCl"]
        )
        for n, _ in blk.unit.effects.items()
    )

    print(f"\n{'SUM EFFECT MASS FLOW WATER:':<50} {total_mass_flow_water}")
    print(f"\n{'SUM EFFECT MASS FLOW TDS:':<50} {total_mass_flow_tds}")

    print()

    display_mec_dof(m, blk)


if __name__ == "__main__":
    solver = get_solver()
    m = build_system()
    blk = m.fs.MEC

    set_system_operating_conditions(m, blk)
    init_system(m, blk)

    results = solver.solve(m)
    assert_optimal_termination(results)
    display_mec_streams(m, blk)
    m.fs.costing.LCOW.display()

# def set_mec_operating_conditions(
#     m,
#     blk,
#     Qin=2.8,  # MGD
#     tds=12,  # g/L
#     operating_pressures=[0.45, 0.25, 0.208, 0.095],
#     upstream_recovery=0.9,  # assumed
#     crystallizer_yield=0.5,
#     saturated_steam_pressure_gage=3,
#     heat_transfer_coefficient=0.1,
#     eps=1e-8,
#     **kwargs
# ):
#     rho = 1000 * pyunits.kg / pyunits.m**3
#     Qin = Qin * pyunits.Mgallons / pyunits.day
#     tds = tds * pyunits.gram / pyunits.liter
#     tds = tds / (1 - upstream_recovery)

#     mec = blk.unit
#     assert len(operating_pressures) == mec.config.number_effects

#     set_mec_initial_scaling(m, blk)

#     # m.flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
#     # m.flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)

#     """
#     Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
#     which is align to the default value, in order to guarantee a solution in the initial solve.
#     """

#     display_mec_dof(m, blk, where="Beginning set_mec_operating_conditions")
#     flow_mass_phase_water_per = rho / (rho + tds) * 1 * pyunits.kg / pyunits.s
#     flow_mass_phase_salt_per = tds / (rho + tds) * 1 * pyunits.kg / pyunits.s

#     saturated_steam_pressure = atm_pressure + pyunits.convert(
#         saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
#     )

#     ### FIX UNIT MODEL PARAMETERS
#     for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):

#         eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
#             flow_mass_phase_water_per
#         )
#         eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
#             flow_mass_phase_salt_per
#         )

#         eff.effect.properties_in[0].pressure.fix(feed_pressure)
#         eff.effect.properties_in[0].temperature.fix(feed_temperature)

#         eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
#         eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
#         eff.effect.properties_in[0].conc_mass_phase_comp[...]

#         eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
#         eff.effect.crystal_growth_rate.fix()
#         eff.effect.souders_brown_constant.fix()
#         eff.effect.crystal_median_length.fix()

#         eff.effect.pressure_operating.fix(
#             pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
#         )
#         eff.effect.overall_heat_transfer_coefficient.fix(heat_transfer_coefficient)

#     first_effect = mec.effects[1].effect

#     first_effect.overall_heat_transfer_coefficient.fix(heat_transfer_coefficient)
#     first_effect.heating_steam[0].pressure_sat
#     first_effect.heating_steam[0].dh_vap_mass
#     first_effect.heating_steam.calculate_state(
#         var_args={
#             ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
#             ("pressure", None): saturated_steam_pressure,
#             ("pressure_sat", None): saturated_steam_pressure,
#         },
#         hold_state=True,
#     )
#     first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

#     ### FIX CV PROPERTIES EXCEPT FOR THE LIQUID FLOW RATES
#     mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
#     mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
#     mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)


#     """
#     Check DOF
#     By this point, each effect is fully specified, so their DOF should be 0.
#     However, the multi-effect flowsheet is over-constrianted by energy_flow_constr,
#     which connect energy flow between effects, and the DOF should be negative (n_effects - 1)
#     """

#     display_mec_dof(m, blk, where="After 1 kg/s flow rate ")
#     # assert degrees_of_freedom(blk) == -3

#     ### INITIALIZE FOR EACH EFFECT
#     """
#     Note: this is essentially to have an initial guess of the crysts,
#     and to populate the feed concentration to all effects
#     """
#     for n, eff in mec.effects.items():
#         eff.effect.initialize()

#     display_mec_dof(m, blk, where="After effect initialize")

#     ### UNFIX THE INLET FLOW RATES OF EACH EFFECT
#     """
#     Note: this is to release the volumetric inlet flow entering different effects
#     with the same salinity
#     """
#     for n, eff in mec.effects.items():
#         if n > 1:
#             eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
#             eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
#             eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

#         eff.effect.properties_in[0].conc_mass_phase_comp.display()

#     display_mec_dof(m, blk, where="After unfixing effect flow mass, fixing conc mass")

#     """
#     Note: by this point, the multi-effect flow sheet should be fully specified (DOF=0),
#     while inlet flows to effect 2-4 should be subject to the previous effect (DOF=1)
#     """


#     display_mec_dof(m, blk, where="before first MEC solve")
#     # assert False
#     results = solver.solve(mec)
#     assert_optimal_termination(results)
#     display_mec_streams(m, blk)
#     results = solver.solve(blk)
#     assert_optimal_termination(results)
#     display_mec_streams(m, blk)
#     # assert False

#     # mec.control_volume.properties_in[0].pressure.unfix()
#     # mec.control_volume.properties_in[0].temperature.unfix()


#     first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
#     first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
#     mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
#         m.flow_mass_water
#     )
#     mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
#         m.flow_mass_tds
#     )
#     display_mec_dof(m, blk, where="After unfixing first effect flow_mass, fixing CV flow mass")

#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp",
#         1 / value(m.flow_mass_water),
#         index=("Liq", "H2O"),
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp",
#         1 / value(m.flow_mass_tds),
#         index=("Liq", "NaCl"),
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 10, index=("Vap", "H2O")
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-2, index=("Sol", "NaCl")
#     )
#     m.fs.vapor_properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-2, index=("Vap", "H2O")
#     )
#     m.fs.vapor_properties.set_default_scaling(
#         "flow_mass_phase_comp", 1, index=("Liq", "H2O")
#     )

#     calculate_scaling_factors(m)
#     results = solver.solve(mec)
#     assert_optimal_termination(results)
#     display_mec_streams(m, blk)
#     results = solver.solve(blk)
#     assert_optimal_termination(results)
#     display_mec_streams(m, blk)
#     # assert False
#     prop_in_state_dict = blk.unit.control_volume.properties_in[0].define_port_members()
#     for k, v in prop_in_state_dict.items():
#         # print(k, v)
#         if k == "flow_mass_phase_comp":
#             for i, vv in v.items():
#                 # print(k, i, vv)
#                 vv.unfix()
#         else:
#             v.unfix()

#     display_mec_dof(m, blk, where="After unfixing blk.steam state vars")
#     results = solver.solve(m)
#     assert_optimal_termination(results)
#     display_mec_streams(m, blk)
#     assert False

#     first_effect.properties_in[0].display()
#     assert False
#     # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
#     # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
#     # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].setlb(m.flow_mass_water / 4)
#     # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].setlb(m.flow_mass_tds/ 4)
#     blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].set_value(m.flow_mass_water)
#     blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(m.flow_mass_tds)


#     # display_mec_dof(m, blk, where="After unfixing first effect flow_mass, fixing CV flow mass")
#     # first_effect.heating_steam[0].display()

#     display_mec_streams(m, blk)
#     # display_mec_dof(m, blk)


#     # mec.control_volume.initialize()
