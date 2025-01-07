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
from idaes.core.util.initialization import propagate_state as _prop_state
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


def build_system():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.steam = Feed(property_package=m.fs.vapor_properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.solids = Product(property_package=m.fs.properties)

    m.fs.MEC = mec = FlowsheetBlock(dynamic=False)

    build_mec(m, m.fs.MEC)

    m.fs.feed_to_unit = Arc(source=m.fs.feed.outlet, destination=mec.feed.inlet)

    m.fs.steam_to_unit = Arc(source=m.fs.steam.outlet, destination=mec.steam.inlet)

    m.fs.mec_to_product = Arc(source=mec.product.outlet, destination=m.fs.product.inlet)

    m.fs.mec_to_solids = Arc(source=mec.solids.outlet, destination=m.fs.solids.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mec(
    m,
    blk,
    property_package=None,
    property_package_vapor=None,
    number_effects=4,
):

    if property_package is None:
        property_package = m.fs.properties
    if property_package_vapor is None:
        property_package_vapor = m.fs.vapor_properties

    blk.feed = StateJunction(property_package=property_package)
    blk.steam = StateJunction(property_package=property_package_vapor)
    blk.product = StateJunction(property_package=property_package)
    blk.solids = StateJunction(property_package=property_package)

    blk.unit = MultiEffectCrystallizer(
        property_package=m.fs.properties,
        property_package_vapor=m.fs.vapor_properties,
        number_effects=number_effects,
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.steam_to_unit = Arc(source=blk.steam.outlet, destination=blk.unit.steam)

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    blk.unit_to_solids = Arc(
        source=blk.unit.solids,
        destination=blk.solids.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def init_mec(m, blk, verbose=True, solver=None):
    pass


def set_mec_operating_conditions(
    m,
    blk,
    Qin=4,  # MGD
    tds=12,  # g/L
    operating_pressures = [0.45, 0.25, 0.208, 0.095], 
    upstream_recovery=0.9,  # assumed
    crystallizer_yield=0.5,
    saturated_steam_pressure_gage=3,
    heat_transfer_coefficient=0.1,
    eps=1e-8,
):
    rho = 1000 * pyunits.kg / pyunits.m**3
    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.gram / pyunits.liter
    tds = tds / (1 - upstream_recovery)

    mec = blk.unit
    assert len(operating_pressures) == mec.config.number_effects

    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)
    """
    Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """
    flow_mass_phase_water_per = rho / (rho + tds) * 1 * pyunits.kg / pyunits.s
    flow_mass_phase_salt_per = tds / (rho + tds) * 1 * pyunits.kg / pyunits.s
    
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
    )

    ### FIX UNIT MODEL PARAMETERS
    for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):

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

        eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(heat_transfer_coefficient)

    first_effect = mec.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(0.1)
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

    # blk.steam.properties[0].pressure_sat
    # blk.steam.properties[0].dh_vap_mass
    # blk.steam.properties.calculate_state(
    #     var_args={
    #         ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
    #         ("pressure", None): saturated_steam_pressure,
    #         ("pressure_sat", None): saturated_steam_pressure,
    #     },
    #     hold_state=True,
    # )
    # blk.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    


    ### FIX CV PROPERTIES EXCEPT FOR THE LIQUID FLOW RATES
    mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
    mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    """
    Check DOF
    By this point, each effect is fully specified, so their DOF should be 0.
    However, the multi-effect flowsheet is over-constrianted by energy_flow_constr, 
    which connect energy flow between effects, and the DOF should be negative (n_effects - 1)
    """
    for n, eff in mec.effects.items():
        print(degrees_of_freedom(eff.effect))
        if n != 1:
            eff.display()
            assert False
        # assert degrees_of_freedom(eff.effect) == 0
    assert degrees_of_freedom(m) == -3

if __name__ == "__main__":
    m = build_system()
    set_mec_operating_conditions(m, m.fs.MEC)
    print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
    # m.fs.MEC.unit.inlet.display()
