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


def build_system():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # m.fs.costing = TreatmentCosting()
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
    display_dof(m, blk, where="After build")


def init_mec(m, blk, verbose=True, solver=None):
    pass


def set_system_operating_conditions(
    m,
    blk,
    Qin=2.8,  # MGD
    tds=12,  # g/L
    upstream_recovery=0.9,  # assumed
    saturated_steam_pressure_gage=3,
):
    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.gram / pyunits.liter
    tds = tds / (1 - upstream_recovery)
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
    )

    m.flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    m.flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(m.flow_mass_water)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(m.flow_mass_tds)
    m.fs.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)

    # m.fs.feed.properties.calculate_state(
    #     var_args={
    #         ("flow_mass_phase_comp", ("Liq", "H2O")): m.flow_mass_water,
    #         ("flow_mass_phase_comp", ("Liq", "NaCl")): m.flow_mass_tds,
    #         ("temperature", None): 298.15,
    #         ("pressure", None): 101325,
    #     },
    #     hold_state=True,
    # )

    m.fs.steam.properties[0].pressure_sat
    m.fs.steam.properties[0].dh_vap_mass
    m.fs.steam.properties.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    m.fs.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    display_dof(m, blk, where="After set_system_operating_conditions")


def set_mec_operating_conditions(
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
):
    rho = 1000 * pyunits.kg / pyunits.m**3
    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.gram / pyunits.liter
    tds = tds / (1 - upstream_recovery)

    mec = blk.unit
    assert len(operating_pressures) == mec.config.number_effects

    # m.flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    # m.flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)

    """
    Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """

    display_dof(m, blk, where="Beginning set_mec_operating_conditions")
    flow_mass_phase_water_per = rho / (rho + tds) * 1 * pyunits.kg / pyunits.s
    flow_mass_phase_salt_per = tds / (rho + tds) * 1 * pyunits.kg / pyunits.s

    # print(flow_mass_phase_water_per(), flow_mass_phase_salt_per())
    # assert False

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

    # first_effect.overall_heat_transfer_coefficient.fix(0.1)
    blk.steam.properties[0].pressure_sat
    blk.steam.properties[0].dh_vap_mass
    blk.steam.properties.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=False,
    )
    blk.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    # first_effect.heating_steam[0].display()
    # assert False
    # blk.steam.initialize()
    # propagate_state(blk.steam_to_unit)
    # blk.unit.steam.display()
    # blk.steam.outlet.display()
    # assert False


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
    # for n, eff in mec.effects.items():
    #     print(f"DOF effect {n}: {degrees_of_freedom(eff.effect)}")
    #     assert degrees_of_freedom(eff.effect) == 0
    display_dof(m, blk, where="After 1 kg/s flow rate ")
    # assert degrees_of_freedom(blk) == -3

    ### INITIALIZE FOR EACH EFFECT
    """
    Note: this is essentially to have an initial guess of the crysts,
    and to populate the feed concentration to all effects
    """
    for n, eff in mec.effects.items():
        eff.effect.initialize()

    display_dof(m, blk, where="After effect initialize")

    ### UNFIX THE INLET FLOW RATES OF EACH EFFECT
    """
    Note: this is to release the volumetric inlet flow entering different effects
    with the same salinity
    """
    for n, eff in mec.effects.items():
        if n > 1:
            conc_in = value(eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"])
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            # eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].set_value(conc_in)
            # eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(conc_in)
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()
            # if n != 1:
            #     eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    display_dof(m, blk, where="After unfixing effect flow mass, fixing conc mass")

    """
    Note: by this point, the multi-effect flow sheet should be fully specified (DOF=0),
    while inlet flows to effect 2-4 should be subject to the previous effect (DOF=1)
    """

    # mec.control_volume.properties_in[0].pressure.unfix()
    # mec.control_volume.properties_in[0].temperature.unfix()
    # display_dof(m, blk, where="After unfixing CV temp and pressure, before first MEC solve")
    display_dof(m, blk, where="before first MEC solve")

    results = solver.solve(mec)
    assert_optimal_termination(results)
    results = solver.solve(blk)
    assert_optimal_termination(results)
    display_mec_streams(m, blk)
    blk.feed.properties[0].display()
    # blk.unit.steam.display()
    blk.steam.properties[0].display()
    first_effect.heating_steam[0].display()

    # blk.steam.outlet.display()
    # assert False
    # assert False
    mec.control_volume.properties_in[0].pressure.unfix()
    mec.control_volume.properties_in[0].temperature.unfix()

    steam_state_dict = first_effect.heating_steam[0].define_port_members()
    display_dof(m, blk, where="After unfixing CV temp and pressure, after first MEC solve")
    # for k, v in steam_state_dict.items():
    #     # print(k, v)
    #     if k == "flow_mass_phase_comp":
    #         for i, vv in v.items():
    #             # print(k, i, vv, value(vv))
    #             vv.unfix()
    #             m.fs.steam.properties[0].flow_mass_phase_comp[i].set_value(value(vv))
    #     else:
    #         v.unfix()
    #         xx = getattr(m.fs.steam.properties[0], k)
    #         xx.set_value(value(v))
    # display_dof(m, blk, where="After unfixing first effect heating_steam state vars")
    # assert False
    # blk.heating_steam_flow_mass_water_vap = Constraint(
    #     expr=first_effect.heating_steam[0].flow_mass_phase_comp["Vap","H2O"] == blk.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"]
    # )
    # blk.heating_steam_flow_mass_water_liq = Constraint(
    #     expr=first_effect.heating_steam[0].flow_mass_phase_comp["Liq","H2O"] == blk.steam.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    # )
    # blk.heating_steam_press = Constraint(
    #     expr=first_effect.heating_steam[0].pressure == blk.steam.properties[0].pressure
    # )
    # blk.heating_steam_temp = Constraint(
    #     expr=first_effect.heating_steam[0].temperature == blk.steam.properties[0].temperature
    # )
    prop_in_state_dict = first_effect.properties_in[0].define_port_members()
    
    for k, v in prop_in_state_dict.items():
        # print(k, v)
        if k == "flow_mass_phase_comp":
            for i, vv in v.items():
                # print(k, i, vv)
                vv.unfix()
        else:
            v.unfix()
    steam_state_vars = blk.steam.properties[0].define_port_members()
    display_dof(m, blk, where="After unfixing first effect prop_in state vars")
    
    for k, v in steam_state_vars.items():
        # print(k, v)
        if k == "flow_mass_phase_comp":
            for i, vv in v.items():
                # print(k, i, vv)
                vv.unfix()
        else:
            v.unfix()
    display_dof(m, blk, where="After unfixing blk.steam state vars")
    # assert False

    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].setlb(m.flow_mass_water / 4)
    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].setlb(m.flow_mass_tds/ 4)
    blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].set_value(m.flow_mass_water)
    blk.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(m.flow_mass_tds)
    # display_dof(m, blk, where="After unfixing first effect flow_mass, fixing CV flow mass")
    # first_effect.heating_steam[0].display()
    
    # display_mec_streams(m, blk)
    # display_dof(m, blk)
    # assert False


    # mec.control_volume.initialize()

def set_mec_scaling(m, blk):
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



    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp",
    #     1 / value(m.flow_mass_water),
    #     index=("Liq", "H2O"),
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp",
    #     1 / value(m.flow_mass_tds),
    #     index=("Liq", "NaCl"),
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 10, index=("Vap", "H2O")
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Sol", "NaCl")
    # )
    # m.fs.vapor_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Vap", "H2O")
    # )
    # m.fs.vapor_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    # )

    # calculate_scaling_factors(blk)



def add_mec_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )

def display_dof(m, blk, where=None):
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

    fe = blk.unit.effects[1].effect # first effect

    print("\nm.fs.feed")
    print(f"\t{'Feed Mass Flow In Liq, H2O:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'Feed Mass Flow In Liq, TDS:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
    print(f"\t{'Feed Mass Flow In Sol, TDS:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
    print(f"\t{'Feed Mass Flow In Vap, H2O:':<50} {value(m.fs.feed.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nMEC Feed blk")
    print(f"\t{'MEC Feed Mass Flow In Liq, H2O:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'MEC Feed Mass Flow In Liq, TDS:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
    print(f"\t{'MEC Feed Mass Flow In Sol, TDS:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
    print(f"\t{'MEC Feed Mass Flow In Vap, H2O:':<50} {value(blk.feed.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    for n, eff in blk.unit.effects.items():
        print(f"\nMEC Effect {n}")
        print(f"\t{f'Effect {n} TDS Conc:':<50} {value(blk.unit.effects[n].effect.properties_in[0].conc_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/m3")
        print(f"\t{f'Effect {n} Mass Flow In Liq, H2O:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
        print(f"\t{f'Effect {n} Mass Flow In Liq, TDS:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
        print(f"\t{f'Effect {n} Mass Flow In Sol, TDS:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
        print(f"\t{f'Effect {n} Mass Flow In Vap, H2O:':<50} {value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
        if n == 1:
            print(f"\t{f'Effect {n} Mass Flow In Vap, H2O:':<50} {value(blk.unit.effects[n].effect.heating_steam[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nm.fs.steam")
    print(f"\t{'Steam Mass Flow In Liq, H2O:':<50} {value(m.fs.steam.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'Steam Mass Flow In Vap, H2O:':<50} {value(m.fs.steam.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nMEC Steam blk")
    print(f"\t{'MEC Steam Mass Flow In Liq, H2O:':<50} {value(blk.steam.properties[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'MEC Steam Mass Flow In Vap, H2O:':<50} {value(blk.steam.properties[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nMEC CV IN")
    print(f"\t{'MEC Control Volume Mass Flow In Liq, H2O:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow In Liq, TDS:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow In Sol, TDS:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow In Vap, H2O:':<50} {value(blk.unit.control_volume.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")
    print("\nMEC CV OUT")
    print(f"\t{'MEC Control Volume Mass Flow Out Liq, H2O:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']):<20.2f} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow Out Liq, TDS:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'NaCl']):<20.2f} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow Out Sol, TDS:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Sol', 'NaCl']):<20.2e} kg/s")
    print(f"\t{'MEC Control Volume Mass Flow Out Vap, H2O:':<50} {value(blk.unit.control_volume.properties_out[0].flow_mass_phase_comp['Vap', 'H2O']):<20.2f} kg/s")

    total_mass_flow_water = sum(value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]) for n, _ in blk.unit.effects.items())
    total_mass_flow_tds = sum(value(blk.unit.effects[n].effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]) for n, _ in blk.unit.effects.items())

    print(f"\n{'SUM EFFECT MASS FLOW WATER:':<50} {total_mass_flow_water}")
    print(f"\n{'SUM EFFECT MASS FLOW TDS:':<50} {total_mass_flow_tds}")

    print()

    display_dof(m, blk)

if __name__ == "__main__":
    solver = get_solver()
    m = build_system()
    blk = m.fs.MEC
    display_dof(m, m.fs.MEC)
    set_system_operating_conditions(m, m.fs.MEC)
    # display_dof(m, m.fs.MEC)
    set_mec_scaling(m, m.fs.MEC)
    m.fs.feed.initialize()
    # display_dof(m, m.fs.MEC)
    propagate_state(m.fs.feed_to_unit)
    # display_dof(m, m.fs.MEC)
    m.fs.MEC.feed.initialize()
    propagate_state(m.fs.MEC.feed_to_unit)
    m.fs.steam.initialize()
    
    propagate_state(m.fs.steam_to_unit)
    # display_dof(m, m.fs.MEC)
    m.fs.MEC.steam.initialize()
    propagate_state(m.fs.MEC.steam_to_unit)

    display_dof(m, m.fs.MEC, where="After initializing/propagating feed/steam blks and Arcs\nbefore set_mec_operating_conditions")

    set_mec_operating_conditions(m, m.fs.MEC)

    display_dof(m, blk, where="After set MEC operating conditions")

    # m.fs.steam.properties.display()
    # m.fs.MEC.steam.properties.display()
    # m.fs.MEC.unit.effects[1].effect.heating_steam.display()
    # m.fs.steam.initialize()
    # m.fs.steam.outlet.display()
    # m.fs.MEC.steam.initialize()
    # m.fs.MEC.steam.inlet.display()
    # m.fs.MEC.steam.outlet.display()
    # m.fs.MEC.unit.steam.display()
    # m.fs.steam.display()

    # m.fs.MEC.steam.properties[0].display()
    # m.fs.MEC.unit.effects[1].effect.heating_steam[0].display()
    # m.fs.MEC.steam.initialize()
    # m.fs.MEC.steam.properties[0].display()
    # m.fs.steam_to_unit.display()
    # m.fs.steam_to_unit_expanded.pprint()
    # m.fs.MEC.steam_to_unit_expanded.pprint()
    # assert False
    # display_dof(m, m.fs.MEC)
    
    propagate_state(m.fs.MEC.unit_to_product)

    # display_dof(m, m.fs.MEC)

    m.fs.MEC.product.initialize()
    propagate_state(m.fs.MEC.unit_to_solids)
    # display_dof(m, m.fs.MEC)
    m.fs.MEC.solids.initialize()

    # display_dof(m, blk, where="After unfixing first effect flow_mass, fixing CV flow mass")
    # assert False
    # m.fs.steam.display()
    # m.fs.MEC.unit.effects[1].effect.display()
    m.fs.feed.properties[0].display()
    m.fs.steam.properties[0].display()
    # blk.feed.initializecalculat
    blk.feed.initialize(state_args=m.fs.feed.properties[0].define_port_members()) 
    blk.feed.properties[0].flow_mass_phase_comp.display()
    # assert False
    blk_feed_state_vars = m.fs.feed.properties[0].define_port_members()
    for k, v in blk_feed_state_vars.items():
        # print(k, v)
        if k == "flow_mass_phase_comp":
            for i, vv in v.items():
                # print(k, i, vv, value(vv))
                vv.fix()
                blk.feed.properties[0].flow_mass_phase_comp[i].fix(value(vv))
                blk.unit.control_volume.properties_in[0].flow_mass_phase_comp[i].set_value(value(vv))
        else:
            v.fix()
            xx = getattr(blk.feed.properties[0], k)
            xx.fix(value(v))
            # xx.set_value(value(v))
    

    blk.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    blk.feed.initialize()
    blk.feed.properties.display()
    blk.steam.properties[0].flow_mass_phase_comp["Liq", "H2O"].setub(1e-5)
    blk.steam.initialize()
    blk.steam.properties.display()

    propagate_state(m.fs.MEC.feed_to_unit)
    propagate_state(m.fs.MEC.steam_to_unit)

    display_dof(m, m.fs.MEC)
    results = solver.solve(blk)
    assert_optimal_termination(results)
    blk.steam.properties[0].flow_mass_phase_comp["Liq", "H2O"].setub(None)
    m.fs.MEC.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].unfix()
    # blk.unit.effects[1].effect.properties_in[0].display()
    # print_infeasible_constraints(blk)
    display_mec_streams(m, blk)
    # assert False
    blk.feed.properties.display()
    blk.steam.properties.display()
    blk.unit.control_volume.properties_in[0].display()
    # blk.solids.display()
    # print(f"termination {results.solver.termination_condition}")
    # assert False

    # m.fs.MEC.unit.inlet.flow_mass_phase_comp["Sol", "NaCl"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].unfix()
    m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
    m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].set_value(value(blk.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"]))
    m.fs.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"].set_value(value(blk.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"]))
    blk.unit.effects[1].effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].set_value(value(blk.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"]))
    blk.feed.properties[0].temperature.unfix()
    blk.feed.properties[0].pressure.unfix()

    # m.fs.constr = Constraint(
    #     expr=m.fs.steam.properties[0].flow_mass_phase_comp["Vap", "H2O"] == m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"]
    # )
    m.fs.feed.properties.display()
    blk.feed.properties.display()
    m.fs.steam.properties.display()
    blk.steam.properties.display()
    blk.unit.control_volume.properties_in[0].display()
    display_dof(m, m.fs.MEC)
    # assert False
    # m.fs.steam.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0)
    display_dof(m, m.fs.MEC, where="before new scaling")

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(m.flow_mass_water),
        1e-2,
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(m.flow_mass_tds),
        1,
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

    calculate_scaling_factors(m)
    solver.options["max_iter"] = 5000
    solver.options["halt_on_ampl_error"] = "yes"
    # blk.heating_steam_flow_mass_water_vap.deactivate()
    display_dof(m, blk, where="BEFORE SOLVE")
    results = solver.solve(m, tee=False)
    display_mec_streams(m, blk)
    print(f"termination {results.solver.termination_condition}")
    # blk.heating_steam_flow_mass_water_vap.activate()
    # blk.unit.effects[1].effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
    #     value(blk.unit.effects[2].effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"])
    # )
    blk.unit.effects[1].effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    display_dof(m, blk, where="BEFORE SOLVE")
    # blk.unit.effects[1].effect.initialize()
    blk.unit.control_volume.properties_in[0].conc_mass_phase_comp
    m.fs.feed.properties[0].conc_mass_phase_comp
    results = solver.solve(m, tee=False)
    
    # results = solver.solve(m, tee=False)
    
    display_mec_streams(m, blk)
    blk.unit.control_volume.properties_in[0].conc_mass_phase_comp.display() 
    m.fs.feed.properties[0].conc_mass_phase_comp.display()
    # print_infeasible_constraints(m)
    print(f"termination {results.solver.termination_condition}")
    # blk.
    # assert_optimal_termination(results)
    # assert False
    # display_mec_streams(m, m.fs.MEC)
    # m.fs.feed.properties[0].flow_mass_phase_comp.display()
    # m.fs.MEC.feed.properties[0].flow_mass_phase_comp.display()
    # m.fs.MEC.unit.control_volume.properties_in[0].flow_mass_phase_comp.display()
    # m.fs.MEC.unit.effects[1].effect.properties_in[0].flow_mass_phase_comp.display()
    # m.fs.MEC.product.properties[0].flow_mass_phase_comp.display()
    # m.fs.MEC.unit.effects[1].effect.properties_solids[0].flow_mass_phase_comp.display()
    # # m.fs.MEC.solids.properties[0].flow_mass_phase_comp.display()
    # m.fs.MEC.unit.solids.display()
    # m.fs.MEC.solids.properties[0].display()
    # m.fs.MEC.unit.effects[1].effect.properties_in.display()

    # m.fs.MEC.unit.initialize()

    #####

    # set_mec_operating_conditions(m, m.fs.MEC)

    # print(f"Degrees of Freedom: {degrees_of_freedom(m)}")

    # set_mec_scaling(m, m.fs.MEC)
    # results = solver.solve(m)
    # assert_optimal_termination(results)
    # first_effect = m.fs.MEC.unit.effects[1].effect
    # mec = m.fs.MEC.unit
    # # first_effect.display()
    # # assert False

    # ### Release 1st effect flow rate and fix total flow rate instead
    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    # first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    # mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
    #     m.flow_mass_water
    # )
    # mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
    #     m.flow_mass_tds
    # )

    # """
    # Note: Rescaling is probably needed for extremely large feed flow,
    # """
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp",
    #     1 / value(m.flow_mass_water),
    #     index=("Liq", "H2O"),
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp",
    #     1 / value(m.flow_mass_tds),
    #     index=("Liq", "NaCl"),
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 10, index=("Vap", "H2O")
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Sol", "NaCl")
    # )
    # m.fs.vapor_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Vap", "H2O")
    # )
    # m.fs.vapor_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    # )

    # calculate_scaling_factors(m)
    # print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
    # results = solver.solve(m)
    # assert_optimal_termination(results)
    # print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
    # mec.effects[1].effect.properties_vapor.display()
