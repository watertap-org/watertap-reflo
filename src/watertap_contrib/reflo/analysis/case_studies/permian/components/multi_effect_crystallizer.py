import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from watertap.core.util.model_diagnostics import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver

from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
    NaClStateBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock,
    WaterStateBlock,
)
from watertap_contrib.reflo.unit_models.multi_effect_crystallizer import (
    MultiEffectCrystallizer,
)
from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

__all__ = [
    "build_system",
    "build_mec",
    "set_mec_op_conditions",
    "add_mec_costing",
    "init_mec",
    "unfix_mec",
    "mec_rescaling",
]

def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = TreatmentCosting()
    build_mec(m, m.fs)

    return m

def build_mec(m, blk) -> None:

    blk.properties = NaClParameterBlock()
    blk.vapor_properties = WaterParameterBlock()

    blk.unit = MultiEffectCrystallizer(
        property_package=blk.properties, property_package_vapor=blk.vapor_properties
    )

def set_mec_op_conditions(m, 
                          blk,
                          operating_pressures = [0.4455, 0.2758, 0.1651, 0.095],
                          nacl_yield = 0.8,
                          ) -> None :
    
    mec = blk.unit

    # Guessed values for initialization
    flow_in = 3.5
    rho = 1000 * pyunits.kg / pyunits.m**3
    conc_in = 160 * pyunits.g / pyunits.L
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    ### TOTAL GOING INTO MEC
    flow_vol_in = pyunits.convert(
        flow_in * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    flow_mass_phase_water_total = pyunits.convert(
        flow_vol_in * rho, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_salt_total = pyunits.convert(
        flow_vol_in * conc_in, to_units=pyunits.kg / pyunits.s
    )
    ### TOTAL INTO EACH EFFECT INITIAL
    """
    Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """
    # flow_mass_phase_water_per = 116.2473764168908 / 100 * pyunits.kg / pyunits.s
    # flow_mass_phase_salt_per = 28.478213652777765 / 100 * pyunits.kg / pyunits.s
    flow_mass_phase_water_per = 116 / 100 * pyunits.kg / pyunits.s
    flow_mass_phase_salt_per = 28 / 100 * pyunits.kg / pyunits.s

    saturated_steam_pressure = 101325 * pyunits.Pa + pyunits.convert(
        3 * pyunits.bar, to_units=pyunits.Pa
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

        eff.effect.crystallization_yield["NaCl"].fix(nacl_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(0.1)

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

    ### FIX CV PROPERTIES EXCEPT FOR THE LIQUID FLOW RATES
    mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
    mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    """
    Check DOF
    """
    for n, eff in mec.effects.items():
        assert degrees_of_freedom(eff.effect) == 0

def init_mec(blk):
    mec = blk.unit

    ### INITIALIZE FOR EACH EFFECT
    for n, eff in mec.effects.items():
        eff.effect.initialize()

    ### UNFIX THE INLET FLOW RATES OF EACH EFFECT
    for n, eff in mec.effects.items():
        if n > 1:
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    for n, eff in mec.effects.items():
        if n == 1:
            assert degrees_of_freedom(eff.effect) == 0
        else:
            assert degrees_of_freedom(eff.effect) == 1

    ### FULLY SOLVE THE MODEL
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    calculate_scaling_factors(blk)

    solver = get_solver()
    results = solver.solve(blk)
    assert_optimal_termination(results)

def unfix_mec(blk):
    mec = blk.unit

    first_effect = mec.effects[1].effect
    ### Release 1st effect flow rate and fix total flow rate instead
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    mec.inlet.temperature[0].unfix()
    mec.inlet.pressure[0].unfix()

def mec_rescaling(blk,
                  flow_mass_phase_water_total = 116.247/10,
                  flow_mass_phase_salt_total = 28.478/10):

    """
    Note: Rescaling is probably needed for extremely large feed flow,
    """
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
        index=("Liq", "NaCl"),
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 10, index=("Vap", "H2O")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Sol", "NaCl")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Vap", "H2O")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )


def add_mec_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=flowsheet_costing_block,
        )


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results

if __name__ == "__main__":

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = TreatmentCosting()
    build_mec(m, m.fs)

    set_mec_op_conditions(m, m.fs)
    init_mec(m.fs)
    unfix_mec(m.fs)

    flow_mass_phase_water_total = 11.6
    flow_mass_phase_salt_total = 2.8

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        flow_mass_phase_water_total
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        flow_mass_phase_salt_total
    )

    m.fs.unit.inlet.temperature[0].fix(273.15 + 30.51)
    m.fs.unit.inlet.pressure[0].fix(101325)
    mec_rescaling(m.fs)
    add_mec_costing(m, m.fs)

    print('')
    print('here')
    print('')
    m.fs.unit.inlet.display()

    solve(m)
