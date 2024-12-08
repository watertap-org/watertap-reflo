from pprint import pprint

from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc


from idaes.core.util.scaling import *
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
)
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.multi_effect_crystallizer import (
    MultiEffectCrystallizer,
)
from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect


solver = get_solver()
rho = 1000 * pyunits.kg / pyunits.m**3
feed_pressure = 101325
feed_temperature = 273.15 + 20


def build_kbhdp_mec(
    flow_in=4,  # MGD
    kbhdp_salinity=12,  # g/L
    assumed_lssro_recovery=0.9,
    number_effects=4,
    crystallizer_yield=0.5,
    saturated_steam_pressure_gage=3,
    heat_transfer_coefficient=0.1,
    eps=1e-8,
):
    """
    Build MultiEffectCrystallizer for KBHDP case study.
        flow_in: volumetric flow rate in MGD
        kbhdp_salinity: salinity of KBHDP brine
        assumed_lssro_recovery: recovery of LSRRO process that is assumed to be before MEC; used to approximate influent concentration
        number_effects: number of effects for MEC model
    """

    global flow_mass_phase_water_total, flow_mass_phase_salt_total, conc_in

    atm_pressure = 101325 * pyunits.Pa

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.mec = mec = MultiEffectCrystallizer(
        property_package=m.fs.properties,
        property_package_vapor=m.fs.vapor_properties,
        number_effects=number_effects,
    )

    operating_pressures = [0.45, 0.25, 0.208, 0.095]

    conc_in = pyunits.convert(
        (kbhdp_salinity * pyunits.g / pyunits.liter) / (1 - assumed_lssro_recovery),
        to_units=pyunits.kg / pyunits.m**3,
    )

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
    flow_mass_phase_water_per = rho / (rho + conc_in) * 1 * pyunits.kg / pyunits.s
    flow_mass_phase_salt_per = conc_in / (rho + conc_in) * 1 * pyunits.kg / pyunits.s

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
        eff.effect.overall_heat_transfer_coefficient.fix(0.1)

    first_effect = m.fs.mec.effects[1].effect

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
    m.fs.mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
    m.fs.mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
    m.fs.mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    """
    Check DOF
    By this point, each effect is fully specified, so their DOF should be 0.
    However, the multi-effect flowsheet is over-constrianted by energy_flow_constr, 
    which connect energy flow between effects, and the DOF should be negative (n_effects - 1)
    """
    for n, eff in m.fs.mec.effects.items():
        assert degrees_of_freedom(eff.effect) == 0
    assert degrees_of_freedom(m) == -3

    ### INITIALIZE FOR EACH EFFECT
    """
    Note: this is essentially to have an initial guess of the crysts,
    and to populate the feed concentration to all effects
    """
    for n, eff in m.fs.mec.effects.items():
        eff.effect.initialize()

    ### UNFIX THE INLET FLOW RATES OF EACH EFFECT
    """
    Note: this is to release the volumetric inlet flow entering different effects
    with the same salinity
    """
    for n, eff in m.fs.mec.effects.items():
        if n > 1:
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    """
    Note: by this point, the multi-effect flow sheet should be fully specified (DOF=0),
    while inlet flows to effect 2-4 should be subject to the previous effect (DOF=1)
    """
    assert degrees_of_freedom(m) == 0
    for n, eff in m.fs.mec.effects.items():
        if n == 1:
            assert degrees_of_freedom(eff.effect) == 0
        else:
            assert degrees_of_freedom(eff.effect) == 1

    ### FULLY SOLVE THE MODEL
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

    calculate_scaling_factors(m)
    results = solver.solve(m)
    assert_optimal_termination(results)

    ### Release 1st effect flow rate and fix total flow rate instead
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    m.fs.mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        flow_mass_phase_water_total
    )
    m.fs.mec.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        flow_mass_phase_salt_total
    )

    """
    Note: Rescaling is probably needed for extremely large feed flow,
    """
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
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
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


if __name__ == "__main__":

    m = build_kbhdp_mec(
        # m=m
        flow_in=2.8,  # MGD
        # kbhdp_salinity=12,  # g/L
        # assumed_lssro_recovery=0.95,
        # number_effects=4,
        crystallizer_yield=0.8,
        # eps=1e-12,
        # saturated_steam_pressure_gage=3,
        # heat_transfer_coefficient=0.1
    )
    calculate_scaling_factors(m)

    assert degrees_of_freedom(m) == 0
    try:
        results = solver.solve(m)
        assert_optimal_termination(results)
    except:
        print("Failed")
