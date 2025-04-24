from pprint import pprint
import numpy as np
import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    Objective,
    Constraint,
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
    kbhdp_salinity=229,  # g/L
    assumed_lssro_recovery=0.9,
    number_effects=4,
    crystallizer_yield=0.8,
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
        (kbhdp_salinity * pyunits.g / pyunits.liter),
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
        eff.effect.properties_out[0].conc_mass_phase_comp[...]

        eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(1.3)

    first_effect = m.fs.mec.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(1.3)
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

    print("")
    print("solve for normalized volume")

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

    print("")
    print("solve for after init")
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


def get_model_performance(m):
    # Print result
    effs = [
        m.fs.mec.effects[1].effect,
        m.fs.mec.effects[2].effect,
        m.fs.mec.effects[3].effect,
        m.fs.mec.effects[4].effect,
    ]
    effect_names = ["Effect 1", "Effect 2", "Effect 3", "Effect 4"]

    feed_salinities = [
        i.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value for i in effs
    ]
    feed_flow_rates = [
        sum(
            i.properties_in[0].flow_mass_phase_comp["Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        for i in effs
    ]
    feed_vol_flow_rates = [
        i.properties_in[0].flow_vol_phase["Liq"].value * 1000 for i in effs
    ]
    temp_operating = [i.temperature_operating.value - 273.15 for i in effs]
    temp_vapor_cond = [
        i.properties_pure_water[0].temperature.value - 273.15 for i in effs
    ]
    p_operating = [i.pressure_operating.value / 1e5 for i in effs]
    water_prod = [
        i.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value for i in effs
    ]
    solid_prod = [
        i.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"].value for i in effs
    ]
    liquid_prod = [
        sum(
            i.properties_out[0].flow_mass_phase_comp["Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        for i in effs
    ]
    liquid_flow_rate = [
        i.properties_out[0].flow_vol_phase["Liq"].value * 1000 for i in effs
    ]

    liquid_salinity = [
        i.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"].value for i in effs
    ]
    power_required = [i.work_mechanical[0].value for i in effs]
    power_provided = [i.energy_flow_superheated_vapor.value for i in effs]
    vapor_enth = [
        i.properties_vapor[0].dh_vap_mass_solvent.value
        * i.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
        for i in effs
    ]
    STEC = [
        i.work_mechanical[0].value
        / i.properties_in[0].flow_vol_phase["Liq"].value
        / 3600
        for i in effs
    ]

    overall_STEC = (
        effs[1].work_mechanical[0].value
        / sum(i.properties_in[0].flow_vol_phase["Liq"].value for i in effs)
        / 3600
    )

    area = [i.heat_exchanger_area.value for i in effs]

    model_output = np.array(
        [
            feed_flow_rates,
            feed_vol_flow_rates,
            feed_salinities,
            temp_operating,
            temp_vapor_cond,
            p_operating,
            water_prod,
            solid_prod,
            liquid_prod,
            liquid_flow_rate,
            liquid_salinity,
            power_required,
            power_provided,
            vapor_enth,
            STEC,
            area,
        ]
    )

    data_table = pd.DataFrame(
        data=model_output,
        columns=effect_names,
        index=[
            "Feed mass flow rate (kg/s)",
            "Feed volumetric flow rate (L/s)",
            "Feed salinities (g/L)",
            "Operating temperature (C)",
            "Vapor condensation temperature (C)",
            "Operating pressure (bar)",
            "Water production (kg/s)",
            "Solid production (kg/s)",
            "Liquid waste (kg/s)",
            "Liquid waste volumetric flow rate (L/s)",
            "Liquid waste salinity (g/L)",
            "Thermal energy requirement (kW)",
            "Thermal energy available from vapor (kW)",
            "Vapor enthalpy (kJ)",
            "STEC (kWh/m3 feed)",
            "Heat transfer area (m2)",
        ],
    )

    overall_performance = {
        "Capacity (m3/day)": sum(feed_vol_flow_rates) * 86400 / 1000,
        "Feed brine salinity (g/L)": effs[1]
        .properties_in[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
        .value,
        "Total brine disposed (kg/s)": sum(feed_flow_rates),
        "Total water production (kg/s)": sum(water_prod),
        "Total solids collected (kg/s)": sum(solid_prod),
        "Total waste water remained (kg/s)": sum(liquid_prod),
        "Initial thermal energy consumption (kW)": effs[1].work_mechanical[0].value,
        "Overall STEC (kWh/m3 feed)": overall_STEC,
        "Total heat transfer area (m2)": sum(i.heat_exchanger_area.value for i in effs),
    }

    return data_table, overall_performance


if __name__ == "__main__":
    m = build_kbhdp_mec(
        # m=m
        flow_in=1,  # MGD
        kbhdp_salinity=230,  # g/L
        # assumed_lssro_recovery=0.95,
        # number_effects=4,
        crystallizer_yield=0.9,
        # eps=1e-12,
        # saturated_steam_pressure_gage=3,
        # heat_transfer_coefficient=0.1
    )
    calculate_scaling_factors(m)

    print("")
    print("here")
    print("")
    assert degrees_of_freedom(m) == 0
    try:
        results = solver.solve(m)
        assert_optimal_termination(results)
    except:
        print("Failed")
    overall_STEC = (
        m.fs.mec.effects[1].effect.work_mechanical[0].value
        / sum(
            m.fs.mec.effects[i].effect.properties_in[0].flow_vol_phase["Liq"].value
            for i in m.fs.mec.Effects
        )
        / 3600
    )
    total_area = sum(
        m.fs.mec.effects[i].effect.heat_exchanger_area for i in m.fs.mec.Effects
    )

    print(value(overall_STEC))
    print(value(total_area))

    data_table, performance = get_model_performance(m)

    m.fs.mec.effects[1].effect.pressure_operating.unfix()
    m.fs.mec.effects[1].effect.pressure_operating.setub(0.5 * 1e5)
    m.fs.mec.effects[2].effect.pressure_operating.unfix()
    m.fs.mec.effects[3].effect.pressure_operating.unfix()
    m.fs.mec.effects[4].effect.pressure_operating.unfix()
    m.fs.mec.effects[4].effect.pressure_operating.setlb(0.022 * 1e5)

    @m.Constraint(m.fs.mec.Effects, doc="Pressure decreasing")
    def pressure_bound1(b, j):
        if j < 4:
            return (
                b.fs.mec.effects[j + 1].effect.pressure_operating
                <= b.fs.mec.effects[j].effect.pressure_operating
            )
        else:
            return Constraint.Skip

    @m.Constraint(m.fs.mec.Effects, doc="Temperature difference")
    def temp_bound1(b, j):
        if j < 4:
            return (
                b.fs.mec.effects[j + 1].effect.temperature_operating
                >= b.fs.mec.effects[j].effect.temperature_operating - 12
            )
        else:
            return Constraint.Skip

    overall_STEC = (
        m.fs.mec.effects[1].effect.work_mechanical[0].value
        / sum(
            m.fs.mec.effects[i].effect.properties_in[0].flow_vol_phase["Liq"].value
            for i in m.fs.mec.Effects
        )
        / 3600
    )
    total_area = sum(
        m.fs.mec.effects[i].effect.heat_exchanger_area for i in m.fs.mec.Effects
    )
    m.fs.objective = Objective(expr=total_area)

    optimization_results = solver.solve(m, tee=False)
    assert (
        optimization_results.solver.termination_condition
        == TerminationCondition.optimal
    )

    data_table2, performance2 = get_model_performance(m)
