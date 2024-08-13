import pandas as pd
import numpy as np
import pytest
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    SolverStatus,
    Objective,
    Expression,
    maximize,
    value,
    Set,
    Var,
    log,
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.models.unit_models import HeatExchanger
from idaes.models.unit_models.heat_exchanger import (
    delta_temperature_lmtd_callback,
    delta_temperature_underwood_callback,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.unit_models.zero_order.crystallizer_zo_watertap import (
    Crystallization,
)
import watertap_contrib.reflo.property_models.cryst_prop_pack as props
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    CrystallizerCostType,
    _compute_steam_properties,
)

solver = get_solver()


def build_fs_multi_effect_crystallizer(
    m=None,
    # num_effect=3,
    operating_pressure_eff1=0.78,  # bar
    operating_pressure_eff2=0.25,  # bar
    operating_pressure_eff3=0.208,  # bar
    operating_pressure_eff4=0.095,  # bar
    feed_flow_mass=1,  # kg/s
    feed_mass_frac_NaCl=0.3,
    feed_pressure=101325,  # Pa
    feed_temperature=273.15 + 20,  # K
    crystallizer_yield=0.5,
    steam_pressure=1.5,  # bar (gauge pressure)
):
    """
    This flowsheet depicts a 4-effect crystallizer, with brine fed in parallel
    to each effect, and the operating pressure is specfied individually.
    """

    if m is None:
        m = ConcreteModel()

    mfs = m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = props.NaClParameterBlock()

    # Create 4 effects of crystallizer
    eff_1 = m.fs.eff_1 = Crystallization(property_package=m.fs.props)
    eff_2 = m.fs.eff_2 = Crystallization(property_package=m.fs.props)
    eff_3 = m.fs.eff_3 = Crystallization(property_package=m.fs.props)
    eff_4 = m.fs.eff_4 = Crystallization(property_package=m.fs.props)

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    for eff in [eff_1, eff_2, eff_3, eff_4]:
        #  Define feed for all effects
        eff.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
        eff.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

        eff.inlet.pressure[0].fix(feed_pressure)
        eff.inlet.temperature[0].fix(feed_temperature)

        # Fix growth rate, crystal length and Sounders brown constant to default values
        eff.crystal_growth_rate.fix()
        eff.souders_brown_constant.fix()
        eff.crystal_median_length.fix()

        # Fix yield
        eff.crystallization_yield["NaCl"].fix(crystallizer_yield)

    # Define operating conditions
    m.fs.eff_1.pressure_operating.fix(operating_pressure_eff1 * pyunits.bar)
    m.fs.eff_2.pressure_operating.fix(operating_pressure_eff2 * pyunits.bar)
    m.fs.eff_3.pressure_operating.fix(operating_pressure_eff3 * pyunits.bar)
    m.fs.eff_4.pressure_operating.fix(operating_pressure_eff4 * pyunits.bar)

    steam_temp = add_heat_exchanger_eff1(m, steam_pressure)
    add_heat_exchanger_eff2(m)
    add_heat_exchanger_eff3(m)
    add_heat_exchanger_eff4(m)

    return m


def add_heat_exchanger_eff1(m, steam_pressure):
    eff_1 = m.fs.eff_1

    eff_1.delta_temperature_in = Var(
        eff_1.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the inlet side",
    )
    eff_1.delta_temperature_out = Var(
        eff_1.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the outlet side",
    )
    delta_temperature_chen_callback(eff_1)

    eff_1.area = Var(
        bounds=(0, None),
        initialize=1000.0,
        doc="Heat exchange area",
        units=pyunits.m**2,
    )

    eff_1.overall_heat_transfer_coefficient = Var(
        eff_1.flowsheet().time,
        bounds=(0, None),
        initialize=100.0,
        doc="Overall heat transfer coefficient",
        units=pyunits.W / pyunits.m**2 / pyunits.K,
    )

    eff_1.overall_heat_transfer_coefficient[0].fix(100)

    # Compute saturation temperature of steam: computed from El-Dessouky expression
    steam_pressure_sat = steam_pressure * pyunits.bar

    tsat_constants = [
        42.6776 * pyunits.K,
        -3892.7 * pyunits.K,
        1000 * pyunits.kPa,
        -9.48654 * pyunits.dimensionless,
    ]
    psat = (
        pyunits.convert(steam_pressure_sat, to_units=pyunits.kPa)
        + 101.325 * pyunits.kPa
    )

    temperature_sat = tsat_constants[0] + tsat_constants[1] / (
        log(psat / tsat_constants[2]) + tsat_constants[3]
    )

    @m.Constraint(eff_1.flowsheet().time, doc="delta_temperature_in at the 1st effect")
    def delta_temperature_in_eff1(b, t):
        return (
            b.fs.eff_1.delta_temperature_in[t]
            == temperature_sat - b.fs.eff_1.temperature_operating
        )

    @m.Constraint(eff_1.flowsheet().time, doc="delta_temperature_out at the 1st effect")
    def delta_temperature_out_eff1(b, t):
        return (
            b.fs.eff_1.delta_temperature_out[t]
            == temperature_sat - b.fs.eff_1.properties_in[0].temperature
        )

    @m.Constraint(eff_1.flowsheet().time)
    def heat_transfer_equation_eff_1(b, t):
        return b.fs.eff_1.work_mechanical[0] == (
            b.fs.eff_1.overall_heat_transfer_coefficient[t]
            * b.fs.eff_1.area
            * b.fs.eff_1.delta_temperature[0]
        )

    iscale.set_scaling_factor(eff_1.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(eff_1.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(eff_1.area, 1e-1)
    iscale.set_scaling_factor(eff_1.overall_heat_transfer_coefficient, 1e-1)

    return temperature_sat()


def add_heat_exchanger_eff2(m):
    eff_1 = m.fs.eff_1
    eff_2 = m.fs.eff_2

    eff_2.delta_temperature_in = Var(
        eff_2.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the inlet side",
    )
    eff_2.delta_temperature_out = Var(
        eff_2.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the outlet side",
    )
    delta_temperature_chen_callback(eff_2)

    eff_2.area = Var(
        bounds=(0, None),
        initialize=1000.0,
        doc="Heat exchange area",
        units=pyunits.m**2,
    )

    eff_2.overall_heat_transfer_coefficient = Var(
        eff_2.flowsheet().time,
        bounds=(0, None),
        initialize=100.0,
        doc="Overall heat transfer coefficient",
        units=pyunits.W / pyunits.m**2 / pyunits.K,
    )

    eff_2.overall_heat_transfer_coefficient[0].fix(100)

    @m.Constraint(eff_2.flowsheet().time, doc="delta_temperature_in at the 2nd effect")
    def delta_temperature_in_eff2(b, t):
        return (
            b.fs.eff_2.delta_temperature_in[t]
            == b.fs.eff_1.properties_vapor[0].temperature
            - b.fs.eff_2.temperature_operating
        )

    @m.Constraint(eff_2.flowsheet().time, doc="delta_temperature_out at the 2nd effect")
    def delta_temperature_out_eff2(b, t):
        return (
            b.fs.eff_2.delta_temperature_out[t]
            == b.fs.eff_1.properties_pure_water[0].temperature
            - b.fs.eff_2.properties_in[0].temperature
        )

    @m.Constraint(eff_2.flowsheet().time)
    def heat_transfer_equation_eff_2(b, t):
        return b.fs.eff_1.energy_flow_superheated_vapor == (
            b.fs.eff_2.overall_heat_transfer_coefficient[t]
            * b.fs.eff_2.area
            * b.fs.eff_2.delta_temperature[0]
        )

    iscale.set_scaling_factor(eff_2.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(eff_2.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(eff_2.area, 1e-1)
    iscale.set_scaling_factor(eff_2.overall_heat_transfer_coefficient, 1e-1)


def add_heat_exchanger_eff3(m):
    eff_2 = m.fs.eff_2
    eff_3 = m.fs.eff_3

    eff_3.delta_temperature_in = Var(
        eff_3.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the inlet side",
    )
    eff_3.delta_temperature_out = Var(
        eff_3.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the outlet side",
    )
    delta_temperature_chen_callback(eff_3)

    eff_3.area = Var(
        bounds=(0, None),
        initialize=1000.0,
        doc="Heat exchange area",
        units=pyunits.m**2,
    )

    eff_3.overall_heat_transfer_coefficient = Var(
        eff_3.flowsheet().time,
        bounds=(0, None),
        initialize=100.0,
        doc="Overall heat transfer coefficient",
        units=pyunits.W / pyunits.m**2 / pyunits.K,
    )

    eff_3.overall_heat_transfer_coefficient[0].fix(100)

    @m.Constraint(eff_3.flowsheet().time, doc="delta_temperature_in at the 2nd effect")
    def delta_temperature_in_eff3(b, t):
        return (
            eff_3.delta_temperature_in[t]
            == eff_2.properties_vapor[0].temperature - eff_3.temperature_operating
        )

    @m.Constraint(eff_3.flowsheet().time, doc="delta_temperature_out at the 2nd effect")
    def delta_temperature_out_eff3(b, t):
        return (
            eff_3.delta_temperature_out[t]
            == eff_2.properties_pure_water[0].temperature
            - eff_3.properties_in[0].temperature
        )

    @m.Constraint(eff_3.flowsheet().time)
    def heat_transfer_equation_eff3(b, t):
        return eff_2.energy_flow_superheated_vapor == (
            eff_3.overall_heat_transfer_coefficient[t]
            * eff_3.area
            * eff_3.delta_temperature[0]
        )

    iscale.set_scaling_factor(eff_3.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(eff_3.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(eff_3.area, 1e-1)
    iscale.set_scaling_factor(eff_3.overall_heat_transfer_coefficient, 1e-1)


def add_heat_exchanger_eff4(m):
    eff_3 = m.fs.eff_3
    eff_4 = m.fs.eff_4

    eff_4.delta_temperature_in = Var(
        eff_4.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the inlet side",
    )
    eff_4.delta_temperature_out = Var(
        eff_4.flowsheet().time,
        initialize=35,
        bounds=(None, None),
        units=pyunits.K,
        doc="Temperature differnce at the outlet side",
    )
    delta_temperature_chen_callback(eff_4)

    eff_4.area = Var(
        bounds=(0, None),
        initialize=1000.0,
        doc="Heat exchange area",
        units=pyunits.m**2,
    )

    eff_4.overall_heat_transfer_coefficient = Var(
        eff_4.flowsheet().time,
        bounds=(0, None),
        initialize=100.0,
        doc="Overall heat transfer coefficient",
        units=pyunits.W / pyunits.m**2 / pyunits.K,
    )

    eff_4.overall_heat_transfer_coefficient[0].fix(100)

    @m.Constraint(eff_4.flowsheet().time, doc="delta_temperature_in at the 2nd effect")
    def delta_temperature_in_eff4(b, t):
        return (
            eff_4.delta_temperature_in[t]
            == eff_3.properties_vapor[0].temperature - eff_4.temperature_operating
        )

    @m.Constraint(eff_4.flowsheet().time, doc="delta_temperature_out at the 2nd effect")
    def delta_temperature_out_eff4(b, t):
        return (
            eff_4.delta_temperature_out[t]
            == eff_3.properties_pure_water[0].temperature
            - eff_4.properties_in[0].temperature
        )

    @m.Constraint(eff_4.flowsheet().time)
    def heat_transfer_equation_eff4(b, t):
        return eff_3.energy_flow_superheated_vapor == (
            eff_4.overall_heat_transfer_coefficient[t]
            * eff_4.area
            * eff_4.delta_temperature[0]
        )

    iscale.set_scaling_factor(eff_4.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(eff_4.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(eff_4.area, 1e-1)
    iscale.set_scaling_factor(eff_4.overall_heat_transfer_coefficient, 1e-1)


def add_costings(m):
    effs = [m.fs.eff_1, m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]

    m.fs.capex_heat_exchanger = Expression(
        expr=(420 * sum(i.area for i in effs)), doc="Capital cost of heat exchangers"
    )

    m.fs.capex_end_plates = Expression(
        expr=(1020 * (sum(i.area for i in effs) / 10) ** 0.6),
        doc="Capital cost of heat exchanger endplates",
    )

    m.fs.costing = TreatmentCosting()
    m.fs.eff_1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )
    m.fs.eff_2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )
    m.fs.eff_3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )
    m.fs.eff_4.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
    )

    # Effect 2-4 doesn't need additional heating steam and the flows are removed
    m.fs.eff_2.costing.costing_package.cost_flow(
        -pyunits.convert(
            (
                m.fs.eff_2.work_mechanical[0]
                / _compute_steam_properties(m.fs.eff_2.costing)
            ),
            to_units=pyunits.m**3 / pyunits.s,
        ),
        "steam",
    )
    m.fs.eff_3.costing.costing_package.cost_flow(
        -pyunits.convert(
            (
                m.fs.eff_3.work_mechanical[0]
                / _compute_steam_properties(m.fs.eff_3.costing)
            ),
            to_units=pyunits.m**3 / pyunits.s,
        ),
        "steam",
    )

    m.fs.eff_4.costing.costing_package.cost_flow(
        -pyunits.convert(
            (
                m.fs.eff_4.work_mechanical[0]
                / _compute_steam_properties(m.fs.eff_4.costing)
            ),
            to_units=pyunits.m**3 / pyunits.s,
        ),
        "steam",
    )

    m.fs.costing.cost_process()

    feed_vol_flow_rates = sum(i.properties_in[0].flow_vol_phase["Liq"] for i in effs)

    # Add a term for the unit cost of treated brine ($/m3)
    m.fs.levelized_cost_of_feed_brine = Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.capital_recovery_factor
                * (m.fs.capex_heat_exchanger + m.fs.capex_end_plates)
            )
            / pyunits.convert(
                feed_vol_flow_rates, to_units=pyunits.m**3 / pyunits.year
            )
        ),
        doc="Levelized cost of feed brine",
    )


def multi_effect_crystallizer_initialization(m):
    # Set scaling factors
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Vap", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl"))

    calculate_scaling_factors(m)

    m.fs.eff_1.initialize()
    m.fs.eff_2.initialize()
    m.fs.eff_3.initialize()
    m.fs.eff_4.initialize()

    # Unfix dof
    brine_salinity = (
        m.fs.eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value
    )

    for eff in [m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]:
        eff.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
        eff.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
        eff.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(brine_salinity)

    # Energy is provided from the previous effect
    @m.Constraint(doc="Energy supplied to the 2nd effect")
    def eqn_energy_from_eff1(b):
        return b.fs.eff_2.work_mechanical[0] == b.fs.eff_1.energy_flow_superheated_vapor

    @m.Constraint(doc="Energy supplied to the 3rd effect")
    def eqn_energy_from_eff2(b):
        return b.fs.eff_3.work_mechanical[0] == b.fs.eff_2.energy_flow_superheated_vapor

    @m.Constraint(doc="Energy supplied to the 4th effect")
    def eqn_energy_from_eff3(b):
        return b.fs.eff_4.work_mechanical[0] == b.fs.eff_3.energy_flow_superheated_vapor


def get_model_performance(m):
    # Print result
    effs = [m.fs.eff_1, m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]
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
        m.fs.eff_1.work_mechanical[0].value
        / sum(i.properties_in[0].flow_vol_phase["Liq"].value for i in effs)
        / 3600
    )

    area = [i.area.value for i in effs]

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
        "Feed brine salinity (g/L)": m.fs.eff_1.properties_in[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
        .value,
        "Total brine disposed (kg/s)": sum(feed_flow_rates),
        "Total water production (kg/s)": sum(water_prod),
        "Total solids collected (kg/s)": sum(solid_prod),
        "Total waste water remained (kg/s)": sum(liquid_prod),
        "Initial thermal energy consumption (kW)": m.fs.eff_1.work_mechanical[0].value,
        "Overall STEC (kWh/m3 feed)": overall_STEC,
        "Total heat transfer area (m2)": sum(i.area.value for i in effs),
        "Levelized cost of feed brine ($/m3)": value(m.fs.levelized_cost_of_feed_brine),
    }

    return data_table, overall_performance


if __name__ == "__main__":
    m = build_fs_multi_effect_crystallizer(
        operating_pressure_eff1=0.45,  # bar
        operating_pressure_eff2=0.25,  # bar
        operating_pressure_eff3=0.208,  # bar
        operating_pressure_eff4=0.095,  # bar
        feed_flow_mass=1,  # kg/s
        feed_mass_frac_NaCl=0.15,
        feed_pressure=101325,  # Pa
        feed_temperature=273.15 + 20,  # K
        crystallizer_yield=0.5,
        steam_pressure=1.5,  # bar (gauge pressure)
    )
    add_costings(m)

    # Negative value for salt recovery value ($/kg)
    m.fs.costing.crystallizer.NaCl_recovery_value.fix(-0.024)

    multi_effect_crystallizer_initialization(m)

    results = solver.solve(m)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    data_table, overall_performance = get_model_performance(m)
