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
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
from watertap_contrib.reflo.unit_models.zero_order.crystallizer_zo_watertap import Crystallization
import watertap_contrib.reflo.property_models.cryst_prop_pack as props
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from idaes.core.solvers import get_solver
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

from watertap.costing import WaterTAPCosting, CrystallizerCostType

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

    feed_mass_frac_H2O = 1- feed_mass_frac_NaCl
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

    add_heat_exchanger_eff2(m)
    add_heat_exchanger_eff3(m)
    add_heat_exchanger_eff4(m)
    return m


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
            == eff_2.properties_vapor[0].temperature
            - eff_3.temperature_operating
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
            == eff_3.properties_vapor[0].temperature
            - eff_4.temperature_operating
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

def multi_effect_crystallizer_initialization(m):
    # Set scaling factors
    m.fs.props.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.props.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    m.fs.props.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    m.fs.props.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )

    calculate_scaling_factors(m)

    m.fs.eff_1.initialize()
    m.fs.eff_2.initialize()
    m.fs.eff_3.initialize()
    m.fs.eff_4.initialize()

    # Unfix dof
    brine_salinity = m.fs.eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value

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


if __name__ == "__main__":
    m = build_fs_multi_effect_crystallizer(
        operating_pressure_eff1=0.45,  # bar
        operating_pressure_eff2=0.25,  # bar
        operating_pressure_eff3=0.208,  # bar
        operating_pressure_eff4=0.095,  # bar
        feed_flow_mass=1,  # kg/s
        feed_mass_frac_NaCl=0.3,
        feed_pressure=101325,  # Pa
        feed_temperature=273.15 + 20,  # K
        crystallizer_yield=0.5,
        steam_pressure=1.5,  # bar (gauge pressure)
    )

    multi_effect_crystallizer_initialization(m)

    print(degrees_of_freedom(m))
    solver = get_solver()

    results = solver.solve(m)

    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    print(m.fs.eff_1.energy_flow_superheated_vapor.value)
    print(m.fs.eff_2.energy_flow_superheated_vapor.value)
    # print(m.fs.eff_3.energy_flow_superheated_vapor.value)
    # print(m.fs.eff_4.energy_flow_superheated_vapor.value)
    print(m.fs.eff_1.temperature_operating.value - 273.15)
    print(m.fs.eff_2.temperature_operating.value - 273.15)
    # print(m.fs.eff_3.temperature_operating.value - 273.15)
    # print(m.fs.eff_4.temperature_operating.value - 273.15)
