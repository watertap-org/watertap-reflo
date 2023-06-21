import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.seto.unit_models.surrogate import VAGMDSurrogate
from watertap_contrib.seto.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.util.model_diagnostics.infeasible import *

# Get default solver for testing
solver = get_solver()


@pytest.fixture(scope="module")
def build_vagmd_fs():
    # System specification (Input variables)
    feed_flow_rate = 600  # 400 - 1100 L/h
    evap_inlet_temp = 80  # 60 - 80 deg C
    cond_inlet_temp = 25  # 20 - 30 deg C
    feed_temp = 25  # 20 - 30 deg C
    feed_salinity = 35  # 35 - 292 g/L
    initial_batch_volume = 50  # > 50 L
    recovery_ratio = 0.5  # -
    module_type = "AS7C1.5L"
    cooling_system_type = "closed"
    cooling_inlet_temp = 25  # deg C, not required when cooling system type is "closed"

    # Create the flowsheet the the system specifications
    m = build_vagmd_flowsheet(
        feed_flow_rate=feed_flow_rate,
        evap_inlet_temp=evap_inlet_temp,
        cond_inlet_temp=cond_inlet_temp,
        feed_temp=feed_temp,
        feed_salinity=feed_salinity,
        recovery_ratio=recovery_ratio,
        initial_batch_volume=initial_batch_volume,
        module_type=module_type,
        cooling_system_type=cooling_system_type,
    )

    # Set scaling factors and initialize the flowsheet
    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,
    )

    return m


@pytest.mark.unit
def test_dof(build_vagmd_fs):
    m = build_vagmd_fs
    assert degrees_of_freedom(m) == 0


@pytest.mark.unit
def test_calculate_scaling(build_vagmd_fs):
    m = build_vagmd_fs
    # check that all variables have scaling factors
    unscaled_var_list = list(unscaled_variables_generator(m))
    assert len(unscaled_var_list) == 0

    # check that all constraints have been scaled
    unscaled_constraint_list = list(unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0


@pytest.mark.component
def test_var_scaling(build_vagmd_fs):
    m = build_vagmd_fs
    badly_scaled_var_lst = list(badly_scaled_var_generator(m))
    assert badly_scaled_var_lst == []


@pytest.mark.unit
def test_fs_solution(build_vagmd_fs):

    m = build_vagmd_fs
    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(69.4664, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(34.5742, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.3404, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(7.0701, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(33.921, abs=1e-3)


@pytest.mark.unit
def test_fs_solution(build_vagmd_fs):

    m = build_vagmd_fs
    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(69.4664, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(34.5742, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.3404, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(7.0701, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(33.921, abs=1e-3)


@pytest.mark.unit
def test_fs2_solution():
    # Create model, flowsheet for configuration of module AS7C1.5L,
    # with high brine salinity (> 173.5 g/L) and closed cooling system
    feed_flow_rate = 600  # 400 - 1100 L/h
    evap_inlet_temp = 80  # 60 - 80 deg C
    cond_inlet_temp = 25  # 20 - 30 deg C
    feed_temp = 25  # 20 - 30 deg C
    feed_salinity = 100  # 35 - 292 g/L
    initial_batch_volume = 50  # > 50 L
    recovery_ratio = 0.5  # -
    module_type = "AS7C1.5L"
    cooling_system_type = "closed"
    cooling_inlet_temp = 25  # deg C, not required when cooling system type is "closed"

    # Create the flowsheet the the system specifications
    m = build_vagmd_flowsheet(
        feed_flow_rate=feed_flow_rate,
        evap_inlet_temp=evap_inlet_temp,
        cond_inlet_temp=cond_inlet_temp,
        feed_temp=feed_temp,
        feed_salinity=feed_salinity,
        recovery_ratio=recovery_ratio,
        initial_batch_volume=initial_batch_volume,
        module_type=module_type,
        cooling_system_type=cooling_system_type,
    )

    # Set scaling factors and initialize the flowsheet
    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,
    )

    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(65.6459, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(37.4374, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(7.6420, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(17.2473, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(18.502, abs=1e-3)


@pytest.mark.unit
def test_fs3_solution():
    # Create model, flowsheet for configuration of module AS7C1.5L,
    # with low brine salinity (< 173.5 g/L) and open cooling system
    feed_flow_rate = 600  # 400 - 1100 L/h
    evap_inlet_temp = 80  # 60 - 80 deg C
    cond_inlet_temp = 25  # 20 - 30 deg C
    feed_temp = 25  # 20 - 30 deg C
    feed_salinity = 50  # 35 - 292 g/L
    initial_batch_volume = 50  # > 50 L
    recovery_ratio = 0.5  # -
    module_type = "AS7C1.5L"
    cooling_system_type = "open"
    cooling_inlet_temp = 25  # deg C, not required when cooling system type is "closed"

    # Create the flowsheet the the system specifications
    m = build_vagmd_flowsheet(
        feed_flow_rate=feed_flow_rate,
        evap_inlet_temp=evap_inlet_temp,
        cond_inlet_temp=cond_inlet_temp,
        feed_temp=feed_temp,
        feed_salinity=feed_salinity,
        recovery_ratio=recovery_ratio,
        initial_batch_volume=initial_batch_volume,
        module_type=module_type,
        cooling_system_type=cooling_system_type,
    )

    # Set scaling factors and initialize the flowsheet
    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,
    )

    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(69.3523, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(34.5056, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.2033, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(7.1045, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(33.921, abs=1e-3)


@pytest.mark.unit
def test_fs4_solution():
    # Create model, flowsheet for configuration of module AS26C7.2L,
    # with closed cooling system
    feed_flow_rate = 600  # 400 - 1100 L/h
    evap_inlet_temp = 80  # 60 - 80 deg C
    cond_inlet_temp = 25  # 20 - 30 deg C
    feed_temp = 25  # 20 - 30 deg C
    feed_salinity = 35  # 35 - 292 g/L
    initial_batch_volume = 50  # > 50 L
    recovery_ratio = 0.5  # -
    module_type = "AS26C7.2L"
    cooling_system_type = "closed"
    cooling_inlet_temp = 25  # deg C, not required when cooling system type is "closed"

    # Create the flowsheet the the system specifications
    m = build_vagmd_flowsheet(
        feed_flow_rate=feed_flow_rate,
        evap_inlet_temp=evap_inlet_temp,
        cond_inlet_temp=cond_inlet_temp,
        feed_temp=feed_temp,
        feed_salinity=feed_salinity,
        recovery_ratio=recovery_ratio,
        initial_batch_volume=initial_batch_volume,
        module_type=module_type,
        cooling_system_type=cooling_system_type,
    )

    # Set scaling factors and initialize the flowsheet
    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,
    )

    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(75.7306, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(27.8301, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(1.6197, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(2.8616, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(122.115, abs=1e-3)


@pytest.mark.unit
def test_fs5_solution():
    # Create model, flowsheet for configuration of module AS26C7.2L,
    # with open cooling system
    feed_flow_rate = 600  # 400 - 1100 L/h
    evap_inlet_temp = 80  # 60 - 80 deg C
    cond_inlet_temp = 25  # 20 - 30 deg C
    feed_temp = 25  # 20 - 30 deg C
    feed_salinity = 50  # 35 - 292 g/L
    initial_batch_volume = 50  # > 50 L
    recovery_ratio = 0.5  # -
    module_type = "AS26C7.2L"
    cooling_system_type = "open"
    cooling_inlet_temp = 25  # deg C, not required when cooling system type is "closed"

    # Create the flowsheet the the system specifications
    m = build_vagmd_flowsheet(
        feed_flow_rate=feed_flow_rate,
        evap_inlet_temp=evap_inlet_temp,
        cond_inlet_temp=cond_inlet_temp,
        feed_temp=feed_temp,
        feed_salinity=feed_salinity,
        recovery_ratio=recovery_ratio,
        initial_batch_volume=initial_batch_volume,
        module_type=module_type,
        cooling_system_type=cooling_system_type,
    )

    # Set scaling factors and initialize the flowsheet
    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,
    )

    results = solver.solve(m)

    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(75.5880, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(28.0441, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(1.4949, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(2.9398, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(122.115, abs=1e-3)
