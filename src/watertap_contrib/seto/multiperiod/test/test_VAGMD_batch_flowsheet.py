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
from watertap_contrib.seto.multiperiod.VAGMD_batch_flowsheet import (
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
    feed_flow_rate=600
    evap_inlet_temp=80
    cond_inlet_temp=25
    feed_temp=25
    feed_salinity=100
    recovery_ratio=0.5
    initial_batch_volume=50
    module_type="AS7C1.5L"
    cooling_system_type="closed"

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

    fix_dof_and_initialize(
        m,
        feed_flow_rate=feed_flow_rate,
        feed_salinity=feed_salinity,
        feed_temp=feed_temp,)

    return m


# @pytest.mark.unit
# def test_dof(build_vagmd_fs):
#     m = build_vagmd_fs
#     assert degrees_of_freedom(m) == 0


@pytest.mark.unit
def test_fs(build_vagmd_fs):

    m = build_vagmd_fs
    assert degrees_of_freedom(m) == 0

    unscaled_var_list = list(unscaled_variables_generator(m))
    assert len(unscaled_var_list) == 0

    # check that all constraints have been scaled
    unscaled_constraint_list = list(unscaled_constraints_generator(m))
    assert len(unscaled_constraint_list) == 0

    results = solver.solve(m)


    print('TCO: ', m.fs.vagmd.condenser_out_props[0].temperature.value -273.15)
    print('TEO: ', m.fs.vagmd.evaporator_out_props[0].temperature.value-273.15)
    print('ThPower: ', m.fs.vagmd.thermal_power.value)
    print('FFR',pyunits.convert(m.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].value * pyunits.m**3/ pyunits.s, to_units=pyunits.L / pyunits.h))
    print('CpF',m.fs.vagmd.avg_feed_props[0].cp_mass_phase["Liq"].value)
    print('RhoF',m.fs.vagmd.avg_feed_props[0].dens_mass_phase["Liq"].value)
    assert m.fs.vagmd.condenser_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(65.6459, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[
        0
    ].temperature.value - 273.15 == pytest.approx(37.4374, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(
        25, abs=1e-3
    )
    assert m.fs.vagmd.feed_props[0].conc_mass_phase_comp[
        "Liq", "TDS"
    ].value == pytest.approx(100, abs=1e-3)
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(7.6420, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(7.0701, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(33.921, abs=1e-3)



    assert False