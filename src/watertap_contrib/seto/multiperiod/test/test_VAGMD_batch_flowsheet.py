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
from watertap_contrib.seto.multiperiod.VAGMD_batch_flowsheet import build_vagmd_flowsheet, fix_dof_and_initialize

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
    m = build_vagmd_flowsheet(  feed_flow_rate=600,
                                evap_inlet_temp=80,
                                cond_inlet_temp=25,
                                feed_temp=25,
                                feed_salinity=35,
                                initial_batch_volume=50,
                                module_type="AS7C1.5L",
                                high_brine_salinity=False,
                                cooling_system_type="closed",)   

    fix_dof_and_initialize(m)

    return m

# @pytest.mark.unit
# def test_dof(build_vagmd_fs):
#     m = build_vagmd_fs
#     assert degrees_of_freedom(m) == 0
        

@pytest.mark.unit
def test_fs(build_vagmd_fs):
    m = build_vagmd_flowsheet(  feed_flow_rate=600,
                                evap_inlet_temp=80,
                                cond_inlet_temp=25,
                                feed_temp=25,
                                feed_salinity=35,
                                initial_batch_volume=50,
                                module_type="AS7C1.5L",
                                high_brine_salinity=False,
                                cooling_system_type="closed",)   

    fix_dof_and_initialize(m)
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert m.fs.vagmd.condenser_out_props[0].temperature.value - 273.15 == pytest.approx(69.4664, abs=1e-3)
    assert m.fs.vagmd.evaporator_out_props[0].temperature.value - 273.15 == pytest.approx(34.5742, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
    assert m.fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].value == pytest.approx(35, abs=1e-3)
    assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.3404, abs=1e-3)
    assert m.fs.vagmd.thermal_power.value == pytest.approx(7.0701, abs=1e-3)
    assert m.fs.acc_thermal_energy.value == pytest.approx(0, abs=1e-3)
    assert m.fs.acc_distillate_volume.value == pytest.approx(0, abs=1e-3)
    assert m.fs.dt.value == pytest.approx(33.921, abs=1e-3)