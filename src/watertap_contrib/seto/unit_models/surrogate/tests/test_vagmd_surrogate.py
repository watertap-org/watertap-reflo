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


class TestVAGMD:
    @pytest.fixture(scope="class")
    def VAGMD_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        feed_flow_rate=600
        evap_inlet_temp=80
        cond_inlet_temp=25
        feed_temp=25
        feed_salinity=35
        initial_batch_volume=50
        module_type="AS7C1.5L"
        high_brine_salinity=False
        cooling_system_type="closed"

        m.fs.vagmd = VAGMDSurrogate(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # System specification
        # Fix the model inputs

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow_rate / 1000 / 3600,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)
        # m.fs.vagmd.dt.fix(20352.55 / feed_flow_rate)

        if cooling_system_type == "closed":  # TODO: update closed cooling
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_frame):
        m = VAGMD_frame
        assert degrees_of_freedom(m) == 0
        
    @pytest.mark.unit
    def test_calculate_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, VAGMD_frame):
        m = VAGMD_frame
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)

        # assert False

    @pytest.mark.component
    def test_solve(self, VAGMD_frame):
        m = VAGMD_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

