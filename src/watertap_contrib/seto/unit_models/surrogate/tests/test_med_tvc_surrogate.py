import pytest
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    Constraint,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock
from watertap_contrib.seto.unit_models.surrogate import MEDTVCSurrogate
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from watertap.property_models.seawater_ion_generic import configuration
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap.core.util.initialization import assert_no_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestMEDTVC:
    @pytest.fixture(scope="class")
    def MED_TVC_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties1 = SeawaterParameterBlock()
        m.fs.properties2 = WaterParameterBlock()
        m.fs.unit = MEDTVCSurrogate(
            property_package=m.fs.properties1, property_package2=m.fs.properties2
        )

        # System specification

        feed_salinity = 35  # g/L
        feed_temperature = 25  # degC
        motive_pressure = 24  # bar
        sys_capacity = 2000  # m3/day
        recovery_rate = 0.3  # dimensionless

        m.fs.unit.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(feed_salinity)
        m.fs.unit.feed_props[0].temperature.fix(feed_temperature + 273.15)
        m.fs.unit.motive_props[0].pressure.fix(motive_pressure * 1e5)
        m.fs.unit.motive_props[0].pressure_sat.fix(motive_pressure * 1e5)
        m.fs.unit.Capacity.fix(sys_capacity)
        m.fs.unit.RR.fix(recovery_rate)

        return m

    @pytest.mark.unit
    def test_config(self, MED_TVC_frame):
        m = MED_TVC_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 5

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties1
        assert m.fs.unit.config.property_package2 is m.fs.properties2

    @pytest.mark.unit
    def test_build(self, MED_TVC_frame):
        m = MED_TVC_frame

        # test ports
        port_lst = ["feed", "distillate", "brine"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert isinstance(port, Port)
            # number of state variables for seawater property package
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 191
        assert number_total_constraints(m) == 61
        assert number_unused_variables(m) == 79  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, MED_TVC_frame):
        m = MED_TVC_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, MED_TVC_frame):
        m = MED_TVC_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, MED_TVC_frame):
        m = MED_TVC_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, MED_TVC_frame):
        m = MED_TVC_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, MED_TVC_frame):
        m = MED_TVC_frame
        assert pytest.approx(12.9102, rel=1e-3) == value(m.fs.unit.GOR)
        assert pytest.approx(5.1664, rel=1e-3) == value(m.fs.unit.sA)
        assert pytest.approx(5.3603e1, rel=1e-3) == value(m.fs.unit.STEC)
        assert pytest.approx(4.46695e3, rel=1e-3) == value(m.fs.unit.P_req)
        assert pytest.approx(2.6766, rel=1e-3) == value(
            m.fs.unit.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(1.2175, rel=1e-3) == value(
            m.fs.unit.motive_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(1.2621e2, rel=1e-3) == value(m.fs.unit.q_cooling)
