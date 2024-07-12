#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.unit_models.zero_order.battery import BatteryStorage

# Get default solver for testing
solver = get_solver()


class TestBattery:
    @pytest.fixture(scope="class")
    def battery_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.battery = BatteryStorage()

        m.fs.battery.nameplate_power.fix(100)
        m.fs.battery.nameplate_energy.fix(1000)
        m.fs.battery.elec_in.fix(100)
        m.fs.battery.elec_out.fix(0)
        m.fs.battery.initial_state_of_charge.fix(0)
        m.fs.battery.initial_energy_throughput.fix(0)

        return m

    @pytest.mark.unit
    def test_config(self, battery_frame):
        m = battery_frame
        # check Chemical softening config arguments

        assert len(m.fs.battery.config) == 2

        assert not m.fs.battery.config.dynamic
        assert not m.fs.battery.config.has_holdup

    @pytest.mark.unit
    def test_build(self, battery_frame):
        # Check battery model
        m = battery_frame

        # Test ports
        port_list = ["power_in", "power_out"]
        for port_str in port_list:
            port = getattr(m.fs.battery, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 1

        # test statistics
        assert number_variables(m) == 8
        assert number_total_constraints(m) == 5
        assert number_unused_variables(m) == 0

    @pytest.mark.unit
    def test_dof(self, battery_frame):
        m = battery_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, battery_frame):
        m = battery_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        if len(unscaled_var_list) > 0:
            for v in unscaled_var_list:
                print(v.name)
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, battery_frame):
        m = battery_frame
        initialization_tester(m, unit=m.fs.battery, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, battery_frame):
        m = battery_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        if len(badly_scaled_var_lst) > 0:
            for i in badly_scaled_var_generator(m):
                print(i[0].name, i[1])
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, battery_frame):
        m = battery_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, battery_frame):
        m = battery_frame

        assert pytest.approx(value(m.fs.battery.state_of_charge[0]), rel=1e-3) == 95
