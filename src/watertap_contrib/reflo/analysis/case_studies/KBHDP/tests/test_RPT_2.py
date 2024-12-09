#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.testing import initialization_tester
from idaes.core.util.exceptions import InitializationError, ConfigurationError
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
from watertap.costing.watertap_costing_package import (
    WaterTAPCostingData,
    WaterTAPCostingBlockData,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
)
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_2 import *
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


class TestRPT2:
    @pytest.fixture(scope="class")
    def RPT_2_frame(self):
        m = build_system(RE=True)
        add_connections(m)
        add_constraints(m)
        set_operating_conditions(m)
        # # apply_scaling(m)
        # # init_system(m)
        # # add_costing(m)
        # # scale_costing(m)

        return m

    @pytest.mark.unit
    def test_config(self, RPT_2_frame):
        m = RPT_2_frame
        treatment = m.fs.treatment
        energy = m.fs.energy

        assert (
            len(
                [
                    v
                    for v in treatment.component_data_objects(
                        ctype=Block, active=True, descend_into=False
                    )
                ]
            )
            == 21
        )
        assert (
            len(
                [
                    v
                    for v in energy.component_data_objects(
                        ctype=Block, active=True, descend_into=False
                    )
                ]
            )
            == 1
        )

        # assert isinstance(treatment.costing, TreatmentCosting)
        # assert isinstance(energy.costing, EnergyCosting)
        # assert isinstance(m.fs.costing, WaterTAPCostingBlockData)

    @pytest.mark.unit
    def test_dof(self, RPT_2_frame):
        m = RPT_2_frame

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_build(self, RPT_2_frame):
        m = RPT_2_frame

        assert number_variables(m) == 401
        assert number_total_constraints(m) == 198
        assert number_unused_variables(m) == 120

    @pytest.mark.component
    def test_initialize(self, RPT_2_frame):
        m = RPT_2_frame
        apply_scaling(m)
        init_system(m)

    @pytest.mark.component
    def test_scaling(self, RPT_2_frame):
        m = RPT_2_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, RPT_2_frame):
        m = RPT_2_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, RPT_2_frame):
        m = RPT_2_frame

        assert pytest.approx(value(m.fs.water_recovery), rel=1e-3) == 0.486

    #     assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.935

    @pytest.mark.component
    def test_costing(self, RPT_2_frame):
        m = RPT_2_frame
        add_costing(m)
        results = solver.solve(m)
        assert isinstance(m.fs.costing, WaterTAPCostingBlockData)

        assert value(m.fs.costing.utilization_factor) == 1
        assert value(m.fs.costing.plant_lifetime) == 20
        assert value(m.fs.costing.TIC) == 2.0
        assert value(m.fs.costing.plant_lifetime) == 20
        assert value(m.fs.costing.heat_cost_buy) == 0.01

        cost_results = {
            "total_capital_cost": 16275469.798,
            "total_operating_cost": 4880943.681138365,
            "aggregate_flow_electricity": 940,
            "aggregate_flow_heat": 19577,
            "total_electric_operating_cost": 577104,
            "total_heat_operating_cost": 1716205,
            "aggregate_flow_electricity_purchased": 940.49,
            "aggregate_flow_electricity_sold": 0,
            "aggregate_flow_heat_purchased": 19577.977,
            "aggregate_flow_heat_sold": 0,
            "frac_heat_from_grid": 0.983,
            "aggregate_capital_cost": 16275469.798,
            "aggregate_fixed_operating_cost": 691261.7801448428,
            "aggregate_variable_operating_cost": 159727,
        }

        for v, r in cost_results.items():
            cost_item = getattr(m.fs.costing, v)
            assert pytest.approx(value(cost_item), rel=1e-3) == r
