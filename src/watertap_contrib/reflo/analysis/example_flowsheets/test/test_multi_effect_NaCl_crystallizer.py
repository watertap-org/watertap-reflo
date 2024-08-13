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
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    Objective,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap_contrib.reflo.analysis.example_flowsheets.multi_effect_NaCl_crystallizer import (
    build_fs_multi_effect_crystallizer,
    add_costings,
    multi_effect_crystallizer_initialization,
    get_model_performance,
)
import watertap_contrib.reflo.property_models.cryst_prop_pack as props
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing import TreatmentCosting, CrystallizerCostType

solver = get_solver()


class TestMultiEffectCrystallization:
    @pytest.fixture(scope="class")
    def MultiEffectCrystallizer_frame(self):
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
        # Negative value to represent salt recovery value ($/kg)
        m.fs.costing.crystallizer.NaCl_recovery_value.fix(-0.024)

        multi_effect_crystallizer_initialization(m)

        return m

    @pytest.mark.unit
    def test_dof(self, MultiEffectCrystallizer_frame):
        m = MultiEffectCrystallizer_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, MultiEffectCrystallizer_frame):
        m = MultiEffectCrystallizer_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, MultiEffectCrystallizer_frame):
        m = MultiEffectCrystallizer_frame

        data_table, overall_performance = get_model_performance(m)

        # Check solid mass in solids stream
        assert (
            pytest.approx(overall_performance["Capacity (m3/day)"], rel=1e-3) == 277.39
        )
        assert (
            pytest.approx(overall_performance["Feed brine salinity (g/L)"], rel=1e-3)
            == 166.3
        )
        assert (
            pytest.approx(overall_performance["Total brine disposed (kg/s)"], rel=1e-3)
            == 3.559
        )
        assert (
            pytest.approx(
                overall_performance["Total water production (kg/s)"], rel=1e-3
            )
            == 2.313
        )
        assert (
            pytest.approx(
                overall_performance["Total solids collected (kg/s)"], rel=1e-3
            )
            == 0.267
        )
        assert (
            pytest.approx(
                overall_performance["Total waste water remained (kg/s)"], rel=1e-3
            )
            == 0.979
        )
        assert (
            pytest.approx(
                overall_performance["Initial thermal energy consumption (kW)"], rel=1e-3
            )
            == 1704.47
        )
        assert (
            pytest.approx(overall_performance["Overall STEC (kWh/m3 feed)"], rel=1e-3)
            == 147.47
        )
        assert (
            pytest.approx(
                overall_performance["Total heat transfer area (m2)"], rel=1e-3
            )
            == 2.03
        )
        assert (
            pytest.approx(
                overall_performance["Levelized cost of feed brine ($/m3)"], rel=1e-3
            )
            == 1.719
        )

    @pytest.mark.component
    def test_optimization(self, MultiEffectCrystallizer_frame):
        m = MultiEffectCrystallizer_frame

        # optimization scenario
        m.fs.eff_1.pressure_operating.unfix()
        m.fs.eff_1.pressure_operating.setub(1.2 * 1e5)
        m.fs.eff_2.pressure_operating.unfix()
        m.fs.eff_3.pressure_operating.unfix()
        m.fs.eff_4.pressure_operating.unfix()
        m.fs.eff_4.pressure_operating.setlb(0.02 * 1e5)

        effs = [m.fs.eff_1, m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]
        total_area = sum(i.area for i in effs)

        # Optimize the heat transfer area
        m.fs.objective = Objective(expr=total_area)

        @m.Constraint(doc="Pressure decreasing")
        def pressure_bound1(b):
            return b.fs.eff_2.pressure_operating <= b.fs.eff_1.pressure_operating

        @m.Constraint(doc="Pressure decreasing")
        def pressure_bound2(b):
            return b.fs.eff_3.pressure_operating <= b.fs.eff_2.pressure_operating

        @m.Constraint(doc="Temperature difference")
        def temp_bound1(b):
            return (
                b.fs.eff_2.temperature_operating
                >= b.fs.eff_1.temperature_operating - 12
            )

        @m.Constraint(doc="Temperature difference")
        def temp_bound2(b):
            return (
                b.fs.eff_3.temperature_operating
                >= b.fs.eff_2.temperature_operating - 12
            )

        @m.Constraint(doc="Pressure decreasing")
        def pressure_bound3(b):
            return b.fs.eff_4.pressure_operating <= b.fs.eff_3.pressure_operating

        @m.Constraint(doc="Temperature difference")
        def temp_bound3(b):
            return (
                b.fs.eff_4.temperature_operating
                >= b.fs.eff_3.temperature_operating - 12
            )

        optimization_results = solver.solve(m, tee=False)
        assert (
            optimization_results.solver.termination_condition
            == TerminationCondition.optimal
        )
        data_table2, overall_performance2 = get_model_performance(m)

        assert (
            pytest.approx(
                data_table2["Effect 1"]["Operating temperature (C)"], rel=1e-3
            )
            == 113.86
        )
        assert (
            pytest.approx(
                data_table2["Effect 2"]["Operating temperature (C)"], rel=1e-3
            )
            == 101.86
        )
        assert (
            pytest.approx(
                data_table2["Effect 3"]["Operating temperature (C)"], rel=1e-3
            )
            == 89.86
        )
        assert (
            pytest.approx(
                data_table2["Effect 4"]["Operating temperature (C)"], rel=1e-3
            )
            == 77.86
        )
