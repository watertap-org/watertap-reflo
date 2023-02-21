import pytest
from pyomo.environ import Var, value, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core.solvers import get_solver
import watertap_contrib.seto.analysis.net_metering.PV_RO as PV_RO
from idaes.core.util.testing import initialization_tester
from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from watertap.unit_models.reverse_osmosis_0D import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

import idaes.core.util.model_statistics as stats

import idaes.logger as idaeslog


class TestPVRO:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = PV_RO.build_ro_pv()
        return m

    @pytest.mark.unit
    def test_config(self, system_frame):
        m = system_frame
        assert hasattr(m.fs, "energy")
        assert hasattr(m.fs, "treatment")

        for component in ["feed", "product", "disposal", "p1", "ro", "erd"]:
            assert hasattr(m.fs.treatment, component)
            unit = getattr(m.fs.treatment, component)
            assert not unit.config.dynamic
            assert not unit.config.has_holdup
            if component == "feed":
                assert unit.properties[0].flow_vol_phase
            if component == "ro":
                assert unit.config.property_package is m.fs.properties
                assert (
                    unit.config.material_balance_type == MaterialBalanceType.useDefault
                )
                assert unit.config.energy_balance_type == EnergyBalanceType.useDefault
                assert (
                    unit.config.momentum_balance_type
                    == MomentumBalanceType.pressureTotal
                )
                assert unit.config.pressure_change_type == PressureChangeType.calculated
                assert (
                    unit.config.concentration_polarization_type
                    == ConcentrationPolarizationType.calculated
                )
                assert (
                    unit.config.mass_transfer_coefficient
                    == MassTransferCoefficient.calculated
                )
                assert unit.config.has_pressure_change
                assert isinstance(unit.feed_side.deltaP, Var)
                assert isinstance(unit.deltaP, Var)

    @pytest.mark.unit
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        PV_RO.set_operating_conditions(m)

        assert stats.number_unused_variables(m) == 2
        assert (
            stats.number_fixed_variables(
                m.fs.treatment.p1.control_volume.properties_out[0]
            )
            == 1
        )
        assert stats.number_total_blocks(m.fs.treatment) == 29
        assert stats.number_total_blocks(m.fs.energy) == 2

    @pytest.mark.unit
    def test_initialization(self, system_frame):
        m = system_frame
        # PV_RO.set_operating_conditions(m)
        solver = get_solver()
        optarg = solver.options
        m.fs.treatment.feed.initialize(optarg=optarg)
        PV_RO.initialize_treatment(m)
        PV_RO.initialize_energy(m)

        assert stats.degrees_of_freedom(m) == 0
        assert stats.number_unused_variables(m) == 2

        for component in ["feed", "product", "disposal", "p1", "ro", "erd"]:
            assert hasattr(m.fs.treatment, component)
            unit = getattr(m.fs.treatment, component)
            assert_units_consistent(unit)

        initialization_tester(m, unit=m.fs.treatment.ro, dof=0, outlvl=idaeslog.DEBUG)

    @pytest.mark.unit
    def test_costing(self, system_frame):
        m = system_frame
        PV_RO.add_costing(m)
        PV_RO.fix_pv_costing(m)
        PV_RO.fix_treatment_global_params(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        PV_RO.optimize_setup(m, m.fs.sys_costing.LCOW)
        results = PV_RO.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, system_frame):
        m = system_frame

        solver = get_solver()
        results = solver.solve(m)
        m.fs.treatment.ro.area.fix()
        m.fs.treatment.ro.feed_side.velocity[0, 0].fix()

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert stats.degrees_of_freedom(m) == 0
        assert pytest.approx(80853.054, rel=1e-2) == value(
            m.fs.sys_costing.annual_energy_generated
        )
        assert pytest.approx(218541.469, rel=1e-2) == value(
            m.fs.sys_costing.total_capital_cost
        )
        assert pytest.approx(31770.751, rel=1e-2) == value(
            m.fs.sys_costing.total_operating_cost
        )
        assert pytest.approx(-4212.056, rel=1e-2) == value(
            m.fs.energy.costing.total_operating_cost
        )
        assert pytest.approx(28510.665, rel=1e-2) == value(
            m.fs.energy.costing.total_capital_cost
        )
        assert pytest.approx(35982.808, rel=1e-2) == value(
            m.fs.treatment.costing.total_operating_cost
        )
        assert pytest.approx(190030.804, rel=1e-2) == value(
            m.fs.treatment.costing.total_capital_cost
        )
        assert pytest.approx(31771, rel=1e-2) == value(
            m.fs.sys_costing.total_operating_cost
        )
        assert pytest.approx(0.31493, rel=1e-2) == value(m.fs.sys_costing.LCOW)
        assert pytest.approx(0.0552, rel=1e-2) == value(m.fs.sys_costing.LCOE)
        assert pytest.approx(1.5062, rel=1e-1) == value(
            m.fs.sys_costing.specific_electric_energy_consumption
        )
