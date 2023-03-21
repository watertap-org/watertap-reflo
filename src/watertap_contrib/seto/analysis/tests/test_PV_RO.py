import pytest
from pyomo.environ import Var, value, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
import watertap_contrib.seto.analysis.net_metering.PV_RO as PV_RO
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

        PV_RO.initialize_system(m)

        # TODO: fix organization of flowsheet and test; e.g., dof should be 0 before initialization
        assert stats.degrees_of_freedom(m) == 0
        assert stats.number_unused_variables(m) == 2

        for component in ["feed", "product", "disposal", "p1", "ro", "erd"]:
            assert hasattr(m.fs.treatment, component)
            unit = getattr(m.fs.treatment, component)
            assert_units_consistent(unit)

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

        assert pytest.approx(0.330, rel=1e-2) == value(m.fs.sys_costing.LCOW)
        assert pytest.approx(0.083, rel=1e-2) == value(m.fs.sys_costing.LCOE)
        assert pytest.approx(1.51, rel=1e-1) == value(
            m.fs.sys_costing.specific_electric_energy_consumption
        )
        assert pytest.approx(247052, rel=1e2) == value(
            m.fs.sys_costing.total_capital_cost
        )
        assert pytest.approx(31771, rel=1e2) == value(
            m.fs.sys_costing.total_operating_cost
        )
