import pytest

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.solar_models.zero_order.flat_plate_physical import (
    FlatPlatePhysical,
)
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.costing import EnergyCosting

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

from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)

# Get default solver for testing
solver = get_solver()


class TestFlatPlatePhysical:
    @pytest.fixture(scope="class")
    def flat_plate_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = WaterParameterBlock()
        m.fs.flatplate = FlatPlatePhysical(property_package=m.fs.properties)

        # Define the model inputs
        m.fs.flatplate.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.045528)
        m.fs.flatplate.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.flatplate.inlet.temperature.fix(38.2 + 273.15)
        m.fs.flatplate.inlet.pressure.fix(101325)
        m.fs.flatplate.outlet.pressure.fix()

        m.fs.flatplate.G_total.fix(540)
        m.fs.flatplate.collector_area.fix(2.98)
        m.fs.flatplate.number_collectors.set_value(2)
        m.fs.flatplate.mdot_test.set_value(0.045528)

        m.fs.flatplate.T_amb.set_value(12 + 273.15)

        m.fs.flatplate.pump_power.set_value(45)
        m.fs.flatplate.pump_eff.set_value(0.85)
        m.fs.flatplate.ta.set_value(1)

        m.fs.flatplate.FR.set_value(0.689)
        m.fs.flatplate.UL.set_value(3.85 / 0.689)

        return m

    @pytest.mark.unit
    def test_config(self, flat_plate_frame):
        m = flat_plate_frame

        assert len(m.fs.flatplate.config) == 7
        assert not m.fs.flatplate.config.dynamic
        assert not m.fs.flatplate.config.has_holdup
        assert m.fs.flatplate.config.property_package is m.fs.properties
        assert m.fs.flatplate.config.energy_balance_type is EnergyBalanceType.useDefault
        assert (
            m.fs.flatplate.config.momentum_balance_type
            is MomentumBalanceType.pressureTotal
        )

        assert m.fs.flatplate._tech_type == "flat_plate"

    @pytest.mark.unit
    def test_build(self, flat_plate_frame):
        m = flat_plate_frame

        # Test ports
        port_list = ["inlet", "outlet"]
        for port_str in port_list:
            port = getattr(m.fs.flatplate, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # Test statistics
        assert number_variables(m.fs.flatplate) == 20
        assert number_unused_variables(m.fs.flatplate) == 5
        assert number_total_constraints(m.fs.flatplate) == 11

    @pytest.mark.unit
    def test_dof(self, flat_plate_frame):
        m = flat_plate_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, flat_plate_frame):
        m = flat_plate_frame
        calculate_scaling_factors(m)
        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_flat_plate_physical_model(self, flat_plate_frame):

        m = flat_plate_frame

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, flat_plate_frame):
        m = flat_plate_frame

        assert pytest.approx(3.97, rel=1e-3) == value(m.fs.flatplate.Fprime_UL)
        assert pytest.approx(1.00, rel=1e-3) == value(m.fs.flatplate.r)
        assert pytest.approx(1616.29, rel=1e-3) == value(
            m.fs.flatplate.heat
        )  # [W], with no pipe or heat exchanger losses
        assert pytest.approx(52.94, rel=1e-3) == value(
            m.fs.flatplate.electricity
        )  # [W]

    @pytest.mark.component
    def test_costing(self, flat_plate_frame):
        m = flat_plate_frame

        m.fs.test_flow = 50 * pyunits.Mgallons / pyunits.day

        m.fs.costing = EnergyCosting()
        m.fs.flatplate.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )

        m.fs.costing.factor_maintenance_labor_chemical.fix(0)
        m.fs.costing.factor_total_investment.fix(1)

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(flow_rate=m.fs.test_flow)

        results = solver.solve(m)
        assert_optimal_termination(results)
