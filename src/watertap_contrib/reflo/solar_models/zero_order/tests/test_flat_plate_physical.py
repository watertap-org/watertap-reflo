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
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
    EnergyBalanceType,
    MomentumBalanceType,
)
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
)

from watertap.core.solvers import get_solver
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.core import SolarModelType
from watertap_contrib.reflo.solar_models.zero_order.flat_plate_physical import (
    FlatPlatePhysical,
)

# Get default solver for testing
solver = get_solver()


class TestFlatPlatePhysical:
    @pytest.fixture(scope="class")
    def flat_plate_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = WaterParameterBlock()
        m.fs.flatplate = FlatPlatePhysical(
            property_package=m.fs.properties, solar_model_type=SolarModelType.physical
        )

        # Define the model inputs
        m.fs.flatplate.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.045528)
        m.fs.flatplate.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.flatplate.inlet.temperature.fix(38.2 + 273.15)
        m.fs.flatplate.inlet.pressure.fix(101325)
        m.fs.flatplate.outlet.pressure.fix(101325)

        m.fs.flatplate.total_irradiance.fix(540)
        m.fs.flatplate.collector_area.fix(2.98)
        m.fs.flatplate.number_collectors.set_value(2)
        m.fs.flatplate.mdot_test.set_value(0.045528)

        m.fs.flatplate.temperature_ambient.set_value(12 + 273.15)

        m.fs.flatplate.pump_power.set_value(45)
        m.fs.flatplate.pump_eff.set_value(0.85)
        m.fs.flatplate.trans_absorb_prod.set_value(1)

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / 0.045528, index=("Liq", "H2O")
        )

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
        )

        return m

    @pytest.mark.unit
    def test_config(self, flat_plate_frame):
        m = flat_plate_frame

        assert len(m.fs.flatplate.config) == 15
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
        assert number_variables(m.fs.flatplate) == 18
        assert number_unused_variables(m.fs.flatplate) == 5
        assert number_total_constraints(m.fs.flatplate) == 9

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
    def test_initialization(self, flat_plate_frame):
        m = flat_plate_frame
        initialization_tester(m, unit=m.fs.flatplate)

    @pytest.mark.component
    def test_flat_plate_physical_model(self, flat_plate_frame):
        m = flat_plate_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, flat_plate_frame):
        m = flat_plate_frame
        fpc_results = {
            "electricity": 52.9411,
            "heat": 1616.26,
            "collector_area": 2.98,
            "collector_area_total": 5.96,
            "total_irradiance": 540.0,
            "Fprime_UL": 3.97081,
            "ratio_FRta": 0.99998311,
            "net_heat_gain": {0.0: 1616.26},
            "heat_load": 3.41806,
            "heat_annual": 14168179.94,
        }
        for v, r in fpc_results.items():
            cv = getattr(m.fs.flatplate, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(cv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)

    @pytest.mark.component
    def test_costing(self, flat_plate_frame):

        fpc_costing_results = {
            "capital_cost": 3760.6909,
            "fixed_operating_cost": 54.68,
            "direct_capital_cost": 3576.0,
            "indirect_capital_cost": 5.89,
            "sales_tax": 178.7999,
            "land_area": 0.00147274,
        }

        sys_costing_results = {
            "aggregate_capital_cost": 3760.69,
            "aggregate_fixed_operating_cost": 54.68,
            "aggregate_flow_heat": -1616.2651,
            "aggregate_flow_electricity": 52.9411,
            "aggregate_flow_costs": {"heat": 0, "electricity": 38136.16},
            "total_capital_cost": 3760.69,
            "total_operating_cost": 38190.85,
            "aggregate_direct_capital_cost": 3576.0,
            "LCOW": 2.792657,
        }

        m = flat_plate_frame

        m.fs.test_flow = 0.01 * pyunits.Mgallons / pyunits.day

        m.fs.costing = EnergyCosting()
        m.fs.costing.heat_cost.set_value(0)
        m.fs.flatplate.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )

        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.total_investment_factor.fix(1)

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(flow_rate=m.fs.test_flow)

        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_costing_results.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(cv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)

        for v, r in fpc_costing_results.items():
            cv = getattr(m.fs.flatplate.costing, v)
            assert pytest.approx(r, rel=1e-3) == value(cv)
