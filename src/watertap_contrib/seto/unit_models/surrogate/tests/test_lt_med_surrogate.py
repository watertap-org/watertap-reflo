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
from watertap_contrib.seto.unit_models.surrogate import LTMEDSurrogate

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.solvers import get_solver
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

# Get default solver for testing
solver = get_solver()


class TestLTMED:
    @pytest.fixture(scope="class")
    def LT_MED_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.water_prop = SeawaterParameterBlock()
        m.fs.steam_prop = WaterParameterBlock()
        m.fs.lt_med = LTMEDSurrogate(
            property_package_water=m.fs.water_prop,
            property_package_steam=m.fs.steam_prop,
        )
        lt_med = m.fs.lt_med
        feed = lt_med.feed_props[0]
        dist = lt_med.distillate_props[0]
        steam = lt_med.steam_props[0]

        # System specification
        feed_salinity = 35  # g/L = kg/m3
        feed_temperature = 25  # degC
        steam_temperature = 80  # degC
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day  # m3/day
        recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless
        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )

        feed.conc_mass_phase_comp["Liq", "TDS"].fix(feed_salinity)
        feed.flow_vol_phase["Liq"].fix(feed_flow)
        feed.temperature.fix(feed_temperature + 273.15)
        steam.temperature.fix(steam_temperature + 273.15)
        # flow rate of liquid steam is zero
        steam.flow_mass_phase_comp["Liq", "H2O"].fix(0)
        dist.flow_mass_phase_comp["Liq", "TDS"].fix(0)  # salinity in distillate is zero

        lt_med.recovery_ratio.fix(recovery_ratio)

        return m

    @pytest.mark.unit
    def test_config(self, LT_MED_frame):
        m = LT_MED_frame
        # check LT-MED config arguments
        assert len(m.fs.lt_med.config) == 5

        assert not m.fs.lt_med.config.dynamic
        assert not m.fs.lt_med.config.has_holdup
        assert m.fs.lt_med.config.property_package_water is m.fs.water_prop
        assert m.fs.lt_med.config.property_package_steam is m.fs.steam_prop

    @pytest.mark.unit
    def test_num_effects_domain(self, LT_MED_frame):
        m = LT_MED_frame
        error_msg = re.escape(
            "Invalid parameter value: fs.lt_med.number_effects[None] = '100', value type=<class 'int'>.\n\tValue not in parameter domain fs.lt_med.number_effects_domain"
        )
        with pytest.raises(ValueError, match=error_msg):
            m.fs.lt_med.number_effects.set_value(100)

    @pytest.mark.unit
    def test_build(self, LT_MED_frame):
        m = LT_MED_frame

        # test ports
        port_lst = ["feed", "distillate", "brine", "steam"]
        for port_str in port_lst:
            port = getattr(m.fs.lt_med, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 178
        assert number_total_constraints(m) == 47
        assert number_unused_variables(m) == 90  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, LT_MED_frame):
        m = LT_MED_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, LT_MED_frame):
        m = LT_MED_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, LT_MED_frame):
        m = LT_MED_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, LT_MED_frame):
        m = LT_MED_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_mass_balance(self, LT_MED_frame):
        m = LT_MED_frame

        lt_med = m.fs.lt_med

        feed_flow_m3_hr = 166.66
        dist_flow_m3_hr = 83.33
        brine_flow_m3_hr = 83.33
        cool_flow_m3_hr = 394.56

        feed_mass_flow_tot = 47.359
        cool_mass_flow_tot = 111.87
        feed_mass_flow_tds = 1.62037
        brine_mass_flow_tds = 1.62037
        recovery = dist_flow_m3_hr / feed_flow_m3_hr

        assert value(lt_med.recovery_ratio) == pytest.approx(recovery, rel=1e-3)
        assert value(
            pyunits.convert(
                lt_med.feed_props[0].flow_vol_phase["Liq"]
                - lt_med.distillate_props[0].flow_vol_phase["Liq"]
                - lt_med.brine_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )
        ) == pytest.approx(
            feed_flow_m3_hr - dist_flow_m3_hr - brine_flow_m3_hr, rel=1e-3
        )
        assert value(lt_med.feed_cool_mass_flow) == pytest.approx(
            feed_mass_flow_tot + cool_mass_flow_tot, rel=1e-2
        )  # mass flow calculated two different ways
        assert value(lt_med.feed_cool_vol_flow) == pytest.approx(
            (feed_flow_m3_hr + cool_flow_m3_hr), rel=1e-3
        )
        assert value(
            lt_med.brine_props[0].flow_mass_phase_comp["Liq", "TDS"]
            - lt_med.feed_props[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(feed_mass_flow_tds - brine_mass_flow_tds, rel=1e-6)

    @pytest.mark.component
    def test_solution(self, LT_MED_frame):
        m = LT_MED_frame
        assert pytest.approx(9.9127, rel=1e-3) == value(m.fs.lt_med.gain_output_ratio)
        assert pytest.approx(3.9592, rel=1e-3) == value(m.fs.lt_med.specific_area)
        assert pytest.approx(6.4290e1, rel=1e-3) == value(
            m.fs.lt_med.specific_thermal_energy_consumption
        )
        assert pytest.approx(5.3575e3, rel=1e-3) == value(
            m.fs.lt_med.thermal_power_requirement
        )
        assert pytest.approx(2.3211, rel=1e-3) == value(
            m.fs.lt_med.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )

    @pytest.mark.component
    def test_costing(self, LT_MED_frame):
        m = LT_MED_frame
        lt_med = m.fs.lt_med
        dist = lt_med.distillate_props[0]
        m.fs.costing = SETOWaterTAPCosting()
        lt_med.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        # Fix some global costing params for better comparison to Pyomo model
        # (This is not necessary in general)
        m.fs.costing.factor_total_investment.fix(1)
        m.fs.costing.factor_maintenance_labor_chemical.fix(0)
        m.fs.costing.factor_capital_annualization.fix(0.08764)

        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(dist.flow_vol_phase["Liq"])
        m.fs.costing.add_LCOW(dist.flow_vol_phase["Liq"])

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(2254.6578, rel=1e-3) == value(
            m.fs.lt_med.costing.med_specific_cost
        )
        assert pytest.approx(4662455.767, rel=1e-3) == value(
            m.fs.lt_med.costing.capital_cost
        )
        assert pytest.approx(2705589.356, rel=1e-3) == value(
            m.fs.lt_med.costing.membrane_system_cost
        )
        assert pytest.approx(1956866.4109, rel=1e-3) == value(
            m.fs.lt_med.costing.evaporator_system_cost
        )
        assert pytest.approx(208604.394, rel=1e-3) == value(
            m.fs.lt_med.costing.fixed_operating_cost
        )

        assert pytest.approx(1.58427, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(748697.4468, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(4662455.7674, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
