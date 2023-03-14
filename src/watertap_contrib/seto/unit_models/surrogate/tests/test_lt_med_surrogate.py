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

from idaes.core.util.testing import initialization_tester
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
import idaes.logger as idaeslog

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
        feed_salinity = 35 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        feed_dens = 1000 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        feed_temperature = 25  # degC
        steam_temperature = 80  # degC
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day  # m3/day
        recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless
        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )

        feed.flow_mass_phase_comp["Liq", "TDS"].fix(feed_salinity * feed_flow)
        feed.flow_mass_phase_comp["Liq", "H2O"].fix(feed_dens * feed_flow)
        feed.temperature.fix(feed_temperature + 273.15)
        feed.pressure.fix(101325)
        steam.temperature.fix(steam_temperature + 273.15)
        # flow rate of liquid steam is zero
        steam.flow_mass_phase_comp["Liq", "H2O"].fix(0)
        dist.flow_mass_phase_comp["Liq", "TDS"].fix(0)  # salinity in distillate is zero

        lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
        m.fs.water_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.water_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "TDS")
        )
        m.fs.steam_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.steam_prop.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
        )

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
        assert number_variables(m) == 190
        assert number_total_constraints(m) == 50
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
    def test_initialize(self, LT_MED_frame):
        m = LT_MED_frame
        initialization_tester(m, unit=m.fs.lt_med, outlvl=idaeslog.DEBUG)

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

        feed_flow_m3_hr = 168.6778
        dist_flow_m3_hr = 84.3389
        brine_flow_m3_hr = 84.3389
        cool_flow_m3_hr = 410.7059

        feed_mass_flow_tot = 47.91666
        cool_mass_flow_tot = 116.4225
        feed_mass_flow_tds = 1.62037
        brine_mass_flow_tds = 1.62037
        recovery = dist_flow_m3_hr / feed_flow_m3_hr

        assert value(lt_med.recovery_vol_phase[0, "Liq"]) == pytest.approx(
            recovery, rel=1e-3
        )
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
            m.fs.lt_med.specific_energy_consumption_thermal
        )
        assert pytest.approx(5.4222e3, rel=1e-3) == value(
            m.fs.lt_med.thermal_power_requirement
        )
        assert pytest.approx(2.3490, rel=1e-3) == value(
            m.fs.lt_med.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )

    @pytest.mark.component
    def test_costing(self, LT_MED_frame):
        m = LT_MED_frame
        lt_med = m.fs.lt_med
        dist = lt_med.distillate_props[0]
        m.fs.costing = SETOWaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020
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

        assert pytest.approx(2251.009, rel=1e-3) == value(
            m.fs.lt_med.costing.med_specific_cost
        )
        assert pytest.approx(4710164.3678, rel=1e-3) == value(
            m.fs.lt_med.costing.capital_cost
        )
        assert pytest.approx(2733806.960, rel=1e-3) == value(
            m.fs.lt_med.costing.membrane_system_cost
        )
        assert pytest.approx(1976357.40761, rel=1e-3) == value(
            m.fs.lt_med.costing.evaporator_system_cost
        )
        assert pytest.approx(210907.779, rel=1e-3) == value(
            m.fs.lt_med.costing.fixed_operating_cost
        )

        assert pytest.approx(1.58295, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(757501.748, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(4710164.367, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
