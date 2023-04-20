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
        m.fs.liquid_prop = SeawaterParameterBlock()
        m.fs.vapor_prop = WaterParameterBlock()
        m.fs.lt_med = LTMEDSurrogate(
            property_package_liquid=m.fs.liquid_prop,
            property_package_vapor=m.fs.vapor_prop,
            number_effects=12,  # assuming 12 effects by default
        )

        lt_med = m.fs.lt_med
        feed = lt_med.feed_props[0]
        dist = lt_med.distillate_props[0]
        steam = lt_med.steam_props[0]

        # System specification
        # Input variable 1: Feed salinity (30-60 g/L = kg/m3)
        feed_salinity = 35 * pyunits.kg / pyunits.m**3  # g/L = kg/m3

        # Input variable 2: Feed temperature (15-35 deg C)
        feed_temperature = 25  # degC

        # Input variable 3: Heating steam temperature (60-85 deg C)
        steam_temperature = 80  # degC

        # Input variable 4: System capacity (> 2000 m3/day)
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day  # m3/day

        # Input variable 5: Recovery ratio (30%- 50%)
        recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless

        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )

        """
        Specify feed flow state properties
        """
        # Specify feed flow state properties
        lt_med.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temperature + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        # Fix input steam temperature
        steam.temperature.fix(steam_temperature + 273.15)

        # Fix target recovery rate
        lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)

        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "TDS")
        )
        m.fs.vapor_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.vapor_prop.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
        )

        return m

    @pytest.mark.unit
    def test_num_effects_domain(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.liquid_prop = SeawaterParameterBlock()
        m.fs.vapor_prop = WaterParameterBlock()

        error_msg = "invalid value for configuration 'number_effects'"
        with pytest.raises(ValueError, match=error_msg):
            m.fs.lt_med = LTMEDSurrogate(
                property_package_liquid=m.fs.liquid_prop,
                property_package_vapor=m.fs.vapor_prop,
                number_effects=15,
            )

    @pytest.mark.unit
    def test_config(self, LT_MED_frame):
        m = LT_MED_frame
        # check LT-MED config arguments
        assert len(m.fs.lt_med.config) == 6

        assert not m.fs.lt_med.config.dynamic
        assert not m.fs.lt_med.config.has_holdup
        assert m.fs.lt_med.config.property_package_liquid is m.fs.liquid_prop
        assert m.fs.lt_med.config.property_package_vapor is m.fs.vapor_prop
        assert m.fs.lt_med.config.number_effects in range(3, 15)

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
        assert number_variables(m) == 193
        assert number_total_constraints(m) == 52
        assert number_unused_variables(m) == 89  # vars from property package parameters

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

        feed_flow_m3_hr = 166.66
        dist_flow_m3_hr = 83.33
        brine_flow_m3_hr = 83.33
        cool_flow_m3_hr = 406.00

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
        assert value(
            lt_med.brine_props[0].flow_mass_phase_comp["Liq", "TDS"]
            - lt_med.feed_props[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(0, abs=1e-3)

    @pytest.mark.component
    def test_solution(self, LT_MED_frame):
        m = LT_MED_frame
        assert pytest.approx(9.9127, rel=1e-3) == value(m.fs.lt_med.gain_output_ratio)
        assert pytest.approx(3.9592, rel=1e-3) == value(
            m.fs.lt_med.specific_area_per_m3_day
        )
        assert pytest.approx(6.4290e1, rel=1e-3) == value(
            m.fs.lt_med.specific_energy_consumption_thermal
        )
        assert pytest.approx(5.3575e3, rel=1e-3) == value(
            m.fs.lt_med.thermal_power_requirement
        )
        assert pytest.approx(2.321, rel=1e-3) == value(
            m.fs.lt_med.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(406.00, rel=1e-3) == value(
            pyunits.convert(
                m.fs.lt_med.cooling_out_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )
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

        assert pytest.approx(2254.658, rel=1e-3) == value(
            m.fs.lt_med.costing.med_specific_cost
        )
        assert pytest.approx(4662455.768, rel=1e-3) == value(
            m.fs.lt_med.costing.capital_cost
        )
        assert pytest.approx(2705589.357, rel=1e-3) == value(
            m.fs.lt_med.costing.membrane_system_cost
        )
        assert pytest.approx(1956866.411, rel=1e-3) == value(
            m.fs.lt_med.costing.evaporator_system_cost
        )
        assert pytest.approx(208604.394, rel=1e-3) == value(
            m.fs.lt_med.costing.fixed_operating_cost
        )

        assert pytest.approx(1.58295, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(748697.447, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(4662455.768, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )

    @pytest.mark.parametrize("number_effects", [4, 5, 7, 8, 10, 11, 13])
    def test_interp_values(self, number_effects):
        # create flowsheet with interpolated values of number of effects
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.liquid_prop = SeawaterParameterBlock()
        m.fs.vapor_prop = WaterParameterBlock()
        m.fs.lt_med = LTMEDSurrogate(
            property_package_liquid=m.fs.liquid_prop,
            property_package_vapor=m.fs.vapor_prop,
            number_effects=number_effects,
        )

        lt_med = m.fs.lt_med
        feed = lt_med.feed_props[0]
        dist = lt_med.distillate_props[0]
        steam = lt_med.steam_props[0]

        # System specification
        # Input variable 1: Feed salinity (30-60 g/L = kg/m3)
        feed_salinity = 35 * pyunits.kg / pyunits.m**3  # g/L = kg/m3

        # Input variable 2: Feed temperature (15-35 deg C)
        feed_temperature = 25  # degC

        # Input variable 3: Heating steam temperature (60-85 deg C)
        steam_temperature = 80  # degC

        # Input variable 4: System capacity (> 2000 m3/day)
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day  # m3/day

        # Input variable 5: Recovery ratio (30%- 50%)
        recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless

        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )

        """
        Specify feed flow state properties
        """
        # Specify feed flow state properties
        lt_med.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temperature + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        # Fix input steam temperature
        steam.temperature.fix(steam_temperature + 273.15)

        # Fix target recovery rate
        lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)

        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "TDS")
        )
        m.fs.vapor_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
        )
        m.fs.vapor_prop.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
        )

        calculate_scaling_factors(m)
        # For 10 and 11 effects, specific area may need a new initialization value,
        # because the surrogate eqn for 9 and 12 effects are largely different
        if number_effects == 11:
            m.fs.lt_med.specific_area_per_m3_day.value = 3
        initialization_tester(m, unit=m.fs.lt_med, outlvl=idaeslog.DEBUG)
        results = solver.solve(m)

        # Check interpolated results for different number of effects: {number_effects: result}
        gain_output_ratios = {
            4: 3.584,
            5: 4.431,
            7: 6.073,
            8: 6.869,
            10: 8.414,
            11: 9.163,
            13: 10.626,
        }
        specific_areas = {
            4: 2.757,
            5: 2.507,
            7: 2.507,
            8: 2.758,
            10: 3.325,
            11: 3.642,
            13: 4.159,
        }

        assert pytest.approx(gain_output_ratios[number_effects], rel=1e-3) == value(
            m.fs.lt_med.gain_output_ratio
        )
        assert pytest.approx(specific_areas[number_effects], rel=1e-3) == value(
            m.fs.lt_med.specific_area_per_m3_day
        )
