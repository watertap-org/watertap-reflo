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

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
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
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

from watertap.core.solvers import get_solver
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.unit_models.surrogate import MEDTVCSurrogate
from watertap_contrib.reflo.costing import REFLOCosting

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestMEDTVC:
    @pytest.fixture(scope="class")
    def MED_TVC_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.liquid_prop = SeawaterParameterBlock()
        m.fs.vapor_prop = WaterParameterBlock()
        m.fs.med_tvc = MEDTVCSurrogate(
            property_package_liquid=m.fs.liquid_prop,
            property_package_vapor=m.fs.vapor_prop,
            number_effects=12,  # assuming 12 effects by default
        )

        med_tvc = m.fs.med_tvc
        steam = med_tvc.heating_steam_props[0]
        motive = med_tvc.motive_steam_props[0]

        # System specification
        # Input variable 1: Feed salinity (30-60 g/L = kg/m3)
        feed_salinity = 35 * pyunits.kg / pyunits.m**3

        # Input variable 2: Feed temperature (25-35 deg C)
        feed_temperature = 25

        # Input variable 3: Motive steam pressure (4-45 bar)
        motive_pressure = 24

        # Input variable 4: System capacity (2,000 - 100,000 m3/day)
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day

        # Input variable 5: Recovery ratio (30%- 40%)
        recovery_ratio = 0.3 * pyunits.dimensionless

        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )  # feed volumetric flow rate [m3/s]

        """
        Specify feed flow state properties
        """
        # Specify feed flow state properties
        med_tvc.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temperature + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        """
        Specify heating steam state properties
        """
        # Heating steam temperature (saturated) is fixed at 70 C in this configuration
        steam.temperature.fix(70 + 273.15)

        # Calculate heating steam pressure (saturated)
        med_tvc.heating_steam_props.calculate_state(
            var_args={
                ("pressure_sat", None): value(steam.pressure),
            },
            hold_state=True,
        )
        # Release mass flow rate
        steam.flow_mass_phase_comp["Vap", "H2O"].unfix()
        steam.flow_mass_phase_comp["Liq", "H2O"].unfix()

        """
        Specify motive steam state properties
        """
        # Calculate temperature of the motive steam (saturated)
        med_tvc.motive_steam_props.calculate_state(
            var_args={
                ("pressure", None): motive_pressure * 1e5,
                ("pressure_sat", None): motive_pressure * 1e5,
            },
            hold_state=True,
        )
        # Release mass flow rate
        motive.flow_mass_phase_comp["Vap", "H2O"].unfix()
        motive.flow_mass_phase_comp["Liq", "H2O"].unfix()

        # Fix target recovery rate
        med_tvc.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)

        # Set scaling factors for mass flow rates
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
        m.fs.water_prop = SeawaterParameterBlock()
        m.fs.steam_prop = WaterParameterBlock()

        tested_number_effects = 50
        error_msg = f"The number of effects was specified as {tested_number_effects}. The number of effects should be specified as an integer between 8 to 16."
        with pytest.raises(ConfigurationError, match=error_msg):
            m.fs.lt_med = MEDTVCSurrogate(
                property_package_liquid=m.fs.water_prop,
                property_package_vapor=m.fs.steam_prop,
                number_effects=tested_number_effects,
            )

    @pytest.mark.unit
    def test_config(self, MED_TVC_frame):
        m = MED_TVC_frame
        # check unit config arguments
        assert len(m.fs.med_tvc.config) == 6

        assert not m.fs.med_tvc.config.dynamic
        assert not m.fs.med_tvc.config.has_holdup
        assert m.fs.med_tvc.config.property_package_liquid is m.fs.liquid_prop
        assert m.fs.med_tvc.config.property_package_vapor is m.fs.vapor_prop
        assert m.fs.med_tvc.config.number_effects in range(8, 17)

    @pytest.mark.unit
    def test_build(self, MED_TVC_frame):
        m = MED_TVC_frame

        # test ports
        port_lst = ["feed", "distillate", "brine", "steam", "motive"]
        for port_str in port_lst:
            port = getattr(m.fs.med_tvc, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 216
        assert number_total_constraints(m) == 61
        assert number_unused_variables(m) == 86  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, MED_TVC_frame):
        m = MED_TVC_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, MED_TVC_frame):
        m = MED_TVC_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, MED_TVC_frame):
        m = MED_TVC_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, MED_TVC_frame):
        m = MED_TVC_frame
        initialization_tester(m, unit=m.fs.med_tvc, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, MED_TVC_frame):
        m = MED_TVC_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_mass_balance(self, MED_TVC_frame):
        m = MED_TVC_frame

        med_tvc = m.fs.med_tvc

        feed_flow_m3_hr = 277.777
        dist_flow_m3_hr = 83.333
        brine_flow_m3_hr = 194.444
        cool_flow_m3_hr = 131.593

        feed_mass_flow_tot = 78.93
        cool_mass_flow_tot = 37.31
        feed_mass_flow_tds = 2.70
        brine_mass_flow_tds = 2.70
        recovery = dist_flow_m3_hr / feed_flow_m3_hr

        assert value(med_tvc.recovery_vol_phase[0, "Liq"]) == pytest.approx(
            recovery, rel=1e-3
        )
        assert value(
            pyunits.convert(
                med_tvc.feed_props[0].flow_vol_phase["Liq"]
                - med_tvc.distillate_props[0].flow_vol_phase["Liq"]
                - med_tvc.brine_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )
        ) == pytest.approx(
            feed_flow_m3_hr - dist_flow_m3_hr - brine_flow_m3_hr, rel=1e-3
        )
        assert value(med_tvc.feed_cool_mass_flow) == pytest.approx(
            feed_mass_flow_tot + cool_mass_flow_tot, rel=1e-2
        )  # mass flow calculated two different ways
        assert value(med_tvc.feed_cool_vol_flow) == pytest.approx(
            (feed_flow_m3_hr + cool_flow_m3_hr), rel=1e-3
        )
        assert value(
            med_tvc.brine_props[0].flow_mass_phase_comp["Liq", "TDS"]
            - med_tvc.feed_props[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(feed_mass_flow_tds - brine_mass_flow_tds, rel=1e-6)

        assert value(
            med_tvc.brine_props[0].flow_mass_phase_comp["Liq", "TDS"]
            - med_tvc.feed_props[0].flow_mass_phase_comp["Liq", "TDS"]
        ) == pytest.approx(0, abs=1e-3)

    @pytest.mark.component
    def test_solution(self, MED_TVC_frame):
        m = MED_TVC_frame

        assert pytest.approx(12.9102, rel=1e-3) == value(m.fs.med_tvc.gain_output_ratio)
        assert pytest.approx(5.1664, rel=1e-3) == value(
            m.fs.med_tvc.specific_area_per_m3_day
        )
        assert pytest.approx(53.1622, rel=1e-3) == value(
            m.fs.med_tvc.specific_energy_consumption_thermal
        )
        assert pytest.approx(4430.19, rel=1e-3) == value(
            m.fs.med_tvc.thermal_power_requirement
        )
        assert pytest.approx(2.6766, rel=1e-3) == value(
            m.fs.med_tvc.heating_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(1.2175, rel=1e-3) == value(
            m.fs.med_tvc.motive_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(131.593, rel=1e-3) == value(
            pyunits.convert(
                m.fs.med_tvc.cooling_out_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )

    @pytest.mark.component
    def test_costing(self, MED_TVC_frame):
        m = MED_TVC_frame
        med_tvc = m.fs.med_tvc
        dist = med_tvc.distillate_props[0]
        m.fs.costing = REFLOCosting()
        med_tvc.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.total_investment_factor.fix(1)
        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.wacc.unfix()
        m.fs.costing.capital_recovery_factor.fix(0.08764)

        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(dist.flow_vol_phase["Liq"])
        m.fs.costing.add_LCOW(dist.flow_vol_phase["Liq"])

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(2254.658, rel=1e-3) == value(
            m.fs.med_tvc.costing.med_specific_cost
        )
        assert pytest.approx(5126761.859, rel=1e-3) == value(
            m.fs.med_tvc.costing.capital_cost
        )
        assert pytest.approx(2705589.357, rel=1e-3) == value(
            m.fs.med_tvc.costing.membrane_system_cost
        )
        assert pytest.approx(2421172.502, rel=1e-3) == value(
            m.fs.med_tvc.costing.evaporator_system_cost
        )
        assert pytest.approx(239692.046, rel=1e-3) == value(
            m.fs.med_tvc.costing.fixed_operating_cost
        )

        assert pytest.approx(1.6905, rel=1e-3) == value(m.fs.costing.LCOW)

        assert pytest.approx(785633.993, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(5126761.859, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )

    @pytest.mark.parametrize("number_effects", [9, 11, 13, 15])
    def test_interp_values(self, number_effects):
        # create flowsheet for an interpolated number of effects
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.liquid_prop = SeawaterParameterBlock()
        m.fs.vapor_prop = WaterParameterBlock()
        m.fs.med_tvc = MEDTVCSurrogate(
            property_package_liquid=m.fs.liquid_prop,
            property_package_vapor=m.fs.vapor_prop,
            number_effects=number_effects,  # Interpolated values include [9, 11, 13, 15]
        )

        med_tvc = m.fs.med_tvc
        feed = med_tvc.feed_props[0]
        cool = med_tvc.cooling_out_props[0]
        dist = med_tvc.distillate_props[0]
        steam = med_tvc.heating_steam_props[0]
        motive = med_tvc.motive_steam_props[0]

        # System specification
        # Input variable 1: Feed salinity (30-60 g/L = kg/m3)
        feed_salinity = 35 * pyunits.kg / pyunits.m**3

        # Input variable 2: Feed temperature (25-35 deg C)
        feed_temperature = 25

        # Input variable 3: Motive steam pressure (4-45 bar)
        motive_pressure = 24

        # Input variable 4: System capacity (2,000 - 100,000 m3/day)
        sys_capacity = 2000 * pyunits.m**3 / pyunits.day

        # Input variable 5: Recovery ratio (30%- 40%)
        recovery_ratio = 0.3 * pyunits.dimensionless

        feed_flow = pyunits.convert(
            (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
        )  # feed volumetric flow rate [m3/s]

        """
        Specify feed flow state properties
        """
        # Specify feed flow state properties
        med_tvc.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temperature + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        """
        Specify heating steam state properties
        """
        # Heating steam temperature (saturated) is fixed at 70 C in this configuration
        steam.temperature.fix(70 + 273.15)

        # Calculate heating steam pressure (saturated)
        med_tvc.heating_steam_props.calculate_state(
            var_args={
                ("pressure_sat", None): value(steam.pressure),
            },
            hold_state=True,
        )
        # Release mass flow rate
        steam.flow_mass_phase_comp["Vap", "H2O"].unfix()
        steam.flow_mass_phase_comp["Liq", "H2O"].unfix()

        """
        Specify motive steam state properties
        """
        # Calculate temperature of the motive steam (saturated)
        med_tvc.motive_steam_props.calculate_state(
            var_args={
                ("pressure", None): motive_pressure * 1e5,
                ("pressure_sat", None): motive_pressure * 1e5,
            },
            hold_state=True,
        )
        # Release mass flow rate
        motive.flow_mass_phase_comp["Vap", "H2O"].unfix()
        motive.flow_mass_phase_comp["Liq", "H2O"].unfix()

        # Fix target recovery rate
        med_tvc.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)

        # Set scaling factors for mass flow rates
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
        initialization_tester(m, unit=m.fs.med_tvc, outlvl=idaeslog.DEBUG)
        results = solver.solve(m)

        # Check interpolated results for different number of effects: {number_effects: result}
        gain_output_ratios = {9: 10.299, 11: 12.067, 13: 13.709, 15: 15.218}
        specific_areas = {9: 3.514, 11: 4.682, 13: 5.797, 15: 7.689}
        heating_steam_mass_flow_rates = {9: 2.642, 11: 2.834, 13: 2.518, 15: 2.786}
        motive_steam_mass_flow_rates = {9: 1.162, 11: 1.334, 13: 1.152, 15: 1.407}

        assert pytest.approx(gain_output_ratios[number_effects], rel=1e-3) == value(
            m.fs.med_tvc.gain_output_ratio
        )
        assert pytest.approx(specific_areas[number_effects], rel=1e-3) == value(
            m.fs.med_tvc.specific_area_per_m3_day
        )
        assert pytest.approx(
            heating_steam_mass_flow_rates[number_effects], rel=1e-3
        ) == value(
            m.fs.med_tvc.heating_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(
            motive_steam_mass_flow_rates[number_effects], rel=1e-3
        ) == value(
            m.fs.med_tvc.motive_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
