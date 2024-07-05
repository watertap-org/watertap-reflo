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
    assert_optimal_termination,
    units as pyunits,
)

from idaes.core import FlowsheetBlock
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
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.unit_models.surrogate import VAGMDSurrogateBase

# Get default solver for testing
solver = get_solver()


class TestVAGMD_unit_model:
    @pytest.fixture(scope="class")
    def VAGMD_frame(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with low brine salinity (< 173.5 g/L) and closed cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600  # 400 - 1100 L/h
        evap_inlet_temp = 80  # 60 - 80 deg C
        cond_inlet_temp = 25  # 20 - 30 deg C
        feed_temp = 25  # 20 - 30 deg C
        feed_salinity = 35  # 35 - 292 g/L
        module_type = "AS7C1.5L"
        cooling_system_type = "closed"
        high_brine_salinity = False  # True if brine salinity > 175.3 g/L
        cooling_inlet_temp = (
            25  # deg C, not required when cooling system type is "closed"
        )

        m.fs.vagmd = VAGMDSurrogateBase(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # Run helper function to determine the salinity mode
        (
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            cooling_system_type,
        ) = m.fs.vagmd._determine_salinity_mode(
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            module_type,
            high_brine_salinity,
            cooling_system_type,
        )

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        # Specify evaporator inlet temperature
        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        # Closed circuit, in which TCI is forced to be constant and the cooling water temperature can be adjusted.
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        # Open circuit, in which cooling is available at a constant water temperature and condenser inlet temperature varies.
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        return m

    @pytest.mark.unit
    def test_config(self, VAGMD_frame):
        m = VAGMD_frame
        # check VAGMD config arguments
        assert len(m.fs.vagmd.config) == 8

        assert not m.fs.vagmd.config.dynamic
        assert not m.fs.vagmd.config.has_holdup
        assert m.fs.vagmd.config.property_package_seawater is m.fs.seawater_properties
        assert m.fs.vagmd.config.property_package_water is m.fs.water_properties
        assert m.fs.vagmd.config.module_type in ["AS7C1.5L", "AS26C7.2L"]
        assert m.fs.vagmd.config.cooling_system_type in ["open", "closed"]

    @pytest.mark.unit
    def test_build(self, VAGMD_frame):
        m = VAGMD_frame

        # test statistics
        assert number_variables(m) == 243
        assert number_total_constraints(m) == 72
        assert number_unused_variables(m) == 139

    @pytest.mark.unit
    def test_dof(self, VAGMD_frame):
        m = VAGMD_frame

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, VAGMD_frame):
        m = VAGMD_frame
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, VAGMD_frame):
        m = VAGMD_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_frame):
        m = VAGMD_frame
        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(69.4664, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(34.5742, abs=1e-3)
        assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.3404, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(7.0701, abs=1e-3)

    @pytest.mark.component
    def test_solution_2(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with high brine salinity (> 173.5 g/L) and closed cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C
        feed_temp = 25  # deg C
        feed_salinity = 100  # g/L
        module_type = "AS7C1.5L"
        cooling_system_type = "closed"
        cooling_inlet_temp = (
            25  # deg C, not required when cooling system type is "closed"
        )
        high_brine_salinity = True

        # If the brine salinity is high (>175.3 g/L) for module AS7C1.5L,
        # the operational parameters need to be fixed at at certain values,
        # and cooling circuit is closed to maintain condenser inlet temperature constant
        feed_flow_rate = 1100  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C

        m.fs.vagmd = VAGMDSurrogateBase(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # Run helper function to determine salinity mode
        (
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            cooling_system_type,
        ) = m.fs.vagmd._determine_salinity_mode(
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            module_type,
            high_brine_salinity,
            cooling_system_type,
        )

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        # Specify evaporator inlet temperature
        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        # Calculate scaling factor and initialize the model
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)
        # Solve the model
        results = solver.solve(m)

        # Test solution
        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(65.6459, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(37.4374, abs=1e-3)
        assert m.fs.vagmd.permeate_flux.value == pytest.approx(7.6420, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(17.2473, abs=1e-3)

    @pytest.mark.component
    def test_solution_3(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with low brine salinity (< 173.5 g/L) and open cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C
        feed_temp = 25  # deg C
        feed_salinity = 50  # g/L
        module_type = "AS7C1.5L"
        cooling_system_type = "open"
        cooling_inlet_temp = (
            25  # deg C, not required when cooling system type is "closed"
        )
        high_brine_salinity = False

        m.fs.vagmd = VAGMDSurrogateBase(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # Run helper function to determine salinity mode
        (
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            cooling_system_type,
        ) = m.fs.vagmd._determine_salinity_mode(
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            module_type,
            high_brine_salinity,
            cooling_system_type,
        )

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        # Specify evaporator inlet temperature
        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        # Calculate scaling factor and initialize the model
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)
        # Solve the model
        results = solver.solve(m)

        # Test solution
        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(69.3523, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(34.5056, abs=1e-3)
        assert m.fs.vagmd.permeate_flux.value == pytest.approx(5.2033, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(7.1045, abs=1e-3)

    @pytest.mark.component
    def test_solution_4(self):
        # Create model, flowsheet for configuration of module AS26C7.2L,
        # with closed cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C
        feed_temp = 25  # deg C
        feed_salinity = 35  # g/L
        initial_batch_volume = 50  # L
        recovery_ratio = 0.5  # -
        module_type = "AS26C7.2L"
        cooling_system_type = "closed"
        cooling_inlet_temp = (
            25  # deg C, not required when cooling system type is "closed"
        )
        high_brine_salinity = False

        m.fs.vagmd = VAGMDSurrogateBase(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # No need to run _determine_salinity_mode for module AS26C7.2L

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        # Specify evaporator inlet temperature
        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        # Calculate scaling factor and initialize the model
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)
        # Solve the model
        results = solver.solve(m)

        # Test solution
        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(75.7306, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(27.8301, abs=1e-3)
        assert m.fs.vagmd.permeate_flux.value == pytest.approx(1.6197, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(2.8616, abs=1e-3)

    @pytest.mark.component
    def test_solution_5(self):
        # Create model, flowsheet for configuration of module AS26C7.2L,
        # with open cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C
        feed_temp = 25  # deg C
        feed_salinity = 50  # g/L
        initial_batch_volume = 50  # L
        recovery_ratio = 0.5  # -
        module_type = "AS26C7.2L"
        cooling_system_type = "open"
        cooling_inlet_temp = (
            25  # deg C, not required when cooling system type is "closed"
        )
        high_brine_salinity = False

        m.fs.vagmd = VAGMDSurrogateBase(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # No need to run _determine_salinity_mode for module AS26C7.2L

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        # Specify evaporator inlet temperature
        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        # Calculate scaling factor and initialize the model
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)
        # Solve the model
        results = solver.solve(m)

        # Test solution
        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(75.5880, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[
            0
        ].temperature.value - 273.15 == pytest.approx(28.0441, abs=1e-3)
        assert m.fs.vagmd.permeate_flux.value == pytest.approx(1.4949, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(2.9398, abs=1e-3)
