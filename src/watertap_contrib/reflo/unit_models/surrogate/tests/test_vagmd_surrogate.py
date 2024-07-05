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

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
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

from watertap_contrib.reflo.unit_models.surrogate import VAGMDSurrogate
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


class TestVAGMD_unit_model:
    @pytest.fixture(scope="class")
    def VAGMD_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        system_capacity = 2000  # m3/day
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

        m.fs.vagmd = VAGMDSurrogate(
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

        # Specify system capacity
        m.fs.vagmd.system_capacity.fix(system_capacity)

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
        else:  # cooling_system_type == "open"
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
        assert number_variables(m) == 246
        assert number_total_constraints(m) == 74
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
        assert value(m.fs.vagmd.num_modules) == pytest.approx(2167.25, abs=1e-3)
        assert m.fs.vagmd.recovery_ratio.value == pytest.approx(0.064085, abs=1e-5)

    @pytest.mark.component
    def test_costing(self, VAGMD_frame):
        m = VAGMD_frame
        vagmd = m.fs.vagmd

        m.fs.costing = TreatmentCosting()
        m.fs.costing.base_currency = pyunits.USD_2020
        vagmd.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        # Fix some global costing params for better comparison to Pyomo model
        m.fs.costing.total_investment_factor.fix(1)
        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.capital_recovery_factor.fix(0.08764)
        m.fs.costing.wacc.unfix()

        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(vagmd.system_capacity)
        m.fs.costing.add_LCOW(vagmd.system_capacity)

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(2167.25, rel=1e-3) == value(vagmd.num_modules)
        assert pytest.approx(2763840.418, rel=1e-3) == value(vagmd.costing.module_cost)
        assert pytest.approx(602433.068, rel=1e-3) == value(
            vagmd.costing.other_capital_cost
        )
        assert pytest.approx(3366273.443, rel=1e-3) == value(vagmd.costing.capital_cost)
        assert pytest.approx(294352.922, rel=1e-3) == value(
            vagmd.costing.fixed_operating_cost
        )
        assert pytest.approx(2.648, rel=1e-3) == value(m.fs.costing.LCOW)
