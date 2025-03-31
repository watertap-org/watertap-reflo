#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
    Set,
    Var,
    Param,
    Expression,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
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
    set_scaling_factor,
)
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.property_models import AirWaterEq
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.air_stripping_0D import (
    AirStripping0D,
    PackingMaterial,
)

# Get default solver for testing
solver = get_solver()


def build_ax_1():

    target = "TCA"
    props = {
        "volatile_solute_list": [target],
        "mw_data": {target: 0.1334},
        "dynamic_viscosity_data": {"Liq": 0.00115, "Vap": 1.75e-5},
        "henry_constant_data": {target: 0.725},  # salinity adjusted
        "standard_enthalpy_change_data": {target: 28.7e3},
        "temperature_boiling_data": {target: 347},
        "molar_volume_data": {target: 9.81e-5},
        "critical_molar_volume_data": {target: 2.94e-4},
        "density_data": {"Liq": 999.15, "Vap": 1.22},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    ax_config = {"property_package": m.fs.properties, "target": target}

    m.fs.ax = ax = AirStripping0D(**ax_config)
    prop_in = ax.control_volume.properties_in[0]
    prop_out = ax.control_volume.properties_out[0]

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.0063345, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 38358.266, index=("Liq", target)
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", target)
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.741114, index=("Vap", "Air")
    )
    set_scaling_factor(prop_out.flow_mass_phase_comp["Vap", target], 1e6)

    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(157.8657)
    prop_in.flow_mass_phase_comp["Liq", target].fix(2.61e-5)
    prop_in.flow_mass_phase_comp["Vap", "Air"].fix(1.34932)
    prop_in.flow_mass_phase_comp["Vap", target].fix(0)  # assume pure air into unit
    prop_in.flow_mass_phase_comp["Vap", "H2O"].fix(0.10)
    prop_in.temperature["Liq"].fix(283)
    prop_in.temperature["Vap"].fix(283)
    prop_in.pressure.fix(101325)

    ax.pressure_drop_gradient.fix(75)
    ax.packing_surf_tension.fix(0.033)
    ax.packing_diam_nominal.fix(0.0889)
    ax.packing_surface_area_total.fix(242)
    ax.packing_factor.fix(33)
    ax.surf_tension_water.fix(0.0735)
    ax.target_reduction_frac[target].set_value(0.97)

    return m


def build_ax_2():

    target = "DCP"
    props = {
        "volatile_solute_list": [target],
        "mw_data": {target: 0.11298},
        "dynamic_viscosity_data": {"Liq": 1.307e-3, "Vap": 0.0000179},
        "henry_constant_data": {target: 0.146},  # salinity adjusted
        "standard_enthalpy_change_data": {target: 31.1e3},
        "temperature_boiling_data": {target: 369.1},
        "molar_volume_data": {target: 1.10071e-4},
        "critical_molar_volume_data": {target: 2.26e-4},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    ax_config = {"property_package": m.fs.properties, "target": target}

    m.fs.ax = ax = AirStripping0D(**ax_config)
    prop_in = ax.control_volume.properties_in[0]
    prop_out = ax.control_volume.properties_out[0]

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.01, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e4, index=("Liq", target)
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", target)
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.133654, index=("Vap", "Air")
    )
    set_scaling_factor(prop_out.flow_mass_phase_comp["Vap", target], 1e6)

    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    prop_in.flow_mass_phase_comp["Liq", target].fix(1e-4)
    prop_in.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    prop_in.flow_mass_phase_comp["Vap", target].fix(0)  # assume pure air into unit
    prop_in.flow_mass_phase_comp["Vap", "H2O"].fix(0.10)
    prop_in.temperature["Liq"].fix(283)
    prop_in.temperature["Vap"].fix(283)
    prop_in.pressure.fix(101325)

    ax.pressure_drop_gradient.fix(50)
    ax.packing_surf_tension.fix(0.033)
    ax.packing_diam_nominal.fix(0.0889)
    ax.packing_surface_area_total.fix(125)
    ax.packing_factor.fix(39)
    ax.surf_tension_water.fix(0.0742)
    ax.target_reduction_frac[target].set_value(0.9)

    return m


class TestAirStripping0D:
    @pytest.fixture(scope="class")
    def ax_frame1(self):

        m = build_ax_1()

        return m

    @pytest.mark.unit
    def test_config1(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax

        assert len(ax.config) == 9

        assert not ax.config.dynamic
        assert not ax.config.has_holdup
        assert ax.config.property_package is m.fs.properties
        assert ax.config.material_balance_type == MaterialBalanceType.componentPhase
        assert ax.config.energy_balance_type is EnergyBalanceType.none
        assert ax.config.momentum_balance_type is MomentumBalanceType.pressureTotal
        assert ax.config.target == "TCA"
        assert ax.config.packing_material == PackingMaterial.PVC

    @pytest.mark.unit
    def test_build1(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax

        # test ports
        port_lst = ["inlet", "outlet"]
        for port_str in port_lst:
            port = getattr(ax, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert isinstance(ax.control_volume, ControlVolume0DBlock)
        assert hasattr(ax.control_volume, "mass_transfer_term")
        assert hasattr(ax.control_volume, "deltaP")

        assert isinstance(ax.target_set, Set)
        assert len(ax.target_set) == 1
        assert ax.target_set.at(1) == ax.config.target

        assert isinstance(ax.liq_target_set, Set)
        assert len(ax.liq_target_set) == 1
        assert ax.liq_target_set.at(1) == ("Liq", ax.config.target)

        assert isinstance(ax.phase_target_set, Set)
        assert len(ax.phase_target_set) == 2
        assert ax.phase_target_set.at(1) == ("Liq", ax.config.target)
        assert ax.phase_target_set.at(2) == ("Vap", ax.config.target)

        assert hasattr(ax, "build_oto")

        # test statistics
        assert number_variables(m) == 91
        assert number_total_constraints(m) == 58
        assert number_unused_variables(m) == 8

        ax_params = [
            "air_water_ratio_param",
            "pressure_drop_tower_param",
            "tower_height_factor",
            "tower_port_diameter",
            "tower_pipe_diameter",
            "target_reduction_frac",
            "overall_mass_transfer_coeff_sf",
            "blower_efficiency",
            "pump_efficiency",
            "power_blower_denom_coeff",
            "power_blower_exponent",
            "oto_a0_param1",
            "oto_a0_param2",
            "oto_a0_param3",
            "oto_a0_param4",
            "oto_a1_param1",
            "oto_a1_param2",
            "oto_a1_param3",
            "oto_a1_param4",
            "oto_a2_param1",
            "oto_a2_param2",
            "oto_a2_param3",
            "oto_a2_param4",
            "oto_aw_param",
            "oto_aw_exp1",
            "oto_aw_exp2",
            "oto_aw_exp3",
            "oto_aw_exp4",
            "oto_liq_mass_xfr_param",
            "oto_liq_mass_xfr_exp1",
            "oto_liq_mass_xfr_exp2",
            "oto_liq_mass_xfr_exp3",
            "oto_liq_mass_xfr_exp4",
            "oto_gas_mass_xfr_param",
            "oto_gas_mass_xfr_exp1",
            "oto_gas_mass_xfr_exp2",
            "oto_gas_mass_xfr_exp3",
        ]

        for pname in ax_params:
            assert hasattr(ax, pname)
            assert isinstance(getattr(ax, pname), Param)

        ax_vars = [
            "packing_surface_area_total",
            "packing_surface_area_wetted",
            "packing_diam_nominal",
            "packing_factor",
            "packing_surf_tension",
            "surf_tension_water",
            "stripping_factor",
            "air_water_ratio_min",
            "packing_height",
            "mass_loading_rate",
            "height_transfer_unit",
            "number_transfer_unit",
            "pressure_drop_gradient",
            "overall_mass_transfer_coeff",
            "oto_E",
            "oto_F",
            "oto_a0",
            "oto_a1",
            "oto_a2",
            "oto_M",
            "oto_mass_transfer_coeff",
        ]

        for vname in ax_vars:
            assert hasattr(ax, vname)
            assert isinstance(getattr(ax, vname), Var)

        ax_expr = [
            "air_water_ratio_op",
            "packing_efficiency_number",
            "tower_area",
            "tower_diam",
            "tower_height",
            "tower_volume",
            "packing_volume",
            "target_remaining_frac",
            "pressure_drop",
            "pressure_drop_tower",
        ]

        for ename in ax_expr:
            assert hasattr(ax, ename)
            assert isinstance(getattr(ax, ename), Expression)

    @pytest.mark.unit
    def test_dof1(self, ax_frame1):
        m = ax_frame1
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling1(self, ax_frame1):
        m = ax_frame1

        calculate_scaling_factors(m)
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize1(self, ax_frame1):
        m = ax_frame1
        initialization_tester(m, unit=m.fs.ax, outlvl=idaeslog.INFO_LOW)

    @pytest.mark.component
    def test_solve1(self, ax_frame1):
        m = ax_frame1
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_mass_balance1(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax
        prop_in = ax.control_volume.properties_in[0]
        prop_out = ax.control_volume.properties_out[0]

        assert pytest.approx(
            sum(
                value(prop_in.flow_mass_phase_comp[p, ax.config.target])
                for p in m.fs.properties.phase_list
            ),
            rel=1e-5,
        ) == sum(
            value(prop_out.flow_mass_phase_comp[p, ax.config.target])
            for p in m.fs.properties.phase_list
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", "H2O"]) == value(
            prop_out.flow_mass_phase_comp["Liq", "H2O"]
        )

        assert value(prop_in.flow_mass_phase_comp["Vap", "Air"]) == value(
            prop_out.flow_mass_phase_comp["Vap", "Air"]
        )

    @pytest.mark.component
    def test_solution1(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax

        ax_results = {
            "blower_power": 2.7578,
            "pump_power": 22.62,
            "packing_surface_area_total": 242.0,
            "packing_surface_area_wetted": 146.29,
            "packing_diam_nominal": 0.0889,
            "packing_factor": 33.0,
            "packing_surf_tension": 0.033,
            "surf_tension_water": 0.0735,
            "stripping_factor": {"TCA": 2.9501},
            "air_water_ratio_min": 2.4721,
            "packing_height": 10.34,
            "mass_loading_rate": {"Liq": 37.77, "Vap": 0.346832},
            "height_transfer_unit": {"TCA": 2.2014},
            "number_transfer_unit": {"TCA": 4.7015},
            "pressure_drop_gradient": 75.0,
            "overall_mass_transfer_coeff": {"TCA": 0.017175},
            "air_water_ratio_op": 7.5187,
            "packing_efficiency_number": 21.51,
            "tower_area": 4.1787,
            "tower_diam": 2.3066,
            "tower_height": 12.41,
            "tower_volume": 51.89,
            "packing_volume": 43.24,
            "target_remaining_frac": {"TCA": 0.03},
            "pressure_drop": 931.49,
            "pressure_drop_tower": 22.22,
            "N_Sc": {("Liq", "TCA"): 1622.79, ("Vap", "TCA"): 1.8693},
            "N_Re": 135.74,
            "N_Fr": 0.035279,
            "N_We": 0.080306,
            "oto_E": 0.580753,
            "oto_F": 1.875,
            "oto_a0": -2.2799,
            "oto_a1": -0.729372,
            "oto_a2": -0.228704,
            "oto_M": 0.001657,
            "oto_mass_transfer_coeff": {
                ("Liq", "TCA"): 0.000358,
                ("Vap", "TCA"): 0.000804,
            },
            "oto_kl_term": 88595.6,
        }

        for v, r in ax_results.items():
            axv = getattr(ax, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(axv[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(axv) == pytest.approx(r, rel=1e-3)

    @pytest.mark.component
    def test_costing1(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax
        prop_out = ax.control_volume.properties_out[0]

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021

        ax.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_out.flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            prop_out.flow_vol_phase["Liq"], name="SEC"
        )
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        # return m

        ax_costing_results = {
            "capital_cost": 455940.87,
            "tower_cost": 23918.06,
            "port_cost": 643.05,
            "piping_liq_cost": 2189.25,
            "piping_air_cost": 2298.71,
            "tray_ring_cost": 1209.42,
            "distributor_cost": 3688.19,
            "plate_cost": 1803.67,
            "tower_internals_cost": 5491.87,
            "packing_cost": 305762.33,
            "mist_eliminator_cost": 4008.36,
            "pump_cost": 45842.5,
            "blower_cost": 64577.26,
            "cost_factor": 1.0,
            "direct_capital_cost": 455940.87,
            "electricity_flow": 25.37,
        }

        for v, r in ax_costing_results.items():
            axc = getattr(ax.costing, v)
            assert value(axc) == pytest.approx(r, rel=1e-3)

        sys_costing_results = {
            "aggregate_capital_cost": 455940.87,
            "aggregate_flow_electricity": 25.37,
            "aggregate_flow_costs": {"electricity": 18281.51},
            "aggregate_direct_capital_cost": 455940.87,
            "total_capital_cost": 455940.87,
            "total_operating_cost": 31959.73,
            "maintenance_labor_chemical_operating_cost": 13678.22,
            "total_fixed_operating_cost": 13678.22,
            "total_variable_operating_cost": 18281.51,
            "total_annualized_cost": 83005.02,
            "LCOW": 0.016647,
            "SEC": 0.044617,
        }

        for v, r in sys_costing_results.items():
            mc = getattr(m.fs.costing, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(mc[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(mc) == pytest.approx(r, rel=1e-3)

    @pytest.fixture(scope="class")
    def ax_frame2(self):

        m = build_ax_2()

        return m

    @pytest.mark.unit
    def test_config2(self, ax_frame2):
        m = ax_frame2
        ax = m.fs.ax

        assert len(ax.config) == 9

        assert not ax.config.dynamic
        assert not ax.config.has_holdup
        assert ax.config.property_package is m.fs.properties
        assert ax.config.material_balance_type == MaterialBalanceType.componentPhase
        assert ax.config.energy_balance_type is EnergyBalanceType.none
        assert ax.config.momentum_balance_type is MomentumBalanceType.pressureTotal
        assert ax.config.target == "DCP"
        assert ax.config.packing_material == PackingMaterial.PVC

    @pytest.mark.unit
    def test_build2(self, ax_frame2):
        m = ax_frame2
        ax = m.fs.ax

        # test ports
        port_lst = ["inlet", "outlet"]
        for port_str in port_lst:
            port = getattr(ax, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert isinstance(ax.control_volume, ControlVolume0DBlock)
        assert hasattr(ax.control_volume, "mass_transfer_term")
        assert hasattr(ax.control_volume, "deltaP")

        assert isinstance(ax.target_set, Set)
        assert len(ax.target_set) == 1
        assert ax.target_set.at(1) == ax.config.target

        assert isinstance(ax.liq_target_set, Set)
        assert len(ax.liq_target_set) == 1
        assert ax.liq_target_set.at(1) == ("Liq", ax.config.target)

        assert isinstance(ax.phase_target_set, Set)
        assert len(ax.phase_target_set) == 2
        assert ax.phase_target_set.at(1) == ("Liq", ax.config.target)
        assert ax.phase_target_set.at(2) == ("Vap", ax.config.target)

        assert hasattr(ax, "build_oto")

        # test statistics
        assert number_variables(m) == 91
        assert number_total_constraints(m) == 58
        assert number_unused_variables(m) == 8

        ax_params = [
            "air_water_ratio_param",
            "pressure_drop_tower_param",
            "tower_height_factor",
            "tower_port_diameter",
            "tower_pipe_diameter",
            "target_reduction_frac",
            "overall_mass_transfer_coeff_sf",
            "blower_efficiency",
            "pump_efficiency",
            "power_blower_denom_coeff",
            "power_blower_exponent",
            "oto_a0_param1",
            "oto_a0_param2",
            "oto_a0_param3",
            "oto_a0_param4",
            "oto_a1_param1",
            "oto_a1_param2",
            "oto_a1_param3",
            "oto_a1_param4",
            "oto_a2_param1",
            "oto_a2_param2",
            "oto_a2_param3",
            "oto_a2_param4",
            "oto_aw_param",
            "oto_aw_exp1",
            "oto_aw_exp2",
            "oto_aw_exp3",
            "oto_aw_exp4",
            "oto_liq_mass_xfr_param",
            "oto_liq_mass_xfr_exp1",
            "oto_liq_mass_xfr_exp2",
            "oto_liq_mass_xfr_exp3",
            "oto_liq_mass_xfr_exp4",
            "oto_gas_mass_xfr_param",
            "oto_gas_mass_xfr_exp1",
            "oto_gas_mass_xfr_exp2",
            "oto_gas_mass_xfr_exp3",
        ]

        for pname in ax_params:
            assert hasattr(ax, pname)
            assert isinstance(getattr(ax, pname), Param)

        ax_vars = [
            "blower_power",
            "pump_power",
            "packing_surface_area_total",
            "packing_surface_area_wetted",
            "packing_diam_nominal",
            "packing_factor",
            "packing_surf_tension",
            "surf_tension_water",
            "stripping_factor",
            "air_water_ratio_min",
            "packing_height",
            "mass_loading_rate",
            "height_transfer_unit",
            "number_transfer_unit",
            "pressure_drop_gradient",
            "overall_mass_transfer_coeff",
            "oto_E",
            "oto_F",
            "oto_a0",
            "oto_a1",
            "oto_a2",
            "oto_M",
            "oto_mass_transfer_coeff",
        ]

        for vname in ax_vars:
            assert hasattr(ax, vname)
            assert isinstance(getattr(ax, vname), Var)

        ax_expr = [
            "air_water_ratio_op",
            "packing_efficiency_number",
            "tower_area",
            "tower_diam",
            "tower_height",
            "tower_volume",
            "packing_volume",
            "target_remaining_frac",
            "pressure_drop",
            "pressure_drop_tower",
        ]

        for ename in ax_expr:
            assert hasattr(ax, ename)
            assert isinstance(getattr(ax, ename), Expression)

    @pytest.mark.unit
    def test_dof2(self, ax_frame2):
        m = ax_frame2
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling2(self, ax_frame2):
        m = ax_frame2

        calculate_scaling_factors(m)
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize2(self, ax_frame2):
        m = ax_frame2
        initialization_tester(m, unit=m.fs.ax, outlvl=idaeslog.INFO_LOW)

    @pytest.mark.component
    def test_solve2(self, ax_frame2):
        m = ax_frame2
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_mass_balance2(self, ax_frame2):
        m = ax_frame2
        ax = m.fs.ax
        prop_in = ax.control_volume.properties_in[0]
        prop_out = ax.control_volume.properties_out[0]

        assert pytest.approx(
            sum(
                value(prop_in.flow_mass_phase_comp[p, ax.config.target])
                for p in m.fs.properties.phase_list
            ),
            rel=1e-5,
        ) == sum(
            value(prop_out.flow_mass_phase_comp[p, ax.config.target])
            for p in m.fs.properties.phase_list
        )

        assert value(prop_in.flow_mass_phase_comp["Liq", "H2O"]) == value(
            prop_out.flow_mass_phase_comp["Liq", "H2O"]
        )

        assert value(prop_in.flow_mass_phase_comp["Vap", "Air"]) == value(
            prop_out.flow_mass_phase_comp["Vap", "Air"]
        )

    @pytest.mark.component
    def test_solution2(self, ax_frame2):
        m = ax_frame2
        ax = m.fs.ax

        ax_results = {
            "blower_power": 11.603,
            "pump_power": 16.318,
            "packing_surface_area_total": 125.0,
            "packing_surface_area_wetted": 56.6,
            "packing_diam_nominal": 0.0889,
            "packing_factor": 39.0,
            "packing_surf_tension": 0.033,
            "surf_tension_water": 0.0742,
            "stripping_factor": {"DCP": 4.56389},
            "air_water_ratio_min": 11.99,
            "packing_height": 11.79,
            "mass_loading_rate": {"Liq": 7.625, "Vap": 0.5783},
            "height_transfer_unit": {"DCP": 4.420},
            "number_transfer_unit": {"DCP": 2.6673},
            "pressure_drop_gradient": 50.0,
            "overall_mass_transfer_coeff": {"DCP": 0.0017256},
            "air_water_ratio_op": 60.801,
            "packing_efficiency_number": 11.112,
            "tower_area": 13.11,
            "tower_diam": 4.085733,
            "tower_height": 14.148,
            "tower_volume": 185.497,
            "packing_volume": 154.581,
            "target_remaining_frac": {"DCP": 0.1},
            "pressure_drop": 707.421,
            "pressure_drop_tower": 59.143,
            "N_Sc": {("Liq", "DCP"): 1813.198, ("Vap", "DCP"): 1.673158},
            "N_Re": 46.671,
            "N_Fr": 0.000741,
            "N_We": 0.006270,
            "oto_E": -0.331644,
            "oto_F": 1.69897,
            "oto_a0": -2.457618,
            "oto_a1": -0.633128,
            "oto_a2": -0.186834,
            "oto_M": 0.005392,
            "oto_mass_transfer_coeff": {
                ("Liq", "DCP"): 0.0001614,
                ("Vap", "DCP"): 0.000794489,
            },
            "oto_kl_term": 77996.197,
        }
        for v, r in ax_results.items():
            axv = getattr(ax, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(axv[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(axv) == pytest.approx(r, rel=1e-3)

    @pytest.mark.component
    def test_costing2(self, ax_frame2):
        m = ax_frame2
        ax = m.fs.ax
        prop_out = ax.control_volume.properties_out[0]

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021

        ax.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_out.flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            prop_out.flow_vol_phase["Liq"], name="SEC"
        )
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        ax_costing_results = {
            "capital_cost": 1424757.1,
            "tower_cost": 37199.7,
            "port_cost": 643.0,
            "piping_liq_cost": 2189.2,
            "piping_air_cost": 2298.7,
            "tray_ring_cost": 2417.8,
            "distributor_cost": 10395.5,
            "plate_cost": 5305.3,
            "tower_internals_cost": 15700.8,
            "packing_cost": 1092845.9,
            "mist_eliminator_cost": 10120.8,
            "pump_cost": 38305.6,
            "blower_cost": 223035.2,
            "cost_factor": 1.0,
            "direct_capital_cost": 1424757.1,
            "electricity_flow": 27.9,
        }

        for v, r in ax_costing_results.items():
            axc = getattr(ax.costing, v)
            assert value(axc) == pytest.approx(r, rel=1e-3)

        sys_costing_results = {
            "aggregate_capital_cost": 1424757.1,
            "aggregate_flow_electricity": 27.9,
            "aggregate_flow_costs": {"electricity": 20113.8},
            "aggregate_direct_capital_cost": 1424757.1,
            "total_capital_cost": 1424757.1,
            "total_operating_cost": 62856.6,
            "maintenance_labor_chemical_operating_cost": 42742.7,
            "total_fixed_operating_cost": 42742.7,
            "total_variable_operating_cost": 20113.8,
            "total_annualized_cost": 222366.6,
            "LCOW": 0.070463,
            "SEC": 0.07756,
        }

        for v, r in sys_costing_results.items():
            mc = getattr(m.fs.costing, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(mc[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(mc) == pytest.approx(r, rel=1e-3)
