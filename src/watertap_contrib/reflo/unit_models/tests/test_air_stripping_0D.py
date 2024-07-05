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
    Set,
    Var,
    Param,
    Expression,
    value,
    assert_optimal_termination,
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
from watertap_contrib.reflo.costing import REFLOCosting
from watertap_contrib.reflo.unit_models.air_stripping_0D import (
    AirStripping0D,
    PackingMaterial,
)

# Get default solver for testing
solver = get_solver()


class TestAirStripping0D:
    @pytest.fixture(scope="class")
    def ax_frame1(self):
        target = "TCA"
        props = {
            "solute_list": [target],
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
        prop_in = ax.process_flow.properties_in[0]
        prop_out = ax.process_flow.properties_out[0]

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
        prop_in.temperature["Liq"].fix(288)
        prop_in.temperature["Vap"].fix(288)
        prop_in.pressure.fix(101325)

        ax.pressure_drop_gradient.fix(75)
        ax.packing_surf_tension.fix(0.033)
        ax.packing_diam_nominal.fix(0.0889)
        ax.packing_surface_area_total.fix(242)
        ax.packing_factor.fix(33)
        ax.surf_tension_water.fix(0.0735)
        ax.target_reduction_frac[target].set_value(0.97)

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

        assert isinstance(ax.process_flow, ControlVolume0DBlock)
        assert hasattr(ax.process_flow, "mass_transfer_term")
        assert hasattr(ax.process_flow, "deltaP")

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
        assert number_variables(m) == 86
        assert number_total_constraints(m) == 64
        assert number_unused_variables(m) == 1

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
        prop_in = ax.process_flow.properties_in[0]
        prop_out = ax.process_flow.properties_out[0]

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
            "blower_power": 2.239159,
            "pump_power": 19.33556,
            "packing_surface_area_total": 242.0,
            "packing_surface_area_wetted": 147.561,
            "packing_diam_nominal": 0.0889,
            "packing_factor": 33.0,
            "packing_surf_tension": 0.033,
            "surf_tension_water": 0.0735,
            "stripping_factor": {"TCA": 3.394407},
            "air_water_ratio_min": 2.000348,
            "packing_height": 8.841,
            "mass_loading_rate": {"Liq": 39.15, "Vap": 0.334627822},
            "height_transfer_unit": {"TCA": 1.967472},
            "number_transfer_unit": {"TCA": 4.49394},
            "pressure_drop_gradient": 75.0,
            "overall_mass_transfer_coeff": {"TCA": 0.019915699},
            "oto_E": 0.611802769,
            "oto_F": 1.875061,
            "oto_a0": -2.279915,
            "oto_a1": -0.729372463,
            "oto_a2": -0.228704898,
            "oto_M": 0.00154258,
            "oto_mass_transfer_coeff": {
                ("Liq", "TCA"): 0.00036455,
                ("Vap", "TCA"): 0.000842984,
            },
            "air_water_ratio_op": 6.999,
            "packing_efficiency_number": 21.513,
            "tower_area": 4.0323,
            "tower_diam": 2.265851,
            "tower_height": 10.61,
            "tower_volume": 42.782,
            "packing_volume": 35.652,
            "target_remaining_frac": {"TCA": 0.03},
            "pressure_drop": 795.753,
            "pressure_drop_tower": 20.688,
            "N_Sc": {("Liq", "TCA"): 1622.798, ("Vap", "TCA"): 1.8083},
            "N_Re": 140.676,
            "N_Fr": 0.037888131,
            "N_We": 0.086245507,
            "oto_kl_term": 88595.604,
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
        prop_out = ax.process_flow.properties_out[0]

        m.fs.costing = REFLOCosting()
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
            "capital_cost": 391313.521,
            "tower_cost": 20201.182,
            "port_cost": 643.059,
            "piping_liq_cost": 2189.254,
            "piping_air_cost": 2298.717,
            "tray_ring_cost": 1185.632,
            "distributor_cost": 3584.056,
            "plate_cost": 1745.286,
            "tower_internals_cost": 5329.343,
            "packing_cost": 252197.627,
            "mist_eliminator_cost": 3899.822,
            "pump_cost": 42038.547,
            "blower_cost": 61305.385,
            "electricity_flow": 21.562,
        }

        for v, r in ax_costing_results.items():
            axc = getattr(ax.costing, v)
            assert value(axc) == pytest.approx(r, rel=1e-3)

        m_costing_results = {
            "aggregate_capital_cost": 391313.521,
            "aggregate_fixed_operating_cost": 0.0,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 21.56,
            "aggregate_flow_costs": {"electricity": 15532.43},
            "total_capital_cost": 391313.521,
            "maintenance_labor_chemical_operating_cost": 11739.41,
            "total_operating_cost": 27280.755,
            "LCOW": 0.0142577,
            "SEC": 0.03793025,
        }

        for v, r in m_costing_results.items():
            mc = getattr(m.fs.costing, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(mc[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(mc) == pytest.approx(r, rel=1e-3)

    @pytest.fixture(scope="class")
    def ax_frame2(self):
        target = "DCP"
        props = {
            "solute_list": [target],
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
        prop_in = ax.process_flow.properties_in[0]
        prop_out = ax.process_flow.properties_out[0]

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

        assert isinstance(ax.process_flow, ControlVolume0DBlock)
        assert hasattr(ax.process_flow, "mass_transfer_term")
        assert hasattr(ax.process_flow, "deltaP")

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
        assert number_variables(m) == 86
        assert number_total_constraints(m) == 64
        assert number_unused_variables(m) == 1

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
        prop_in = ax.process_flow.properties_in[0]
        prop_out = ax.process_flow.properties_out[0]

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
            "blower_power": 11.475,
            "pump_power": 16.365,
            "packing_surface_area_total": 125.0,
            "packing_surface_area_wetted": 56.766,
            "packing_diam_nominal": 0.0889,
            "packing_factor": 39.0,
            "packing_surf_tension": 0.033,
            "surf_tension_water": 0.0742,
            "stripping_factor": {"DCP": 4.503699},
            "air_water_ratio_min": 11.99,
            "packing_height": 11.824,
            "mass_loading_rate": {"Liq": 7.7, "Vap": 0.5763451},
            "height_transfer_unit": {"DCP": 4.423251},
            "number_transfer_unit": {"DCP": 2.673204},
            "pressure_drop_gradient": 50.0,
            "overall_mass_transfer_coeff": {"DCP": 0.0017415},
            "oto_E": -0.325878163,
            "oto_F": 1.69897,
            "oto_a0": -2.457618,
            "oto_a1": -0.633128697,
            "oto_a2": -0.186834171,
            "oto_M": 0.00535628,
            "oto_mass_transfer_coeff": {
                ("Liq", "DCP"): 0.000162215,
                ("Vap", "DCP"): 0.000800008,
            },
            "air_water_ratio_op": 59.999,
            "packing_efficiency_number": 11.112,
            "tower_area": 12.981,
            "tower_diam": 4.06558,
            "tower_height": 14.189,
            "tower_volume": 184.2,
            "packing_volume": 153.5,
            "target_remaining_frac": {"DCP": 0.099999999},
            "pressure_drop": 709.455,
            "pressure_drop_tower": 58.744,
            "N_Sc": {("Liq", "DCP"): 1813.198, ("Vap", "DCP"): 1.673158},
            "N_Re": 47.135,
            "N_Fr": 0.000756345,
            "N_We": 0.006395675,
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
        prop_out = ax.process_flow.properties_out[0]

        m.fs.costing = REFLOCosting()
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
            "capital_cost": 1414661.385,
            "tower_cost": 37232.248,
            "port_cost": 643.059,
            "piping_liq_cost": 2189.254,
            "piping_air_cost": 2298.717,
            "tray_ring_cost": 2402.264,
            "distributor_cost": 10295.84,
            "plate_cost": 5255.23,
            "tower_internals_cost": 15551.071,
            "packing_cost": 1085202.793,
            "mist_eliminator_cost": 10036.566,
            "pump_cost": 38366.194,
            "blower_cost": 220739.214,
            "electricity_flow": 27.84,
        }

        for v, r in ax_costing_results.items():
            axc = getattr(ax.costing, v)
            assert value(axc) == pytest.approx(r, rel=1e-3)

        m_costing_results = {
            "aggregate_capital_cost": 1414661.385,
            "aggregate_fixed_operating_cost": 0.0,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 27.84,
            "aggregate_flow_costs": {"electricity": 20054.989},
            "total_capital_cost": 1414661.385,
            "maintenance_labor_chemical_operating_cost": 42439.841,
            "total_operating_cost": 62494.830,
            "LCOW": 0.0699909,
            "SEC": 0.077335064,
        }

        for v, r in m_costing_results.items():
            mc = getattr(m.fs.costing, v)
            if isinstance(r, dict):
                for i, s in r.items():
                    assert value(mc[i]) == pytest.approx(s, rel=1e-3)
            else:
                assert value(mc) == pytest.approx(r, rel=1e-3)
