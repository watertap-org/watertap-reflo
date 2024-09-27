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

from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import (
    ConcreteModel,
    Var,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock

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
from idaes.core import UnitModelCostingBlock

from watertap.core.solvers import get_solver
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
    NaClStateBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock,
    WaterStateBlock,
)

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.multi_effect_crystallizer import (
    MultiEffectCrystallizer,
)
from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect

__author__ = "Kurban Sitterley"

solver = get_solver()

feed_pressure = 101325
feed_temperature = 273.15 + 20
eps = 1e-6


def build_mec4():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.unit = mec = MultiEffectCrystallizer(
        property_package=m.fs.properties, property_package_vapor=m.fs.vapor_properties
    )

    operating_pressures = [0.45, 0.25, 0.208, 0.095]

    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.15
    crystallizer_yield = 0.5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    atm_pressure = 101325 * pyunits.Pa
    saturated_steam_pressure_gage = 3 * pyunits.bar
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage, to_units=pyunits.Pa
    )

    for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        eff.effect.properties_in[0].pressure.fix(feed_pressure)
        eff.effect.properties_in[0].temperature.fix(feed_temperature)

        eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(eps)
        eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(eps)
        eff.effect.properties_in[0].conc_mass_phase_comp[...]
        eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()
        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.set_value(0.1)

    first_effect = m.fs.unit.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(0.1)
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    return m


def build_mec3():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.unit = mec = MultiEffectCrystallizer(
        property_package=m.fs.properties,
        property_package_vapor=m.fs.vapor_properties,
        number_effects=3,
    )

    feed_flow_mass = 10
    feed_mass_frac_NaCl = 0.25
    crystallizer_yield = 0.55
    operating_pressures = [0.45, 0.25, 0.208]
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    atm_pressure = 101325 * pyunits.Pa
    saturated_steam_pressure_gage = 3.5 * pyunits.bar
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage, to_units=pyunits.Pa
    )

    for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        eff.effect.properties_in[0].pressure.fix(feed_pressure)
        eff.effect.properties_in[0].temperature.fix(feed_temperature)

        eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(eps)
        eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(eps)
        eff.effect.properties_in[0].conc_mass_phase_comp[...]
        eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()
        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.set_value(0.1)

    first_effect = m.fs.unit.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(0.1)
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    return m


class TestMultiEffectCrystallizer_4Effects:
    @pytest.fixture(scope="class")
    def MEC4_frame(self):
        """
        Test crystallizer with 4 effects
        """
        m = build_mec4()
        return m

    @pytest.mark.unit
    def test_config(self, MEC4_frame):
        m = MEC4_frame
        assert len(m.fs.unit.config) == 6

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.property_package_vapor is m.fs.vapor_properties
        assert m.fs.unit.config.number_effects == 4
        assert m.fs.unit.config.number_effects == len(m.fs.unit.effects)

        for _, eff in m.fs.unit.effects.items():
            assert isinstance(eff.effect, CrystallizerEffect)
            assert not eff.effect.config.standalone

    @pytest.mark.unit
    def test_build(self, MEC4_frame):
        m = MEC4_frame

        # test ports and variables
        port_lst = ["inlet", "outlet", "solids", "vapor", "steam", "pure_water"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        state_blks = [
            "properties_in",
            "properties_out",
            "properties_pure_water",
            "properties_solids",
            "properties_vapor",
        ]

        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        for n, eff in m.fs.unit.effects.items():
            for b in state_blks:
                assert hasattr(eff.effect, b)
                sb = getattr(eff.effect, b)
                assert isinstance(sb, NaClStateBlock)
                assert sb[0].temperature.ub == 1000

        effect_params = [
            "approach_temperature_heat_exchanger",
            "dimensionless_crystal_length",
        ]
        effect_vars = [
            "crystal_growth_rate",
            "crystal_median_length",
            "crystallization_yield",
            "dens_mass_magma",
            "dens_mass_slurry",
            "diameter_crystallizer",
            "height_crystallizer",
            "height_slurry",
            "magma_circulation_flow_vol",
            "pressure_operating",
            "product_volumetric_solids_fraction",
            "relative_supersaturation",
            "souders_brown_constant",
            "t_res",
            "temperature_operating",
            "volume_suspension",
            "work_mechanical",
        ]

        effect_exprs = [
            "delta_temperature",
            "eq_max_allowable_velocity",
            "eq_minimum_height_diameter_ratio",
            "eq_vapor_space_height",
        ]

        effect_constr = [
            "eq_T_con1",
            "eq_T_con2",
            "eq_T_con3",
            "eq_crystallizer_height_constraint",
            "eq_dens_magma",
            "eq_dens_mass_slurry",
            "eq_enthalpy_balance",
            "eq_mass_balance_constraints",
            "eq_minimum_hex_circulation_rate_constraint",
            "eq_operating_pressure_constraint",
            "eq_p_con1",
            "eq_p_con2",
            "eq_p_con3",
            "eq_p_con4",
            "eq_p_con5",
            "eq_pure_vapor_flow_rate",
            "eq_relative_supersaturation",
            "eq_removal_balance",
            "eq_residence_time",
            "eq_slurry_height_constraint",
            "eq_solubility_massfrac_equality_constraint",
            "eq_suspension_volume",
            "eq_vapor_energy_constraint",
            "eq_vapor_head_diameter_constraint",
            "eq_vol_fraction_solids",
        ]

        for n, eff in m.fs.unit.effects.items():
            for p in effect_params:
                assert hasattr(eff.effect, p)
            for v in effect_vars:
                assert hasattr(eff.effect, v)
            for e in effect_exprs:
                assert hasattr(eff.effect, e)
            for c in effect_constr:
                assert hasattr(eff.effect, c)
            assert hasattr(eff.effect, f"eq_delta_temperature_inlet_effect_{n}")
            assert hasattr(eff.effect, f"eq_delta_temperature_outlet_effect_{n}")
            assert hasattr(eff.effect, f"eq_heat_transfer_effect_{n}")

            if n == 1:
                assert number_variables(eff.effect) == 146
                assert number_total_constraints(eff.effect) == 120
                assert number_unused_variables(eff.effect) == 3
                assert hasattr(eff.effect, "heating_steam")
                assert isinstance(eff.effect.heating_steam, WaterStateBlock)
                assert eff.effect.heating_steam[0].temperature.ub == 1000
                assert hasattr(eff.effect, "eq_heating_steam_flow_rate")

            if n != 1:
                assert number_variables(eff.effect) == 140
                assert number_total_constraints(eff.effect) == 118
                assert number_unused_variables(eff.effect) == 1
                assert hasattr(
                    eff.effect, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )

        assert number_variables(m) == 708
        assert number_total_constraints(m) == 474
        assert number_unused_variables(m) == 46

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, MEC4_frame):
        m = MEC4_frame
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 0
            else:
                assert degrees_of_freedom(eff.effect) == 3

    @pytest.mark.unit
    def test_calculate_scaling(self, MEC4_frame):
        m = MEC4_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
        )
        m.fs.vapor_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
        )
        m.fs.vapor_properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )

        calculate_scaling_factors(m)
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        badly_scaled_var_list = list(badly_scaled_var_generator(m))
        assert len(badly_scaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, MEC4_frame):
        m = MEC4_frame
        m.fs.unit.initialize()
        for n, eff in m.fs.unit.effects.items():
            if n == 1:
                htc = value(eff.effect.overall_heat_transfer_coefficient)
                c0 = value(
                    eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                )
                assert degrees_of_freedom(eff.effect) == 0
            if n != 1:
                assert (
                    not eff.effect.properties_in[0]
                    .flow_mass_phase_comp["Liq", "H2O"]
                    .is_fixed()
                )
                assert (
                    not eff.effect.properties_in[0]
                    .flow_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                assert (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                assert (
                    value(
                        eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                    )
                    == c0
                )
                c = getattr(eff.effect, f"eq_energy_for_effect_{n}_from_effect_{n - 1}")
                assert c.active
                assert eff.effect.overall_heat_transfer_coefficient.is_fixed()
                assert value(eff.effect.overall_heat_transfer_coefficient) == htc
                assert degrees_of_freedom(eff.effect) == 3

    @pytest.mark.component
    def test_solve(self, MEC4_frame):

        m = MEC4_frame

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, MEC4_frame):

        m = MEC4_frame

        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]

        for _, e in m.fs.unit.effects.items():
            eff = e.effect

            phase_comp_list = [
                (p, j)
                for j in comp_lst
                for p in phase_lst
                if (p, j) in eff.properties_in[0].phase_component_set
            ]
            flow_mass_in = sum(
                eff.properties_in[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_out = sum(
                eff.properties_out[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_solids = sum(
                eff.properties_solids[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_vapor = sum(
                eff.properties_vapor[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )

            assert (
                abs(
                    value(
                        flow_mass_in
                        - flow_mass_out
                        - flow_mass_solids
                        - flow_mass_vapor
                    )
                )
                <= 1e-6
            )

            assert (
                abs(
                    value(
                        flow_mass_in * eff.properties_in[0].enth_mass_phase["Liq"]
                        - flow_mass_out * eff.properties_out[0].enth_mass_phase["Liq"]
                        - flow_mass_vapor
                        * eff.properties_vapor[0].enth_mass_solvent["Vap"]
                        - flow_mass_solids
                        * eff.properties_solids[0].enth_mass_solute["Sol"]
                        - flow_mass_solids
                        * eff.properties_solids[0].dh_crystallization_mass_comp["NaCl"]
                        + eff.work_mechanical[0]
                    )
                )
                <= 1e-2
            )

        assert (
            pytest.approx(value(m.fs.unit.total_flow_vol_in), rel=1e-3) == 0.00321056811
        )

        assert (
            pytest.approx(
                sum(
                    value(eff.effect.properties_in[0].flow_vol_phase["Liq"])
                    for _, eff in m.fs.unit.effects.items()
                ),
                rel=1e-3,
            )
            == 0.00321056811
        )

    @pytest.mark.component
    def test_solution(self, MEC4_frame):
        m = MEC4_frame

        unit_results_dict = {
            1: {
                "delta_temperature": {0.0: 86.31},
                "delta_temperature_in": {0.0: 57.39},
                "delta_temperature_out": {0.0: 123.72},
                "dens_mass_magma": 281.24,
                "dens_mass_slurry": 1298.23,
                "diameter_crystallizer": 1.0799,
                "energy_flow_superheated_vapor": 1518.73,
                "eq_max_allowable_velocity": 2.6304,
                "eq_minimum_height_diameter_ratio": 1.6199,
                "eq_vapor_space_height": 0.809971,
                "heat_exchanger_area": 197.47,
                "height_crystallizer": 1.883,
                "height_slurry": 1.073,
                "magma_circulation_flow_vol": 0.119395,
                "product_volumetric_solids_fraction": 0.132975,
                "relative_supersaturation": {"NaCl": 0.567038},
                "t_res": 1.0228,
                "temperature_operating": 359.48,
                "volume_suspension": 0.982965,
                "work_mechanical": {0.0: 1704.47},
            },
            2: {
                "delta_temperature": {0.0: 31.59},
                "delta_temperature_in": {0.0: 14.61},
                "delta_temperature_out": {0.0: 58.77},
                "dens_mass_magma": 279.46,
                "dens_mass_slurry": 1301.86,
                "diameter_crystallizer": 1.1773,
                "energy_flow_superheated_vapor": 1396.71,
                "eq_max_allowable_velocity": 3.4641,
                "eq_minimum_height_diameter_ratio": 1.766,
                "eq_vapor_space_height": 0.883027,
                "heat_exchanger_area": 480.72,
                "height_crystallizer": 1.766,
                "height_slurry": 0.82868,
                "magma_circulation_flow_vol": 0.105985,
                "product_volumetric_solids_fraction": 0.132132,
                "relative_supersaturation": {"NaCl": 0.571251},
                "t_res": 1.0228,
                "temperature_operating": 344.86,
                "volume_suspension": 0.9022,
                "work_mechanical": {0.0: 1518.73},
            },
            3: {
                "delta_temperature": {0.0: 16.85},
                "delta_temperature_in": {0.0: 4.3421},
                "delta_temperature_out": {0.0: 44.83},
                "dens_mass_magma": 279.07,
                "dens_mass_slurry": 1303.08,
                "diameter_crystallizer": 1.1811,
                "energy_flow_superheated_vapor": 1296.7,
                "eq_max_allowable_velocity": 3.7763,
                "eq_minimum_height_diameter_ratio": 1.7717,
                "eq_vapor_space_height": 0.885888,
                "heat_exchanger_area": 828.74,
                "height_crystallizer": 1.7717,
                "height_slurry": 0.763759,
                "magma_circulation_flow_vol": 0.09735,
                "product_volumetric_solids_fraction": 0.131951,
                "relative_supersaturation": {"NaCl": 0.572396},
                "t_res": 1.0228,
                "temperature_operating": 340.52,
                "volume_suspension": 0.836915,
                "work_mechanical": {0.0: 1396.71},
            },
            4: {
                "delta_temperature": {0.0: 27.32},
                "delta_temperature_in": {0.0: 17.29},
                "delta_temperature_out": {0.0: 40.69},
                "dens_mass_magma": 278.12,
                "dens_mass_slurry": 1308.54,
                "diameter_crystallizer": 1.3795,
                "energy_flow_superheated_vapor": 1250.34,
                "eq_max_allowable_velocity": 5.4596,
                "eq_minimum_height_diameter_ratio": 2.0692,
                "eq_vapor_space_height": 1.0346,
                "heat_exchanger_area": 474.46,
                "height_crystallizer": 2.0692,
                "height_slurry": 0.537503,
                "magma_circulation_flow_vol": 0.089893,
                "product_volumetric_solids_fraction": 0.131503,
                "relative_supersaturation": {"NaCl": 0.576471},
                "t_res": 1.0228,
                "temperature_operating": 323.22,
                "volume_suspension": 0.803391,
                "work_mechanical": {0.0: 1296.7},
            },
        }

        for n, d in unit_results_dict.items():
            eff = m.fs.unit.effects[n].effect
            for v, r in d.items():
                effv = getattr(eff, v)
                if effv.is_indexed():
                    for i, s in r.items():
                        assert pytest.approx(value(effv[i]), rel=1e-3) == s
                else:
                    assert pytest.approx(value(effv), rel=1e-3) == r

        steam_results_dict = {
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 0.799053},
            "temperature": 416.87,
            "pressure": 401325.0,
            "dh_vap_mass": 2133119.0,
            "pressure_sat": 401325.0,
        }

        for v, r in steam_results_dict.items():
            sv = getattr(m.fs.unit.effects[1].effect.heating_steam[0], v)
            if sv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, MEC4_frame):
        m = MEC4_frame
        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

        m.fs.costing.nacl_recovered.cost.set_value(-0.024)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.total_flow_vol_in)
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.total_flow_vol_in, name="SEC"
        )
        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_dict = {
            "aggregate_capital_cost": 4527295.03,
            "aggregate_flow_electricity": 7.5296,
            "aggregate_flow_NaCl_recovered": 0.266962,
            "aggregate_flow_steam": 0.383078,
            "aggregate_flow_costs": {
                "electricity": 5423.99,
                "NaCl_recovered": -240108.02,
                "steam": 56766.9,
            },
            "aggregate_direct_capital_cost": 2263647.51,
            "total_capital_cost": 4527295.03,
            "total_operating_cost": -42098.28,
            "maintenance_labor_chemical_operating_cost": 135818.85,
            "total_fixed_operating_cost": 135818.85,
            "total_variable_operating_cost": -177917.13,
            "total_annualized_cost": 464759.33,
            "LCOW": 4.5871,
            "SEC": 0.651465,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4527295.03,
            "direct_capital_cost": 2263647.51,
            "capital_cost_crystallizer_effect_1": 659171.48,
            "capital_cost_heat_exchanger_effect_1": 209071.19,
            "capital_cost_effect_1": 434121.33,
            "capital_cost_crystallizer_effect_2": 627459.32,
            "capital_cost_heat_exchanger_effect_2": 498502.0,
            "capital_cost_effect_2": 562980.66,
            "capital_cost_crystallizer_effect_3": 602356.42,
            "capital_cost_heat_exchanger_effect_3": 851138.46,
            "capital_cost_effect_3": 726747.44,
            "capital_cost_crystallizer_effect_4": 587458.93,
            "capital_cost_heat_exchanger_effect_4": 492137.2,
            "capital_cost_effect_4": 539798.07,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing_by_volume(self):

        m = build_mec4()

        for _, eff in m.fs.unit.effects.items():
            eff.effect.crystal_median_length.fix(0.6e-3)
            eff.effect.crystal_growth_rate.fix(5e-9)

        m.fs.unit.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": "volume_basis"},
        )

        m.fs.costing.nacl_recovered.cost.set_value(-0.024)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.total_flow_vol_in)
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.total_flow_vol_in, name="SEC"
        )
        results = solver.solve(m)

        assert_optimal_termination(results)

        sys_costing_dict = {
            "aggregate_capital_cost": 4672521.02,
            "aggregate_fixed_operating_cost": 0.0,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 7.5296,
            "aggregate_flow_NaCl_recovered": 0.266962,
            "aggregate_flow_steam": 0.383078,
            "aggregate_flow_costs": {
                "electricity": 5423.99,
                "NaCl_recovered": -240108.02,
                "steam": 56766.9,
            },
            "aggregate_direct_capital_cost": 2336260.51,
            "total_capital_cost": 4672521.02,
            "total_operating_cost": -37741.5,
            "maintenance_labor_chemical_operating_cost": 140175.63,
            "total_fixed_operating_cost": 140175.63,
            "total_variable_operating_cost": -177917.13,
            "total_annualized_cost": 485375.02,
            "LCOW": 4.7906,
            "SEC": 0.651465,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4672521.02,
            "direct_capital_cost": 2336260.51,
            "capital_cost_crystallizer_effect_1": 675665.09,
            "capital_cost_heat_exchanger_effect_1": 209071.19,
            "capital_cost_effect_1": 442368.14,
            "capital_cost_crystallizer_effect_2": 658735.08,
            "capital_cost_heat_exchanger_effect_2": 498502.0,
            "capital_cost_effect_2": 578618.54,
            "capital_cost_crystallizer_effect_3": 638712.62,
            "capital_cost_heat_exchanger_effect_3": 851138.46,
            "capital_cost_effect_3": 744925.54,
            "capital_cost_crystallizer_effect_4": 648559.35,
            "capital_cost_heat_exchanger_effect_4": 492137.2,
            "capital_cost_effect_4": 570348.28,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r


class TestMultiEffectCrystallizer_3Effects:
    @pytest.fixture(scope="class")
    def MEC3_frame(self):
        """
        Test 3 effect crystallizer
        """
        m = build_mec3()
        return m

    @pytest.mark.unit
    def test_config(self, MEC3_frame):
        m = MEC3_frame
        assert len(m.fs.unit.config) == 6

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.property_package_vapor is m.fs.vapor_properties
        assert m.fs.unit.config.number_effects == 3
        assert m.fs.unit.config.number_effects == len(m.fs.unit.effects)

        assert isinstance(m.fs.unit.effects, FlowsheetBlock)

        for _, eff in m.fs.unit.effects.items():
            assert isinstance(eff.effect, CrystallizerEffect)
            assert not eff.effect.config.standalone

    @pytest.mark.unit
    def test_dof(self, MEC3_frame):
        m = MEC3_frame
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 0
            else:
                assert degrees_of_freedom(eff.effect) == 3

    @pytest.mark.component
    def test_initialize(self, MEC3_frame):
        m = MEC3_frame
        m.fs.unit.initialize()
        for n, eff in m.fs.unit.effects.items():
            if n == 1:
                htc = value(eff.effect.overall_heat_transfer_coefficient)
                c0 = value(
                    eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                )
                assert degrees_of_freedom(eff.effect) == 0
            if n != 1:
                assert (
                    not eff.effect.properties_in[0]
                    .flow_mass_phase_comp["Liq", "H2O"]
                    .is_fixed()
                )
                assert (
                    not eff.effect.properties_in[0]
                    .flow_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                assert (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                assert (
                    pytest.approx(
                        value(
                            eff.effect.properties_in[0].conc_mass_phase_comp[
                                "Liq", "NaCl"
                            ]
                        ),
                        rel=1e-3,
                    )
                    == c0
                )

                c = getattr(eff.effect, f"eq_energy_for_effect_{n}_from_effect_{n - 1}")
                assert c.active
                assert eff.effect.overall_heat_transfer_coefficient.is_fixed()
                assert value(eff.effect.overall_heat_transfer_coefficient) == htc
                assert degrees_of_freedom(eff.effect) == 3

    @pytest.mark.component
    def test_solve(self, MEC3_frame):

        m = MEC3_frame

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, MEC3_frame):

        m = MEC3_frame

        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]

        for _, e in m.fs.unit.effects.items():
            eff = e.effect

            phase_comp_list = [
                (p, j)
                for j in comp_lst
                for p in phase_lst
                if (p, j) in eff.properties_in[0].phase_component_set
            ]
            flow_mass_in = sum(
                eff.properties_in[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_out = sum(
                eff.properties_out[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_solids = sum(
                eff.properties_solids[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )
            flow_mass_vapor = sum(
                eff.properties_vapor[0].flow_mass_phase_comp[p, j]
                for p in phase_lst
                for j in comp_lst
                if (p, j) in phase_comp_list
            )

            assert (
                abs(
                    value(
                        flow_mass_in
                        - flow_mass_out
                        - flow_mass_solids
                        - flow_mass_vapor
                    )
                )
                <= 1e-6
            )

            assert (
                abs(
                    value(
                        flow_mass_in * eff.properties_in[0].enth_mass_phase["Liq"]
                        - flow_mass_out * eff.properties_out[0].enth_mass_phase["Liq"]
                        - flow_mass_vapor
                        * eff.properties_vapor[0].enth_mass_solvent["Vap"]
                        - flow_mass_solids
                        * eff.properties_solids[0].enth_mass_solute["Sol"]
                        - flow_mass_solids
                        * eff.properties_solids[0].dh_crystallization_mass_comp["NaCl"]
                        + eff.work_mechanical[0]
                    )
                )
                <= 1e-2
            )

        assert pytest.approx(value(m.fs.unit.total_flow_vol_in), rel=1e-3) == 0.0232595

        assert (
            pytest.approx(
                sum(
                    value(eff.effect.properties_in[0].flow_vol_phase["Liq"])
                    for _, eff in m.fs.unit.effects.items()
                ),
                rel=1e-3,
            )
            == 0.0232595
        )

    @pytest.mark.component
    def test_solution(self, MEC3_frame):
        m = MEC3_frame

        unit_results_dict = {
            1: {
                "dens_mass_magma": 333.63,
                "dens_mass_slurry": 1321.57,
                "diameter_crystallizer": 2.8505,
                "energy_flow_superheated_vapor": 10580.97,
                "eq_max_allowable_velocity": 2.6304,
                "eq_minimum_height_diameter_ratio": 4.2758,
                "eq_vapor_space_height": 2.1379,
                "heat_exchanger_area": 1322.35,
                "height_crystallizer": 4.5158,
                "height_slurry": 2.3779,
                "magma_circulation_flow_vol": 0.852729,
                "product_volumetric_solids_fraction": 0.157748,
                "relative_supersaturation": {"NaCl": 0.661231},
                "t_res": 1.0228,
                "temperature_operating": 359.48,
                "volume_suspension": 15.17,
                "work_mechanical": {0.0: 12008.03},
            },
            2: {
                "dens_mass_magma": 331.37,
                "dens_mass_slurry": 1324.86,
                "diameter_crystallizer": 3.1042,
                "energy_flow_superheated_vapor": 9709.33,
                "eq_max_allowable_velocity": 3.4641,
                "eq_minimum_height_diameter_ratio": 4.6563,
                "eq_vapor_space_height": 2.3281,
                "heat_exchanger_area": 3349.19,
                "height_crystallizer": 4.6563,
                "height_slurry": 1.8467,
                "magma_circulation_flow_vol": 0.748561,
                "product_volumetric_solids_fraction": 0.156677,
                "relative_supersaturation": {"NaCl": 0.666441},
                "t_res": 1.0228,
                "temperature_operating": 344.86,
                "volume_suspension": 13.97,
                "work_mechanical": {0.0: 10580.97},
            },
            3: {
                "dens_mass_magma": 330.8,
                "dens_mass_slurry": 1325.96,
                "diameter_crystallizer": 3.1152,
                "energy_flow_superheated_vapor": 9019.89,
                "eq_max_allowable_velocity": 3.7763,
                "eq_minimum_height_diameter_ratio": 4.6729,
                "eq_vapor_space_height": 2.3364,
                "heat_exchanger_area": 5761.07,
                "height_crystallizer": 4.6729,
                "height_slurry": 1.7045,
                "magma_circulation_flow_vol": 0.686049,
                "product_volumetric_solids_fraction": 0.15641,
                "relative_supersaturation": {"NaCl": 0.667858},
                "t_res": 1.0228,
                "temperature_operating": 340.52,
                "volume_suspension": 12.99,
                "work_mechanical": {0.0: 9709.33},
            },
        }

        for n, d in unit_results_dict.items():
            eff = m.fs.unit.effects[n].effect
            for v, r in d.items():
                effv = getattr(eff, v)
                if effv.is_indexed():
                    for i, s in r.items():
                        assert pytest.approx(value(effv[i]), rel=1e-3) == s
                else:
                    assert pytest.approx(value(effv), rel=1e-3) == r

        steam_results_dict = {
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 5.66421262},
            "temperature": 421.1592,
            "pressure": 451325.0,
            "dh_vap_mass": 2119982.2092,
            "pressure_sat": 451324.9999,
        }

        for v, r in steam_results_dict.items():
            sv = getattr(m.fs.unit.effects[1].effect.heating_steam[0], v)
            if sv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, MEC3_frame):
        m = MEC3_frame
        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

        m.fs.costing.nacl_recovered.cost.set_value(-0.004)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.total_flow_vol_in)
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.total_flow_vol_in, name="SEC"
        )
        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_dict = {
            "LCOW": 3.502332,
            "aggregate_capital_cost": 19360250.738,
            "aggregate_direct_capital_cost": 9680125.369,
            "aggregate_flow_NaCl_recovered": 3.799812,
            "aggregate_flow_costs": {
                "NaCl_recovered": -569596.699,
                "electricity": 30561.497,
                "steam": 361498.197,
            },
            "aggregate_flow_electricity": 42.425,
            "aggregate_flow_steam": 2.439485,
            "capital_recovery_factor": 0.111955949,
            "maintenance_labor_chemical_operating_cost": 580807.522,
            "total_annualized_cost": 2570765.766,
            "total_capital_cost": 19360250.738,
            "total_fixed_operating_cost": 580807.522,
            "total_operating_cost": 403270.518,
            "total_variable_operating_cost": -177537.004,
            "wacc": 0.093073397,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 19360250.738,
            "capital_cost_crystallizer_effect_1": 3079736.682,
            "capital_cost_crystallizer_effect_2": 2937547.941,
            "capital_cost_crystallizer_effect_3": 2823493.637,
            "capital_cost_effect_1": 2214301.874,
            "capital_cost_effect_2": 3159294.994,
            "capital_cost_effect_3": 4306528.5,
            "capital_cost_heat_exchanger_effect_1": 1348867.065,
            "capital_cost_heat_exchanger_effect_2": 3381042.046,
            "capital_cost_heat_exchanger_effect_3": 5789563.363,
            "direct_capital_cost": 9680125.369,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r
