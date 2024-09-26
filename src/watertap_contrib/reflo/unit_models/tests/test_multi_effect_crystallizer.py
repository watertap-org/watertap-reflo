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

    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.15
    crystallizer_yield = 0.5
    operating_pressures = [0.45, 0.25, 0.208, 0.095]
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

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
    first_effect.heating_steam[0].pressure_sat
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): 101325,
            ("temperature", None): 393,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_vol_phase
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
    operating_pressures = [0.85, 0.25, 0.208]
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

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
    first_effect.heating_steam[0].pressure_sat
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): 101325,
            ("temperature", None): 500,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_vol_phase
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
                assert number_variables(eff.effect) == 150
                assert number_total_constraints(eff.effect) == 124
                assert number_unused_variables(eff.effect) == 1
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

        assert number_variables(m) == 712
        assert number_total_constraints(m) == 478
        assert number_unused_variables(m) == 37

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
                "delta_temperature": {0.0: 60.65},
                "delta_temperature_in": {0.0: 33.51},
                "delta_temperature_out": {0.0: 99.85},
                "dens_mass_magma": 281.24,
                "dens_mass_slurry": 1298.23,
                "diameter_crystallizer": 1.0799,
                "energy_flow_superheated_vapor": 1518.73,
                "eq_max_allowable_velocity": 2.6304,
                "eq_minimum_height_diameter_ratio": 1.6199,
                "eq_vapor_space_height": 0.809971,
                "heat_exchanger_area": 281,
                "height_crystallizer": 1.88304,
                "height_slurry": 1.0730,
                "magma_circulation_flow_vol": 0.119395,
                "product_volumetric_solids_fraction": 0.1329756,
                "relative_supersaturation": {"NaCl": 0.56703},
                "temperature_operating": 359.48,
                "volume_suspension": 0.98296,
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
                "eq_max_allowable_velocity": 3.46414453,
                "eq_minimum_height_diameter_ratio": 1.7660549,
                "eq_vapor_space_height": 0.883027,
                "heat_exchanger_area": 480.72,
                "height_crystallizer": 1.7660549,
                "height_slurry": 0.828680,
                "magma_circulation_flow_vol": 0.105985,
                "product_volumetric_solids_fraction": 0.132132,
                "relative_supersaturation": {"NaCl": 0.571251},
                "t_res": 1.0228,
                "temperature_operating": 344.86,
                "volume_suspension": 0.9022,
                "work_mechanical": {0.0: 1518.73},
            },
            3: {
                "delta_temperature": {0.0: 16.853},
                "delta_temperature_in": {0.0: 4.3421},
                "delta_temperature_out": {0.0: 44.835},
                "dens_mass_magma": 279.07,
                "dens_mass_slurry": 1303.08,
                "diameter_crystallizer": 1.1811,
                "energy_flow_superheated_vapor": 1296.7,
                "eq_max_allowable_velocity": 3.776,
                "eq_minimum_height_diameter_ratio": 1.7717,
                "eq_vapor_space_height": 0.88588,
                "heat_exchanger_area": 828.74,
                "height_crystallizer": 1.7717,
                "height_slurry": 0.763759,
                "magma_circulation_flow_vol": 0.097350,
                "product_volumetric_solids_fraction": 0.13195,
                "relative_supersaturation": {"NaCl": 0.57239},
                "temperature_operating": 340.5214,
                "volume_suspension": 0.836915,
                "work_mechanical": {0.0: 1396.71},
            },
            4: {
                "delta_temperature": {0.0: 27.3299},
                "delta_temperature_in": {0.0: 17.299},
                "delta_temperature_out": {0.0: 40.695},
                "dens_mass_magma": 278.12,
                "dens_mass_slurry": 1308.54,
                "diameter_crystallizer": 1.3795,
                "energy_flow_superheated_vapor": 1250.3,
                "eq_max_allowable_velocity": 5.4596,
                "eq_minimum_height_diameter_ratio": 2.0692,
                "eq_vapor_space_height": 1.03464,
                "heat_exchanger_area": 474.4,
                "height_crystallizer": 2.0692,
                "height_slurry": 0.53750,
                "magma_circulation_flow_vol": 0.089893,
                "product_volumetric_solids_fraction": 0.131503,
                "relative_supersaturation": {"NaCl": 0.576471},
                "t_res": 1.0228,
                "temperature_operating": 323.2,
                "volume_suspension": 0.803391,
                "work_mechanical": {0.0: 1296.70},
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
            "dens_mass_phase": {"Liq": 943.14, "Vap": 0.55862},
            "dh_vap_mass": 2202682,
            "flow_mass_phase_comp": {
                ("Liq", "H2O"): 0.0,
                ("Vap", "H2O"): 0.77381,
            },
            "flow_vol_phase": {"Liq": 0.0, "Vap": 1.3852},
            "pressure_sat": 197743,
            "temperature": 393,
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
            "LCOW": 5.1605,
            "SEC": 0.651465,
            "NaCl_recovered_cost": -0.024,
            "aggregate_capital_cost": 4613044.1,
            "aggregate_direct_capital_cost": 2306522,
            "aggregate_flow_NaCl_recovered": 0.26696,
            "aggregate_flow_costs": {
                "NaCl_recovered": -240108.0,
                "electricity": 5423.99,
                "steam": 102690.74,
            },
            "aggregate_flow_electricity": 7.5296,
            "aggregate_flow_steam": 0.69298,
            "capital_recovery_factor": 0.111955,
            "maintenance_labor_chemical_operating_cost": 138391.3,
            "total_annualized_cost": 522855.76,
            "total_capital_cost": 4613044.1,
            "total_fixed_operating_cost": 138391.32,
            "total_operating_cost": 6398.03,
            "total_variable_operating_cost": -131993.29,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4613044.11,
            "direct_capital_cost": 2306522.05,
            "capital_cost_crystallizer_effect_1": 659171.48,
            "capital_cost_heat_exchanger_effect_1": 294820.27,
            "capital_cost_effect_1": 476995.87,
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
            "LCOW": 5.364,
            "SEC": 0.65146,
            "aggregate_capital_cost": 4758270.1,
            "aggregate_direct_capital_cost": 2379135.0,
            "aggregate_flow_NaCl_recovered": 0.26696,
            "aggregate_flow_costs": {
                "NaCl_recovered": -240108.02,
                "electricity": 5423.9,
                "steam": 102690.74,
            },
            "aggregate_flow_electricity": 7.529,
            "aggregate_flow_steam": 0.69298,
            "capital_recovery_factor": 0.11195,
            "maintenance_labor_chemical_operating_cost": 142748.1,
            "total_annualized_cost": 543471.45,
            "total_capital_cost": 4758270.10,
            "total_fixed_operating_cost": 142748.1,
            "total_operating_cost": 10754.8,
            "total_variable_operating_cost": -131993.29,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4758270.1,
            "capital_cost_effect_1": 485242.6,
            "capital_cost_effect_2": 578618.5,
            "capital_cost_effect_3": 744925.5,
            "capital_cost_effect_4": 570348.2,
            "direct_capital_cost": 2379135.0,
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

        assert pytest.approx(value(m.fs.unit.total_flow_vol_in), rel=1e-3) == 0.02321134

        assert (
            pytest.approx(
                sum(
                    value(eff.effect.properties_in[0].flow_vol_phase["Liq"])
                    for _, eff in m.fs.unit.effects.items()
                ),
                rel=1e-3,
            )
            == 0.02321134
        )

    @pytest.mark.component
    def test_solution(self, MEC3_frame):
        m = MEC3_frame

        unit_results_dict = {
            1: {
                "delta_temperature": {0.0: 161.4093},
                "delta_temperature_in": {0.0: 123.1938},
                "delta_temperature_out": {0.0: 206.85},
                "dens_mass_magma": 337.00,
                "dens_mass_slurry": 1318.39,
                "diameter_crystallizer": 2.4814,
                "energy_flow_superheated_vapor": 10546.62,
                "eq_max_allowable_velocity": 1.95486,
                "eq_minimum_height_diameter_ratio": 3.72212,
                "eq_vapor_space_height": 1.8610,
                "heat_exchanger_area": 776.59,
                "height_crystallizer": 4.96774,
                "height_slurry": 3.10668,
                "magma_circulation_flow_vol": 0.89371,
                "product_volumetric_solids_fraction": 0.15933,
                "relative_supersaturation": {"NaCl": 0.654196},
                "temperature_operating": 376.80,
                "volume_suspension": 15.0239,
                "work_mechanical": {0.0: 12535},
            },
            2: {
                "delta_temperature": {0.0: 50.4936},
                "delta_temperature_in": {0.0: 31.9425},
                "delta_temperature_out": {0.0: 75.2194},
                "dens_mass_magma": 331.37,
                "dens_mass_slurry": 1324.86,
                "diameter_crystallizer": 3.0991,
                "energy_flow_superheated_vapor": 9677.8,
                "eq_max_allowable_velocity": 3.4641,
                "eq_minimum_height_diameter_ratio": 4.6487,
                "eq_vapor_space_height": 2.3243,
                "heat_exchanger_area": 2088.70,
                "height_crystallizer": 4.6487,
                "height_slurry": 1.8467,
                "magma_circulation_flow_vol": 0.7461,
                "product_volumetric_solids_fraction": 0.15667,
                "relative_supersaturation": {"NaCl": 0.6664},
                "temperature_operating": 344.86,
                "volume_suspension": 13.9312,
                "work_mechanical": {0.0: 10546.6},
            },
            3: {
                "delta_temperature": {0.0: 16.8533},
                "delta_temperature_in": {0.0: 4.3421},
                "delta_temperature_out": {0.0: 44.8351},
                "dens_mass_magma": 330.8,
                "dens_mass_slurry": 1325.9,
                "diameter_crystallizer": 3.1102,
                "energy_flow_superheated_vapor": 8990.6,
                "eq_max_allowable_velocity": 3.7763,
                "eq_minimum_height_diameter_ratio": 4.6653,
                "eq_vapor_space_height": 2.3326,
                "heat_exchanger_area": 5742.36,
                "height_crystallizer": 4.6653,
                "height_slurry": 1.7045,
                "magma_circulation_flow_vol": 0.68382,
                "product_volumetric_solids_fraction": 0.156410,
                "relative_supersaturation": {"NaCl": 0.667858},
                "temperature_operating": 340.5,
                "volume_suspension": 12.950,
                "work_mechanical": {0.0: 9677.8},
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
            "dens_mass_phase": {"Liq": 828.03, "Vap": 0.439081},
            "dh_vap_mass": 1827723.26,
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 6.858},
            "flow_vol_phase": {"Liq": 0.0, "Vap": 15.619},
            "pressure_sat": 2639870.3,
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
            "LCOW": 2.7665,
            "SEC": 0.51529,
            "aggregate_capital_cost": 17518907.2,
            "aggregate_direct_capital_cost": 8759453.6,
            "aggregate_flow_NaCl_recovered": 3.7919,
            "aggregate_flow_costs": {
                "NaCl_recovered": -568416.4,
                "electricity": 31017.1,
                "steam": 76976.0,
            },
            "aggregate_flow_electricity": 43.05,
            "aggregate_flow_steam": 0.51945,
            "capital_recovery_factor": 0.11195,
            "total_annualized_cost": 2026489.8,
            "total_capital_cost": 17518907.2,
            "total_fixed_operating_cost": 525567.2,
            "total_operating_cost": 65143.9,
            "total_variable_operating_cost": -460423.2,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 17518907.2,
            "capital_cost_crystallizer_effect_1": 3079736.6,
            "capital_cost_crystallizer_effect_2": 2932488.6,
            "capital_cost_crystallizer_effect_3": 2818630.8,
            "capital_cost_effect_1": 1939078.0,
            "capital_cost_effect_2": 2525607.6,
            "capital_cost_effect_3": 4294767.8,
            "capital_cost_heat_exchanger_effect_1": 798419.4,
            "capital_cost_heat_exchanger_effect_2": 2118726.6,
            "capital_cost_heat_exchanger_effect_3": 5770904.9,
            "direct_capital_cost": 8759453.6,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r
