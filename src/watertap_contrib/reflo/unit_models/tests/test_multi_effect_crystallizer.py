#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
    Objective,
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
from idaes.core.util.exceptions import ConfigurationError

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
eps = 1e-12


def build_mec2():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.unit = mec = MultiEffectCrystallizer(
        property_package=m.fs.properties,
        property_package_vapor=m.fs.vapor_properties,
        number_effects=2,
    )

    num_effects = 2
    feed_flow_mass = 2
    total_feed_flow_mass = num_effects * feed_flow_mass
    feed_mass_frac_NaCl = 0.45
    crystallizer_yield = 0.6
    operating_pressures = [0.45, 0.25]
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    atm_pressure = 101325 * pyunits.Pa
    saturated_steam_pressure_gage = 5 * pyunits.bar
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage, to_units=pyunits.Pa
    )

    for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].set_value(
            feed_flow_mass * feed_mass_frac_H2O
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(
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
        eff.effect.overall_heat_transfer_coefficient.set_value(100)

    first_effect = m.fs.unit.effects[1].effect
    first_effect.overall_heat_transfer_coefficient.fix(100)
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

    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        total_feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        total_feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.unit.control_volume.properties_in[0].pressure.fix(feed_pressure)
    m.fs.unit.control_volume.properties_in[0].temperature.fix(feed_temperature)

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

    num_effects = 3
    feed_flow_mass = 1
    total_feed_flow_mass = num_effects * feed_flow_mass
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
        eff.effect.overall_heat_transfer_coefficient.set_value(100)

    first_effect = m.fs.unit.effects[1].effect
    first_effect.overall_heat_transfer_coefficient.fix(100)
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

    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        total_feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        total_feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.unit.control_volume.properties_in[0].pressure.fix(feed_pressure)
    m.fs.unit.control_volume.properties_in[0].temperature.fix(feed_temperature)

    return m


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

    total_feed_flow_mass = 4

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
        eff.effect.overall_heat_transfer_coefficient.set_value(100)

    first_effect = m.fs.unit.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(100)
    first_effect.heating_steam[0].pressure_sat
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

    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        total_feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        total_feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.unit.control_volume.properties_in[0].pressure.fix(feed_pressure)
    m.fs.unit.control_volume.properties_in[0].temperature.fix(feed_temperature)

    return m


@pytest.mark.component
def test_single_effect_error():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    error_msg = (
        "The MultiEffectCrystallizer model requires more than 1 effect."
        "To model a crystallizer with one effect, use the CrystallizerEffect model with 'standalone=True'."
    )

    with pytest.raises(ConfigurationError, match=error_msg):
        m.fs.unit = MultiEffectCrystallizer(
            property_package=m.fs.properties,
            property_package_vapor=m.fs.vapor_properties,
            number_effects=1,
        )


class TestMultiEffectCrystallizer_2Effects:
    @pytest.fixture(scope="class")
    def MEC2_frame(self):
        """
        Test crystallizer with 2 effects
        """
        m = build_mec2()
        return m

    @pytest.mark.unit
    def test_config(self, MEC2_frame):
        m = MEC2_frame
        assert len(m.fs.unit.config) == 6

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.property_package_vapor is m.fs.vapor_properties
        assert m.fs.unit.config.number_effects == 2
        assert m.fs.unit.config.number_effects == len(m.fs.unit.effects)

        for _, eff in m.fs.unit.effects.items():
            assert isinstance(eff.effect, CrystallizerEffect)
            assert not eff.effect.config.standalone

    @pytest.mark.unit
    def test_build(self, MEC2_frame):
        m = MEC2_frame

        # test ports and variables
        port_lst = ["inlet", "outlet", "solids", "vapor", "steam"]
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
            "eq_mass_balance_constraints",
            "eq_solubility_massfrac_equality_constraint",
            "eq_removal_balance",
            "eq_vol_fraction_solids",
            "eq_dens_magma",
            "eq_operating_pressure_constraint",
            "eq_relative_supersaturation",
            "eq_enthalpy_balance",
            "eq_p_con3",
            "eq_T_con1",
            "eq_T_con2",
            "eq_T_con3",
            "eq_minimum_hex_circulation_rate_constraint",
            "eq_dens_mass_slurry",
            "eq_residence_time",
            "eq_suspension_volume",
            "eq_vapor_head_diameter_constraint",
            "eq_slurry_height_constraint",
            "eq_crystallizer_height_constraint",
            "eq_pure_water_mass_flow_rate",
            "eq_vapor_energy_constraint",
            "eq_p_con1",
            "eq_p_con2",
            "eq_p_con4",
            "eq_p_con5",
        ]

        assert hasattr(m.fs.unit, "recovery_vol_phase")
        assert isinstance(m.fs.unit.recovery_vol_phase, Var)

        for n, eff in m.fs.unit.effects.items():
            for p in effect_params:
                assert hasattr(eff.effect, p)
            for v in effect_vars:
                assert hasattr(eff.effect, v)
            for e in effect_exprs:
                assert hasattr(eff.effect, e)
            for c in effect_constr:
                assert hasattr(eff.effect, c)
            if n == 1:
                assert number_variables(eff.effect) == 154
                assert number_total_constraints(eff.effect) == 128
                assert number_unused_variables(eff.effect) == 2
                assert hasattr(eff.effect, "heating_steam")
                assert isinstance(eff.effect.heating_steam, WaterStateBlock)
                assert eff.effect.heating_steam[0].temperature.ub == 1000
                assert hasattr(eff.effect, "eq_heating_steam_flow_rate")
            if n != 1:
                assert number_variables(eff.effect) == 148
                assert number_total_constraints(eff.effect) == 122
                assert number_unused_variables(eff.effect) == 4
                assert hasattr(m.fs.unit, f"eq_delta_temperature_inlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_delta_temperature_outlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_heat_transfer_effect_{n}")
                assert hasattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )

        assert number_variables(m) == 461
        assert number_total_constraints(m) == 268
        assert number_unused_variables(m) == 43

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, MEC2_frame):
        m = MEC2_frame
        # With mass flow rates into each of the effects unfixed,
        # the model is under specified
        assert degrees_of_freedom(m) == 2

        # Fixing flow rates into individual effects will reduce DOF...
        for n, eff in m.fs.unit.effects.items():
            eff.effect.properties_in[0].flow_mass_phase_comp.fix()
            assert degrees_of_freedom(eff.effect) == 0
        # ... and result in an overspecified model.
        assert degrees_of_freedom(m) == -2

        for n, eff in m.fs.unit.effects.items():
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
        assert degrees_of_freedom(m) == 2

    @pytest.mark.unit
    def test_calculate_scaling(self, MEC2_frame):
        m = MEC2_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 10, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Sol", "NaCl")
        )
        m.fs.vapor_properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
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
    def test_initialize(self, MEC2_frame):
        m = MEC2_frame
        m.fs.unit.initialize()
        c0 = value(
            m.fs.unit.effects[1]
            .effect.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
        )
        htc = value(m.fs.unit.effects[1].effect.overall_heat_transfer_coefficient)
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
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
                pytest.approx(
                    value(
                        eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                    ),
                    rel=1e-6,
                )
                == c0
            )
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 2
            if n != 1:
                assert (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                linking_constr = getattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )
                assert linking_constr.active
                assert eff.effect.overall_heat_transfer_coefficient.is_fixed()
                assert value(eff.effect.overall_heat_transfer_coefficient) == htc
                assert degrees_of_freedom(eff.effect) == 1

    @pytest.mark.component
    def test_solve(self, MEC2_frame):

        m = MEC2_frame

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, MEC2_frame):

        m = MEC2_frame

        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]

        total_mass_flow_water_in = 0
        total_mass_flow_salt_in = 0
        total_mass_flow_water_out = 0

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
                        + pyunits.convert(
                            eff.work_mechanical[0], to_units=pyunits.J * pyunits.s**-1
                        )
                    )
                )
                <= 1e-2
            )

            total_mass_flow_water_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            total_mass_flow_salt_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
            )
            total_mass_flow_water_out += value(
                eff.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        # Test control volume mass balance
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "NaCl"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_salt_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_out[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_out
        )

    @pytest.mark.component
    def test_solution(self, MEC2_frame):
        m = MEC2_frame

        assert (
            pytest.approx(value(m.fs.unit.recovery_vol_phase["Liq"]), rel=1e-3)
            == 0.10289
        )

        unit_results_dict = {
            1: {
                "temperature_operating": 359.4825,
                "pressure_operating": 45000.0,
                "dens_mass_magma": 395.3173,
                "dens_mass_slurry": 1349.0471,
                "work_mechanical": {0.0: 370.3557},
                "diameter_crystallizer": 0.50292237869,
                "height_slurry": 22.8532,
                "height_crystallizer": 23.2304,
                "magma_circulation_flow_vol": 0.026733701982,
                "relative_supersaturation": {"NaCl": 0.767477499944},
                "t_res": 1.02282118,
                "volume_suspension": 4.53982873,
                "eq_max_allowable_velocity": 2.63046746,
                "eq_vapor_space_height": 0.377191784018,
                "eq_minimum_height_diameter_ratio": 0.754383568036,
                "energy_flow_superheated_vapor": 329356.7611,
                "delta_temperature_in": {0.0: 72.5775},
                "delta_temperature_out": {0.0: 138.9101},
                "delta_temperature": {0.0: 102.1559},
                "heat_exchanger_area": 36.2539,
                "overall_heat_transfer_coefficient": 100.0,
            },
            2: {
                "temperature_operating": 344.8635,
                "pressure_operating": 25000.0,
                "dens_mass_magma": 392.7294,
                "dens_mass_slurry": 1352.0474,
                "work_mechanical": {0.0: 329.3567},
                "diameter_crystallizer": 0.601158860753,
                "height_slurry": 19.5821,
                "height_crystallizer": 20.0329,
                "magma_circulation_flow_vol": 0.023686095554,
                "relative_supersaturation": {"NaCl": 0.773914820448},
                "t_res": 1.02282118,
                "volume_suspension": 5.5581223,
                "eq_max_allowable_velocity": 3.46414453,
                "eq_vapor_space_height": 0.450869145565,
                "eq_minimum_height_diameter_ratio": 0.90173829113,
                "energy_flow_superheated_vapor": 364132.1326,
                "delta_temperature_in": {0.0: 14.6189},
                "delta_temperature_out": {0.0: 58.7762},
                "delta_temperature": {0.0: 31.5926},
                "heat_exchanger_area": 104.2511,
                "overall_heat_transfer_coefficient": 100.0,
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
            "flow_mass_phase_comp": {("Vap", "H2O"): 0.177583},
            "temperature": 432.06,
            "pressure": 601325.0,
            "dh_vap_mass": 2085530.7,
            "pressure_sat": 601325.0,
        }

        for v, r in steam_results_dict.items():
            sv = getattr(m.fs.unit.effects[1].effect.heating_steam[0], v)
            if sv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, MEC2_frame):
        m = MEC2_frame
        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": "mass_basis"},
        )

        m.fs.costing.nacl_recovered.cost.set_value(-0.011)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.total_flow_vol_in)
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.total_flow_vol_in, name="SEC"
        )
        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_dict = {
            "aggregate_capital_cost": 3902338.7,
            "aggregate_flow_electricity": 0.95390,
            "aggregate_flow_NaCl_recovered": 1.0799,
            "aggregate_flow_steam": 0.05888,
            "aggregate_flow_costs": {
                "electricity": 687.144,
                "NaCl_recovered": -445206.6,
                "steam": 8726.6,
            },
            "aggregate_direct_capital_cost": 1951169.3,
            "total_capital_cost": 3902338.7,
            "total_operating_cost": -318722.72,
            "maintenance_labor_chemical_operating_cost": 117070.16,
            "total_fixed_operating_cost": 117070.16,
            "total_variable_operating_cost": -435792.88,
            "total_annualized_cost": 118167.30,
            "LCOW": 1.27058636,
            "SEC": 0.08991097,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 3902338.7,
            "direct_capital_cost": 1951169.3,
            "capital_cost_crystallizer_effect_1": 1777273.3,
            "capital_cost_heat_exchanger_effect_1": 40936.8,
            "capital_cost_effect_1": 909105.0,
            "capital_cost_crystallizer_effect_2": 1971550.7,
            "capital_cost_heat_exchanger_effect_2": 112577.7,
            "capital_cost_effect_2": 1042064.2,
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
    def test_build(self, MEC3_frame):
        m = MEC3_frame

        # test ports and variables
        port_lst = ["inlet", "outlet", "solids", "vapor", "steam"]
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
            "eq_mass_balance_constraints",
            "eq_solubility_massfrac_equality_constraint",
            "eq_removal_balance",
            "eq_vol_fraction_solids",
            "eq_dens_magma",
            "eq_operating_pressure_constraint",
            "eq_relative_supersaturation",
            "eq_enthalpy_balance",
            "eq_p_con3",
            "eq_T_con1",
            "eq_T_con2",
            "eq_T_con3",
            "eq_minimum_hex_circulation_rate_constraint",
            "eq_dens_mass_slurry",
            "eq_residence_time",
            "eq_suspension_volume",
            "eq_vapor_head_diameter_constraint",
            "eq_slurry_height_constraint",
            "eq_crystallizer_height_constraint",
            "eq_pure_water_mass_flow_rate",
            "eq_vapor_energy_constraint",
            "eq_p_con1",
            "eq_p_con2",
            "eq_p_con4",
            "eq_p_con5",
        ]

        assert hasattr(m.fs.unit, "recovery_vol_phase")
        assert isinstance(m.fs.unit.recovery_vol_phase, Var)

        for n, eff in m.fs.unit.effects.items():
            for p in effect_params:
                assert hasattr(eff.effect, p)
            for v in effect_vars:
                assert hasattr(eff.effect, v)
            for e in effect_exprs:
                assert hasattr(eff.effect, e)
            for c in effect_constr:
                assert hasattr(eff.effect, c)
            if n == 1:
                assert number_variables(eff.effect) == 154
                assert number_total_constraints(eff.effect) == 128
                assert number_unused_variables(eff.effect) == 2
                assert hasattr(eff.effect, "heating_steam")
                assert isinstance(eff.effect.heating_steam, WaterStateBlock)
                assert eff.effect.heating_steam[0].temperature.ub == 1000
                assert hasattr(eff.effect, "eq_heating_steam_flow_rate")
            if n != 1:
                assert number_variables(eff.effect) == 148
                assert number_total_constraints(eff.effect) == 122
                assert number_unused_variables(eff.effect) == 4
                assert hasattr(m.fs.unit, f"eq_delta_temperature_inlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_delta_temperature_outlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_heat_transfer_effect_{n}")
                assert hasattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )

        assert number_variables(m) == 609
        assert number_total_constraints(m) == 394
        assert number_unused_variables(m) == 43

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, MEC3_frame):
        m = MEC3_frame
        # With mass flow rates into each of the effects fixed,
        # the model is over specified
        assert degrees_of_freedom(m) == -2
        # Unfixing the mass flow rates for CV will result in 0 DOF
        # (This is only done for testing purposes)
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].unfix()
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "NaCl"
        ].unfix()
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
            assert degrees_of_freedom(eff.effect) == 0
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].fix()
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "NaCl"
        ].fix()
        assert degrees_of_freedom(m) == -2
        # Alternatively, one can .set_value() each of the mass flows for the individual effects,
        # resulting in 4 DOF:
        for n, eff in m.fs.unit.effects.items():
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
        assert degrees_of_freedom(m) == 4

    @pytest.mark.component
    def test_calculate_scaling(self, MEC3_frame):
        m = MEC3_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Sol", "NaCl")
        )
        m.fs.vapor_properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", "H2O")
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
    def test_initialize(self, MEC3_frame):
        m = MEC3_frame
        m.fs.unit.initialize()
        c0 = value(
            m.fs.unit.effects[1]
            .effect.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
        )
        htc = value(m.fs.unit.effects[1].effect.overall_heat_transfer_coefficient)
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
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
                pytest.approx(
                    value(
                        eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                    ),
                    rel=1e-6,
                )
                == c0
            )
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 2
            if n != 1:
                assert (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                linking_constr = getattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )
                assert linking_constr.active
                assert eff.effect.overall_heat_transfer_coefficient.is_fixed()
                assert value(eff.effect.overall_heat_transfer_coefficient) == htc
                assert degrees_of_freedom(eff.effect) == 1

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

        total_mass_flow_water_in = 0
        total_mass_flow_salt_in = 0
        total_mass_flow_water_out = 0

        # Test mass balance for each effect
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

            total_mass_flow_water_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            total_mass_flow_salt_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
            )
            total_mass_flow_water_out += value(
                eff.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        # Test control volume mass balance
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "NaCl"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_salt_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_out[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_out
        )

    @pytest.mark.component
    def test_solution(self, MEC3_frame):
        m = MEC3_frame

        assert (
            pytest.approx(value(m.fs.unit.recovery_vol_phase["Liq"]), rel=1e-3)
            == 0.5485
        )

        unit_results_dict = {
            1: {
                "product_volumetric_solids_fraction": 0.15779,
                "temperature_operating": 359.4825,
                "dens_mass_magma": 333.7,
                "dens_mass_slurry": 1321.6,
                "work_mechanical": {0.0: 1303.5},
                "diameter_crystallizer": 0.93920,
                "height_slurry": 2.3783,
                "height_crystallizer": 3.082,
                "magma_circulation_flow_vol": 0.09257,
                "relative_supersaturation": {"NaCl": 0.66123},
                "volume_suspension": 1.647,
                "eq_max_allowable_velocity": 2.6304,
                "eq_vapor_space_height": 0.704405,
                "eq_minimum_height_diameter_ratio": 1.4088,
                "energy_flow_superheated_vapor": 1148650.3773,
                "delta_temperature_in": {0.0: 61.676},
                "delta_temperature_out": {0.0: 128.00},
                "delta_temperature": {0.0: 90.8077},
                "heat_exchanger_area": 143.552,
            },
            2: {
                "product_volumetric_solids_fraction": 0.15677,
                "temperature_operating": 344.86,
                "dens_mass_magma": 331.57,
                "dens_mass_slurry": 1324.95,
                "work_mechanical": {0.0: 1148.65},
                "diameter_crystallizer": 1.02278,
                "height_slurry": 1.8472,
                "height_crystallizer": 2.6143,
                "magma_circulation_flow_vol": 0.081266,
                "relative_supersaturation": {"NaCl": 0.66644},
                "volume_suspension": 1.5176,
                "eq_max_allowable_velocity": 3.4641,
                "eq_vapor_space_height": 0.76709,
                "eq_minimum_height_diameter_ratio": 1.5341,
                "energy_flow_superheated_vapor": 1054026.64,
                "delta_temperature_in": {0.0: 14.6189},
                "delta_temperature_out": {0.0: 58.7762},
                "delta_temperature": {0.0: 31.5926},
                "heat_exchanger_area": 363.5818,
            },
            3: {
                "product_volumetric_solids_fraction": 0.15653,
                "temperature_operating": 340.52,
                "dens_mass_magma": 331.0627,
                "dens_mass_slurry": 1326.0758,
                "work_mechanical": {0.0: 1054.0266},
                "diameter_crystallizer": 1.0264,
                "height_slurry": 1.7051,
                "height_crystallizer": 2.47501,
                "magma_circulation_flow_vol": 0.07448,
                "relative_supersaturation": {"NaCl": 0.66785},
                "volume_suspension": 1.41098144,
                "eq_max_allowable_velocity": 3.77639066,
                "eq_vapor_space_height": 0.76982,
                "eq_minimum_height_diameter_ratio": 1.5396,
                "energy_flow_superheated_vapor": 979181.8,
                "delta_temperature_in": {0.0: 4.3421},
                "delta_temperature_out": {0.0: 44.835},
                "delta_temperature": {0.0: 16.85334},
                "heat_exchanger_area": 625.411,
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
            "flow_mass_phase_comp": {
                ("Liq", "H2O"): 0.0,
                ("Vap", "H2O"): 0.614896112148,
            },
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
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": "mass_basis"},
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
            "aggregate_capital_cost": 3903131.16,
            "aggregate_flow_electricity": 4.6061,
            "aggregate_flow_NaCl_recovered": 0.412499,
            "aggregate_flow_steam": 0.264825,
            "aggregate_flow_costs": {
                "electricity": 3318.05,
                "NaCl_recovered": -371005.47,
                "steam": 39243.55,
            },
            "aggregate_direct_capital_cost": 1951565.58,
            "total_capital_cost": 3903131.16,
            "total_operating_cost": -211349.92,
            "maintenance_labor_chemical_operating_cost": 117093.93,
            "total_fixed_operating_cost": 117093.93,
            "total_variable_operating_cost": -328443.86,
            "total_annualized_cost": 225628.82,
            "LCOW": 2.8315,
            "LCOW_component_direct_capex": {"fs.unit": 2.7419},
            "LCOW_component_indirect_capex": {"fs.unit": 2.7419},
            "LCOW_component_fixed_opex": {"fs.unit": 1.4694},
            "LCOW_component_variable_opex": {
                "fs.unit": 0.0,
                "fs.unit.effects[1].effect": -1.1768,
                "fs.unit.effects[2].effect": -1.5274,
                "fs.unit.effects[3].effect": -1.4176,
            },
            "LCOW_aggregate_direct_capex": {"MultiEffectCrystallizer": 2.7419},
            "LCOW_aggregate_indirect_capex": {"MultiEffectCrystallizer": 2.7419},
            "LCOW_aggregate_fixed_opex": {"MultiEffectCrystallizer": 1.4694},
            "LCOW_aggregate_variable_opex": {
                "MultiEffectCrystallizer": 0.0,
                "electricity": 0.04164,
                "CrystallizerEffect": -4.1218,
                "NaCl_recovered": -4.656,
                "steam": 0.492494,
            },
            "SEC": 0.506727,
            "SEC_component": {
                "fs.unit.effects[1].effect": 0.18856,
                "fs.unit.effects[2].effect": 0.165947,
                "fs.unit.effects[3].effect": 0.15222,
            },
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost_crystallizer_effect_1": 949323.7,
            "capital_cost_heat_exchanger_effect_1": 153401.52,
            "capital_cost_effect_1": 551362.61,
            "capital_cost_crystallizer_effect_2": 905494.07,
            "capital_cost_heat_exchanger_effect_2": 379213.97,
            "capital_cost_effect_2": 642354.02,
            "capital_cost_crystallizer_effect_3": 870336.83,
            "capital_cost_heat_exchanger_effect_3": 645361.05,
            "capital_cost_effect_3": 757848.94,
        }
        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r


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
        port_lst = ["inlet", "outlet", "solids", "vapor", "steam"]
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
            "eq_mass_balance_constraints",
            "eq_solubility_massfrac_equality_constraint",
            "eq_removal_balance",
            "eq_vol_fraction_solids",
            "eq_dens_magma",
            "eq_operating_pressure_constraint",
            "eq_relative_supersaturation",
            "eq_enthalpy_balance",
            "eq_p_con3",
            "eq_T_con1",
            "eq_T_con2",
            "eq_T_con3",
            "eq_minimum_hex_circulation_rate_constraint",
            "eq_dens_mass_slurry",
            "eq_residence_time",
            "eq_suspension_volume",
            "eq_vapor_head_diameter_constraint",
            "eq_slurry_height_constraint",
            "eq_crystallizer_height_constraint",
            "eq_pure_water_mass_flow_rate",
            "eq_vapor_energy_constraint",
            "eq_p_con1",
            "eq_p_con2",
            "eq_p_con4",
            "eq_p_con5",
        ]

        assert hasattr(m.fs.unit, "recovery_vol_phase")
        assert isinstance(m.fs.unit.recovery_vol_phase, Var)

        for n, eff in m.fs.unit.effects.items():
            for p in effect_params:
                assert hasattr(eff.effect, p)
            for v in effect_vars:
                assert hasattr(eff.effect, v)
            for e in effect_exprs:
                assert hasattr(eff.effect, e)
            for c in effect_constr:
                assert hasattr(eff.effect, c)
            if n == 1:
                assert number_variables(eff.effect) == 154
                assert number_total_constraints(eff.effect) == 128
                assert number_unused_variables(eff.effect) == 2
                assert hasattr(eff.effect, "heating_steam")
                assert isinstance(eff.effect.heating_steam, WaterStateBlock)
                assert eff.effect.heating_steam[0].temperature.ub == 1000
                assert hasattr(eff.effect, "eq_heating_steam_flow_rate")
            if n != 1:
                assert number_variables(eff.effect) == 148
                assert number_total_constraints(eff.effect) == 122
                assert number_unused_variables(eff.effect) == 4
                assert hasattr(m.fs.unit, f"eq_delta_temperature_inlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_delta_temperature_outlet_effect_{n}")
                assert hasattr(m.fs.unit, f"eq_heat_transfer_effect_{n}")
                assert hasattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )

        assert number_variables(m) == 757
        assert number_total_constraints(m) == 520
        assert number_unused_variables(m) == 43

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, MEC4_frame):

        m = MEC4_frame
        # With mass flow rates into each of the effects fixed,
        # and the control_volume flows fixed, the model is over specified
        assert degrees_of_freedom(m) == -2
        # Unfixing the mass flow rates for CV will result in 0 DOF
        # (This is only done for testing purposes)
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].unfix()
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "NaCl"
        ].unfix()
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
            assert degrees_of_freedom(eff.effect) == 0
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ].fix()
        m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
            "Liq", "NaCl"
        ].fix()
        assert degrees_of_freedom(m) == -2
        # Alternatively, one can .set_value() each of the mass flows for the individual effects,
        # resulting in 6 DOF:
        for n, eff in m.fs.unit.effects.items():
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
        assert degrees_of_freedom(m) == 6

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
        c0 = value(
            m.fs.unit.effects[1]
            .effect.properties_in[0]
            .conc_mass_phase_comp["Liq", "NaCl"]
        )
        htc = value(m.fs.unit.effects[1].effect.overall_heat_transfer_coefficient)
        assert degrees_of_freedom(m) == 0
        for n, eff in m.fs.unit.effects.items():
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
                pytest.approx(
                    value(
                        eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                    ),
                    rel=1e-6,
                )
                == c0
            )
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 2
            if n != 1:
                assert (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .is_fixed()
                )
                linking_constr = getattr(
                    m.fs.unit, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )
                assert linking_constr.active
                assert eff.effect.overall_heat_transfer_coefficient.is_fixed()
                assert value(eff.effect.overall_heat_transfer_coefficient) == htc
                assert degrees_of_freedom(eff.effect) == 1

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

        total_mass_flow_water_in = 0
        total_mass_flow_salt_in = 0
        total_mass_flow_water_out = 0

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
                        + pyunits.convert(
                            eff.work_mechanical[0], to_units=pyunits.J * pyunits.s**-1
                        )
                    )
                )
                <= 1e-2
            )

            total_mass_flow_water_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            total_mass_flow_salt_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
            )
            total_mass_flow_water_out += value(
                eff.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        # Test control volume mass balance
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "NaCl"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_salt_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_out[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_out
        )

    @pytest.mark.component
    def test_solution(self, MEC4_frame):
        m = MEC4_frame

        assert (
            pytest.approx(value(m.fs.unit.recovery_vol_phase["Liq"]), rel=1e-3)
            == 0.734452
        )

        unit_results_dict = {
            1: {
                "crystallization_yield": {"NaCl": 0.5},
                "product_volumetric_solids_fraction": 0.132962,
                "temperature_operating": 359.48,
                "pressure_operating": 45000.0,
                "dens_mass_magma": 281.21,
                "dens_mass_slurry": 1298.22,
                "work_mechanical": {0.0: 1915.44},
                "diameter_crystallizer": 1.1448,
                "height_slurry": 1.073,
                "height_crystallizer": 1.9316,
                "magma_circulation_flow_vol": 0.134172,
                "relative_supersaturation": {"NaCl": 0.567033},
                "t_res": 1.0228,
                "volume_suspension": 1.1045,
                "eq_max_allowable_velocity": 2.6304,
                "eq_vapor_space_height": 0.858635,
                "eq_minimum_height_diameter_ratio": 1.7172,
                "energy_flow_superheated_vapor": 1706708.64,
                "delta_temperature_in": {0.0: 57.39},
                "delta_temperature_out": {0.0: 123.72},
                "delta_temperature": {0.0: 86.31},
                "heat_exchanger_area": 221.91,
            },
            2: {
                "product_volumetric_solids_fraction": 0.132109,
                "temperature_operating": 344.86,
                "pressure_operating": 25000.0,
                "dens_mass_magma": 279.41,
                "dens_mass_slurry": 1301.84,
                "work_mechanical": {0.0: 1706.708},
                "diameter_crystallizer": 1.2481,
                "height_slurry": 0.828632,
                "height_crystallizer": 1.8721,
                "magma_circulation_flow_vol": 0.119101,
                "relative_supersaturation": {"NaCl": 0.571246},
                "t_res": 1.0228,
                "volume_suspension": 1.0138,
                "eq_max_allowable_velocity": 3.4641,
                "eq_vapor_space_height": 0.936079,
                "eq_minimum_height_diameter_ratio": 1.8721,
                "energy_flow_superheated_vapor": 1569579.6,
                "delta_temperature_in": {0.0: 14.61},
                "delta_temperature_out": {0.0: 58.77},
                "delta_temperature": {0.0: 31.59},
                "heat_exchanger_area": 540.22,
            },
            3: {
                "product_volumetric_solids_fraction": 0.131922,
                "temperature_operating": 340.52,
                "pressure_operating": 20800.0,
                "dens_mass_magma": 279.01,
                "dens_mass_slurry": 1303.05,
                "work_mechanical": {0.0: 1569.57},
                "diameter_crystallizer": 1.2521,
                "height_slurry": 0.763702,
                "height_crystallizer": 1.8782,
                "magma_circulation_flow_vol": 0.109397,
                "relative_supersaturation": {"NaCl": 0.572391},
                "t_res": 1.0228,
                "volume_suspension": 0.940428,
                "eq_max_allowable_velocity": 3.7763,
                "eq_vapor_space_height": 0.939111,
                "eq_minimum_height_diameter_ratio": 1.8782,
                "energy_flow_superheated_vapor": 1457196.7,
                "delta_temperature_in": {0.0: 4.3421},
                "delta_temperature_out": {0.0: 44.83},
                "delta_temperature": {0.0: 16.85},
                "heat_exchanger_area": 931.31,
            },
            4: {
                "product_volumetric_solids_fraction": 0.131442,
                "temperature_operating": 323.22,
                "pressure_operating": 9500.0,
                "dens_mass_magma": 278.0,
                "dens_mass_slurry": 1308.49,
                "work_mechanical": {0.0: 1457.19},
                "diameter_crystallizer": 1.4623,
                "height_slurry": 0.537417,
                "height_crystallizer": 2.1935,
                "magma_circulation_flow_vol": 0.101015,
                "relative_supersaturation": {"NaCl": 0.576466},
                "t_res": 1.0228,
                "volume_suspension": 0.902678,
                "eq_max_allowable_velocity": 5.4596,
                "eq_vapor_space_height": 1.0967,
                "eq_minimum_height_diameter_ratio": 2.1935,
                "energy_flow_superheated_vapor": 1405094.4,
                "delta_temperature_in": {0.0: 17.29},
                "delta_temperature_out": {0.0: 40.69},
                "delta_temperature": {0.0: 27.32},
                "heat_exchanger_area": 533.18,
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
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 0.89795},
            "temperature": 416.8,
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
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
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
            "aggregate_capital_cost": 4934127.42,
            "aggregate_flow_electricity": 8.4612,
            "aggregate_flow_NaCl_recovered": 0.299999,
            "aggregate_flow_steam": 0.4304,
            "aggregate_flow_costs": {
                "electricity": 6095.09,
                "NaCl_recovered": -269822.09,
                "steam": 63793.11,
            },
            "aggregate_direct_capital_cost": 2467063.71,
            "total_capital_cost": 4934127.42,
            "total_operating_cost": -51910.07,
            "maintenance_labor_chemical_operating_cost": 148023.82,
            "total_fixed_operating_cost": 148023.82,
            "total_variable_operating_cost": -199933.89,
            "total_annualized_cost": 500494.84,
            "LCOW": 4.3957,
            "SEC": 0.65144,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4934127.42,
            "cost_factor": 2.0,
            "direct_capital_cost": 2467063.71,
            "capital_cost_crystallizer_effect_1": 701221.05,
            "capital_cost_heat_exchanger_effect_1": 234213.59,
            "capital_cost_effect_1": 467717.32,
            "capital_cost_crystallizer_effect_2": 667484.84,
            "capital_cost_heat_exchanger_effect_2": 558948.69,
            "capital_cost_effect_2": 613216.76,
            "capital_cost_crystallizer_effect_3": 640779.56,
            "capital_cost_heat_exchanger_effect_3": 954745.08,
            "capital_cost_effect_3": 797762.32,
            "capital_cost_crystallizer_effect_4": 624930.83,
            "capital_cost_heat_exchanger_effect_4": 551803.75,
            "capital_cost_effect_4": 588367.29,
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

        self.test_calculate_scaling(m)

        m.fs.unit.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
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
            "aggregate_capital_cost": 5078402.93,
            "aggregate_flow_electricity": 8.4612,
            "aggregate_flow_NaCl_recovered": 0.299999,
            "aggregate_flow_steam": 0.430492,
            "aggregate_flow_costs": {
                "electricity": 6095.09,
                "NaCl_recovered": -269822.09,
                "steam": 63793.11,
            },
            "aggregate_direct_capital_cost": 2539201.46,
            "total_capital_cost": 5078402.93,
            "total_operating_cost": -47581.81,
            "maintenance_labor_chemical_operating_cost": 152352.08,
            "total_fixed_operating_cost": 152352.08,
            "total_variable_operating_cost": -199933.89,
            "total_annualized_cost": 520975.61,
            "LCOW": 4.5756,
            "SEC": 0.65144,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 5078402.93,
            "cost_factor": 2.0,
            "direct_capital_cost": 2539201.46,
            "capital_cost_crystallizer_effect_1": 715324.02,
            "capital_cost_heat_exchanger_effect_1": 234213.59,
            "capital_cost_effect_1": 474768.8,
            "capital_cost_crystallizer_effect_2": 697956.34,
            "capital_cost_heat_exchanger_effect_2": 558948.69,
            "capital_cost_effect_2": 628452.52,
            "capital_cost_crystallizer_effect_3": 676895.98,
            "capital_cost_heat_exchanger_effect_3": 954745.08,
            "capital_cost_effect_3": 815820.53,
            "capital_cost_crystallizer_effect_4": 688515.45,
            "capital_cost_heat_exchanger_effect_4": 551803.75,
            "capital_cost_effect_4": 620159.6,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing_work_as_steam(self):

        m = build_mec4()

        for _, eff in m.fs.unit.effects.items():
            eff.effect.crystal_median_length.fix(0.6e-3)
            eff.effect.crystal_growth_rate.fix(5e-9)

        self.test_calculate_scaling(m)

        m.fs.unit.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_work_as": "heat"},
        )

        m.fs.costing.nacl_recovered.cost.set_value(-0.024)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.total_flow_vol_in)
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.total_flow_vol_in, name="SEC"
        )
        results = solver.solve(m)

        assert_optimal_termination(results)

        assert not hasattr(m.fs.costing, "aggregate_flow_steam")

        sys_costing_dict = {
            "aggregate_capital_cost": 4934127.42,
            "aggregate_flow_electricity": 8.4612,
            "aggregate_flow_NaCl_recovered": 0.299999,
            "aggregate_flow_heat": 1915.44,
            "aggregate_flow_costs": {
                "electricity": 6095.09,
                "NaCl_recovered": -269822.09,
                "heat": 148989.56,
            },
            "aggregate_direct_capital_cost": 2467063.71,
            "total_capital_cost": 4934127.42,
            "total_operating_cost": 33286.38,
            "maintenance_labor_chemical_operating_cost": 148023.82,
            "total_fixed_operating_cost": 148023.82,
            "total_variable_operating_cost": -114737.44,
            "total_annualized_cost": 585691.30,
            "LCOW": 5.1441,
            "SEC": 0.65144,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 4934127.42,
            "cost_factor": 2.0,
            "direct_capital_cost": 2467063.71,
            "capital_cost_crystallizer_effect_1": 701221.05,
            "capital_cost_heat_exchanger_effect_1": 234213.59,
            "capital_cost_effect_1": 467717.32,
            "capital_cost_crystallizer_effect_2": 667484.84,
            "capital_cost_heat_exchanger_effect_2": 558948.69,
            "capital_cost_effect_2": 613216.76,
            "capital_cost_crystallizer_effect_3": 640779.56,
            "capital_cost_heat_exchanger_effect_3": 954745.08,
            "capital_cost_effect_3": 797762.32,
            "capital_cost_crystallizer_effect_4": 624930.83,
            "capital_cost_heat_exchanger_effect_4": 551803.75,
            "capital_cost_effect_4": 588367.2,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

    @pytest.mark.component
    def test_optimization(self, MEC4_frame):
        m = MEC4_frame

        # Optimization scenario
        # Looking for the best operating pressures within a typical range

        for n, eff in m.fs.unit.effects.items():
            eff.effect.pressure_operating.unfix()
            if n == m.fs.unit.first_effect:
                eff.effect.pressure_operating.setub(1.2e5)
            if n == m.fs.unit.last_effect:
                eff.effect.pressure_operating.setlb(0.02e5)

        eff_idx = [n for n in m.fs.unit.Effects if n != m.fs.unit.first_effect]

        @m.fs.unit.Constraint(eff_idx, doc="Pressure decreasing across effects")
        def eq_pressure_decreasing_across_eff(b, n):
            return (
                b.effects[n].effect.pressure_operating
                <= b.effects[n - 1].effect.pressure_operating
            )

        m.fs.unit.temperature_diff_typical = Var(
            initialize=12,
            bounds=(0, None),
            units=pyunits.degK,
            doc="Typical temperature difference limit between effects in industry",
        )

        m.fs.unit.temperature_diff_typical.fix()

        @m.fs.unit.Constraint(
            eff_idx,
            doc="Temperature difference limit (based on industrial convention)",
        )
        def eq_temperature_difference_limit(b, n):
            return (
                b.effects[n].effect.temperature_operating
                >= b.effects[n - 1].effect.temperature_operating
                - b.temperature_diff_typical
            )

        m.fs.objective = Objective(expr=m.fs.costing.LCOW)

        opt_results = solver.solve(m)
        assert_optimal_termination(opt_results)

    @pytest.mark.component
    def test_optimization_solution(self, MEC4_frame):
        m = MEC4_frame

        # Check the optimized LCOW
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 3.865

        assert (
            pytest.approx(value(m.fs.unit.recovery_vol_phase["Liq"]), rel=1e-3)
            == 0.7544
        )

        unit_results_dict = {
            1: {
                "product_volumetric_solids_fraction": 0.135236,
                "temperature_operating": 387.01,
                "pressure_operating": 120000.0,
                "work_mechanical": {0.0: 2163.65},
                "heat_exchanger_area": 329.25,
            },
            2: {
                "product_volumetric_solids_fraction": 0.134156,
                "temperature_operating": 375.01,
                "pressure_operating": 79812.1,
                "work_mechanical": {0.0: 1818.52},
                "heat_exchanger_area": 495.94,
            },
            3: {
                "product_volumetric_solids_fraction": 0.13323,
                "temperature_operating": 363.01,
                "pressure_operating": 51490.6,
                "work_mechanical": {0.0: 1567.57},
                "heat_exchanger_area": 467.40,
            },
            4: {
                "product_volumetric_solids_fraction": 0.132473,
                "temperature_operating": 351.01,
                "pressure_operating": 32193.7,
                "work_mechanical": {0.0: 1386.22},
                "heat_exchanger_area": 458.50,
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
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 1.014},
            "temperature": 416.8,
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
    def test_optimization_conservation(self, MEC4_frame):

        m = MEC4_frame

        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]

        total_mass_flow_water_in = 0
        total_mass_flow_salt_in = 0
        total_mass_flow_water_out = 0

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
                        + pyunits.convert(
                            eff.work_mechanical[0], to_units=pyunits.J * pyunits.s**-1
                        )
                    )
                )
                <= 1e-2
            )

            total_mass_flow_water_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            total_mass_flow_salt_in += value(
                eff.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
            )
            total_mass_flow_water_out += value(
                eff.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        # Test control volume mass balance
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_in[0]
                .flow_mass_phase_comp["Liq", "NaCl"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_salt_in
        )
        assert (
            pytest.approx(
                m.fs.unit.control_volume.properties_out[0]
                .flow_mass_phase_comp["Liq", "H2O"]
                .value,
                rel=1e-6,
            )
            == total_mass_flow_water_out
        )
