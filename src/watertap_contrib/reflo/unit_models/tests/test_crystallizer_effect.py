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
    Var,
    Param,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
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
from pyomo.util.check_units import assert_units_consistent
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
    NaClParameterData,
    NaClStateBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock,
    WaterStateBlock,
)

from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


def build_effect():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.unit = eff = CrystallizerEffect(
        property_package=m.fs.properties,
        property_package_vapor=m.fs.vapor_properties,
        standalone=True,
    )

    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.15
    atm_pressure = 101325 * pyunits.Pa
    saturated_steam_pressure_gage = 3 * pyunits.bar
    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage, to_units=pyunits.Pa
    )
    feed_temperature = 273.15 + 20
    crystallizer_yield = 0.5
    operating_pressure_eff = 0.45

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    eff.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    eff.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    eff.inlet.pressure[0].fix(atm_pressure)
    eff.inlet.temperature[0].fix(feed_temperature)

    eff.heating_steam[0].pressure_sat

    eff.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    eff.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    eff.crystallization_yield["NaCl"].fix(crystallizer_yield)
    eff.crystal_growth_rate.fix()
    eff.souders_brown_constant.fix()
    eff.crystal_median_length.fix()
    eff.overall_heat_transfer_coefficient.fix(100)

    eff.pressure_operating.fix(operating_pressure_eff * pyunits.bar)

    return m


class TestCrystallizerEffect:
    @pytest.fixture(scope="class")
    def effect_frame(self):
        m = build_effect()
        return m

    @pytest.mark.unit
    def test_config(self, effect_frame):
        m = effect_frame

        assert len(m.fs.unit.config) == 6

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.property_package_vapor is m.fs.vapor_properties
        assert m.fs.unit.config.standalone
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, effect_frame):
        m = effect_frame

        # test ports and variables
        port_lst = ["inlet", "outlet", "solids", "vapor", "steam", "pure_water"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        state_blks = [
            "properties_in",
            "properties_out",
            "properties_pure_water",
            "properties_solids",
            "properties_vapor",
            "heating_steam",
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

        for b in state_blks:
            assert hasattr(m.fs.unit, b)
            sb = getattr(m.fs.unit, b)
            assert sb[0].temperature.ub == 1000
            if b == "heating_steam":
                assert isinstance(sb, WaterStateBlock)
            else:
                assert isinstance(sb, NaClStateBlock)

        effect_vars = [
            "crystal_median_length",
            "crystal_growth_rate",
            "souders_brown_constant",
            "crystallization_yield",
            "product_volumetric_solids_fraction",
            "temperature_operating",
            "pressure_operating",
            "dens_mass_magma",
            "dens_mass_slurry",
            "work_mechanical",
            "diameter_crystallizer",
            "height_slurry",
            "height_crystallizer",
            "magma_circulation_flow_vol",
            "relative_supersaturation",
            "energy_flow_superheated_vapor",
            "delta_temperature_in",
            "delta_temperature_out",
            "heat_exchanger_area",
            "overall_heat_transfer_coefficient",
        ]

        for ev in effect_vars:
            v = getattr(m.fs.unit, ev)
            assert isinstance(v, Var)

        effect_params = [
            "approach_temperature_heat_exchanger",
            "dimensionless_crystal_length",
            "steam_pressure",
            "efficiency_pump",
            "pump_head_height",
        ]

        for ep in effect_params:
            p = getattr(m.fs.unit, ep)
            assert isinstance(p, Param)

        assert number_variables(m) == 288
        assert number_total_constraints(m) == 120
        assert number_unused_variables(m) == 44

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, effect_frame):
        m = effect_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, effect_frame):
        m = effect_frame

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
    def test_initialize(self, effect_frame):
        m = effect_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solve(self, effect_frame):
        m = effect_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_conservation(self, effect_frame):

        m = effect_frame

        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]

        phase_comp_list = [
            (p, j)
            for j in comp_lst
            for p in phase_lst
            if (p, j) in m.fs.unit.properties_in[0].phase_component_set
        ]
        flow_mass_in = sum(
            m.fs.unit.properties_in[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_out = sum(
            m.fs.unit.properties_out[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_solids = sum(
            m.fs.unit.properties_solids[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_vapor = sum(
            m.fs.unit.properties_vapor[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )

        assert (
            abs(
                value(flow_mass_in - flow_mass_out - flow_mass_solids - flow_mass_vapor)
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    flow_mass_in * m.fs.unit.properties_in[0].enth_mass_phase["Liq"]
                    - flow_mass_out * m.fs.unit.properties_out[0].enth_mass_phase["Liq"]
                    - flow_mass_vapor
                    * m.fs.unit.properties_vapor[0].enth_mass_solvent["Vap"]
                    - flow_mass_solids
                    * m.fs.unit.properties_solids[0].enth_mass_solute["Sol"]
                    - flow_mass_solids
                    * m.fs.unit.properties_solids[0].dh_crystallization_mass_comp[
                        "NaCl"
                    ]
                    + pyunits.convert(
                        m.fs.unit.work_mechanical[0], to_units=pyunits.J * pyunits.s**-1
                    )
                )
            )
            <= 1e-2 * 1e3
        )

    @pytest.mark.unit
    def test_solution(self, effect_frame):
        m = effect_frame

        eff_dict = {
            "product_volumetric_solids_fraction": 0.132975,
            "temperature_operating": 359.4,
            "dens_mass_magma": 281.2434,
            "dens_mass_slurry": 1298.2382,
            "work_mechanical": {0.0: 1704.47},
            "diameter_crystallizer": 1.07996,
            "height_slurry": 1.07307,
            "height_crystallizer": 1.8830,
            "magma_circulation_flow_vol": 0.119395,
            "relative_supersaturation": {"NaCl": 0.56703},
            "t_res": 1.022821,
            "volume_suspension": 0.98296,
            "eq_max_allowable_velocity": 2.6304,
            "eq_vapor_space_height": 0.80997,
            "eq_minimum_height_diameter_ratio": 1.6199,
            "energy_flow_superheated_vapor": 1518733.4,
            "delta_temperature_in": {0.0: 57.39},
            "delta_temperature_out": {0.0: 123.7},
            "delta_temperature": {0.0: 86.3139},
            "heat_exchanger_area": 197.474,
        }

        for v, r in eff_dict.items():
            effv = getattr(m.fs.unit, v)
            if effv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(effv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(effv), rel=1e-3) == r

        steam_dict = {
            "flow_mass_phase_comp": {("Liq", "H2O"): 0.0, ("Vap", "H2O"): 0.799053},
            "temperature": 416.87,
            "pressure": 401325.0,
            "dh_vap_mass": 2133119.0,
            "pressure_sat": 401325.0,
        }

        for v, r in steam_dict.items():
            sv = getattr(m.fs.unit.heating_steam[0], v)
            if sv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, effect_frame):
        m = effect_frame
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
        m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.properties_in[0].flow_vol_phase["Liq"], name="SEC"
        )
        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_dict = {
            "aggregate_capital_cost": 868242.67,
            "aggregate_flow_electricity": 2.1715,
            "aggregate_flow_steam": 0.38307,
            "aggregate_flow_costs": {
                "electricity": 1564.26,
                "NaCl_recovered": -67456.42,
                "steam": 56766.9,
            },
            "aggregate_direct_capital_cost": 434121.33,
            "total_capital_cost": 868242.67,
            "total_operating_cost": 16922.02,
            "maintenance_labor_chemical_operating_cost": 26047.28,
            "total_fixed_operating_cost": 26047.28,
            "total_variable_operating_cost": -9125.25,
            "total_annualized_cost": 114126.95,
            "LCOW": 4.0094,
            "LCOW_component_direct_capex": {"fs.unit": 1.7074},
            "LCOW_component_indirect_capex": {"fs.unit": 1.7074},
            "LCOW_component_fixed_opex": {"fs.unit": 0.91508},
            "LCOW_component_variable_opex": {"fs.unit": -0.320584},
            "LCOW_aggregate_direct_capex": {"CrystallizerEffect": 1.7074},
            "LCOW_aggregate_indirect_capex": {"CrystallizerEffect": 1.7074},
            "LCOW_aggregate_fixed_opex": {"CrystallizerEffect": 0.91508},
            "LCOW_aggregate_variable_opex": {
                "CrystallizerEffect": -0.320584,
                "electricity": 0.054954,
                "NaCl_recovered": -2.3698,
                "steam": 1.9943,
            },
            "SEC": 0.66875,
            "SEC_component": {"fs.unit": 0.66875},
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        eff_costing_dict = {
            "capital_cost": 868242.67,
            "direct_capital_cost": 434121.33,
            "capital_cost_crystallizer": 659171.48,
            "capital_cost_heat_exchanger": 209071.19,
        }

        for v, r in eff_costing_dict.items():
            cv = getattr(m.fs.unit.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        assert (
            pytest.approx(
                value(m.fs.unit.heating_steam[0].flow_vol_phase["Vap"]), rel=1e-3
            )
            == 0.383078
        )

        assert (
            pytest.approx(
                value(m.fs.unit.heating_steam[0].flow_vol_phase["Liq"]), rel=1e-8
            )
            == 0
        )
