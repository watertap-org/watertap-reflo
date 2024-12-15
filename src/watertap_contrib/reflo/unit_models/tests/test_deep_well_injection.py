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
    Param,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.exceptions import ConfigurationError
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
)

from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.costing.units.deep_well_injection import (
    blm_costing_params_dict,
)

# Get default solver for testing
solver = get_solver()

rho = 1000 * pyunits.kg / pyunits.m**3


def build_dwi_default():
    """
    Build for 2,500 ft deep injection well (default)
    at injection pressure of 5 bar (default)
    """

    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    flow_mgd = 5.08 * pyunits.Mgallons / pyunits.day

    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )
    m.fs.unit = dwi = DeepWellInjection(property_package=m.fs.properties)

    prop = dwi.properties[0]
    prop.temperature.fix()
    prop.pressure.fix()

    prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / mass_flow_solute),
            index=("Liq", solute),
        )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_phase_water),
        index=("Liq", "H2O"),
    )

    return m


def build_dwi_10000():
    """
    Build for 10,000 ft deep injection well (max)
    at injection pressure of 82 bar
    """
    inlet_conc = {
        "Ca_2+": 0.1,
        "Mg_2+": 0.2,
        "SiO2": 0.3,
        "Alkalinity_2-": 0.4,
    }

    flow_mgd = 1 * pyunits.Mgallons / pyunits.day

    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )
    m.fs.unit = dwi = DeepWellInjection(
        property_package=m.fs.properties, injection_well_depth=10000
    )

    dwi.injection_pressure.set_value(82)

    prop = dwi.properties[0]
    prop.temperature.fix()
    prop.pressure.fix()

    prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / mass_flow_solute),
            index=("Liq", solute),
        )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_phase_water),
        index=("Liq", "H2O"),
    )

    return m


@pytest.mark.unit
def test_injection_well_depth_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["Ca_2+"], material_flow_basis=MaterialFlowBasis.mass
    )
    error_msg = (
        f"The injection well depth was specified as 1234567890. "
        "The injection well depth must be 2500, 5000, 7500, or 10000."
    )
    with pytest.raises(ConfigurationError, match=error_msg):

        m.fs.unit = DeepWellInjection(
            property_package=m.fs.properties, injection_well_depth=1234567890
        )


@pytest.mark.component()
def test_smooth_bound_lower():
    """
    Test that smooth_bound will give lower bound on pipe_diameter.
    """

    inlet_conc = {
        "Pb_2+": 200,
    }

    flow_mgd = 0.08 * pyunits.Mgallons / pyunits.day
    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )
    m.fs.unit = dwi = DeepWellInjection(
        property_package=m.fs.properties, injection_well_depth=5000
    )
    prop = dwi.properties[0]
    prop.temperature.fix()
    prop.pressure.fix()

    prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)

    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert pytest.approx(value(m.fs.unit.pipe_diameter), rel=1e-3) == 2


@pytest.mark.component()
def test_smooth_bound_upper():
    """
    Test that smooth_bound will give upper bound on pipe_diameter.
    """

    inlet_conc = {
        "Au_1+": 42,
    }

    flow_mgd = 200 * pyunits.Mgallons / pyunits.day
    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )
    m.fs.unit = dwi = DeepWellInjection(
        property_package=m.fs.properties, injection_well_depth=5000
    )
    prop = dwi.properties[0]
    prop.temperature.fix()
    prop.pressure.fix()

    prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)

    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert pytest.approx(value(m.fs.unit.pipe_diameter), rel=1e-3) == 24


class TestDeepWellInjection_BLMCosting:
    @pytest.fixture(scope="class")
    def dwi_frame(self):
        m = build_dwi_default()
        return m

    @pytest.mark.unit
    def test_config(self, dwi_frame):
        m = dwi_frame
        assert len(m.fs.unit.config) == 5
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.injection_well_depth == 2500
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, dwi_frame):
        m = dwi_frame
        port = getattr(m.fs.unit, "inlet")
        assert isinstance(port, Port)
        assert len(port.vars) == 3

        assert number_variables(m) == 28
        assert number_total_constraints(m) == 7
        assert number_unused_variables(m) == 16

        dwi_params = [
            "pipe_diameter_coeff",
            "pipe_diameter_exponent",
            "injection_pressure",
            "injection_well_depth",
            "monitoring_well_depth",
        ]

        for pr in dwi_params:
            p = getattr(m.fs.unit, pr)
            assert isinstance(p, Param)

        assert (
            value(m.fs.unit.injection_well_depth)
            == m.fs.unit.config.injection_well_depth
        )

    @pytest.mark.unit
    def test_dof(self, dwi_frame):
        m = dwi_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, dwi_frame):
        m = dwi_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, dwi_frame):
        m = dwi_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solve_and_solution(self, dwi_frame):
        m = dwi_frame
        results = solver.solve(m)

        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.unit.pipe_diameter), rel=1e-3) == 12.01

    @pytest.mark.component
    def test_costing(self, dwi_frame):
        m = dwi_frame
        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
        # m.fs.costing.base_currency = pyunits.kUSD_2001 # for comparison to original BLM reference
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol_phase["Liq"])
        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 5775541.91,
            "aggregate_flow_electricity": 111.51,
            "aggregate_flow_costs": {"electricity": 80330.72},
            "total_capital_cost": 5775541.91,
            "total_operating_cost": 253596.98,
            "aggregate_direct_capital_cost": 5775541.91,
            "maintenance_labor_chemical_operating_cost": 173266.25,
            "total_fixed_operating_cost": 173266.25,
            "total_variable_operating_cost": 80330.72,
            "total_annualized_cost": 900203.26,
            "LCOW": 0.127899,
        }

        for v, r in sys_cost_results.items():
            dwiv = getattr(m.fs.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r

        dwi_cost_results = {
            "capital_cost": 5775541.91,
            "drilling_capital_cost": 987081.97,
            "piping_capital_cost": 536954.05,
            "packing_capital_cost": 161437.09,
            "casing_capital_cost": 1133557.72,
            "grouting_capital_cost": 609793.31,
            "monitoring_well_capital_cost": 746214.86,
            "mobilization_capital_cost": 1051412.29,
            "logging_testing_capital_cost": 549090.58,
            "pumping_power_required": 111.51,
            "direct_capital_cost": 5775541.91,
        }
        for v, r in dwi_cost_results.items():
            dwiv = getattr(m.fs.unit.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r


class TestDeepWellInjection_SimpleCosting:

    @pytest.fixture(scope="class")
    def dwi_frame(self):

        m = build_dwi_default()
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)
        return m

    @pytest.mark.unit
    def test_reporting(self, dwi_frame):
        m = dwi_frame
        m.fs.unit.report()
        _ = m.fs.unit._get_stream_table_contents()

    @pytest.mark.component
    def test_costing_as_capex(self, dwi_frame):
        m = dwi_frame

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_method": "as_capex"},
        )
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol_phase["Liq"])

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 3690312.62,
            "aggregate_fixed_operating_cost": -110709.37,
            "total_capital_cost": 3690312.62,
            "aggregate_direct_capital_cost": 3690312.62,
            "maintenance_labor_chemical_operating_cost": 110709.37,
            "total_annualized_cost": 413152.45,
            "LCOW": 0.0587,
        }
        for v, r in sys_cost_results.items():
            dwiv = getattr(m.fs.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r

        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == value(
            m.fs.costing.deep_well_injection.dwi_lcow
        )

    @pytest.mark.component
    def test_costing_as_opex(self, dwi_frame):
        m = dwi_frame

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_method": "as_opex"},
        )
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol_phase["Liq"])

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 0.0,
            "aggregate_variable_operating_cost": 413152.45,
            "total_capital_cost": 0.0,
            "total_operating_cost": 413152.45,
            "aggregate_direct_capital_cost": 0.0,
            "maintenance_labor_chemical_operating_cost": 0.0,
            "total_fixed_operating_cost": 0.0,
            "total_variable_operating_cost": 413152.45,
            "total_annualized_cost": 413152.45,
            "LCOW": 0.0587,
        }

        for v, r in sys_cost_results.items():
            dwiv = getattr(m.fs.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r

        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == value(
            m.fs.costing.deep_well_injection.dwi_lcow
        )


class TestDeepWellInjection_10000ft:

    @pytest.fixture(scope="class")
    def dwi_10000_frame(self):
        m = build_dwi_10000()
        return m

    @pytest.mark.unit
    def test_config(self, dwi_10000_frame):
        m = dwi_10000_frame
        assert len(m.fs.unit.config) == 5
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.injection_well_depth == 10000
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, dwi_10000_frame):
        m = dwi_10000_frame
        port = getattr(m.fs.unit, "inlet")
        assert isinstance(port, Port)
        assert len(port.vars) == 3

        assert number_variables(m) == 28
        assert number_total_constraints(m) == 7
        assert number_unused_variables(m) == 16

        dwi_params = [
            "pipe_diameter_coeff",
            "pipe_diameter_exponent",
            "injection_pressure",
            "injection_well_depth",
            "monitoring_well_depth",
        ]

        for pr in dwi_params:
            p = getattr(m.fs.unit, pr)
            assert isinstance(p, Param)

        assert (
            value(m.fs.unit.injection_well_depth)
            == m.fs.unit.config.injection_well_depth
        )

    @pytest.mark.unit
    def test_dof(self, dwi_10000_frame):
        m = dwi_10000_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, dwi_10000_frame):
        m = dwi_10000_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, dwi_10000_frame):
        m = dwi_10000_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solve_and_solution(self, dwi_10000_frame):
        m = dwi_10000_frame
        results = solver.solve(m)

        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.unit.pipe_diameter), rel=1e-3) == 5.33166

    @pytest.mark.component
    def test_costing(self, dwi_10000_frame):
        m = dwi_10000_frame
        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2021
        # m.fs.costing.base_currency = pyunits.kUSD_2001 # for comparison to original BLM reference
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.properties[0].flow_vol_phase["Liq"], name="SEC"
        )

        results = solver.solve(m)
        assert_optimal_termination(results)

        costing_params_depth_dict = {  # param dict for 10000 ft depth injection well
            "logging_testing": {"intercept": 415.55, "slope": 7.5694},
            "drilling": {"intercept": 1166.5, "slope": 60.407},
            "piping": {"base": 407.38, "exponent": 0.3741},
            "casing": {"intercept": 1228.1, "slope": 65.52},
            "grouting": {"base": 465.73, "exp_coeff": 0.0708},
        }

        assert (
            costing_params_depth_dict
            == blm_costing_params_dict[int(value(m.fs.unit.injection_well_depth))]
        )

        for cv, d in costing_params_depth_dict.items():
            for p, v in d.items():
                cvp = getattr(
                    m.fs.costing.deep_well_injection, f"{cv}_capital_cost_{p}"
                )
                assert value(cvp) == v

        sys_cost_results = {
            "aggregate_capital_cost": 11554864.2,
            "aggregate_flow_electricity": 359.62,
            "aggregate_flow_costs": {"electricity": 259054.26},
            "total_capital_cost": 11554864.2,
            "total_operating_cost": 605700.19,
            "aggregate_direct_capital_cost": 11554864.2,
            "maintenance_labor_chemical_operating_cost": 346645.92,
            "total_fixed_operating_cost": 346645.92,
            "total_variable_operating_cost": 259054.26,
            "total_annualized_cost": 1899335.98,
            "LCOW": 1.3723,
            "SEC": 2.2777,
        }

        for v, r in sys_cost_results.items():
            dwiv = getattr(m.fs.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r

        dwi_cost_results = {
            "capital_cost": 11554864.2,
            "drilling_capital_cost": 2672856.67,
            "piping_capital_cost": 1368124.15,
            "packing_capital_cost": 130562.52,
            "casing_capital_cost": 2832413.96,
            "grouting_capital_cost": 1219772.74,
            "monitoring_well_capital_cost": 746214.86,
            "mobilization_capital_cost": 1766297.67,
            "logging_testing_capital_cost": 818621.61,
            "pumping_power_required": 359.62,
            "direct_capital_cost": 11554864.2,
        }
        for v, r in dwi_cost_results.items():
            dwiv = getattr(m.fs.unit.costing, v)
            if dwiv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(dwiv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(dwiv), rel=1e-3) == r
