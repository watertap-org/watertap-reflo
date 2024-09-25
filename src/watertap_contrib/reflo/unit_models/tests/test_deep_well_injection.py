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
from pyomo.util.check_units import assert_units_consistent

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
)

from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection


# Get default solver for testing
solver = get_solver()


def build_dwi():

    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    flow_mgd = 5.08 * pyunits.Mgallons / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3

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


class TestDeepWellInjection_BLMCosting:
    @pytest.fixture(scope="class")
    def dwi_frame(self):
        m = build_dwi()
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
            "tubing_capital_cost": 536954.05,
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

        m = build_dwi()
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)
        return m

    @pytest.mark.component
    def test_simple_costing(self, dwi_frame):
        m = dwi_frame

        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_method": "simple"},
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
