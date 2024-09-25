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
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)

from watertap_contrib.reflo.unit_models.chemical_softening import (
    ChemicalSofteningZO,
)
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


class TestChemSoft_ExcessLimeSodaSilicaRemoval:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):
        component_list = ["Ca_2+", "Mg_2+", "SiO2", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
        )

        m.fs.soft = soft = ChemicalSofteningZO(
            property_package=m.fs.properties,
            silica_removal=True,
            softening_procedure_type="excess_lime_soda",
        )

        q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3

        inlet_conc = {
            "Ca_2+": 1.43,
            "Mg_2+": 0.1814,
            "SiO2": 0.054,
            "Alkalinity_2-": 0.421,
        }

        prop_in = soft.properties_in[0]
        prop_in.temperature.fix()
        prop_in.pressure.fix()
        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
            pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
        )
        for solute, conc in inlet_conc.items():
            mass_flow_solute = pyunits.convert(
                q_in * conc * pyunits.kg / pyunits.m**3,
                to_units=pyunits.kg / pyunits.s,
            )
            prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
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
        soft.number_mixers.set_value(1)
        soft.number_floc.set_value(2)

        soft.frac_mass_water_recovery.fix()
        soft.removal_efficiency.fix()
        soft.CO2_CaCO3.fix(0.063)

        soft.ca_eff_target.fix(0.020 / 2.5)
        soft.mg_eff_target.fix(0.010 / 4.12)

        soft.retention_time_mixer.fix(0.4)
        soft.retention_time_floc.fix(25)
        soft.retention_time_sed.fix(130)
        soft.retention_time_recarb.fix(20)

        soft.vel_gradient_mix.fix(300)
        soft.vel_gradient_floc.fix(50)

        return m

    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame

        assert len(m.fs.soft.config) == 6

        assert not m.fs.soft.config.dynamic
        assert not m.fs.soft.config.has_holdup
        assert m.fs.soft.config.property_package is m.fs.properties
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):
        m = chem_soft_frame

        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert number_variables(m) == 95
        assert number_total_constraints(m) == 61
        assert number_unused_variables(m) == 17

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "ca_eff_target": 0.008,
            "mg_eff_target": 0.002427,
            "removal_efficiency": {"SiO2": 0.7, "Alkalinity_2-": 0.7},
            "volume_mixer": 1.0535825,
            "volume_floc": 65.8489,
            "volume_sed": 342.4143,
            "volume_recarb": 52.6791,
            "CaO_dosing": 2740.80,
            "Na2CO3_dosing": 15652.67,
            "CO2_first_basin": 119.23,
            "CO2_second_basin": 6547.40,
            "excess_CaO": 0.0614,
            "sludge_prod": 0.2301,
            "Ca_CaCO3": 3.567,
            "Mg_CaCO3": 0.7458,
            "Ca_hardness_CaCO3": 0.420123,
            "Ca_hardness_nonCaCO3": 3.147,
            "Mg_hardness_CaCO3": 0.0,
            "Mg_hardness_nonCaCO3": 0.7458,
            "MgCl2_dosing": 0.0,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        assert value(m.fs.soft.ca_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
            rel=1e-3,
        )

        assert value(m.fs.soft.mg_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
            rel=1e-3,
        )

    @pytest.mark.component
    def test_costing(self, chem_soft_frame):

        m = chem_soft_frame
        prop_in = m.fs.soft.properties_in[0]
        m.fs.costing = TreatmentCosting()
        m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_in.flow_vol)

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 2183265.99,
            "aggregate_fixed_operating_cost": 699005.06,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 0.259444,
            "aggregate_flow_lime": 1086203.02,
            "aggregate_flow_soda_ash": 15652.67,
            "aggregate_flow_co2": 2434990.33,
            "aggregate_flow_costs": {
                "electricity": 186.89,
                "lime": 108620.30,
                "soda_ash": 1372113.8,
                "co2": 1098808.78,
            },
            "total_capital_cost": 2183265.99,
            "total_operating_cost": 3344232.83,
            "aggregate_direct_capital_cost": 2183265.996,
            "maintenance_labor_chemical_operating_cost": 65497.97,
            "total_fixed_operating_cost": 764503.04,
            "total_variable_operating_cost": 2579729.78,
            "total_annualized_cost": 3588662.44,
            "LCOW": 2.59,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 2183265.99,
            "fixed_operating_cost": 699005.06,
            "mixer_power": 94.82,
            "floc_power": 164.62,
            "electricity_flow": 0.259444,
            "mix_tank_capital_cost": 23244.46,
            "floc_tank_capital_cost": 308513.2,
            "sed_basin_capital_cost": 208749.21,
            "recarb_basin_capital_cost": 35051.86,
            "recarb_basin_source_capital_cost": 1042673.54,
            "lime_feed_system_capital_cost": 162549.35,
            "admin_capital_cost": 93971.14,
            "mixer_op_cost": 1775.1,
            "floc_tank_op_cost": 1210.0,
            "sed_basin_op_cost": 11045.92,
            "recarb_basin_op_cost": 65189.35,
            "lime_feed_op_cost": 250916.77,
            "lime_sludge_mngt_op_cost": 280198.55,
            "admin_op_cost": 88669.34,
            "direct_capital_cost": 2183265.99,
            "cao_dosing": 1086203.02,
            "co2_dosing": 2434990.33,
        }

        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r


class TestChemSoft_SingleStageLime:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):

        component_list = ["Ca_2+", "Mg_2+", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
        )

        m.fs.soft = soft = ChemicalSofteningZO(
            property_package=m.fs.properties,
            silica_removal=False,
            softening_procedure_type="single_stage_lime",
        )

        q_in = 50000 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3

        inlet_conc = {"Ca_2+": 0.072, "Mg_2+": 0.0061, "Alkalinity_2-": 0.195}
        prop_in = soft.properties_in[0]
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
            pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
        )
        for solute, conc in inlet_conc.items():
            mass_flow_solute = pyunits.convert(
                q_in * conc * pyunits.kg / pyunits.m**3,
                to_units=pyunits.kg / pyunits.s,
            )
            prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
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

        CO2_in = 0.072 * 1.612 * pyunits.kg / pyunits.m**3

        soft.number_mixers.set_value(1)
        soft.number_floc.set_value(2)

        soft.ca_eff_target.fix(0.020 / 2.5)
        soft.mg_eff_target.fix(0.010 / 4.12)

        soft.retention_time_mixer.fix(0.4)
        soft.retention_time_floc.fix(25)
        soft.retention_time_sed.fix(130)
        soft.retention_time_recarb.fix(20)
        soft.frac_mass_water_recovery.fix()
        soft.removal_efficiency.fix()
        soft.CO2_CaCO3.fix(CO2_in)
        soft.vel_gradient_mix.fix(300)
        soft.vel_gradient_floc.fix(50)
        soft.excess_CaO.fix(0)
        soft.CO2_second_basin.fix(0)
        soft.Na2CO3_dosing.fix(0)
        soft.MgCl2_dosing.fix(0)

        return m

    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame

        assert len(m.fs.soft.config) == 6

        assert not m.fs.soft.config.dynamic
        assert not m.fs.soft.config.has_holdup
        assert m.fs.soft.config.property_package is m.fs.properties

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):
        m = chem_soft_frame

        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert number_variables(m) == 85
        assert number_total_constraints(m) == 50
        assert number_unused_variables(m) == 18

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "ca_eff_target": 0.008,
            "mg_eff_target": 0.002427,
            "removal_efficiency": {"Alkalinity_2-": 0.7},
            "volume_mixer": 13.89,
            "volume_floc": 868.29,
            "volume_sed": 4515.12,
            "volume_recarb": 694.63,
            "vel_gradient_mix": 300.0,
            "vel_gradient_floc": 50.0,
            "frac_mass_water_recovery": 0.99,
            "CaO_dosing": 9414.37,
            "Na2CO3_dosing": 0.0,
            "CO2_first_basin": 770.12,
            "CO2_second_basin": 0.0,
            "excess_CaO": 0.0,
            "CO2_CaCO3": 0.116064,
            "sludge_prod": 0.240284,
            "Ca_CaCO3": 0.17995,
            "Mg_CaCO3": 0.025125,
            "Ca_hardness_CaCO3": 0.17995,
            "Ca_hardness_nonCaCO3": 0.0,
            "Mg_hardness_CaCO3": 0.014995,
            "Mg_hardness_nonCaCO3": 0.010129,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        assert value(m.fs.soft.ca_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
            rel=1e-3,
        )

        assert value(m.fs.soft.mg_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
            rel=1e-3,
        )

    @pytest.mark.component
    def test_costing(self, chem_soft_frame):

        m = chem_soft_frame

        prop_in = m.fs.soft.properties_in[0]
        m.fs.costing = TreatmentCosting()
        m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_in.flow_vol)

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 3739060.59,
            "aggregate_fixed_operating_cost": 1110019.36,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 3.421,
            "aggregate_flow_lime": 3438600.65,
            "aggregate_flow_co2": 281286.38,
            "aggregate_flow_costs": {
                "electricity": 2464.36,
                "lime": 343860.06,
                "soda_ash": 0.0,
                "mgcl2": 0.0,
                "co2": 126932.72,
            },
            "total_capital_cost": 3739060.59,
            "total_operating_cost": 1695448.33,
            "aggregate_direct_capital_cost": 3739060.59,
            "maintenance_labor_chemical_operating_cost": 112171.81,
            "total_fixed_operating_cost": 1222191.17,
            "total_variable_operating_cost": 473257.15,
            "total_annualized_cost": 2114058.41,
            "LCOW": 0.11572,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 3739060.59,
            "fixed_operating_cost": 1110019.36,
            "mixer_power": 1250.34,
            "floc_power": 2170.73,
            "mix_tank_capital_cost": 153068.11,
            "floc_tank_capital_cost": 501015.63,
            "sed_basin_capital_cost": 1464705.96,
            "recarb_basin_capital_cost": 285616.04,
            "recarb_basin_source_capital_cost": 204411.61,
            "lime_feed_system_capital_cost": 238714.65,
            "admin_capital_cost": 390512.93,
            "mixer_op_cost": 13145.02,
            "floc_tank_op_cost": 12104.99,
            "sed_basin_op_cost": 31916.41,
            "recarb_basin_op_cost": 28664.92,
            "lime_feed_op_cost": 442039.08,
            "lime_sludge_mngt_op_cost": 292551.0,
            "admin_op_cost": 289597.92,
            "direct_capital_cost": 3739060.59,
            "cao_dosing": 3438600.65,
            "mgcl2_dosing": 0.0,
            "co2_dosing": 281286.38,
        }

        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r


class TestChemSoft_ExcessLime:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):

        component_list = ["Ca_2+", "Mg_2+", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
        )

        m.fs.soft = soft = ChemicalSofteningZO(
            property_package=m.fs.properties,
            silica_removal=False,
            softening_procedure_type="excess_lime",
        )

        q_in = 50000 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3

        inlet_conc = {"Ca_2+": 0.072, "Mg_2+": 0.0061, "Alkalinity_2-": 0.195}
        prop_in = soft.properties_in[0]
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
            pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
        )
        for solute, conc in inlet_conc.items():
            mass_flow_solute = pyunits.convert(
                q_in * conc * pyunits.kg / pyunits.m**3,
                to_units=pyunits.kg / pyunits.s,
            )
            prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
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

        CO2_in = 0.072 * pyunits.kg / pyunits.m**3

        soft.number_mixers.set_value(1)
        soft.number_floc.set_value(2)

        soft.MgCl2_dosing.fix(0)
        soft.ca_eff_target.fix(0.0007)
        soft.mg_eff_target.fix(0.00006)

        soft.retention_time_mixer.fix(0.4)
        soft.retention_time_floc.fix(25)
        soft.retention_time_sed.fix(130)
        soft.retention_time_recarb.fix(20)
        soft.frac_mass_water_recovery.fix()
        soft.removal_efficiency.fix()
        soft.CO2_CaCO3.fix(CO2_in)
        soft.vel_gradient_mix.fix(300)
        soft.vel_gradient_floc.fix(50)
        soft.CO2_second_basin.fix(0)
        soft.Na2CO3_dosing.fix(0)

        return m

    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame

        assert len(m.fs.soft.config) == 6

        assert not m.fs.soft.config.dynamic
        assert not m.fs.soft.config.has_holdup
        assert m.fs.soft.config.property_package is m.fs.properties
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):

        m = chem_soft_frame

        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert number_variables(m) == 85
        assert number_total_constraints(m) == 51
        assert number_unused_variables(m) == 18

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "mg_eff_target": 6e-05,
            "removal_efficiency": {"Alkalinity_2-": 0.7},
            "volume_mixer": 13.89,
            "volume_floc": 868.29,
            "volume_sed": 4515.12,
            "volume_recarb": 694.63,
            "CaO_dosing": 8589.25,
            "Na2CO3_dosing": 0.0,
            "CO2_first_basin": 365.31,
            "CO2_second_basin": 0.0,
            "excess_CaO": 0.0146035,
            "CO2_CaCO3": 0.072,
            "sludge_prod": 0.24873,
            "MgCl2_dosing": 0.0,
            "Ca_CaCO3": 0.17995,
            "Mg_CaCO3": 0.025125,
            "Ca_hardness_CaCO3": 0.17995,
            "Ca_hardness_nonCaCO3": 0.0,
            "Mg_hardness_CaCO3": 0.014995,
            "Mg_hardness_nonCaCO3": 0.010129,
            "mg_removal_eff": 0.990261,
            "ca_removal_eff": 0.990374,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        assert value(m.fs.soft.ca_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
            rel=1e-3,
        )

        assert value(m.fs.soft.mg_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
            rel=1e-3,
        )

    @pytest.mark.component
    def test_costing(self, chem_soft_frame):

        m = chem_soft_frame
        prop_in = m.fs.soft.properties_in[0]
        m.fs.costing = TreatmentCosting()
        m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_in.flow_vol)

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 3644336.21,
            "aggregate_fixed_operating_cost": 1095006.35,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 3.421,
            "aggregate_flow_lime": 3403997.79,
            "aggregate_flow_soda_ash": 0.0,
            "aggregate_flow_mgcl2": 0.0,
            "aggregate_flow_co2": 133432.11,
            "aggregate_flow_costs": {
                "electricity": 2464.36,
                "lime": 340399.77,
                "soda_ash": 0.0,
                "mgcl2": 0.0,
                "co2": 60212.30,
            },
            "total_capital_cost": 3644336.21,
            "total_operating_cost": 1607412.90,
            "aggregate_direct_capital_cost": 3644336.21,
            "maintenance_labor_chemical_operating_cost": 109330.08,
            "total_fixed_operating_cost": 1204336.44,
            "total_variable_operating_cost": 403076.45,
            "total_annualized_cost": 2015418.02,
            "LCOW": 0.110328,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 3644336.21,
            "fixed_operating_cost": 1095006.35,
            "mixer_power": 1250.34,
            "floc_power": 2170.73,
            "electricity_flow": 3.421,
            "mix_tank_capital_cost": 153068.11,
            "floc_tank_capital_cost": 501015.63,
            "sed_basin_capital_cost": 1464705.96,
            "recarb_basin_capital_cost": 285616.04,
            "recarb_basin_source_capital_cost": 116409.65,
            "lime_feed_system_capital_cost": 231992.23,
            "admin_capital_cost": 390512.93,
            "mixer_op_cost": 13145.02,
            "floc_tank_op_cost": 12104.99,
            "sed_basin_op_cost": 31916.41,
            "recarb_basin_op_cost": 21580.14,
            "lime_feed_op_cost": 423818.60,
            "lime_sludge_mngt_op_cost": 302843.27,
            "admin_op_cost": 289597.92,
            "cost_factor": 1.0,
            "direct_capital_cost": 3644336.21,
            "cao_dosing": 3403997.79,
            "mgcl2_dosing": 0.0,
            "co2_dosing": 133432.11,
        }
        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r


class TestChemSoft_ExcessLimeSodaSilicaRemoval_2:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):

        component_list = ["Ca_2+", "Mg_2+", "SiO2", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(
            solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
        )

        m.fs.soft = soft = ChemicalSofteningZO(
            property_package=m.fs.properties,
            silica_removal=True,
            softening_procedure_type="excess_lime_soda",
        )

        q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3

        inlet_conc = {
            "Ca_2+": 1.43,
            "Mg_2+": 0.1814,
            "Alkalinity_2-": 0.421,
            "SiO2": 0.094,
        }
        prop_in = soft.properties_in[0]
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
            pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
        )
        for solute, conc in inlet_conc.items():
            mass_flow_solute = pyunits.convert(
                q_in * conc * pyunits.kg / pyunits.m**3,
                to_units=pyunits.kg / pyunits.s,
            )
            prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
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
        CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3

        soft.number_mixers.set_value(1)
        soft.number_floc.set_value(2)

        soft.ca_eff_target.fix(0.020 / 2.5)
        soft.mg_eff_target.fix(0.010 / 4.12)

        soft.retention_time_mixer.fix(0.4)
        soft.retention_time_floc.fix(25)
        soft.retention_time_sed.fix(130)
        soft.retention_time_recarb.fix(20)
        soft.frac_mass_water_recovery.fix()
        soft.removal_efficiency.fix()
        soft.CO2_CaCO3.fix(CO2_in)
        soft.vel_gradient_mix.fix(300)
        soft.vel_gradient_floc.fix(50)

        return m

    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame

        assert len(m.fs.soft.config) == 6

        assert not m.fs.soft.config.dynamic
        assert not m.fs.soft.config.has_holdup
        assert m.fs.soft.config.property_package is m.fs.properties
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):
        m = chem_soft_frame

        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert number_variables(m) == 95
        assert number_total_constraints(m) == 61
        assert number_unused_variables(m) == 17

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "ca_eff_target": 0.008,
            "mg_eff_target": 0.002427,
            "volume_mixer": 1.0536,
            "volume_floc": 65.85,
            "volume_sed": 342.42,
            "volume_recarb": 52.68,
            "frac_mass_water_recovery": 0.99,
            "CaO_dosing": 2842.17,
            "Na2CO3_dosing": 15652.67,
            "CO2_first_basin": 123.02,
            "CO2_second_basin": 6547.4,
            "excess_CaO": 0.0637,
            "CO2_CaCO3": 0.108449,
            "sludge_prod": 0.28457,
            "Ca_CaCO3": 3.5674,
            "Mg_CaCO3": 0.745782,
            "Ca_hardness_CaCO3": 0.420106,
            "Ca_hardness_nonCaCO3": 3.1473,
            "Mg_hardness_CaCO3": 0.0,
            "Mg_hardness_nonCaCO3": 0.745782,
            "mg_removal_eff": 0.986751,
            "ca_removal_eff": 0.99446,
            "MgCl2_dosing": 0.515902,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r
        assert value(m.fs.soft.ca_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
            rel=1e-3,
        )

        assert value(m.fs.soft.mg_eff_target) == pytest.approx(
            value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
            rel=1e-3,
        )

    @pytest.mark.component
    def test_costing(self, chem_soft_frame):

        m = chem_soft_frame
        prop_in = m.fs.soft.properties_in[0]
        m.fs.costing = TreatmentCosting()
        m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_in.flow_vol)

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 2185574.46,
            "aggregate_fixed_operating_cost": 769517.24,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 0.259455,
            "aggregate_flow_lime": 1126377.39,
            "aggregate_flow_soda_ash": 15652.67,
            "aggregate_flow_mgcl2": 2800124.40,
            "aggregate_flow_co2": 2436376.63,
            "aggregate_flow_costs": {
                "electricity": 186.89,
                "lime": 112637.73,
                "soda_ash": 1372113.8,
                "mgcl2": 1629351.15,
                "co2": 1099434.36,
            },
            "total_capital_cost": 2185574.46,
            "total_operating_cost": 5048808.43,
            "aggregate_direct_capital_cost": 2185574.46,
            "maintenance_labor_chemical_operating_cost": 65567.233,
            "total_fixed_operating_cost": 835084.47,
            "total_variable_operating_cost": 4213723.96,
            "total_annualized_cost": 5293496.50,
            "LCOW": 3.8208,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 2185574.46,
            "fixed_operating_cost": 769517.24,
            "mixer_power": 94.82,
            "floc_power": 164.62,
            "electricity_flow": 0.259455,
            "mix_tank_capital_cost": 23245.14,
            "floc_tank_capital_cost": 308513.83,
            "sed_basin_capital_cost": 208755.51,
            "recarb_basin_capital_cost": 35052.99,
            "recarb_basin_source_capital_cost": 1043121.664,
            "lime_feed_system_capital_cost": 164398.26,
            "admin_capital_cost": 93973.21,
            "mixer_op_cost": 1775.16,
            "floc_tank_op_cost": 1210.04,
            "sed_basin_op_cost": 11045.99,
            "recarb_basin_op_cost": 65203.48,
            "lime_feed_op_cost": 255133.73,
            "lime_sludge_mngt_op_cost": 346477.85,
            "admin_op_cost": 88670.96,
            "direct_capital_cost": 2185574.46,
            "cao_dosing": 1126377.39,
            "mgcl2_dosing": 2800124.40,
            "co2_dosing": 2436376.63,
        }

        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r
