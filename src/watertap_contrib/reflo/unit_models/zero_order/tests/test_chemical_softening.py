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

from watertap_contrib.reflo.unit_models.zero_order.chemical_softening_zo import (
    ChemicalSofteningZO,
)
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


# class TestChemSoft_ExcessLimeSodaSilicaRemoval:
#     @pytest.fixture(scope="class")
#     def chem_soft_frame(self):
#         component_list = ["Ca_2+", "Mg_2+", "SiO2", "Alkalinity_2-"]
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(
#             solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
#         )

#         m.fs.soft = soft = ChemicalSofteningZO(
#             property_package=m.fs.properties,
#             silica_removal=True,
#             softening_procedure_type="excess_lime_soda",
#         )

#         q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
#         rho = 1000 * pyunits.kg / pyunits.m**3

#         inlet_conc = {
#             "Ca_2+": 1.43,
#             "Mg_2+": 0.1814,
#             "SiO2": 0.054,
#             "Alkalinity_2-": 0.421,
#         }

#         prop_in = soft.properties_in[0]
#         prop_in.temperature.fix()
#         prop_in.pressure.fix()
#         flow_mass_phase_water = pyunits.convert(
#             q_in * rho, to_units=pyunits.kg / pyunits.s
#         )
#         prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
#             pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
#         )
#         for solute, conc in inlet_conc.items():
#             mass_flow_solute = pyunits.convert(
#                 q_in * conc * pyunits.kg / pyunits.m**3, to_units=pyunits.kg / pyunits.s
#             )
#             prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
#             m.fs.properties.set_default_scaling(
#                 "flow_mass_phase_comp",
#                 value(1 / mass_flow_solute),
#                 index=("Liq", solute),
#             )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_water),
#             index=("Liq", "H2O"),
#         )
#         soft.number_mixers.set_value(1)
#         soft.number_floc.set_value(2)

#         soft.frac_vol_recovery.fix()
#         soft.removal_efficiency.fix()
#         soft.CO2_CaCO3.fix(0.063)

#         soft.ca_eff_target.fix(0.020 / 2.5)
#         soft.mg_eff_target.fix(0.010 / 4.12)

#         soft.retention_time_mixer.fix(0.4)
#         soft.retention_time_floc.fix(25)
#         soft.retention_time_sed.fix(130)
#         soft.retention_time_recarb.fix(20)

#         soft.vel_gradient_mix.fix(300)
#         soft.vel_gradient_floc.fix(50)

#         return m

#     @pytest.mark.unit
#     def test_config(self, chem_soft_frame):
#         m = chem_soft_frame
#         # check Chemical softening config arguments

#         assert len(m.fs.soft.config) == 6

#         assert not m.fs.soft.config.dynamic
#         assert not m.fs.soft.config.has_holdup
#         assert m.fs.soft.config.property_package is m.fs.properties
#         assert_units_consistent(m)

#     @pytest.mark.unit
#     def test_build(self, chem_soft_frame):
#         m = chem_soft_frame

#         # Test ports
#         port_list = ["inlet", "outlet", "waste"]
#         for port_str in port_list:
#             port = getattr(m.fs.soft, port_str)
#             assert isinstance(port, Port)
#             assert len(port.vars) == 3

#         # test statistics
#         assert number_variables(m) == 92
#         assert number_total_constraints(m) == 58
#         assert number_unused_variables(m) == 17

#     @pytest.mark.unit
#     def test_dof(self, chem_soft_frame):
#         m = chem_soft_frame
#         assert degrees_of_freedom(m) == 0

#     @pytest.mark.unit
#     def test_calculate_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         calculate_scaling_factors(m)

#         # check that all variables have scaling factors
#         unscaled_var_list = list(unscaled_variables_generator(m))
#         assert len(unscaled_var_list) == 0

#     @pytest.mark.component
#     def test_initialize(self, chem_soft_frame):
#         m = chem_soft_frame
#         initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

#     # @pytest.mark.skip(reason="flow_mol_phase_comp in badly_scaled")
#     @pytest.mark.component
#     def test_var_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         badly_scaled_var_lst = list(badly_scaled_var_generator(m))
#         assert badly_scaled_var_lst == []

#     @pytest.mark.component
#     def test_solve(self, chem_soft_frame):
#         m = chem_soft_frame
#         results = solver.solve(m)

#         # Check for optimal solution
#         assert_optimal_termination(results)

#     @pytest.mark.component
#     def test_solution(self, chem_soft_frame):
#         m = chem_soft_frame

#         soft_results = {
#             "ca_eff_target": 0.008,
#             "mg_eff_target": 0.002427,
#             "removal_efficiency": {"SiO2": 0.7, "Alkalinity_2-": 0.7},
#             "volume_mixer": 1.0535825,
#             "volume_floc": 65.8489,
#             "volume_sed": 342.4143,
#             "volume_recarb": 52.6791,
#             "CaO_dosing": 9760.27,
#             "Na2CO3_dosing": 15652.67,
#             "CO2_first_basin": 381.86,
#             "CO2_second_basin": 6547.40,
#             "excess_CaO": 0.2188184,
#             "sludge_prod": 0.237047,
#             "Ca_CaCO3": 3.567,
#             "Mg_CaCO3": 0.7458,
#             "Ca_hardness_CaCO3": 0.420123,
#             "Ca_hardness_nonCaCO3": 3.147,
#             "Mg_hardness_CaCO3": 0.0,
#             "Mg_hardness_nonCaCO3": 0.7458,
#             "MgCl2_dosing": 0.0,
#         }

#         for v, r in soft_results.items():
#             softv = getattr(m.fs.soft, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#         assert value(m.fs.soft.ca_eff_target) == pytest.approx(
#             value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
#             rel=1e-3,
#         )

#         assert value(m.fs.soft.mg_eff_target) == pytest.approx(
#             value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
#             rel=1e-3,
#         )

#     @pytest.mark.component
#     def test_costing(self, chem_soft_frame):

#         m = chem_soft_frame
#         prop_in = m.fs.soft.properties_in[0]
#         m.fs.costing = TreatmentCosting()
#         m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
#         m.fs.costing.cost_process()
#         m.fs.costing.add_LCOW(prop_in.flow_vol)

#         results = solver.solve(m)
#         assert_optimal_termination(results)

#         sys_cost_results = {
#             "aggregate_capital_cost": 2292991.76,
#             "aggregate_fixed_operating_cost": 906885.05,
#             "aggregate_variable_operating_cost": 0.0,
#             "aggregate_flow_electricity": 0.259444,
#             "aggregate_flow_lime": 3868083.56,
#             "aggregate_flow_soda_ash": 15652.67,
#             "aggregate_flow_co2": 2530917.25,
#             "aggregate_flow_costs": {
#                 "electricity": 186.89,
#                 "lime": 785476.58,
#                 "soda_ash": 4412996.0,
#                 "co2": 1142096.57,
#             },
#             "total_capital_cost": 2292991.76,
#             "total_operating_cost": 7316430.86,
#             "aggregate_direct_capital_cost": 2292991.76,
#             "maintenance_labor_chemical_operating_cost": 68789.75,
#             "total_fixed_operating_cost": 975674.8,
#             "total_variable_operating_cost": 6340756.05,
#             "total_annualized_cost": 7573144.93,
#             "LCOW": 5.4665,
#         }
#         for v, r in sys_cost_results.items():
#             softv = getattr(m.fs.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#         soft_cost_results = {
#             "capital_cost": 2292991.76,
#             "fixed_operating_cost": 906885.05,
#             "mixer_power": 94.82,
#             "floc_power": 164.62,
#             "electricity_flow": 0.259444,
#             "mix_tank_capital_cost": 23244.46,
#             "floc_tank_capital_cost": 308513.2,
#             "sed_basin_capital_cost": 208749.21,
#             "recarb_basin_capital_cost": 35051.86,
#             "recarb_basin_source_capital_cost": 1073536.45,
#             "lime_feed_system_capital_cost": 241412.21,
#             "admin_capital_cost": 93971.14,
#             "mixer_op_cost": 1775.1,
#             "floc_tank_op_cost": 1210.0,
#             "sed_basin_op_cost": 11045.92,
#             "recarb_basin_op_cost": 66155.3,
#             "lime_feed_op_cost": 449419.57,
#             "lime_sludge_mngt_op_cost": 288609.79,
#             "admin_op_cost": 88669.34,
#             "direct_capital_cost": 2292991.76,
#             "cao_dosing": 3868083.56,
#             "co2_dosing": 2530917.25,
#         }

#         for v, r in soft_cost_results.items():
#             softv = getattr(m.fs.soft.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r


# class TestChemSoft_SingleStageLime:
#     @pytest.fixture(scope="class")
#     def chem_soft_frame(self):
#         # create model, flowsheet

#         component_list = ["Ca_2+", "Mg_2+", "Alkalinity_2-"]
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(
#             solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
#         )

#         m.fs.soft = soft = ChemicalSofteningZO(
#             property_package=m.fs.properties,
#             silica_removal=False,
#             softening_procedure_type="single_stage_lime",
#         )

#         q_in = 50000 * pyunits.m**3 / pyunits.day  # m3/d
#         rho = 1000 * pyunits.kg / pyunits.m**3

#         inlet_conc = {"Ca_2+": 0.072, "Mg_2+": 0.0061, "Alkalinity_2-": 0.195}
#         prop_in = soft.properties_in[0]
#         prop_in.temperature.fix()
#         prop_in.pressure.fix()

#         flow_mass_phase_water = pyunits.convert(
#             q_in * rho, to_units=pyunits.kg / pyunits.s
#         )
#         prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(
#             pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
#         )
#         for solute, conc in inlet_conc.items():
#             mass_flow_solute = pyunits.convert(
#                 q_in * conc * pyunits.kg / pyunits.m**3, to_units=pyunits.kg / pyunits.s
#             )
#             prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
#             m.fs.properties.set_default_scaling(
#                 "flow_mass_phase_comp",
#                 value(1 / mass_flow_solute),
#                 index=("Liq", solute),
#             )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_water),
#             index=("Liq", "H2O"),
#         )

#         CO2_in = 0.072 * 1.612 * pyunits.kg / pyunits.m**3

#         soft.number_mixers.set_value(1)
#         soft.number_floc.set_value(2)

#         soft.ca_eff_target.fix(0.020 / 2.5)
#         soft.mg_eff_target.fix(0.010 / 4.12)

#         soft.retention_time_mixer.fix(0.4)
#         soft.retention_time_floc.fix(25)
#         soft.retention_time_sed.fix(130)
#         soft.retention_time_recarb.fix(20)
#         soft.frac_vol_recovery.fix()
#         soft.removal_efficiency.fix()
#         soft.CO2_CaCO3.fix(CO2_in)
#         soft.vel_gradient_mix.fix(300)
#         soft.vel_gradient_floc.fix(50)
#         soft.excess_CaO.fix(0)
#         soft.CO2_second_basin.fix(0)
#         soft.Na2CO3_dosing.fix(0)
#         soft.MgCl2_dosing.fix(0)

#         return m

#     @pytest.mark.unit
#     def test_config(self, chem_soft_frame):
#         m = chem_soft_frame
#         # check Chemical softening config arguments

#         assert len(m.fs.soft.config) == 6

#         assert not m.fs.soft.config.dynamic
#         assert not m.fs.soft.config.has_holdup
#         assert m.fs.soft.config.property_package is m.fs.properties

#         assert_units_consistent(m)

#     @pytest.mark.unit
#     def test_build(self, chem_soft_frame):
#         m = chem_soft_frame

#         # Test ports
#         port_list = ["inlet", "outlet", "waste"]
#         for port_str in port_list:
#             port = getattr(m.fs.soft, port_str)
#             assert isinstance(port, Port)
#             assert len(port.vars) == 3

#         # test statistics
#         assert number_variables(m) == 82
#         assert number_total_constraints(m) == 47
#         assert number_unused_variables(m) == 18

#     @pytest.mark.unit
#     def test_dof(self, chem_soft_frame):
#         m = chem_soft_frame
#         assert degrees_of_freedom(m) == 0

#     @pytest.mark.unit
#     def test_calculate_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         calculate_scaling_factors(m)

#         # check that all variables have scaling factors
#         unscaled_var_list = list(unscaled_variables_generator(m))
#         assert len(unscaled_var_list) == 0

#     @pytest.mark.component
#     def test_initialize(self, chem_soft_frame):
#         m = chem_soft_frame
#         initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

#     @pytest.mark.component
#     def test_var_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         badly_scaled_var_lst = list(badly_scaled_var_generator(m))
#         assert badly_scaled_var_lst == []

#     @pytest.mark.component
#     def test_solve(self, chem_soft_frame):
#         m = chem_soft_frame
#         results = solver.solve(m)

#         # Check for optimal solution
#         assert_optimal_termination(results)

#     @pytest.mark.component
#     def test_solution(self, chem_soft_frame):
#         m = chem_soft_frame

#         soft_results = {
#             "ca_eff_target": 0.008,
#             "mg_eff_target": 0.002427,
#             "removal_efficiency": {"Alkalinity_2-": 0.7},
#             "volume_mixer": 13.89,
#             "volume_floc": 868.29,
#             "volume_sed": 4515.12,
#             "volume_recarb": 694.63,
#             "vel_gradient_mix": 300.0,
#             "vel_gradient_floc": 50.0,
#             "frac_vol_recovery": 0.99,
#             "CaO_dosing": 8290.67,
#             "Na2CO3_dosing": 0.0,
#             "CO2_first_basin": 770.12,
#             "CO2_second_basin": 0.0,
#             "excess_CaO": 0.0,
#             "CO2_CaCO3": 0.116064,
#             "sludge_prod": 0.240284,
#             "Ca_CaCO3": 0.17995,
#             "Mg_CaCO3": 0.025125,
#             "Ca_hardness_CaCO3": 0.17995,
#             "Ca_hardness_nonCaCO3": 0.0,
#             "Mg_hardness_CaCO3": 0.014995,
#             "Mg_hardness_nonCaCO3": 0.010129,
#         }

#         for v, r in soft_results.items():
#             softv = getattr(m.fs.soft, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#         assert value(m.fs.soft.ca_eff_target) == pytest.approx(
#             value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]),
#             rel=1e-3,
#         )

#         assert value(m.fs.soft.mg_eff_target) == pytest.approx(
#             value(m.fs.soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]),
#             rel=1e-3,
#         )

#     @pytest.mark.component
#     def test_costing(self, chem_soft_frame):

#         m = chem_soft_frame
#         prop_in = m.fs.soft.properties_in[0]
#         m.fs.costing = TreatmentCosting()
#         m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
#         m.fs.costing.cost_process()
#         m.fs.costing.add_LCOW(prop_in.flow_vol)

#         results = solver.solve(m)
#         assert_optimal_termination(results)

#         sys_cost_results = {
#             "aggregate_capital_cost": 3729796.07,
#             "aggregate_fixed_operating_cost": 1084973.26,
#             "aggregate_variable_operating_cost": 0.0,
#             "aggregate_flow_electricity": 3.421,
#             "aggregate_flow_lime": 3028170.69,
#             "aggregate_flow_co2": 281286.38,
#             "aggregate_flow_costs": {
#                 "electricity": 2464.36,
#                 "lime": 614918.76,
#                 "soda_ash": 0.0,
#                 "mgcl2": 0.0,
#                 "co2": 126932.72,
#             },
#             "total_capital_cost": 3729796.07,
#             "total_operating_cost": 1941183.0,
#             "aggregate_direct_capital_cost": 3729796.07,
#             "maintenance_labor_chemical_operating_cost": 111893.88,
#             "total_fixed_operating_cost": 1196867.14,
#             "total_variable_operating_cost": 744315.86,
#             "total_annualized_cost": 2358755.86,
#             "LCOW": 0.129123,
#         }
#         for v, r in sys_cost_results.items():
#             softv = getattr(m.fs.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#         soft_cost_results = {
#             "capital_cost": 3729796.07,
#             "fixed_operating_cost": 1084973.26,
#             "mixer_power": 1250.34,
#             "floc_power": 2170.73,
#             "mix_tank_capital_cost": 153068.11,
#             "floc_tank_capital_cost": 501015.63,
#             "sed_basin_capital_cost": 1464705.96,
#             "recarb_basin_capital_cost": 285616.04,
#             "recarb_basin_source_capital_cost": 204411.61,
#             "lime_feed_system_capital_cost": 229450.13,
#             "admin_capital_cost": 390512.93,
#             "mixer_op_cost": 13145.02,
#             "floc_tank_op_cost": 12104.99,
#             "sed_basin_op_cost": 31916.41,
#             "recarb_basin_op_cost": 28664.92,
#             "lime_feed_op_cost": 416992.99,
#             "lime_sludge_mngt_op_cost": 292551.0,
#             "admin_op_cost": 289597.92,
#             "direct_capital_cost": 3729796.07,
#             "cao_dosing": 3028170.69,
#             "mgcl2_dosing": 0.0,
#             "co2_dosing": 281286.38,
#         }

#         for v, r in soft_cost_results.items():
#             softv = getattr(m.fs.soft.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r


class TestChemSoft3:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):
        # create model, flowsheet
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

        prop_in = soft.properties_in[0]

        # System specifications
        ca_in = 0.072 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        mg_in = 0.0061 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        alk_in = 0.195 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        CO2_in = 0.072 * 1.612 * pyunits.kg / pyunits.m**3
        q_in = 50000 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3

        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_ca = pyunits.convert(
            q_in * ca_in, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_mg = pyunits.convert(
            q_in * mg_in, to_units=pyunits.kg / pyunits.s
        )

        flow_mass_phase_alk = pyunits.convert(
            q_in * alk_in, to_units=pyunits.kg / pyunits.s
        )

        soft.ca_eff_target.fix()
        soft.mg_eff_target.fix()

        soft.no_of_mixer.fix(1)
        soft.no_of_floc.fix(2)
        soft.retention_time_mixer.fix(0.4)
        soft.retention_time_floc.fix(25)
        soft.retention_time_sed.fix(130)
        soft.retention_time_recarb.fix(20)
        soft.frac_vol_recovery.fix()
        soft.removal_efficiency.fix()
        soft.CO2_CaCO3.fix(CO2_in)
        soft.vel_gradient_mix.fix(300)
        soft.vel_gradient_floc.fix(50)
        soft.excess_CaO.fix(0)
        soft.CO2_second_basin.fix(0)
        soft.Na2CO3_dosing.fix(0)
        soft.MgCl2_dosing.fix(0)

        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(value(flow_mass_phase_water))
        prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(value(flow_mass_phase_ca))
        prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(value(flow_mass_phase_mg))
        prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(
            value(flow_mass_phase_alk)
        )
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(flow_mass_phase_water),
            index=("Liq", "H2O"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(flow_mass_phase_ca),
            index=("Liq", "Ca_2+"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(flow_mass_phase_mg),
            index=("Liq", "Mg_2+"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(flow_mass_phase_alk),
            index=("Liq", "Alkalinity_2-"),
        )

        return m

    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame
        # check Chemical softening config arguments

        assert len(m.fs.soft.config) == 6

        assert not m.fs.soft.config.dynamic
        assert not m.fs.soft.config.has_holdup
        assert m.fs.soft.config.property_package is m.fs.properties
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):

        m = chem_soft_frame

        # Test ports
        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 76
        assert number_total_constraints(m) == 40
        assert number_unused_variables(m) == 18


# class TestChemSoft4:
#     @pytest.fixture(scope="class")
#     def chem_soft_frame(self):
#         component_list = ["Ca_2+", "Mg_2+", "SiO2", "Alkalinity_2-"]
#         m = ConcreteModel()
#         m.fs = FlowsheetBlock(dynamic=False)
#         m.fs.properties = MCASParameterBlock(
#             solute_list=component_list, material_flow_basis=MaterialFlowBasis.mass
#         )

#         m.fs.soft = soft = ChemicalSofteningZO(
#             property_package=m.fs.properties,
#             silica_removal=True,
#             softening_procedure_type="excess_lime_soda",
#         )

#         ca_in = 1.43 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
#         mg_in = 0.1814 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
#         sio2_in = 0.094 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
#         alk_in = 0.421 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
#         CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
#         q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
#         rho = 1000 * pyunits.kg / pyunits.m**3
#         prop_in = soft.properties_in[0]

#         flow_mass_phase_water = pyunits.convert(
#             q_in * rho, to_units=pyunits.kg / pyunits.s
#         )
#         flow_mass_phase_ca = pyunits.convert(
#             q_in * ca_in, to_units=pyunits.kg / pyunits.s
#         )
#         flow_mass_phase_mg = pyunits.convert(
#             q_in * mg_in, to_units=pyunits.kg / pyunits.s
#         )
#         flow_mass_phase_si = pyunits.convert(
#             q_in * sio2_in, to_units=pyunits.kg / pyunits.s
#         )
#         flow_mass_phase_alk = pyunits.convert(
#             q_in * alk_in, to_units=pyunits.kg / pyunits.s
#         )
#         prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water())

#         prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(value(flow_mass_phase_ca))
#         prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(value(flow_mass_phase_mg))
#         prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(value(flow_mass_phase_si))
#         prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(
#             value(flow_mass_phase_alk)
#         )
#         prop_in.temperature.fix()
#         prop_in.pressure.fix()

#         soft.ca_eff_target.fix(0.020 / 2.5)
#         soft.mg_eff_target.fix(0.010 / 4.12)

#         soft.no_of_mixer.fix(1)
#         soft.no_of_floc.fix(2)
#         soft.retention_time_mixer.fix(0.4)
#         soft.retention_time_floc.fix(25)
#         soft.retention_time_sed.fix(130)
#         soft.retention_time_recarb.fix(20)
#         soft.frac_vol_recovery.fix()
#         soft.removal_efficiency.fix()
#         soft.CO2_CaCO3.fix(CO2_in)
#         soft.vel_gradient_mix.fix(300)
#         soft.vel_gradient_floc.fix(50)

#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_water),
#             index=("Liq", "H2O"),
#         )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_ca),
#             index=("Liq", "Ca_2+"),
#         )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_mg),
#             index=("Liq", "Mg_2+"),
#         )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp", value(1 / flow_mass_phase_si), index=("Liq", "SiO2")
#         )
#         m.fs.properties.set_default_scaling(
#             "flow_mass_phase_comp",
#             value(1 / flow_mass_phase_alk),
#             index=("Liq", "Alkalinity_2-"),
#         )

#         return m

#     @pytest.mark.unit
#     def test_config(self, chem_soft_frame):
#         m = chem_soft_frame
#         # check Chemical softening config arguments

#         assert len(m.fs.soft.config) == 6

#         assert not m.fs.soft.config.dynamic
#         assert not m.fs.soft.config.has_holdup
#         assert m.fs.soft.config.property_package is m.fs.properties
#         assert_units_consistent(m)

#     @pytest.mark.unit
#     def test_build(self, chem_soft_frame):
#         m = chem_soft_frame

#         # Test ports
#         port_list = ["inlet", "outlet", "waste"]
#         for port_str in port_list:
#             port = getattr(m.fs.soft, port_str)
#             assert isinstance(port, Port)
#             assert len(port.vars) == 3

#         # test statistics
#         assert number_variables(m) == 94
#         assert number_total_constraints(m) == 58
#         assert number_unused_variables(m) == 17

#     @pytest.mark.unit
#     def test_dof(self, chem_soft_frame):
#         m = chem_soft_frame
#         assert degrees_of_freedom(m) == 0

#     @pytest.mark.unit
#     def test_calculate_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         calculate_scaling_factors(m)

#         # check that all variables have scaling factors
#         unscaled_var_list = list(unscaled_variables_generator(m))
#         assert len(unscaled_var_list) == 0

#     @pytest.mark.component
#     def test_initialize(self, chem_soft_frame):
#         m = chem_soft_frame
#         initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

#     @pytest.mark.component
#     def test_var_scaling(self, chem_soft_frame):
#         m = chem_soft_frame
#         badly_scaled_var_lst = list(badly_scaled_var_generator(m))
#         assert badly_scaled_var_lst == []

#     @pytest.mark.component
#     def test_solve(self, chem_soft_frame):
#         m = chem_soft_frame
#         results = solver.solve(m)

#         # Check for optimal solution
#         assert_optimal_termination(results)

#     @pytest.mark.component
#     def test_solution(self, chem_soft_frame):
#         m = chem_soft_frame

#         soft_results = {
#             "ca_eff_target": 0.008,
#             "mg_eff_target": 0.002427,
#             "removal_efficiency": {"SiO2": 0.7, "Alkalinity_2-": 0.7},
#             "retention_time_mixer": 0.4,
#             "retention_time_floc": 25.0,
#             "retention_time_sed": 130.0,
#             "retention_time_recarb": 20.0,
#             "sedimentation_overflow": 90.0,
#             "no_of_mixer": 1.0,
#             "no_of_floc": 2.0,
#             "volume_mixer": 1.053582,
#             "volume_floc": 131.697,
#             "volume_sed": 342.414,
#             "volume_recarb": 52.679,
#             "vel_gradient_mix": 300.0,
#             "vel_gradient_floc": 50.0,
#             "frac_vol_recovery": 0.99,
#             "CaO_dosing": 10800.845,
#             "Na2CO3_dosing": 15652.679,
#             "CO2_first_basin": 1123.60,
#             "CO2_second_basin": 6547.40,
#             "excess_CaO": 0.663246,
#             "CO2_CaCO3": 0.10844915,
#             "sludge_prod": 0.2792070,
#             "Ca_CaCO3": 3.567556,
#             "Mg_CaCO3": 0.745811937,
#             "Ca_hardness_CaCO3": 0.420123454,
#             "Ca_hardness_nonCaCO3": 3.147433,
#             "Mg_hardness_CaCO3": 0.0,
#             "Mg_hardness_nonCaCO3": 0.745811937,
#             "MgCl2_dosing": 0.51590,
#         }

#         for v, r in soft_results.items():
#             softv = getattr(m.fs.soft, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#     @pytest.mark.component
#     def test_costing(self, chem_soft_frame):

#         m = chem_soft_frame
#         prop_in = m.fs.soft.properties_in[0]
#         m.fs.costing = TreatmentCosting()
#         m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
#         m.fs.costing.cost_process()
#         m.fs.costing.add_LCOW(prop_in.flow_vol)

#         results = solver.solve(m)
#         assert_optimal_termination(results)

#         sys_cost_results = {
#             "aggregate_capital_cost": 1408330.94,
#             "aggregate_fixed_operating_cost": 986678.53,
#             "aggregate_variable_operating_cost": 0.0,
#             "aggregate_flow_electricity": 0.4240669,
#             "aggregate_flow_lime": 4863882.12,
#             "aggregate_flow_soda_ash": 15652.67,
#             "aggregate_flow_co2": 2801838.57,
#             "aggregate_flow_costs": {
#                 "electricity": 305.476,
#                 "lime": 987689.50,
#                 "soda_ash": 4412996.0,
#                 "co2": 1264351.98,
#             },
#             "total_capital_cost": 1408330.94,
#             "maintenance_labor_chemical_operating_cost": 42249.92,
#             "total_operating_cost": 8161092.69,
#             "aggregate_direct_capital_cost": 1408330.94,
#             "LCOW": 6.00454,
#         }
#         for v, r in sys_cost_results.items():
#             softv = getattr(m.fs.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r

#         soft_cost_results = {
#             "capital_cost": 1408330.94,
#             "fixed_operating_cost": 986678.53,
#             "mixer_power": 94.822,
#             "floc_power": 329.244,
#             "electricity_flow": 0.424066,
#             "mix_tank_capital_cost": 29431.701,
#             "floc_tank_capital_cost": 240667.366,
#             "sed_basin_capital_cost": 253632.68,
#             "recarb_basin_capital_cost": 37274.263,
#             "recarb_basin_source_capital_cost": 564872.72,
#             "lime_feed_system_capital_cost": 213175.66,
#             "admin_capital_cost": 69270.535,
#             "mix_tank_op_cost": 22694.68,
#             "floc_tank_op_cost": 7510.352,
#             "sed_basin_op_cost": 8141.784,
#             "recarb_basin_op_cost": 48914.58,
#             "lime_feed_op_cost": 470805.69,
#             "lime_sludge_mngt_op_cost": 339940.35,
#             "admin_op_cost": 88669.34,
#             "direct_capital_cost": 1408330.94,
#             "cao_dosing": 4863882.12,
#             "co2_dosing": 2801838.57,
#         }

#         for v, r in soft_cost_results.items():
#             softv = getattr(m.fs.soft.costing, v)
#             if softv.is_indexed():
#                 for i, s in r.items():
#                     assert pytest.approx(value(softv[i]), rel=1e-3) == s
#             else:
#                 assert pytest.approx(value(softv), rel=1e-3) == r
