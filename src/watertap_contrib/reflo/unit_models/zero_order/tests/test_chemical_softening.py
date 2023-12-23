import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from idaes.core import MaterialFlowBasis

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.unit_models.zero_order.chemical_softening_zo import (
    ChemicalSofteningZO,
)
from watertap_contrib.reflo.costing import TreatmentCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
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

# Get default solver for testing
solver = get_solver()


class TestChemSoft1:
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
        print(m.fs.soft.properties_in[0].is_property_constructed("flow_mol_phase_comp"))

        ca_in = 1.43 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        mg_in = 0.1814 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        sio2_in = 0.054 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        alk_in = 0.421 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
        q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
        rho = 1000 * pyunits.kg / pyunits.m**3
        prop_in = soft.properties_in[0]

        flow_mass_phase_water = pyunits.convert(
            q_in * rho, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_ca = pyunits.convert(
            q_in * ca_in, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_mg = pyunits.convert(
            q_in * mg_in, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_si = pyunits.convert(
            q_in * sio2_in, to_units=pyunits.kg / pyunits.s
        )
        flow_mass_phase_alk = pyunits.convert(
            q_in * alk_in, to_units=pyunits.kg / pyunits.s
        )
        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water())

        prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(value(flow_mass_phase_ca))
        prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(value(flow_mass_phase_mg))
        prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(value(flow_mass_phase_si))
        prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(
            value(flow_mass_phase_alk)
        )
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        soft.ca_eff_target.fix(3)
        soft.mg_eff_target.fix(0.2)

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

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / flow_mass_phase_water),
            index=("Liq", "H2O"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / flow_mass_phase_ca),
            index=("Liq", "Ca_2+"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / flow_mass_phase_mg),
            index=("Liq", "Mg_2+"),
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", value(1 / flow_mass_phase_si), index=("Liq", "SiO2")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / flow_mass_phase_alk),
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

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):
        # Check electrocoagulation model
        m = chem_soft_frame

        # Test ports
        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 101
        assert number_total_constraints(m) == 61
        assert number_unused_variables(m) == 23

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    # @pytest.mark.skip(reason="flow_mol_phase_comp in badly_scaled")
    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "ca_eff_target": 3.0,
            "mg_eff_target": 0.2,
            "removal_efficiency": {"SiO2": 0.7, "Alkalinity_2-": 0.7},
            "retention_time_mixer": 0.4,
            "retention_time_floc": 25.0,
            "retention_time_sed": 130.0,
            "retention_time_recarb": 20.0,
            "sedimentation_overflow": 90.0,
            "no_of_mixer": 1.0,
            "no_of_floc": 2.0,
            "volume_mixer": 1.053582,
            "volume_floc": 131.697,
            "volume_sed": 342.414,
            "volume_recarb": 52.679,
            "vel_gradient_mix": 300.0,
            "vel_gradient_floc": 50.0,
            "frac_vol_recovery": 0.99,
            "CaO_dosing": 3112.844,
            "Na2CO3_dosing": 15652.679,
            "CO2_first_basin": 652.793,
            "CO2_second_basin": 631.26,
            "excess_CaO": 0.191157681,
            "CO2_CaCO3": 0.10844915,
            "sludge_prod": 0.235832903,
            "Ca_CaCO3": 3.567556,
            "Mg_CaCO3": 0.745811937,
            "Ca_hardness_CaCO3": 0.420123454,
            "Ca_hardness_nonCaCO3": 3.147433,
            "Mg_hardness_CaCO3": 0.0,
            "Mg_hardness_nonCaCO3": 0.745811937,
            "MgCl2_dosing": 0.0,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

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
            "aggregate_capital_cost": 1072612.378,
            "aggregate_fixed_operating_cost": 704966.719,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 0.424066958,
            "aggregate_flow_lime": 1401787.793,
            "aggregate_flow_soda ash": 15652.678,
            "aggregate_flow_co2": 469000.661,
            "aggregate_flow_costs": {
                "electricity": 305.476,
                "lime": 284655.559,
                "soda ash": 4412996.0,
                "co2": 211640.285,
            },
            "total_capital_cost": 1072612.378,
            "maintenance_labor_chemical_operating_cost": 32178.371,
            "total_operating_cost": 5646742.412,
            "aggregate_direct_capital_cost": 1072612.378,
            "LCOW": 4.153448,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 1072612.378,
            "fixed_operating_cost": 704966.719,
            "mixer_power": 94.822,
            "floc_power": 329.244,
            "electricity_flow": 0.424066958,
            "mix_tank_capital_cost": 29431.701,
            "floc_tank_capital_cost": 240667.366,
            "sed_basin_capital_cost": 253632.68,
            "recarb_basin_capital_cost": 37274.263,
            "recarb_basin_source_capital_cost": 243330.372,
            "lime_feed_system_capital_cost": 199005.459,
            "admin_capital_cost": 69270.535,
            "mix_tank_op_cost": 22694.68,
            "floc_tank_op_cost": 7510.352,
            "sed_basin_op_cost": 8141.784,
            "recarb_basin_op_cost": 24809.35,
            "lime_feed_op_cost": 266009.828,
            "lime_sludge_mngt_op_cost": 287131.377,
            "admin_op_cost": 88669.344,
            "cost_factor": 1.0,
            "direct_capital_cost": 1072612.378,
            "cao_dosing": 1401787.793,
            "co2_dosing": 469000.661,
        }

        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r


class TestChemSoft2:
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
            softening_procedure_type="single_stage_lime",
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

        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water())
        prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(value(flow_mass_phase_ca))
        prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(value(flow_mass_phase_mg))
        prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(
            value(flow_mass_phase_alk)
        )
        prop_in.temperature.fix()
        prop_in.pressure.fix()

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / flow_mass_phase_water(), index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / flow_mass_phase_ca(), index=("Liq", "Ca_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / flow_mass_phase_mg(), index=("Liq", "Mg_2+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / flow_mass_phase_alk(),
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

    @pytest.mark.unit
    def test_build(self, chem_soft_frame):
        # Check electrocoagulation model
        m = chem_soft_frame

        # Test ports
        port_list = ["inlet", "outlet", "waste"]
        for port_str in port_list:
            port = getattr(m.fs.soft, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 90
        assert number_total_constraints(m) == 49
        assert number_unused_variables(m) == 24

    @pytest.mark.unit
    def test_dof(self, chem_soft_frame):
        m = chem_soft_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    # @pytest.mark.skip(reason="flow_mol_phase_comp in badly_scaled")
    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        soft_results = {
            "ca_eff_target": 0.02,
            "mg_eff_target": 0.01,
            "removal_efficiency": {"Alkalinity_2-": 0.7},
            "retention_time_mixer": 0.4,
            "retention_time_floc": 25.0,
            "retention_time_sed": 130.0,
            "retention_time_recarb": 20.0,
            "sedimentation_overflow": 90.0,
            "no_of_mixer": 1.0,
            "no_of_floc": 2.0,
            "volume_mixer": 13.892,
            "volume_floc": 1736.585,
            "volume_sed": 4515.121,
            "volume_recarb": 694.634,
            "vel_gradient_mix": 300.0,
            "vel_gradient_floc": 50.0,
            "frac_vol_recovery": 0.99,
            "CaO_dosing": 8290.679,
            "Na2CO3_dosing": 0.0,
            "CO2_first_basin": 770.12,
            "CO2_second_basin": 0.0,
            "excess_CaO": 0.0,
            "CO2_CaCO3": 0.116064,
            "sludge_prod": 0.240284259,
            "MgCl2_dosing": 0.0,
            "Ca_CaCO3": 0.179950855,
            "Mg_CaCO3": 0.025125138,
            "Ca_hardness_CaCO3": 0.179950855,
            "Ca_hardness_nonCaCO3": 0.0,
            "Mg_hardness_CaCO3": 0.014995904,
            "Mg_hardness_nonCaCO3": 0.010129233,
        }

        for v, r in soft_results.items():
            softv = getattr(m.fs.soft, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

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
            "utilization_factor": 1.0,
            "electricity_cost": 0.07,
            "electrical_carbon_intensity": 0.475,
            "TPEC": 4.121212,
            "TIC": 2.0,
            "factor_total_investment": 1.0,
            "factor_maintenance_labor_chemical": 0.03,
            "factor_capital_annualization": 0.1,
            "plant_lifetime": 20.0,
            "aggregate_capital_cost": 2529666.844,
            "aggregate_fixed_operating_cost": 1086190.752,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 5.591,
            "aggregate_flow_lime": 8290.679,
            "aggregate_flow_soda ash": 0.0,
            "aggregate_flow_co2": 770.12,
            "aggregate_flow_costs": {
                "electricity": 4028.055,
                "lime": 1683.555,
                "soda ash": 0.0,
                "co2": 347.522,
            },
            "total_capital_cost": 2529666.844,
            "maintenance_labor_chemical_operating_cost": 75890.005,
            "total_operating_cost": 1168139.892,
            "capital_recovery_factor": 0.1,
            "lime_cost": 0.171,
            "soda ash_cost": 0.65,
            "mgcl2_cost": 1.5,
            "co2_cost": 0.38,
            "aggregate_direct_capital_cost": 2529666.844,
            "LCOW": 0.077794309,
        }
        for v, r in sys_cost_results.items():
            softv = getattr(m.fs.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r

        soft_cost_results = {
            "capital_cost": 2529666.844,
            "fixed_operating_cost": 1086190.752,
            "mixer_power": 1250.341,
            "floc_power": 4341.463,
            "electricity_flow": 5.591,
            "mix_tank_capital_cost": 39806.397,
            "floc_tank_capital_cost": 526375.783,
            "sed_basin_capital_cost": 1062899.167,
            "recarb_basin_capital_cost": 203951.722,
            "recarb_basin_source_capital_cost": 200219.357,
            "lime_feed_system_capital_cost": 208549.018,
            "admin_capital_cost": 287865.396,
            "mix_tank_op_cost": 24169.141,
            "floc_tank_op_cost": 24106.761,
            "sed_basin_op_cost": 19102.499,
            "recarb_basin_op_cost": 19670.431,
            "lime_feed_op_cost": 416992.991,
            "lime_sludge_mngt_op_cost": 292551.006,
            "admin_op_cost": 289597.921,
            "cost_factor": 1.0,
            "direct_capital_cost": 2529666.844,
            "cao_dosing": 3028170.694,
            "mgcl2_dosing": 0.0,
            "co2_dosing": 281286.389,
        }

        for v, r in soft_cost_results.items():
            softv = getattr(m.fs.soft.costing, v)
            if softv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(softv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(softv), rel=1e-3) == r
