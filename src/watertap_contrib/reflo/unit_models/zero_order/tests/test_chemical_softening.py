import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.unit_models.zero_order.chemical_softening_zo import (
    ChemicalSofteningZO,
)

from watertap_contrib.reflo.property_models.basic_water_properties import (
    BasicWaterParameterBlock,
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
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

# Get default solver for testing
solver = get_solver()


class TestChemSoft1:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):
        # create model, flowsheet
        component_list = ["Ca_2+", "Mg_2+", "SiO2", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BasicWaterParameterBlock(solute_list=component_list)

        m.fs.soft = soft = ChemicalSofteningZO(
            property_package=m.fs.properties,
            silica_removal=True,
            softening_procedure_type="excess_lime_soda",
        )

        prop_in = soft.properties_in[0]

        # System specifications
        ca_in = 1.43 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        mg_in = 0.1814 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        sio2_in = 0.054 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        alk_in = 0.421 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
        q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d

        prop_in.flow_vol.fix(q_in)

        prop_in.conc_mass_comp["Ca_2+"].fix(ca_in)
        prop_in.conc_mass_comp["Mg_2+"].fix(mg_in)
        prop_in.conc_mass_comp["SiO2"].fix(sio2_in)
        prop_in.conc_mass_comp["Alkalinity_2-"].fix(alk_in)

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

        m.fs.properties.set_default_scaling("flow_vol", 1)

        for comp in component_list:
            m.fs.properties.set_default_scaling("conc_mass_comp", 1, index=comp)

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
            assert len(port.vars) == 2

        # test statistics
        assert number_variables(m) == 55
        assert number_total_constraints(m) == 35
        assert number_unused_variables(m) == 3

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

        # check that all constraints have been scaled
        # unscaled_constraint_list = list(unscaled_constraints_generator(m))
        # assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        assert pytest.approx(3112.2928969315317, rel=1e-3) == value(
            m.fs.soft.CaO_dosing
        )

        # Components selected do not depend on exit composition
        assert pytest.approx(15652.678552799998, rel=1e-3) == value(
            m.fs.soft.Na2CO3_dosing
        )
        assert pytest.approx(0.1915225725, rel=1e-3) == value(m.fs.soft.excess_CaO)
        assert pytest.approx(0.0, rel=1e-3) == value(m.fs.soft.MgCl2_dosing)
        assert pytest.approx(0.23583140676982073, rel=1e-3) == value(
            m.fs.soft.sludge_prod
        )

        assert pytest.approx(1.051388888888889, rel=1e-3) == value(
            m.fs.soft.volume_mixer
        )
        assert pytest.approx(131.42361111111111, rel=1e-3) == value(
            m.fs.soft.volume_floc
        )
        assert pytest.approx(341.7013888888889, rel=1e-3) == value(m.fs.soft.volume_sed)
        assert pytest.approx(52.56944444444444, rel=1e-3) == value(
            m.fs.soft.volume_recarb
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

        # Capital Costs
        assert pytest.approx(29429.936031519825, rel=1e-3) == value(
            m.fs.soft.costing.mix_tank_capital_cost
        )

        assert pytest.approx(240621.09722222222, rel=1e-3) == value(
            m.fs.soft.costing.floc_tank_capital_cost
        )

        assert pytest.approx(253486.5718140751, rel=1e-3) == value(
            m.fs.soft.costing.sed_basin_capital_cost
        )

        assert pytest.approx(37238.144148918014, rel=1e-3) == value(
            m.fs.soft.costing.recarb_basin_capital_cost
        )

        assert pytest.approx(199004.44248183776, rel=1e-3) == value(
            m.fs.soft.costing.lime_feed_system_capital_cost
        )

        assert pytest.approx(69195.0, rel=1e-3) == value(
            m.fs.soft.costing.admin_capital_cost
        )

        # Operating Costs
        assert pytest.approx(22694.456150427995, rel=1e-3) == value(
            m.fs.soft.costing.mix_tank_op_cost
        )

        assert pytest.approx(7507.309333747968, rel=1e-3) == value(
            m.fs.soft.costing.floc_tank_op_cost
        )

        assert pytest.approx(8139.208862286755, rel=1e-3) == value(
            m.fs.soft.costing.sed_basin_op_cost
        )

        assert pytest.approx(265988.1995885008, rel=1e-3) == value(
            m.fs.soft.costing.lime_feed_op_cost
        )

        assert pytest.approx(88589.0, rel=1e-3) == value(
            m.fs.soft.costing.admin_op_cost
        )

        # Power consumption
        assert pytest.approx(94.62499999999999, rel=1e-3) == value(
            m.fs.soft.costing.mixer_power
        )
        assert pytest.approx(328.5590277777777, rel=1e-3) == value(
            m.fs.soft.costing.floc_power
        )


class TestChemSoft2:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):
        # create model, flowsheet
        component_list = ["Ca_2+", "Mg_2+", "Alkalinity_2-"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BasicWaterParameterBlock(solute_list=component_list)

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

        prop_in.flow_vol.fix(q_in)

        prop_in.conc_mass_comp["Ca_2+"].fix(ca_in)
        prop_in.conc_mass_comp["Mg_2+"].fix(mg_in)
        prop_in.conc_mass_comp["Alkalinity_2-"].fix(alk_in)

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

        m.fs.properties.set_default_scaling("flow_vol", 1)

        for comp in component_list:
            m.fs.properties.set_default_scaling("conc_mass_comp", 1, index=comp)

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
            assert len(port.vars) == 2

        # test statistics
        assert number_variables(m) == 49
        assert number_total_constraints(m) == 27
        assert number_unused_variables(m) == 5

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

        # check that all constraints have been scaled
        # unscaled_constraint_list = list(unscaled_constraints_generator(m))
        # assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, chem_soft_frame):
        m = chem_soft_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, chem_soft_frame):
        m = chem_soft_frame
        initialization_tester(m, unit=m.fs.soft, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, chem_soft_frame):
        m = chem_soft_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, chem_soft_frame):
        m = chem_soft_frame

        assert pytest.approx(8289.792, rel=1e-3) == value(m.fs.soft.CaO_dosing)

        # Components selected do not depend on exit composition
        assert pytest.approx(0.0, rel=1e-3) == value(m.fs.soft.Na2CO3_dosing)
        assert pytest.approx(0.0, rel=1e-3) == value(m.fs.soft.excess_CaO)
        assert pytest.approx(0.0, rel=1e-3) == value(m.fs.soft.MgCl2_dosing)
        assert pytest.approx(0.24028425925925978, rel=1e-3) == value(
            m.fs.soft.sludge_prod
        )

        assert pytest.approx(13.88, rel=1e-3) == value(m.fs.soft.volume_mixer)
        assert pytest.approx(1736.11, rel=1e-3) == value(m.fs.soft.volume_floc)
        assert pytest.approx(4513.89, rel=1e-3) == value(m.fs.soft.volume_sed)
        assert pytest.approx(694.45, rel=1e-3) == value(m.fs.soft.volume_recarb)

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

        # Capital Costs
        assert pytest.approx(39803.3206, rel=1e-3) == value(
            m.fs.soft.costing.mix_tank_capital_cost
        )

        assert pytest.approx(526325.0016, rel=1e-3) == value(
            m.fs.soft.costing.floc_tank_capital_cost
        )

        assert pytest.approx(1062683.65, rel=1e-3) == value(
            m.fs.soft.costing.sed_basin_capital_cost
        )

        assert pytest.approx(203901.925, rel=1e-3) == value(
            m.fs.soft.costing.recarb_basin_capital_cost
        )

        assert pytest.approx(208547.38, rel=1e-3) == value(
            m.fs.soft.costing.lime_feed_system_capital_cost
        )

        assert pytest.approx(287839.2798, rel=1e-3) == value(
            m.fs.soft.costing.admin_capital_cost
        )

        # Operating Costs
        assert pytest.approx(24168.658988826148, rel=1e-3) == value(
            m.fs.soft.costing.mix_tank_op_cost
        )

        assert pytest.approx(24102.201762147426, rel=1e-3) == value(
            m.fs.soft.costing.floc_tank_op_cost
        )

        assert pytest.approx(19100.27097544179, rel=1e-3) == value(
            m.fs.soft.costing.sed_basin_op_cost
        )

        assert pytest.approx(416972.5058115546, rel=1e-3) == value(
            m.fs.soft.costing.lime_feed_op_cost
        )

        assert pytest.approx(289576.0903340209, rel=1e-3) == value(
            m.fs.soft.costing.admin_op_cost
        )

        # Power consumption
        assert pytest.approx(1250, rel=1e-3) == value(m.fs.soft.costing.mixer_power)
        assert pytest.approx(4340.277, rel=1e-3) == value(m.fs.soft.costing.floc_power)
