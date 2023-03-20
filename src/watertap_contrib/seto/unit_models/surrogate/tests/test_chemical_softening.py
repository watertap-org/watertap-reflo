import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.seto.unit_models.chemical_softening_0D import ChemicalSoftening0D

# from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
# from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.seto.property_models.basic_water_properties import (
    BasicWaterParameterBlock,
)

from watertap_contrib.seto.costing import (
    TreatmentCosting,
    EnergyCosting,
    SETOSystemCosting,
)

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


class TestChemSoft:
    @pytest.fixture(scope="class")
    def chem_soft_frame(self):
        # create model, flowsheet
        component_list = ["Ca_2+","Mg_2+","Alkalinity_2-","TSS"]
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties =  BasicWaterParameterBlock(solute_list=component_list)

        m.fs.soft = soft = ChemicalSoftening0D(
        property_package=m.fs.properties, silica_removal= False ,softening_procedure_type= 'single_stage_lime'
        )

        prop_in = soft.properties_in[0]
        # prop_out = soft.properties_out[0]
        # prop_waste = soft.properties_waste[0]

        # System specifications
        ca_in = 0.075  * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        mg_in = 0.0061 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        alk_in = 0.196 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        TSS_in = 0.20 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
        CO2_in  = 0.072*1.612 * pyunits.kg / pyunits.m**3 
        q_in = 50000 * pyunits.m**3 / pyunits.day  # m3/d

        prop_in.flow_vol.fix(q_in)

        prop_in.conc_mass_comp["Ca_2+"].fix(ca_in)
        prop_in.conc_mass_comp["Mg_2+"].fix(mg_in)
        prop_in.conc_mass_comp["Alkalinity_2-"].fix(alk_in)
        prop_in.conc_mass_comp["TSS"].fix(TSS_in)

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

        m.fs.properties.set_default_scaling(
            "flow_vol" , 1
        )

        for comp in component_list:
            m.fs.properties.set_default_scaling(
                "conc_mass_comp" , 1, index = comp
            )

        return m


    @pytest.mark.unit
    def test_config(self, chem_soft_frame):
        m = chem_soft_frame
        # check Chemical softening config arguments

        assert len(m.fs.soft.config) == 7

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
        assert number_variables(m) == 53
        assert number_total_constraints(m) == 32
        assert number_unused_variables(m) == 4


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
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0        


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

        assert pytest.approx(8499.791999999998, rel = 1e-3) == value(m.fs.soft.CaOH_dosing)
        # assert pytest.approx(0.0, rel = 1e-3) == value(m.fs.soft.Na2CO3_dosing)
        assert pytest.approx(847.0, rel = 1e-3) == value(m.fs.soft.CO2_first_basin)
        # assert pytest.approx(0.0, rel = 1e-3) == value(m.fs.soft.CO2_second_basin)
        # assert pytest.approx(0.0, rel = 1e-3) == value(m.fs.soft.excess_CaOH)
        assert pytest.approx(0.0, rel = 1e-3) == value(m.fs.soft.MgCl2_dosing)
        assert pytest.approx(0.36094, rel = 1e-3) == value(m.fs.soft.sludge_prod)

        assert pytest.approx(13.88, rel = 1e-3) == value(m.fs.soft.volume_mixer)
        assert pytest.approx(1736.11, rel = 1e-3) == value(m.fs.soft.volume_floc)
        assert pytest.approx(4513.89, rel = 1e-3) == value(m.fs.soft.volume_sed)
        assert pytest.approx(694.45, rel = 1e-3) == value(m.fs.soft.volume_recarb)


    @pytest.mark.component
    def test_costing(self, chem_soft_frame):
        m = chem_soft_frame
        prop_in = m.fs.soft.properties_in[0]
        m.fs.costing = TreatmentCosting()
        m.fs.soft.costing = UnitModelCostingBlock(flowsheet_costing_block = m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_in.flow_vol)

        
        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(2536319.311176081, rel = 1e-3) == value(m.fs.soft.costing.capital_cost)
        assert pytest.approx(1238684.9810677157, rel = 1e-3) == value(m.fs.soft.costing.fixed_operating_cost)
        assert pytest.approx(30, rel = 1e-3) == value(m.fs.soft.costing.mixer_power)
        assert pytest.approx(104.166, rel = 1e-3) == value(m.fs.soft.costing.floc_power)