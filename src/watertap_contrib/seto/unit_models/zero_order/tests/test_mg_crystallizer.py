import pytest
from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
import re
from pyomo.network import Port
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.seto.unit_models import MgCrystallizerZO

from watertap.property_models.cryst_prop_pack import NaClParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
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

class TestMgCrystallizer:
    @pytest.fixture(scope="class")
    def Mg_crystallizer_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.liquid_prop = NaClParameterBlock()

        m.fs.mg_cryst = MgCrystallizerZO(
                        property_package=m.fs.liquid_prop,
                        )
        cryst = m.fs.mg_cryst

        # Input variables
        feed_flow = 5.04 * pyunits.m**3 / pyunits.day
        feed_temperature = 25 # deg C
        feed_salinity = 319 * pyunits.g / pyunits.L
        
        cryst.Mg_conc_inlet.fix(35.49)
        cryst.Ca_conc_inlet.fix(10.60)
        cryst.conversion_factor_mg.fix(0.9)
        cryst.conversion_factor_ca.fix(1)
        cryst.conc_chemical_dose.fix(0.5)

        # Feed flow properties
        cryst.properties_in.calculate_state(
                var_args={
                    ("flow_vol_phase", "Liq"): pyunits.convert(feed_flow, to_units=pyunits.m**3 / pyunits.s),
                    ("conc_mass_phase_comp", ("Liq", "NaCl")): feed_salinity,
                    ("flow_vol_phase", "Vap"): 0,
                    ("flow_vol_phase", "Sol"): 0, 
                    ("temperature", None): feed_temperature + 273.15,
                    # feed flow is at atmospheric pressure
                    ("pressure", None): 101325,
                },
                hold_state=True,        
        ) 

        # Set scaling factors for mass flow rates
        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.liquid_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )
        
        return m
    
    @pytest.mark.unit
    def test_build(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame

        # test ports
        port_lst = ["inlet", "outlet"]
        for port_str in port_lst:
            port = getattr(m.fs.mg_cryst, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 142
        assert number_total_constraints(m) == 33
        assert number_unused_variables(m) == 85  # vars from property package parameters

    @pytest.mark.unit
    def test_calculate_scaling(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame
        initialization_tester(m, unit=m.fs.mg_cryst, outlvl=idaeslog.DEBUG)

    @pytest.mark.unit
    def test_dof(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame
        results = solver.solve(m)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame

        results = m.fs.mg_cryst._get_performance_contents()['vars']

        assert pytest.approx(78.128, rel=1e-3) == results["NaCl concentration in the waste flow (g/L)"]
        assert pytest.approx(605.388, rel=1e-3) == results["Hydromagnesite precipitated (kg/day)"]
        assert pytest.approx(907.358, rel=1e-3) == results["Calcite precipitated (kg/day)"]
        assert pytest.approx(0.880, rel=1e-3) == results["Mg recovery rate"]
        assert pytest.approx(0.973, rel=1e-3) == results["Ca recovery rate"]
        assert pytest.approx(575.663, rel=1e-3) == results["Ca(OH)2 dose (kg/day)"]
        assert pytest.approx(15.539, rel=1e-3) == results["Ca(OH)2 dose (m3/day)"]
        assert pytest.approx(284.951, rel=1e-3) == results["CO2 dose (kg/day)"]
    
    @pytest.mark.component
    def test_costing(self, Mg_crystallizer_frame):
        m = Mg_crystallizer_frame
        mg_cryst = m.fs.mg_cryst

        m.fs.costing = SETOWaterTAPCosting()
        mg_cryst.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.factor_total_investment.fix(1)
        m.fs.costing.factor_maintenance_labor_chemical.fix(0)

        # Do not include metal recovery revenue
        mg_cryst.costing.include_metal_recovery.fix(0)

        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(mg_cryst.properties_in[0].flow_vol_phase["Liq"])
        m.fs.costing.add_LCOW(mg_cryst.properties_in[0].flow_vol_phase["Liq"])

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(20160, rel=1e-3) == value(mg_cryst.costing.capital_cost)
        assert pytest.approx(10.283, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(16912.904, rel=1e-3) == value(m.fs.costing.total_operating_cost)
        assert pytest.approx(20160, rel=1e-3) == value(m.fs.costing.total_capital_cost)

        # Test to include metal recovery revenue
        mg_cryst.costing.include_metal_recovery.fix(1)
        results = solver.solve(m)
        assert_optimal_termination(results)
        assert pytest.approx(20160, rel=1e-3) == value(mg_cryst.costing.capital_cost)
        assert pytest.approx(-121.788, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(-226211.071, rel=1e-3) == value(m.fs.costing.total_operating_cost)
        assert pytest.approx(20160, rel=1e-3) == value(m.fs.costing.total_capital_cost)