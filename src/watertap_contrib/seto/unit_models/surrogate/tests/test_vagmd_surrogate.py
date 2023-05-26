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
from watertap_contrib.seto.unit_models.surrogate import VAGMDSurrogate

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
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
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.util.model_diagnostics.infeasible import *

# Get default solver for testing
solver = get_solver()


class TestVAGMD:
    @pytest.fixture(scope="class")
    def VAGMD_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.seawater_properties = SeawaterParameterBlock()
        m.fs.water_properties = WaterParameterBlock()

        # System specification (Input variables)
        feed_flow_rate = 600
        evap_inlet_temp = 80
        cond_inlet_temp = 25
        feed_temp = 25
        feed_salinity = 100
        initial_batch_volume = 50
        recovery_ratio = 0.5
        module_type = "AS7C1.5L"
        cooling_system_type = "closed"
        cooling_inlet_temp = 25  # Not required when cooling system type is "close"

        # Identify if the final brine salinity is larger than 175.3 g/L for module "AS7C1.5L"
        # If yes, then operational parameters need to be fixed at a certain value, 
        # and coolying circuit is closed to maintain condenser inlet temperature constant
        final_brine_salinity = feed_salinity / (1- recovery_ratio) # g/L
        if module_type == "AS7C1.5L" and final_brine_salinity > 175.3:
            cooling_system_type = "closed"
            feed_flow_rate = 1100
            evap_inlet_temp = 80
            cond_inlet_temp = 25
            high_brine_salinity = True
        else:
            high_brine_salinity = False

        m.fs.vagmd = VAGMDSurrogate(
            property_package_seawater=m.fs.seawater_properties,
            property_package_water=m.fs.water_properties,
            module_type=module_type,
            high_brine_salinity=high_brine_salinity,
            cooling_system_type=cooling_system_type,
        )

        # Specify feed flow state properties
        m.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_flow_rate / 1000 / 3600,
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

        # Identify cooling system type
        if cooling_system_type == "closed":
            m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
        else:  # "open"
            m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

        m.fs.seawater_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.seawater_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "TDS")
        )
        m.fs.water_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_frame):
        m = VAGMD_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, VAGMD_frame):
        m = VAGMD_frame
        initialization_tester(m, unit=m.fs.vagmd, outlvl=idaeslog.DEBUG)

        # assert False

    @pytest.mark.component
    def test_solve(self, VAGMD_frame):
        m = VAGMD_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert m.fs.vagmd.condenser_in_props[
            0
        ].temperature.value - 273.15 == pytest.approx(25, abs=1e-3)
        assert m.fs.vagmd.condenser_out_props[0].temperature.value - 273.15 == pytest.approx(65.6459, abs=1e-3)
        assert m.fs.vagmd.evaporator_out_props[0].temperature.value - 273.15 == pytest.approx(37.4374, abs=1e-3)
        assert m.fs.vagmd.thermal_power.value == pytest.approx(17.2473, abs=1e-3)

        print(pyunits.convert(m.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].value * pyunits.m**3/ pyunits.s, to_units=pyunits.L / pyunits.h))
        print(m.fs.vagmd.avg_feed_props[0].cp_mass_phase["Liq"].value)
        print(m.fs.vagmd.condenser_out_props[0].temperature.value - 273.15)
        print(m.fs.vagmd.avg_feed_props[0].dens_mass_phase["Liq"].value)


        assert False
