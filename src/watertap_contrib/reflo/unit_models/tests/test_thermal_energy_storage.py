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
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
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
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.thermal_energy_storage import (
    ThermalEnergyStorage,
)

# Get default solver for testing
solver = get_solver()


class TestThermalEnergyStorage:
    @pytest.fixture(scope="class")
    def tes(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = WaterParameterBlock()
        m.fs.tes = ThermalEnergyStorage(property_package=m.fs.properties)

        # Define model inputs

        m.fs.tes.tes_hx_inlet.temperature.fix(65 + 273.15)
        m.fs.tes.tes_hx_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.1)
        m.fs.tes.tes_hx_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.tes.tes_hx_inlet.pressure.fix(101325)

        m.fs.tes.tes_process_inlet.temperature.fix(55 + 273.15)
        m.fs.tes.tes_process_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.05)
        m.fs.tes.tes_process_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.tes.tes_process_inlet.pressure.fix(101325)

        m.fs.tes.tes_hx_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.tes.tes_hx_outlet.pressure.fix(101325)

        TES_outlet_initial = (40 + 273.15) * pyunits.K
        STEC = 68 * pyunits.kWh / pyunits.m**3
        permeate_flow = 0.041 * pyunits.m**3 / pyunits.h
        heat_load = pyunits.convert(STEC * permeate_flow, to_units=pyunits.MW)

        # Fix outlet vapor flow to be 0
        m.fs.tes.tes_hx_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.tes.tes_hx_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.1)
        m.fs.tes.tes_hx_outlet.pressure.fix()
        m.fs.tes.tes_hx_outlet.temperature.fix(TES_outlet_initial)

        m.fs.tes.tes_process_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
        m.fs.tes.tes_process_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.05)
        m.fs.tes.tes_process_outlet.pressure.fix()
        m.fs.tes.tes_process_outlet.temperature.fix(TES_outlet_initial)

        m.fs.tes.dt.fix()
        m.fs.tes.tes_initial_temperature.fix(TES_outlet_initial)

        m.fs.tes.hours_storage.fix(6)
        m.fs.tes.heat_load.fix()

        m.fs.tes.initialize()

        return m

    @pytest.mark.unit
    def test_config(self, tes):
        m = tes

        assert len(m.fs.tes.config) == 4
        assert not m.fs.tes.config.dynamic
        assert not m.fs.tes.config.has_holdup
        assert m.fs.tes.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, tes):
        m = tes

        # Test ports
        port_list = [
            "tes_hx_inlet",
            "tes_hx_outlet",
            "tes_process_inlet",
            "tes_process_outlet",
        ]
        for port_str in port_list:
            port = getattr(m.fs.tes, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # Test statistics
        assert number_variables(m.fs.tes) == 46
        assert number_unused_variables(m.fs.tes) == 0
        assert number_total_constraints(m.fs.tes) == 26

    @pytest.mark.unit
    def test_dof(self, tes):
        m = tes
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, tes):
        m = tes
        calculate_scaling_factors(m)
        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_tes_model(self, tes):

        m = tes
        m.fs.tes.initialize()
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing(self, tes):
        m = tes

        m.fs.test_flow = 50 * pyunits.Mgallons / pyunits.day

        m.fs.costing = TreatmentCosting()
        m.fs.tes.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.total_investment_factor.fix(1)

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(flow_rate=m.fs.test_flow)

        results = solver.solve(m)
        assert_optimal_termination(results)
