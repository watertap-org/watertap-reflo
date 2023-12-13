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
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.unit_models.air_stripping_0D import (
    AirStripping0D,
    PackingMaterial,
)
from watertap_contrib.reflo.property_models import AirWaterEq

from watertap_contrib.reflo.costing import REFLOCosting
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
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

from watertap.core import ControlVolume0DBlock

# Get default solver for testing
solver = get_solver()


class TestAirStripping0D:
    @pytest.fixture(scope="class")
    def ax_frame1(self):
        target = "TCA"
        props = {
            "solute_list": [target],
            "mw_data": {target: 0.1334},
            "dynamic_viscosity_data": {"Liq": 0.00115, "Vap": 1.75e-5},
            "henry_constant_data": {target: 0.725},  # salinity adjusted
            "standard_enthalpy_change_data": {target: 28.7e3},
            "temperature_boiling_data": {target: 347},
            "molar_volume_data": {target: 9.81e-5},
            "critical_molar_volume_data": {target: 2.94e-4},
            "density_data": {"Liq": 999.15, "Vap": 1.22},
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = AirWaterEq(**props)

        ax_config = {"property_package": m.fs.properties, "target": target}

        m.fs.ax = ax = AirStripping0D(**ax_config)
        prop_in = ax.process_flow.properties_in[0]

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 0.0063345, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 38358.266, index=("Liq", target)
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Vap", target)
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 0.741114, index=("Vap", "Air")
        )

        prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(157.8657)
        prop_in.flow_mass_phase_comp["Liq", target].fix(2.61e-5)
        prop_in.flow_mass_phase_comp["Vap", "Air"].fix(1.34932)
        prop_in.flow_mass_phase_comp["Vap", target].fix(0)  # assume pure air into unit
        prop_in.temperature["Liq"].fix(288)
        prop_in.temperature["Vap"].fix(288)
        prop_in.pressure.fix(101325)

        ax.pressure_drop_gradient.fix(75)
        ax.packing_surf_tension.fix(0.033)
        ax.packing_diam_nominal.fix(0.0889)
        ax.packing_surface_area_total.fix(242)
        ax.packing_factor.fix(33)
        ax.surf_tension_water.fix(0.0735)
        ax.target_reduction_frac[target].set_value(0.97)

        return m

    @pytest.mark.unit
    def test_config(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax

        assert len(ax.config) == 9

        assert not ax.config.dynamic
        assert not ax.config.has_holdup
        assert ax.config.property_package is m.fs.properties
        assert ax.config.material_balance_type == MaterialBalanceType.componentPhase
        assert ax.config.energy_balance_type is EnergyBalanceType.none
        assert ax.config.momentum_balance_type is MomentumBalanceType.pressureTotal
        assert ax.config.target == "TCA"
        assert ax.config.packing_material == PackingMaterial.PVC

    @pytest.mark.unit
    def test_build(self, ax_frame1):
        m = ax_frame1
        ax = m.fs.ax

        # test ports
        port_lst = ["inlet", "outlet"]
        for port_str in port_lst:
            port = getattr(ax, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        assert isinstance(ax.process_flow, ControlVolume0DBlock)
        
        assert isinstance(ax.target_set, Set)
        assert len(ax.target_set) == 1
        assert ax.target_set.at(1) == ax.config.target

        assert isinstance(ax.liq_target_set, Set)
        assert len(ax.liq_target_set) == 1
        assert ax.liq_target_set.at(1) == ("Liq", ax.config.target)

        assert isinstance(ax.phase_target_set, Set)
        assert len(ax.phase_target_set) == 2
        assert ax.phase_target_set.at(1) == ("Liq", ax.config.target)
        assert ax.phase_target_set.at(2) == ("Vap", ax.config.target)

    # # test statistics
    # assert number_variables(m) == 193
    # assert number_total_constraints(m) == 51
    # assert number_unused_variables(m) == 90  # vars from property package parameters

    # @pytest.mark.unit
    # def test_dof(self, LT_MED_frame):

    # @pytest.mark.unit
    # def test_calculate_scaling(self, LT_MED_frame):

    # @pytest.mark.component
    # def test_var_scaling(self, LT_MED_frame):

    # @pytest.mark.component
    # def test_initialize(self, LT_MED_frame):

    # @pytest.mark.component
    # def test_solve(self, LT_MED_frame):
    #     pass

    # @pytest.mark.component
    # def test_mass_balance(self, LT_MED_frame):
    #     pass

    # @pytest.mark.component
    # def test_solution(self, LT_MED_frame):
    #     pass

    # @pytest.mark.component
    # def test_costing(self, LT_MED_frame):
    #     pass
