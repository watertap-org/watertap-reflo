###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
    Objective,
    Expression,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    MCASStateBlock,
)

from watertap.core.util.initialization import check_dof
from watertap.core.util.infeasible import *

from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import *
from pyomo.util.check_units import assert_units_consistent
import idaes.logger as idaeslog

from watertap_contrib.seto.unit_models.electrocoagulation import (
    Electrocoagulation,
    ElectrodeMaterial,
    ReactorMaterial,
)

from watertap_contrib.seto.property_models.basic_water_properties import (
    BasicWaterParameterBlock,
)
from watertap_contrib.seto.costing import TreatmentCosting
from idaes.core.util.exceptions import ConfigurationError, InitializationError

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestEC_noTDS:
    @pytest.mark.unit
    def test_no_TDS_in_feed(self):
        error_msg = "TDS must be in feed stream"
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BasicWaterParameterBlock(solute_list=["foo", "bar", "baz"])
        with pytest.raises(ConfigurationError, match=error_msg):
            m.fs.ec = Electrocoagulation(property_package=m.fs.properties)


class TestEC:
    @pytest.fixture(scope="class")
    def EC_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BasicWaterParameterBlock(solute_list=["TDS", "A", "B", "C"])
        m.fs.ec = ec = Electrocoagulation(property_package=m.fs.properties)

        prop_in = ec.properties_in[0]
        prop_in.flow_vol.fix(0.05)
        prop_in.conc_mass_comp["TDS"].fix(10)
        prop_in.conc_mass_comp["A"].fix(0.1)
        prop_in.conc_mass_comp["B"].fix(0.025)
        prop_in.conc_mass_comp["C"].fix(1.5)

        ec.electrode_thick.fix(0.001)
        ec.current_density.fix(500)
        ec.electrolysis_time.fix(60)
        ec.metal_loading.fix(0.0002)
        ec.number_electrode_pairs.fix(10)
        ec.electrode_gap.fix(0.02)
        ec.current_efficiency.fix(1.5)
        ec.potential_balance.fix(2.5)

        return m

    @pytest.mark.unit
    def test_config(self, EC_frame):
        m = EC_frame
        ec = m.fs.ec
        assert len(ec.config) == 6
        assert not ec.config.dynamic
        assert not ec.config.has_holdup
        assert isinstance(m.fs.properties, BasicWaterParameterBlock)
        assert ec.config.property_package is m.fs.properties
        assert ec.config.electrode_material is ElectrodeMaterial.aluminum
        assert ec.config.electrode_material is not ElectrodeMaterial.iron
        assert ec.config.reactor_material is ReactorMaterial.pvc
        assert ec.config.reactor_material is not ReactorMaterial.stainless_steel
        assert len(ec.config.property_package.component_set) == 4
        assert "TDS" in ec.config.property_package.component_set

    @pytest.mark.unit
    def test_build(self, EC_frame):
        m = EC_frame
        ec = m.fs.ec

        assert assert_units_consistent(m) is None
        ports = ["inlet", "outlet", "waste"]
        for port in ports:
            p = getattr(ec, port)
            assert len(p.vars) == 2
            assert isinstance(p, Port)

        var_lst = [
            "electrode_width",
            "electrode_height",
            "electrode_thick",
            "electrode_mass",
            "electrode_area_total",
            "electrode_area_per",
            "electrode_volume_per",
            "electrode_gap",
            "electrolysis_time",
            "number_electrode_pairs",
            "number_cells",
            "applied_current",
            "current_efficiency",
            "cell_voltage",
            "potential_balance",
            "reactor_volume",
            "metal_loading",
            "ohmic_resistance",
            "charge_loading_rate",
            "current_density",
        ]

        for v_str in var_lst:
            v = getattr(ec, v_str)
            assert isinstance(v, Var)

        param_lst = [
            "mw_electrode_material",
            "valence_electrode_material",
            "density_electrode_material",
            "metal_dose_to_toc_ratio",
            "current_per_reactor",
            "tds_to_cond_conversion",
            "removal_efficiency",
            "vol_recovery",
        ]

        for p_str in param_lst:
            p = getattr(ec, p_str)
            assert isinstance(p, Param)
            assert p.mutable  # all Params are mutable

        idx = [*ec.removal_efficiency.index_set()]
        assert idx == ec.config.property_package.component_set

        for i in idx:
            if i == "TDS":
                assert value(ec.removal_efficiency[i]) == 1e-3
            else:
                assert value(ec.removal_efficiency[i]) == 0.7

        assert isinstance(ec.conductivity, Expression)

        assert number_variables(m) == 50
        assert number_total_constraints(m) == 37
        assert number_unused_variables(m) == 0

    @pytest.mark.unit
    def test_dof(self, EC_frame):
        m = EC_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, EC_frame):
        m = EC_frame
        m.fs.properties.set_default_scaling("flow_vol", 1e3)
        m.fs.properties.set_default_scaling("conc_mass_comp", 0.1, index=("TDS"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 10, index=("A"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 100, index=("B"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 1, index=("C"))

        calculate_scaling_factors(m)

        unscaled_var_lst = list(unscaled_variables_generator(m))
        assert len(unscaled_var_lst) == 0

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert len(badly_scaled_var_lst) == 0

    @pytest.mark.component
    def test_initialization(self, EC_frame):
        m = EC_frame
        initialization_tester(m, unit=m.fs.ec)

    @pytest.mark.component
    def test_solve(self, EC_frame):
        m = EC_frame
        results = solver.solve(m)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, EC_frame):
        m = EC_frame
        ec = m.fs.ec
        prop_in = ec.properties_in[0]
        prop_out = ec.properties_out[0]
        prop_waste = ec.properties_waste[0]

        assert pytest.approx(19.3685, rel=1e-3) == value(ec.electrode_mass)
        assert pytest.approx(142.9412, rel=1e-3) == value(ec.electrode_area_total)
        assert pytest.approx(71470.6163, rel=1e-3) == value(ec.applied_current)
        assert pytest.approx(102.4999, rel=1e-3) == value(ec.cell_voltage)
        assert pytest.approx(0.001399176, rel=1e-3) == value(ec.ohmic_resistance)
        assert pytest.approx(1429.4123, rel=1e-3) == value(ec.charge_loading_rate)
        # test volumetric flow balance
        assert pytest.approx(0, rel=1e-8) == value(
            prop_in.flow_vol - prop_out.flow_vol - prop_waste.flow_vol
        )
        for j in ec.config.property_package.component_list:
            # test mass balance of all components including H2O
            assert pytest.approx(0, rel=1e-8) == value(
                prop_in.flow_mass_comp[j]
                - prop_out.flow_mass_comp[j]
                - prop_waste.flow_mass_comp[j]
            )


class TestECCosting:
    @pytest.fixture(scope="class")
    def EC_frame_cost(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BasicWaterParameterBlock(solute_list=["TDS", "A", "B", "C"])
        m.fs.costing = TreatmentCosting()
        m.fs.ec = ec = Electrocoagulation(property_package=m.fs.properties)
        ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ec.properties_out[0].flow_vol)
        m.fs.costing.add_specific_energy_consumption(ec.properties_out[0].flow_vol)

        prop_in = ec.properties_in[0]
        prop_in.flow_vol.fix(0.05)
        prop_in.conc_mass_comp["TDS"].fix(10)
        prop_in.conc_mass_comp["A"].fix(0.1)
        prop_in.conc_mass_comp["B"].fix(0.025)
        prop_in.conc_mass_comp["C"].fix(1.5)

        ec.electrode_thick.fix(0.001)
        ec.current_density.fix(500)
        ec.electrolysis_time.fix(60)
        ec.metal_loading.fix(0.0002)
        ec.number_electrode_pairs.fix(10)
        ec.electrode_gap.fix(0.02)
        ec.current_efficiency.fix(1.5)
        ec.potential_balance.fix(2.5)

        m.fs.properties.set_default_scaling("flow_vol", 1e3)
        m.fs.properties.set_default_scaling("conc_mass_comp", 0.1, index=("TDS"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 10, index=("A"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 100, index=("B"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 1, index=("C"))

        calculate_scaling_factors(m)

        return m

    @pytest.mark.unit
    def test_dof(self, EC_frame_cost):
        m = EC_frame_cost
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization(self, EC_frame_cost):
        m = EC_frame_cost
        initialization_tester(m, unit=m.fs.ec)

    @pytest.mark.component
    def test_solve(self, EC_frame_cost):
        m = EC_frame_cost
        results = solver.solve(m)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, EC_frame_cost):
        m = EC_frame_cost
        ecc = m.fs.ec.costing

        assert pytest.approx(6941480.794, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(6355422.10308, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(4.51287, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(41.109641, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )
        assert pytest.approx(7325.7381, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_electricity
        )
        assert pytest.approx(130621.245, rel=1e-3) == value(ecc.capital_cost_reactor)
        assert pytest.approx(22307.3218, rel=1e-3) == value(ecc.capital_cost_electrodes)
        assert pytest.approx(6714120.086, rel=1e-3) == value(
            ecc.capital_cost_power_supply
        )


class TestEC_NMSUdemo:
    @pytest.fixture(scope="class")
    def EC_frame_NMSU_demo(self):
        # Intended to be as faithful a representation of the conditions presented in the model delivered by NMSU for SETO project
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = TreatmentCosting()
        m.fs.properties = BasicWaterParameterBlock(
            solute_list=["TDS", "TOC", "COD", "Turbidity", "Hardness", "SiO2", "TSS"]
        )

        m.fs.ec = ec = Electrocoagulation(property_package=m.fs.properties)
        m.fs.ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        flow = 0.0438126363888
        prop_in = ec.properties_in[0]
        prop_out = ec.properties_out[0]
        prop_in.flow_vol.fix(flow)
        prop_in.conc_mass_comp["TDS"].fix(130)
        prop_in.conc_mass_comp["TOC"].fix(0.105)
        prop_in.conc_mass_comp["COD"].fix(1.625)
        prop_in.conc_mass_comp["Turbidity"].fix(0.116)
        prop_in.conc_mass_comp["Hardness"].fix(12.621)
        prop_in.conc_mass_comp["SiO2"].fix(0.108)
        prop_in.conc_mass_comp["TSS"].fix(0.345)
        ec.removal_efficiency["TOC"].set_value(0.79)
        ec.removal_efficiency["COD"].set_value(0.88)
        ec.removal_efficiency["Turbidity"].set_value(0.98)
        ec.removal_efficiency["TSS"].set_value(0.7)
        ec.removal_efficiency["Hardness"].set_value(0.75)
        ec.removal_efficiency["SiO2"].set_value(0.95)
        m.fs.costing.plant_lifetime.fix(25)
        m.fs.costing.factor_capital_annualization.fix(0.09367)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(prop_out.flow_vol)
        m.fs.costing.add_specific_energy_consumption(prop_out.flow_vol)
        ec.electrode_thick.fix(0.001)
        ec.current_density.fix(300)
        ec.electrolysis_time.fix(50)
        ec.metal_loading.fix(0.000105)
        ec.number_electrode_pairs.fix(10)
        ec.electrode_gap.fix(0.02)
        ec.current_efficiency.fix(1.6638)
        ec.potential_balance.fix(2.2715091)
        m.fs.properties.set_default_scaling("flow_vol", 1e3)
        m.fs.properties.set_default_scaling("conc_mass_comp", 0.1, index=("TDS"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 100, index=("TOC"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 1, index=("COD"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 10, index=("Turbidity"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 10, index=("TSS"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 10, index=("SiO2"))
        m.fs.properties.set_default_scaling("conc_mass_comp", 0.1, index=("Hardness"))
        calculate_scaling_factors(m)
        return m

    @pytest.mark.unit
    def test_dof(self, EC_frame_NMSU_demo):
        m = EC_frame_NMSU_demo
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization(self, EC_frame_NMSU_demo):
        m = EC_frame_NMSU_demo
        initialization_tester(m, unit=m.fs.ec)

    @pytest.mark.component
    def test_solve(self, EC_frame_NMSU_demo):
        m = EC_frame_NMSU_demo
        results = solver.solve(m)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, EC_frame_NMSU_demo):
        m = EC_frame_NMSU_demo
        ec = m.fs.ec
        prop_in = ec.properties_in[0]
        prop_out = ec.properties_out[0]
        prop_waste = ec.properties_waste[0]
        ecc = ec.costing

        assert pytest.approx(3245208.9485, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        # NMSU: Total CAPEX = $2665869.49
        assert pytest.approx(2245694.351, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        # NMSU: Total OPEX = $2642485.56 / yr
        assert pytest.approx(1.86273, rel=1e-3) == value(m.fs.costing.LCOW)
        # NMSU: LCOW = 2.09 $/m3
        assert pytest.approx(1.307352, rel=1e-3) == value(
            m.fs.costing.specific_energy_consumption
        )
        # NMSU: Specific energy consumption = 1.29 kWh/m^3
        assert pytest.approx(204.14076, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_electricity
        )
        # NMSU: Power = 204.2 kW
        assert pytest.approx(79953.287, rel=1e-3) == value(ecc.capital_cost_reactor)
        # NMSU: Reactor cost (PVC) = $26250.47
        assert pytest.approx(7094.12941, rel=1e-3) == value(ecc.capital_cost_electrodes)
        # NMSU: Electrode cost (Al) = $1097.20
        assert pytest.approx(3088966.531, rel=1e-3) == value(
            ecc.capital_cost_power_supply
        )
        # NMSU: Power supply cost = $2569326.82

        assert pytest.approx(0, rel=1e-8) == value(
            prop_in.flow_vol - prop_out.flow_vol - prop_waste.flow_vol
        )
        for j in ec.config.property_package.component_list:
            # test mass balance of all components including H2O
            assert pytest.approx(0, rel=1e-8) == value(
                prop_in.flow_mass_comp[j]
                - prop_out.flow_mass_comp[j]
                - prop_waste.flow_mass_comp[j]
            )
