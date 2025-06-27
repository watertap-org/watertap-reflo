#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import re

import pytest
from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    Expression,
    Block,
    Reals,
    NonNegativeReals,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap.costing.watertap_costing_package import WaterTAPCostingBlockData
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    REFLOCosting,
    REFLOCostingData,
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap_contrib.reflo.costing.tests.dummy_costing_units import (
    DummyTreatmentUnit,
    DummyTreatmentNoHeatUnit,
    DummyElectricityUnit,
    DummyHeatUnit,
)

solver = get_solver()


def check_proper_aggregation(treat_cost, energy_cost, system_cost):
    """
    Test if capital and operating costs from treatment and energy
    costing blocks are properly aggregated on the system costing block.
    """

    assert pytest.approx(
        value(treat_cost.total_capital_cost + energy_cost.total_capital_cost), rel=1e-3
    ) == value(system_cost.total_capital_cost)

    assert pytest.approx(
        value(
            treat_cost.total_fixed_operating_cost
            + energy_cost.total_fixed_operating_cost
        ),
        rel=1e-3,
    ) == value(system_cost.total_fixed_operating_cost)

    if hasattr(energy_cost, "aggregate_flow_heat"):
        assert pytest.approx(
            value(treat_cost.aggregate_flow_heat + energy_cost.aggregate_flow_heat),
            rel=1e-3,
        ) == value(system_cost.aggregate_flow_heat)

        assert pytest.approx(value(system_cost.aggregate_flow_heat), rel=1e-3) == value(
            system_cost.aggregate_flow_heat_purchased
        )
    assert pytest.approx(
        value(system_cost.aggregate_flow_electricity), rel=1e-3
    ) == value(system_cost.aggregate_flow_electricity_purchased)

    assert pytest.approx(
        value(
            treat_cost.aggregate_flow_electricity
            + energy_cost.aggregate_flow_electricity
        ),
        rel=1e-3,
    ) == value(system_cost.aggregate_flow_electricity)

    assert pytest.approx(
        value(
            treat_cost.total_operating_cost
            - sum(
                treat_cost.aggregate_flow_costs[f]
                for f in ["electricity", "heat"]
                if f in treat_cost.used_flows
            )
            + energy_cost.total_operating_cost
            - sum(
                energy_cost.aggregate_flow_costs[f]
                for f in ["electricity", "heat"]
                if f in energy_cost.used_flows
            )
            + system_cost.total_electric_operating_cost
            + system_cost.total_heat_operating_cost
        ),
        rel=1e-3,
    ) == value(system_cost.total_operating_cost)


def build_electricity_gen_only_with_heat():
    """
    Test flowsheet with only electricity generation units on energy block.
    The treatment unit consumes both heat and electricity.
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    #### TREATMENT BLOCK
    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()

    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.unit.design_var_a.fix()
    m.fs.treatment.unit.design_var_b.fix()
    m.fs.treatment.unit.electricity_consumption.fix(100)
    m.fs.treatment.unit.heat_consumption.fix()
    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
    )

    #### ENERGY BLOCK
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyElectricityUnit()
    m.fs.energy.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.unit.electricity.fix(10)
    m.fs.energy.costing.cost_process()

    #### SYSTEM COSTING
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOE()

    #### SCALING
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "TDS")
    )
    calculate_scaling_factors(m)

    m.fs.treatment.unit.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.04381,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 35,
            ("temperature", None): 293,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    return m


def build_electricity_gen_only_no_heat():
    """
    Test flowsheet with only electricity generation units on energy block.
    The treatment unit consumes only electricity.
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    #### TREATMENT BLOCK
    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()

    m.fs.treatment.unit = DummyTreatmentNoHeatUnit(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.unit.design_var_a.fix()
    m.fs.treatment.unit.design_var_b.fix()
    m.fs.treatment.unit.electricity_consumption.fix(100)
    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
    )

    #### ENERGY BLOCK
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()

    m.fs.energy.unit = DummyElectricityUnit()
    m.fs.energy.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.unit.electricity.fix(75)
    m.fs.energy.costing.cost_process()

    #### SYSTEM COSTING
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOE()
    m.fs.costing.add_LCOT(m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"])

    #### SCALING
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "TDS")
    )
    calculate_scaling_factors(m)

    set_scaling_factor(m.fs.energy.costing.yearly_electricity_production, 1e-7)
    set_scaling_factor(m.fs.energy.costing.lifetime_electricity_production, 1e-9)

    m.fs.treatment.unit.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.04381,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 35,
            ("temperature", None): 293,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    return m


def build_heat_gen_only():
    """
    Test flowsheet with only heat generation unit on energy block.
    The treatment unit consumes both heat and electricity.
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    #### TREATMENT BLOCK
    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.costing.electricity_cost.fix(0.09)

    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.unit.design_var_a.fix()
    m.fs.treatment.unit.design_var_b.fix()
    m.fs.treatment.unit.electricity_consumption.fix(1000)
    m.fs.treatment.unit.heat_consumption.fix(2500)
    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
    )

    #### ENERGY BLOCK
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyHeatUnit()
    m.fs.energy.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.unit.heat.fix(70)
    m.fs.energy.costing.cost_process()

    #### SYSTEM COSTING
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()
    m.fs.costing.add_specific_thermal_energy_consumption(
        m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOH()
    m.fs.costing.add_LCOT(m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"])

    #### SCALING
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "TDS")
    )
    calculate_scaling_factors(m)

    m.fs.treatment.unit.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.4381,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 20,
            ("temperature", None): 293,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    return m


def build_heat_and_elec_gen():
    """
    Test flowsheet with both heat and electricity generation unit on energy block.
    The heat generating unit also consumes electricity.
    The treatment unit consumes both heat and electricity.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    #### TREATMENT BLOCK
    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.costing.base_currency = pyunits.USD_2023

    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.unit.design_var_a.fix()
    m.fs.treatment.unit.design_var_b.fix()
    m.fs.treatment.unit.electricity_consumption.fix(110)
    m.fs.treatment.unit.heat_consumption.fix(250)
    m.fs.treatment.costing.cost_process()

    #### ENERGY BLOCK
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.costing.base_currency = pyunits.USD_2023

    m.fs.energy.heat_unit = DummyHeatUnit()
    m.fs.energy.elec_unit = DummyElectricityUnit()
    m.fs.energy.heat_unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.elec_unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.heat_unit.heat.fix(50)
    m.fs.energy.elec_unit.electricity.fix(100)
    m.fs.energy.costing.cost_process()

    #### SYSTEM COSTING
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
    )

    #### SCALING
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "TDS")
    )
    calculate_scaling_factors(m)

    m.fs.treatment.unit.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.4381,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 20,
            ("temperature", None): 293,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    return m


def build_default():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyHeatUnit()
    m.fs.energy.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )

    return m


class TestCostingPackagesDefault:
    @pytest.fixture(scope="class")
    def default_build(self):

        m = build_default()

        m.fs.energy.costing.cost_process()
        m.fs.treatment.costing.cost_process()
        m.fs.costing = REFLOSystemCosting()
        m.fs.costing.cost_process()

        return m

    def test_default_build(self, default_build):
        m = default_build

        assert_units_consistent(m)

        # check inheritance
        assert isinstance(m.fs.treatment.costing, REFLOCostingData)
        assert isinstance(m.fs.energy.costing, REFLOCostingData)
        assert isinstance(m.fs.costing, WaterTAPCostingBlockData)

        # no case study loaded by default
        assert m.fs.treatment.costing.config.case_study_definition is None
        assert m.fs.energy.costing.config.case_study_definition is None
        assert not hasattr(m.fs.treatment.costing, "case_study_def")
        assert not hasattr(m.fs.energy.costing, "case_study_def")

        # check unique properties of each
        assert hasattr(m.fs.energy.costing, "has_electricity_generation")
        assert not m.fs.energy.costing.has_electricity_generation  # default is False
        assert hasattr(m.fs.energy.costing, "base_energy_units")
        assert hasattr(m.fs.energy.costing, "plant_lifetime_set")
        assert hasattr(m.fs.energy.costing, "annual_electrical_system_degradation")
        assert hasattr(m.fs.energy.costing, "annual_heat_system_degradation")
        assert not hasattr(m.fs.treatment.costing, "has_electricity_generation")

        assert m.fs.treatment.costing.base_currency is pyunits.USD_2023
        assert m.fs.energy.costing.base_currency is pyunits.USD_2023
        assert m.fs.costing.base_currency is pyunits.USD_2023

        assert m.fs.treatment.costing.base_period is pyunits.year
        assert m.fs.energy.costing.base_period is pyunits.year
        assert m.fs.costing.base_period is pyunits.year

        assert hasattr(m.fs.treatment.costing, "sales_tax_frac")
        assert hasattr(m.fs.energy.costing, "sales_tax_frac")
        assert not hasattr(m.fs.costing, "sales_tax_frac")

        # general domain checks
        assert m.fs.costing.total_heat_operating_cost.domain is Reals
        assert m.fs.costing.total_electric_operating_cost.domain is Reals
        assert m.fs.costing.aggregate_flow_electricity.domain is Reals
        assert m.fs.costing.aggregate_flow_heat.domain is Reals
        assert (
            m.fs.costing.aggregate_flow_electricity_purchased.domain is NonNegativeReals
        )
        assert m.fs.costing.aggregate_flow_heat_purchased.domain is NonNegativeReals

        # capital cost is only positive
        assert m.fs.treatment.unit.costing.capital_cost.domain is NonNegativeReals
        # operating costs can be negative
        assert m.fs.treatment.unit.costing.fixed_operating_cost.domain is Reals
        assert m.fs.treatment.unit.costing.variable_operating_cost.domain is Reals

        # default electricity cost is zero
        assert value(m.fs.costing.electricity_cost) == 0
        assert value(m.fs.treatment.costing.electricity_cost) == 0
        assert value(m.fs.energy.costing.electricity_cost) == 0

        # default heat cost is zero and there is no heat cost in system costing block
        assert value(m.fs.treatment.costing.heat_cost) == 0
        assert value(m.fs.energy.costing.heat_cost) == 0
        assert not hasattr(m.fs.costing, "heat_cost")
        assert hasattr(m.fs.costing, "heat_cost_buy")


class TestElectricityGenOnlyWithHeat:

    @pytest.fixture(scope="class")
    def energy_gen_only_with_heat(self):

        m = build_electricity_gen_only_with_heat()

        return m

    @pytest.mark.unit
    def test_build(self, energy_gen_only_with_heat):

        m = energy_gen_only_with_heat

        assert_units_consistent(m)

        assert degrees_of_freedom(m) == 0

        # have heat flows
        assert m.fs.costing.has_heat_flows
        assert not m.fs.costing.aggregate_flow_heat.is_fixed()
        assert not m.fs.costing.aggregate_flow_heat_purchased.is_fixed()
        # no heat generated so nothing to sell
        assert m.fs.energy.costing.has_electricity_generation
        assert hasattr(m.fs.costing, "frac_elec_from_grid_constraint")
        assert not hasattr(m.fs.costing, "frac_heat_from_grid")

    @pytest.mark.component
    def test_init_and_solve(self, energy_gen_only_with_heat):
        m = energy_gen_only_with_heat

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        # check state after initialization
        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.costing.frac_elec_from_grid), rel=1e-3) == 0.9
        assert (
            pytest.approx(
                value(m.fs.costing.aggregate_flow_electricity_purchased), rel=1e-3
            )
            == 90
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_electricity
            + m.fs.energy.costing.aggregate_flow_electricity
        )
        assert pytest.approx(
            value(m.fs.costing.frac_elec_from_grid), rel=1e-3
        ) == 1 - value(m.fs.energy.unit.electricity) / value(
            m.fs.treatment.unit.electricity_consumption
        )

        # no heat is generated
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_heat), rel=1e-3
        ) == value(m.fs.treatment.unit.heat_consumption)

    @pytest.mark.component
    def test_check_proper_aggregation(self, energy_gen_only_with_heat):
        m = energy_gen_only_with_heat
        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )

    @pytest.mark.component
    def test_optimize_frac_from_grid(self):

        m = build_electricity_gen_only_with_heat()

        m.fs.energy.unit.electricity.unfix()
        m.fs.costing.frac_elec_from_grid.fix(0.05)

        assert degrees_of_freedom(m) == 0

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert (
            pytest.approx(
                value(m.fs.costing.aggregate_flow_electricity_purchased), rel=1e-3
            )
            == 5
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_electricity
            + m.fs.energy.costing.aggregate_flow_electricity
        )

        # test aggregation
        assert pytest.approx(
            value(m.fs.costing.aggregate_capital_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_capital_cost
            + m.fs.energy.costing.aggregate_capital_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_fixed_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_fixed_operating_cost
            + m.fs.energy.costing.aggregate_fixed_operating_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_variable_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_variable_operating_cost
            + m.fs.energy.costing.aggregate_variable_operating_cost
        )

        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )


class TestElectricityGenOnlyNoHeat:

    @pytest.fixture(scope="class")
    def elec_gen_only_no_heat(self):

        m = build_electricity_gen_only_no_heat()

        return m

    @pytest.mark.unit
    def test_build(self, elec_gen_only_no_heat):

        m = elec_gen_only_no_heat

        assert_units_consistent(m)

        assert degrees_of_freedom(m) == 0

        # no heat flows
        assert not m.fs.costing.has_heat_flows
        assert m.fs.costing.aggregate_flow_heat_purchased.is_fixed()
        assert m.fs.energy.costing.has_electricity_generation
        assert hasattr(m.fs.costing, "frac_elec_from_grid_constraint")
        assert not hasattr(m.fs.costing, "frac_heat_from_grid")

    @pytest.mark.component
    def test_init_and_solve(self, elec_gen_only_no_heat):
        m = elec_gen_only_no_heat

        # constraints are active before initialization
        assert m.fs.costing.total_heat_operating_cost_constraint.active
        assert m.fs.costing.aggregate_flow_heat_constraint.active

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        # check state after initialization
        assert degrees_of_freedom(m) == 0

        assert (m.fs.costing.total_heat_operating_cost.is_fixed()) and (
            value(m.fs.costing.total_heat_operating_cost) == 0
        )
        assert (m.fs.costing.aggregate_flow_heat.is_fixed()) and (
            value(m.fs.costing.aggregate_flow_heat) == 0
        )
        assert not m.fs.costing.total_heat_operating_cost_constraint.active
        assert not m.fs.costing.aggregate_flow_heat_constraint.active

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.costing.frac_elec_from_grid), rel=1e-3) == 0.25
        assert (
            pytest.approx(
                value(m.fs.costing.aggregate_flow_electricity_purchased), rel=1e-3
            )
            == 25
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_electricity
            + m.fs.energy.costing.aggregate_flow_electricity
        )
        assert pytest.approx(
            value(m.fs.costing.frac_elec_from_grid), rel=1e-3
        ) == 1 - value(m.fs.energy.unit.electricity) / value(
            m.fs.treatment.unit.electricity_consumption
        )

        # no heat is generated or consumed
        assert pytest.approx(value(m.fs.costing.aggregate_flow_heat), rel=1e-3) == 0

        # test aggregation
        assert pytest.approx(
            value(m.fs.costing.aggregate_capital_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_capital_cost
            + m.fs.energy.costing.aggregate_capital_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_fixed_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_fixed_operating_cost
            + m.fs.energy.costing.aggregate_fixed_operating_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_variable_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_variable_operating_cost
            + m.fs.energy.costing.aggregate_variable_operating_cost
        )

    @pytest.mark.component
    def test_check_proper_aggregation(self, elec_gen_only_no_heat):
        m = elec_gen_only_no_heat
        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )

    @pytest.mark.component
    def test_optimize_frac_from_grid(self):

        m = build_electricity_gen_only_no_heat()

        m.fs.energy.unit.electricity.unfix()
        m.fs.costing.frac_elec_from_grid.fix(0.33)

        assert degrees_of_freedom(m) == 0

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert (
            pytest.approx(
                value(m.fs.costing.aggregate_flow_electricity_purchased), rel=1e-3
            )
            == 33
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_electricity
            + m.fs.energy.costing.aggregate_flow_electricity
        )

        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )


class TestHeatGenOnly:
    @pytest.fixture(scope="class")
    def heat_gen_only(self):

        m = build_heat_gen_only()

        return m

    @pytest.mark.unit
    def test_build(self, heat_gen_only):

        m = heat_gen_only

        assert_units_consistent(m)

        assert degrees_of_freedom(m) == 0

        # has heat flows, no electricity generation
        assert m.fs.costing.has_heat_flows
        assert not m.fs.costing.aggregate_flow_heat_purchased.is_fixed()
        assert not m.fs.energy.costing.has_electricity_generation
        assert not hasattr(m.fs.costing, "frac_elec_from_grid_constraint")
        assert m.fs.costing.frac_elec_from_grid.is_fixed()
        assert hasattr(m.fs.costing, "frac_heat_from_grid")
        assert hasattr(m.fs.costing, "frac_heat_from_grid_constraint")

    @pytest.mark.component
    def test_init_and_solve(self, heat_gen_only):

        m = heat_gen_only

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0
        assert not m.fs.costing.aggregate_flow_heat_purchased.is_fixed()

        results = solver.solve(m)
        assert_optimal_termination(results)

        # all electricity comes from grid, none is generated
        assert pytest.approx(value(m.fs.costing.frac_elec_from_grid), rel=1e-3) == 1
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity_purchased), rel=1e-3
        ) == value(
            m.fs.treatment.unit.electricity_consumption + m.fs.energy.unit.electricity
        )
        # both energy and treatment processes are consuming electricity
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.unit.electricity_consumption + m.fs.energy.unit.electricity
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        assert pytest.approx(
            value(m.fs.costing.frac_heat_from_grid), rel=1e-3
        ) == 1 - value(m.fs.energy.unit.heat / m.fs.treatment.unit.heat_consumption)
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_heat), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_heat
            * m.fs.costing.frac_heat_from_grid
        )

    @pytest.mark.component
    def test_check_proper_aggregation(self, heat_gen_only):
        m = heat_gen_only
        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )

    @pytest.mark.component
    def test_optimize_frac_from_grid(self):

        m = build_heat_gen_only()

        m.fs.energy.unit.heat.unfix()
        m.fs.costing.frac_heat_from_grid.fix(0.002)

        assert degrees_of_freedom(m) == 0

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert (
            pytest.approx(value(m.fs.costing.aggregate_flow_heat_purchased), rel=1e-3)
            == 5
        )
        assert pytest.approx(value(m.fs.costing.aggregate_flow_heat), rel=1e-3) == 5
        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(m.fs.costing.aggregate_flow_electricity_purchased)
        # test aggregation
        assert pytest.approx(
            value(m.fs.costing.aggregate_capital_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_capital_cost
            + m.fs.energy.costing.aggregate_capital_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_fixed_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_fixed_operating_cost
            + m.fs.energy.costing.aggregate_fixed_operating_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_variable_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_variable_operating_cost
            + m.fs.energy.costing.aggregate_variable_operating_cost
        )

        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )


class TestElectricityAndHeatGen:

    @pytest.fixture(scope="class")
    def heat_and_elec_gen(self):

        m = build_heat_and_elec_gen()

        assert_units_consistent(m)

        m.fs.costing.add_LCOE()
        assert hasattr(m.fs.energy.costing, "yearly_electricity_production")
        assert hasattr(m.fs.energy.costing, "lifetime_electricity_production")
        m.fs.costing.add_LCOH()
        assert hasattr(m.fs.energy.costing, "yearly_heat_production")
        assert hasattr(m.fs.energy.costing, "lifetime_heat_production")
        m.fs.costing.add_LCOT(m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_electric_energy_consumption(
            m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_thermal_energy_consumption(
            m.fs.treatment.unit.properties[0].flow_vol_phase["Liq"]
        )

        return m

    @pytest.mark.unit
    def test_build(self, heat_and_elec_gen):

        m = heat_and_elec_gen

        assert_units_consistent(m)

        assert degrees_of_freedom(m) == 0

        # has heat and electricity flows
        assert m.fs.costing.has_heat_flows
        assert not m.fs.costing.aggregate_flow_heat_purchased.is_fixed()
        assert not m.fs.costing.aggregate_flow_electricity_purchased.is_fixed()
        assert m.fs.energy.costing.has_electricity_generation
        assert hasattr(m.fs.costing, "frac_elec_from_grid_constraint")
        assert not m.fs.costing.frac_elec_from_grid.is_fixed()
        assert hasattr(m.fs.costing, "frac_heat_from_grid")
        assert hasattr(m.fs.costing, "frac_heat_from_grid_constraint")

        # all metrics end up as references on REFLOSystemCosting
        assert hasattr(m.fs.costing, "LCOT")
        assert hasattr(m.fs.costing, "LCOE")
        assert hasattr(m.fs.costing, "LCOH")
        assert hasattr(m.fs.costing, "SEEC")
        assert hasattr(m.fs.costing, "STEC")
        # energy metrics on EnergyCosting
        assert hasattr(m.fs.energy.costing, "LCOE")
        assert hasattr(m.fs.energy.costing, "LCOH")
        # treatment metrics on TreatmentCosting
        assert not hasattr(m.fs.treatment.costing, "LCOT")
        assert hasattr(m.fs.treatment.costing, "SEEC")
        assert hasattr(m.fs.treatment.costing, "STEC")

    @pytest.mark.component
    def test_init_and_solve(self, heat_and_elec_gen):
        m = heat_and_elec_gen

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(
            value(m.fs.costing.aggregate_flow_electricity), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_flow_electricity
            + m.fs.energy.costing.aggregate_flow_electricity
        )
        # fraction from grid is generated electricity from elec_unit
        # over consumed electricity from heat_unit and treatment unit
        assert pytest.approx(
            value(m.fs.costing.frac_elec_from_grid), rel=1e-3
        ) == 1 - value(
            m.fs.energy.elec_unit.electricity
            / (
                m.fs.treatment.unit.electricity_consumption
                + m.fs.energy.heat_unit.electricity
            )
        )
        # two equivalent ways to calculate fraction heat from grid
        assert pytest.approx(
            value(m.fs.costing.frac_heat_from_grid), rel=1e-3
        ) == 1 - value(
            m.fs.energy.heat_unit.heat / m.fs.treatment.unit.heat_consumption
        )
        assert pytest.approx(
            value(m.fs.costing.frac_heat_from_grid), rel=1e-3
        ) == 1 - value(
            -1
            * m.fs.energy.costing.aggregate_flow_heat
            / m.fs.treatment.costing.aggregate_flow_heat
        )

    @pytest.mark.component
    def test_check_proper_aggregation(self, heat_and_elec_gen):
        m = heat_and_elec_gen
        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )

    @pytest.mark.component
    def test_optimize_frac_from_grid(self):

        m = build_heat_and_elec_gen()

        m.fs.energy.elec_unit.electricity.unfix()
        m.fs.energy.heat_unit.heat.unfix()

        m.fs.costing.frac_elec_from_grid.fix(0.99)
        m.fs.costing.frac_heat_from_grid.fix(0.85)

        assert degrees_of_freedom(m) == 0

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        # fraction from grid is generated electricity from elec_unit
        # over consumed electricity from heat_unit and treatment unit
        assert pytest.approx(
            value(m.fs.costing.frac_elec_from_grid), rel=1e-3
        ) == 1 - value(
            m.fs.energy.elec_unit.electricity
            / (
                m.fs.treatment.unit.electricity_consumption
                + m.fs.energy.heat_unit.electricity
            )
        )
        # two equivalent ways to calculate fraction heat from grid
        assert pytest.approx(
            value(m.fs.costing.frac_heat_from_grid), rel=1e-3
        ) == 1 - value(
            m.fs.energy.heat_unit.heat / m.fs.treatment.unit.heat_consumption
        )
        assert pytest.approx(
            value(m.fs.costing.frac_heat_from_grid), rel=1e-3
        ) == 1 - value(
            -1
            * m.fs.energy.costing.aggregate_flow_heat
            / m.fs.treatment.costing.aggregate_flow_heat
        )
        # test aggregation
        assert pytest.approx(
            value(m.fs.costing.aggregate_capital_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_capital_cost
            + m.fs.energy.costing.aggregate_capital_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_fixed_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_fixed_operating_cost
            + m.fs.energy.costing.aggregate_fixed_operating_cost
        )
        assert pytest.approx(
            value(m.fs.costing.aggregate_variable_operating_cost), rel=1e-3
        ) == value(
            m.fs.treatment.costing.aggregate_variable_operating_cost
            + m.fs.energy.costing.aggregate_variable_operating_cost
        )

        check_proper_aggregation(
            m.fs.treatment.costing, m.fs.energy.costing, m.fs.costing
        )


@pytest.mark.component
def test_no_energy_treatment_block():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)

    with pytest.raises(
        ValueError,
        match="REFLOSystemCosting package requires a EnergyCosting block but one was not found\\.",
    ):
        m.fs.costing = REFLOSystemCosting()


@pytest.mark.component
def test_add_LCOW_to_energy_costing():

    m = build_default()

    m.fs.energy.costing.cost_process()
    m.fs.treatment.costing.cost_process()

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()

    with pytest.raises(ValueError, match="Can't add LCOW to EnergyCosting package\\."):
        m.fs.energy.costing.add_LCOW()


@pytest.mark.component
def test_lazy_flow_costing():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    m.fs.electricity = Var(units=pyunits.kW)
    m.fs.costing.cost_flow(m.fs.electricity, "electricity")

    assert "foo" not in m.fs.costing.flow_types
    with pytest.raises(
        ValueError,
        match="foo is not a recognized flow type. Please check "
        "your spelling and that the flow type has been registered with"
        " the FlowsheetCostingBlock.",
    ):
        m.fs.costing.cost_flow(m.fs.electricity, "foo")

    m.fs.costing.foo_cost = foo_cost = Var(
        initialize=42, doc="foo", units=pyunits.USD_2020 / pyunits.m
    )

    m.fs.costing.register_flow_type("foo", m.fs.costing.foo_cost)

    # make sure the component was not replaced
    # by register_flow_type
    assert foo_cost is m.fs.costing.foo_cost

    assert "foo" in m.fs.costing.flow_types

    # not used until aggregated
    assert "foo" not in m.fs.costing.used_flows

    m.fs.foo = Var(units=pyunits.m / pyunits.year)

    m.fs.costing.cost_flow(m.fs.foo, "foo")
    m.fs.costing.aggregate_costs()

    # now should be used
    assert "foo" in m.fs.costing.used_flows

    m.fs.costing.bar_base_cost = Var(
        initialize=0.42, doc="bar", units=pyunits.USD_2020 / pyunits.g
    )
    m.fs.costing.bar_purity = Param(
        initialize=0.50, doc="bar purity", units=pyunits.dimensionless
    )

    m.fs.costing.register_flow_type(
        "bar", m.fs.costing.bar_base_cost * m.fs.costing.bar_purity
    )

    bar_cost = m.fs.costing.bar_cost
    assert isinstance(bar_cost, Expression)
    assert value(bar_cost) == 0.21

    m.fs.costing.bar_base_cost.value = 1.5
    assert value(bar_cost) == 0.75

    m.fs.costing.baz_cost = Var()

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "Component baz_cost already exists on fs.costing but is not 42*USD_2020/m**2."
        ),
    ):
        m.fs.costing.register_flow_type("baz", 42 * pyunits.USD_2020 / pyunits.m**2)
