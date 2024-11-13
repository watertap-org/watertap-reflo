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
    assert_optimal_termination,
    value,
    units as pyunits,
)

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
    REFLOSystemCosting,
)
from pyomo.environ import ConcreteModel, Var, Param, Expression, value, units as pyunits

from idaes.core import FlowsheetBlock

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.tests.costing_dummy_units import (
    DummyTreatmentUnit,
    DummyElectricityUnit,
    DummyHeatUnit,
)
from watertap_contrib.reflo.costing import (
    REFLOCosting,
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

solver = get_solver()


def build_electricity_gen_only():
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

    #### ENERGY BLOCK
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyElectricityUnit()
    m.fs.energy.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.unit.electricity.fix()
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

    #### INITIALIZE

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


class TestElectricityGenOnly:

    @pytest.fixture(scope="class")
    def energy_gen_only(self):

        m = build_electricity_gen_only()

        return m

    @pytest.mark.unit
    def test_build(slef, energy_gen_only):

        m = energy_gen_only

        assert degrees_of_freedom(m) == 0

        assert m.fs.energy.costing.has_electricity_generation

        m.fs.treatment.unit.initialize()
        m.fs.treatment.costing.initialize()
        m.fs.energy.costing.initialize()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)


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
def test_common_params_not_equivalent():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)

    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyElectricityUnit()

    m.fs.energy.costing.electricity_cost.fix(0.02)

    with pytest.raises(
        ValueError,
        match="The common costing parameter electricity_cost was found to "
        "have a different value on the energy \\(0\\.02\\) and treatment \\(0\\.0\\) costing "
        "blocks\\. Common costing parameters must be equivalent across all"
        " costing blocks to use REFLOSystemCosting\\.",
    ):
        m.fs.costing = REFLOSystemCosting()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.unit = DummyTreatmentUnit(property_package=m.fs.properties)

    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.unit = DummyElectricityUnit()

    m.fs.energy.costing.electricity_cost.fix(0.02)
    m.fs.treatment.costing.electricity_cost.fix(0.02)

    m.fs.costing = REFLOSystemCosting()

    # assert value(m.fs.costing.electricity_cost) == value(m.fs.treatment.electricity_cost)
    # assert value(m.fs.costing.electricity_cost) == value(m.fs.energy.electricity_cost)


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
