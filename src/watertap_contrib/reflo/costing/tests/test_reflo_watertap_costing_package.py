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
import re

import pytest

from pyomo.environ import ConcreteModel, Var, Param, Expression, value, units as pyunits

from idaes.core import FlowsheetBlock

from watertap_contrib.reflo.costing import REFLOCosting


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
