# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Constraint,
    Expression,
    NonNegativeReals,
    units as pyunits,
)

from watertap.costing.util import register_costing_parameter_block

def make_fixed_operating_cost_var(blk):
    blk.fixed_operating_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit fixed operating cost",
    )

def build_disposal_cost_param_block(blk):
    costing = blk.parent_block()
    blk.cost_of_disposal = Param(
        initialize=0.0587,
        mutable=True,
        units=costing.base_currency / pyunits.m**3,
        doc="Assumed LCOW of deep well injection",
    )

@register_costing_parameter_block(
    build_rule=build_disposal_cost_param_block,
    parameter_block_name="brine_disposal",
)
def cost_disposal(blk):

    make_fixed_operating_cost_var(blk)

    # blk.costing_package.display()
    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        == pyunits.convert(
            blk.unit_model.properties[0].flow_vol_phase["Liq"] * blk.costing_package.brine_disposal.cost_of_disposal,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )