import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)

def build_lt_med_surrogate_cost_param_block(blk):

    costing = blk.parent_block()

    
    blk.cost_fraction_evaporator = pyo.Var(
    initialize=0.4,
    units=pyo.units.dimensionless,
    bounds=(0, None),
    doc="Cost fraction of the evaporator",
)

    blk.fix_all_vars()


@register_costing_parameter_block(build_rule=build_lt_med_surrogate_cost_param_block, parameter_block_name="lt_med_surrogate")

def cost_lt_med_surrogate(blk):

    lt_med_params = blk.costing_package.lt_med_surrogate
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    # blk.costing_package.cost_flow(blk.unit_model.electricity, "electricity")

