import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)

def build_tes_cost_param_block(blk):

    costing = blk.parent_block()

    thermal_energy_storage_cost = 10.79 $/kW


@register_costing_parameter_block(
    build_rule=build_tes_cost_param_block,
    parameter_block_name='tes',
)
def cost_tes(blk):

    global_params = blk.costing_package
    tes_params = blk.costing_package.thermal_energy_storage
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = thermal_energy_storage_cost*storage_capacity

    blk.indirect_cost

    blk.sales_tax = 5%

    blk.direct_cost_constraint

    blk.indirect_cost_constraint

    blk.sales_tax_constraint

    blk.capital_cost_constraint

    blk.fixed_operating_cost_constraint

