import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing equations from:
# Kosmadakis G, Papapetrou M, Ortega-Delgado B, Cipollina A, Alarcón-Padilla D-C.
# "Correlations for estimating the specific capital cost of multi-effect distillation plants
#    considering the main design trends and operating conditions"
# doi: 10.1016/j.desal.2018.09.011


def build_electrocoagulation_cost_param_block(blk):

    costing = blk.parent_block()

    blk.fix_all_vars()


@register_costing_parameter_block(
    build_rule=build_electrocoagulation_cost_param_block,
    parameter_block_name="electrocoagulation",
)
def cost_electrocoagulation(blk):

    ec_params = blk.costing_package.electrocoagulation
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    ec = blk.unit_model
    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
