from watertap.core.zero_order_costing import ZeroOrderCostingData
from idaes.core import declare_process_block_class
import pyomo.environ as pyo
from idaes.core.util.misc import StrEnum
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


def build_chem_softening_cost_param_block(blk):
    '''
        Parameters and variables to be used in the costing model
    '''

    blk.unit_cost = pyo.Var(
        initialize = 10,
        doc = 'sample',
        units = pyo.dimensionless,
    )


@register_costing_parameter_block(
    build_rule = build_chem_softening_cost_param_block,
    parameter_block_name = "chemical_softening"
    )

def cost_chem_softening(blk, softeningType = SofteningProcedureType.single_stage_lime):
    '''
        Function for costing chemical softening unit for 
        1. Lime softening
    '''
    if (softeningType == SofteningProcedureType.single_stage_lime):
        cost_chem_softening_single_stage_lime(blk)


def cost_chem_softening_single_stage_lime(blk):
    '''
        Capital cost for single stage lime softening
    '''

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        ==

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        ==