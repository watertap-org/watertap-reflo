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


from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Var,
    Param,
    units as pyunits,
)

from watertap.costing.util import register_costing_parameter_block
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_deep_well_injection_cost_param_block(blk):
    """
    Parameters and variables to be used in the costing model
    References :

    """
    costing = blk.parent_block()



@register_costing_parameter_block(
    build_rule=build_deep_well_injection_cost_param_block,
    parameter_block_name="deep_well_injection",
)
def cost_deep_well_injection(blk):
    """
    Capital and operating costs for deep well injection
    """
    dwi = blk.costing_package.deep_well_injection
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
