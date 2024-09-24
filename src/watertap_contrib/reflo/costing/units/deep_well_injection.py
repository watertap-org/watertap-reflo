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


from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Var,
    Param,
    log,
    exp,
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

    # All default values are for 2,500 ft depth

    blk.logging_testing_capital_cost_intercept = Param(
        initialize=251.48,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Logging + testing capital cost equation - intercept",
    )

    blk.logging_testing_capital_cost_slope = Param(
        initialize=4.5193,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Logging + testing capital cost equation - slope",
    )

    blk.drilling_capital_cost_intercept = Param(
        initialize=338.05,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Drilling capital cost equation - intercept",
    )

    blk.drilling_capital_cost_slope = Param(
        initialize=17.611,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Drilling capital cost equation - slope",
    )

    blk.tubing_capital_cost_base = Param(
        initialize=116.65,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Tubing capital cost equation - base",
    )

    blk.tubing_capital_cost_exponent = Param(
        initialize=0.3786,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Tubing capital cost equation - exponent",
    )

    blk.packing_capital_cost_base = Param(
        initialize=579.85,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Packing capital cost equation - base",
    )

    blk.packing_capital_cost_intercept = Param(
        initialize=596.64,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Packing capital cost equation - intercept",
    )

    blk.casing_capital_cost_intercept = Param(
        initialize=20.545,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Casing capital cost equation - intercept",
    )

    blk.casing_capital_cost_slope = Param(
        initialize=384.36,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Casing capital cost equation - slope",
    )

    blk.grouting_capital_cost_base = Param(
        initialize=145.01,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Grouting capital cost equation - base",
    )

    blk.grouting_capital_cost_exp_coeff = Param(
        initialize=0.0708,
        mutable=True,
        units=pyunits.inches**-1,
        doc="Grouting capital cost equation - exponential coefficient",
    )

    blk.monitoring_well_capital_cost_base = Param(
        initialize=31.337,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Monitoring well capital cost equation - base",
    )

    blk.monitoring_well_capital_cost_exponent = Param(
        initialize=0.3742,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Monitoring well capital cost equation - exponent",
    )

    blk.mobilization_capital_cost_base = Param(
        initialize=31.337,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Mobilization capital cost equation - base",
    )

    blk.mobilization_capital_cost_exponent = Param(
        initialize=0.3742,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Mobilization capital cost equation - exponent",
    )


@register_costing_parameter_block(
    build_rule=build_deep_well_injection_cost_param_block,
    parameter_block_name="deep_well_injection",
)
def cost_deep_well_injection(blk):
    """
    Capital and operating costs for deep well injection
    """
    dwi_params = blk.costing_package.deep_well_injection
    dwi = blk.unit_model
    pipe_diam_dimensionless = pyunits.convert(
        dwi.pipe_diameter * pyunits.inches**-1, to_units=pyunits.dimensionless
    )
    well_depth_dimensionless = pyunits.convert(dwi.well_depth * pyunits.ft**-1, to_units=pyunits.dimensionless)
    make_capital_cost_var(blk)
    # make_fixed_operating_cost_var(blk)

    blk.logging_testing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Logging + testing capital costs",
    )

    blk.drilling_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Drilling capital costs",
    )

    blk.tubing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Tubing capital costs",
    )

    blk.packing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Packing capital costs",
    )

    blk.casing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Casing capital costs",
    )

    blk.grouting_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Grouting capital costs",
    )

    blk.monitoring_well_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Monitoring well capital costs",
    )

    blk.mobilization_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Mobilization capital costs",
    )

    capital_cost_expr = 0

    blk.logging_testing_capital_cost_constraint = Constraint(
        expr=blk.logging_testing_capital_cost
        == pyunits.convert(
            dwi_params.logging_testing_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.logging_testing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.logging_testing_capital_cost

    blk.drilling_capital_cost_constraint = Constraint(
        expr=blk.drilling_capital_cost
        == pyunits.convert(
            dwi_params.drilling_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.drilling_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.drilling_capital_cost

    blk.tubing_capital_cost_constraint = Constraint(
        expr=blk.tubing_capital_cost
        == pyunits.convert(
            dwi_params.tubing_capital_cost_base
            * pipe_diam_dimensionless**dwi_params.tubing_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.tubing_capital_cost

    blk.packing_capital_cost_constraint = Constraint(
        expr=blk.packing_capital_cost
        == pyunits.convert(
            dwi_params.packing_capital_cost_base * log(pipe_diam_dimensionless)
            + dwi_params.packing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.packing_capital_cost

    blk.casing_capital_cost_constraint = Constraint(
        expr=blk.casing_capital_cost
        == pyunits.convert(
            dwi_params.casing_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.casing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.casing_capital_cost

    blk.grouting_capital_cost_constraint = Constraint(
        expr=blk.grouting_capital_cost
        == pyunits.convert(
            dwi_params.grouting_capital_cost_base
            * exp(dwi_params.grouting_capital_cost_exp_coeff * dwi.pipe_diameter),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.grouting_capital_cost

    blk.monitoring_well_capital_cost_constraint = Constraint(
        expr=blk.monitoring_well_capital_cost
        == pyunits.convert(
            dwi_params.monitoring_well_capital_cost_base
            * well_depth_dimensionless**dwi_params.monitoring_well_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.monitoring_well_capital_cost

    blk.mobilization_capital_cost_constraint = Constraint(
        expr=blk.mobilization_capital_cost
        == pyunits.convert(
            dwi_params.mobilization_capital_cost_base
            * well_depth_dimensionless**dwi_params.mobilization_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.mobilization_capital_cost

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = Constraint(expr=blk.capital_cost == capital_cost_expr)
