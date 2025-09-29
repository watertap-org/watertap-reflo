#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)


# Costs are defaults from SAM 2024.12.12
# Wind + Single Owner


def build_wind_cost_param_block(blk):

    costing = blk.parent_block()

    blk.cost_turbine_per_kw = pyo.Var(
        initialize=1112.4,
        units=costing.base_currency / pyo.units.kilowatt,
        bounds=(0, None),
        doc="Cost per kW for wind turbine",
    )

    blk.cost_turbine_per = pyo.Var(
        initialize=0,
        units=costing.base_currency,  # per turbine
        bounds=(0, None),
        doc="Unit cost per turbine",
    )

    blk.cost_balance_system_per_kw = pyo.Var(
        initialize=0,
        units=costing.base_currency / pyo.units.kilowatt,
        bounds=(0, None),
        doc="Balance of system costs per kW",
    )

    blk.contingency_frac_direct_cost = pyo.Var(
        initialize=0,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs for contingency",
    )

    blk.indirect_frac_direct_cost = pyo.Var(
        initialize=0,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs, including contingency, for indirect costs",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=40,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of wind plant per kW capacity",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0,
        units=costing.base_currency / pyo.units.MWh,
        bounds=(0, None),
        doc="Variable operating cost of wind plant per MWh generated",
    )


@register_costing_parameter_block(
    build_rule=build_wind_cost_param_block,
    parameter_block_name="wind",
)
def cost_wind(blk):

    global_params = blk.costing_package
    wind_params = blk.costing_package.wind
    wind = blk.unit_model
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_fixed_operating_cost_var(blk)
    make_variable_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=1e4,
        units=global_params.base_currency,
        bounds=(0, None),
        doc="Direct cost of flat plate system",
    )

    blk.indirect_cost = pyo.Var(
        initialize=1e4,
        units=global_params.base_currency,
        bounds=(0, None),
        doc="Indirect costs of flat plate system",
    )

    blk.sales_tax = pyo.Var(
        initialize=1e2,
        units=global_params.base_currency,
        bounds=(0, None),
        doc="Sales tax for flat plate system",
    )

    capital_cost_expr = 0

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            capital_cost_expr,
            to_units=global_params.base_currency,
        )
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            wind_params.fixed_operating_by_capacity
            * pyo.units.convert(wind.system_capacity, to_units=pyo.units.kW),
            to_units=global_params.base_currency / global_params.base_period,
        )
    )

    blk.costing_package.cost_flow(
        wind.electricity,
        "electricity",
    )
