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
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import ConfigurationError
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)

# Costs are defaults from SAM 2024.12.12
# Photovoltaic > Detailed PV Model > Single Owner


class PVSurrogateCostMethod(StrEnum):
    detailed = "detailed"
    simple = "simple"


def cost_pv(blk, cost_method=PVSurrogateCostMethod.simple):
    blk.cost_method = cost_method
    blk.costing_package.has_electricity_generation = True
    if cost_method == PVSurrogateCostMethod.detailed:
        cost_pv_detailed(blk)
    elif cost_method == PVSurrogateCostMethod.simple:
        cost_pv_simple(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for cost_method:"
            f" {cost_method}. Argument must be a member of the PVSurrogateCostMethod Enum."
        )


def build_cost_pv_simple_param_block(blk):
    costing = blk.parent_block()

    blk.cost_per_watt_installed = pyo.Var(
        initialize=1.6,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for solar module installed",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=31,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of PV system per kW generated",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0,
        units=costing.base_currency / pyo.units.MWh,
        bounds=(0, None),
        doc="Annual operating cost of PV system per MWh generated",
    )

    blk.fix_all_vars()


def build_cost_pv_detailed_param_block(blk):

    costing = blk.parent_block()

    blk.cost_per_watt_module = pyo.Var(
        initialize=0.34,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for solar module",
    )

    blk.cost_per_watt_inverter = pyo.Var(
        initialize=0.03,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for inverter",
    )

    blk.cost_per_watt_other_direct = pyo.Var(
        initialize=0.62,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for balance of system equipment, installation labor, and margin/overhead",
    )

    blk.cost_per_watt_indirect = pyo.Var(
        initialize=0.05,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for permitting, environmental studies, engineering, land prep, and grid interconnection",
    )

    blk.contingency_frac_direct_cost = pyo.Var(
        initialize=0.03,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs to apply contingency",
    )

    blk.tax_frac_direct_cost = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs to apply sales tax",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=31,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of PV system per kW generated",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0,
        units=costing.base_currency / pyo.units.MWh,
        bounds=(0, None),
        doc="Annual operating cost of PV system per MWh generated",
    )

    blk.fix_all_vars()


@register_costing_parameter_block(
    build_rule=build_cost_pv_detailed_param_block,
    parameter_block_name="pv",
)
def cost_pv_detailed(blk):

    global_params = blk.costing_package
    pv_params = blk.costing_package.pv

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Direct costs of PV system",
    )

    blk.indirect_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Indirect costs of PV system",
    )

    blk.land_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Land costs of PV system",
    )

    blk.sales_tax = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Sales tax for PV system",
    )

    system_capacity_watt = pyo.units.convert(
        blk.unit_model.system_capacity, to_units=pyo.units.watt
    )

    inverter_capacity_watt = pyo.units.convert(
        blk.unit_model.inverter_capacity, to_units=pyo.units.watt
    )

    capital_cost_expr = 0

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == pyo.units.convert(
            (
                system_capacity_watt
                * (
                    pv_params.cost_per_watt_module
                    + pv_params.cost_per_watt_other_direct
                )
                + inverter_capacity_watt * pv_params.cost_per_watt_inverter
            )
            * (1 + pv_params.contingency_frac_direct_cost),  # BUG Check this
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.direct_cost

    blk.indirect_capital_cost_constraint = pyo.Constraint(
        expr=blk.indirect_cost
        == pyo.units.convert(
            system_capacity_watt * pv_params.cost_per_watt_indirect,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.indirect_cost

    blk.land_cost_constraint = pyo.Constraint(
        expr=blk.land_cost
        == pyo.units.convert(
            blk.unit_model.land_req * global_params.land_cost,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.land_cost

    blk.sales_tax_constraint = pyo.Constraint(
        expr=blk.sales_tax
        == pyo.units.convert(
            blk.direct_cost
            * pv_params.tax_frac_direct_cost
            * global_params.sales_tax_frac,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.sales_tax

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            capital_cost_expr, to_units=blk.costing_package.base_currency
        )
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            pv_params.fixed_operating_by_capacity * system_capacity_watt,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == pyo.units.convert(
            pv_params.variable_operating_by_generation
            * blk.unit_model.electricity_annual,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.total_installed_cost_per_capacity = pyo.Expression(
        expr=blk.capital_cost / system_capacity_watt,
        doc="Total installed cost per capacity of PV system",
    )

    blk.costing_package.cost_flow(-1 * blk.unit_model.electricity, "electricity")


@register_costing_parameter_block(
    build_rule=build_cost_pv_simple_param_block,
    parameter_block_name="pv",
)
def cost_pv_simple(blk):

    global_params = blk.costing_package
    pv_params = blk.costing_package.pv
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Direct costs of PV system",
    )

    blk.land_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Land costs of PV system",
    )

    capital_cost_expr = 0

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == pyo.units.convert(
            pyo.units.convert(blk.unit_model.system_capacity, to_units=pyo.units.watt)
            * pv_params.cost_per_watt_installed,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.direct_cost

    blk.land_cost_constraint = pyo.Constraint(
        expr=blk.land_cost
        == pyo.units.convert(
            (blk.unit_model.land_req * global_params.land_cost),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.land_cost

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            capital_cost_expr, to_units=blk.costing_package.base_currency
        )
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            pv_params.fixed_operating_by_capacity
            * pyo.units.convert(blk.unit_model.system_capacity, to_units=pyo.units.kW),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == pyo.units.convert(
            pv_params.variable_operating_by_generation
            * blk.unit_model.electricity_annual,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.costing_package.cost_flow(-1 * blk.unit_model.electricity, "electricity")
