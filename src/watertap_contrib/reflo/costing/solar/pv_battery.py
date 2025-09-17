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
from watertap_contrib.reflo.costing.solar.pv import (
    build_cost_pv_simple_param_block,
    build_cost_pv_detailed_param_block,
)

# Costs are defaults from SAM 2024.12.12
# Energy Storage > Detailed PV-Battery > Single Owner


def build_pv_battery_surrogate_cost_param_block(blk):

    build_cost_pv_detailed_param_block(blk)

    costing = blk.parent_block()

    blk.cost_per_kw_battery_power = pyo.Var(
        initialize=233,
        units=costing.base_currency / pyo.units.kilowatt,
        bounds=(0, None),
        doc="Direct cost per kW of battery power",
    )

    blk.cost_per_kwh_battery_storage = pyo.Var(
        initialize=252,
        units=costing.base_currency / pyo.units.kWh,
        bounds=(0, None),
        doc="Direct cost per kWh of battery storage",
    )

    blk.battery_fixed_operating_by_capacity = pyo.Var(
        initialize=7.25,
        units=costing.base_currency / (pyo.units.kWh * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of battery by capacity",
    )

    # blk.battery_variable_operating_by_discharged = pyo.Var(
    #     initialize=0,
    #     units=costing.base_currency / (pyo.units.MWh * costing.base_period),
    #     bounds=(0, None),
    #     doc="Variable operating cost of battery system per MWh discharged",
    # )

    blk.battery_replacement_frequency = pyo.Var(
        initialize=20,
        units=pyo.units.year,
        bounds=(0, None),
        doc="Replacement frequency of battery",
    )

    blk.battery_replacement_cost_by_capacity = pyo.Var(
        initialize=252,
        units=costing.base_currency / pyo.units.kWh,
        bounds=(0, None),
        doc="Replacement cost of battery by capacity",
    )


@register_costing_parameter_block(
    build_rule=build_pv_battery_surrogate_cost_param_block,
    parameter_block_name="pv_battery",
)
def cost_pv_battery(blk):

    blk.costing_package.has_electricity_generation = True
    global_params = blk.costing_package
    pv_batt_params = blk.costing_package.pv_battery

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Direct costs of PV + battery system",
    )

    blk.indirect_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Indirect costs of PV + battery system",
    )

    blk.land_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Land costs of PV + battery system",
    )

    blk.sales_tax = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Sales tax for PV + battery system",
    )

    system_capacity_watt = pyo.units.convert(
        blk.unit_model.system_capacity, to_units=pyo.units.watt
    )

    inverter_capacity_watt = pyo.units.convert(
        blk.unit_model.inverter_capacity, to_units=pyo.units.watt
    )

    batt_size_kwh = pyo.units.convert(
        blk.unit_model.battery_power * blk.unit_model.hours_storage,
        to_units=pyo.units.kWh,
    )

    blk.pv_direct_capital_cost = pyo.Expression(
        expr=pyo.units.convert(
            (
                system_capacity_watt
                * (
                    pv_batt_params.cost_per_watt_module
                    + pv_batt_params.cost_per_watt_other_direct
                )
                + inverter_capacity_watt * pv_batt_params.cost_per_watt_inverter
            ),
            to_units=blk.costing_package.base_currency,
        ),
        doc="Direct capital cost of PV system",
    )

    blk.battery_direct_capital_cost = pyo.Expression(
        expr=pyo.units.convert(
            (
                batt_size_kwh * pv_batt_params.cost_per_kwh_battery_storage
                + blk.unit_model.battery_power
                * pv_batt_params.cost_per_kw_battery_power
            ),
            to_units=blk.costing_package.base_currency,
        ),
        doc="Direct capital cost of battery system",
    )

    capital_cost_expr = 0

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == pyo.units.convert(
            (blk.pv_direct_capital_cost + blk.battery_direct_capital_cost)
            * (1 + pv_batt_params.contingency_frac_direct_cost),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.direct_cost

    blk.land_cost_constraint = pyo.Constraint(
        expr=blk.land_cost
        == pyo.units.convert(
            blk.unit_model.land_req * global_params.land_cost,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.land_cost

    blk.indirect_capital_cost_constraint = pyo.Constraint(
        expr=blk.indirect_cost
        == pyo.units.convert(
            system_capacity_watt * pv_batt_params.cost_per_watt_indirect,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.indirect_cost

    blk.sales_tax_constraint = pyo.Constraint(
        expr=blk.sales_tax
        == pyo.units.convert(
            blk.direct_cost
            * pv_batt_params.tax_frac_direct_cost
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

    fixed_op_cost_expr = 0

    blk.pv_fixed_operating_cost = pyo.Expression(
        expr=pyo.units.convert(
            pv_batt_params.fixed_operating_by_capacity * system_capacity_watt,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        ),
        doc="Fixed operating cost of PV system",
    )

    fixed_op_cost_expr += blk.pv_fixed_operating_cost

    blk.battery_fixed_operating_cost = pyo.Expression(
        expr=pyo.units.convert(
            pv_batt_params.battery_fixed_operating_by_capacity * batt_size_kwh,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        ),
        doc="Fixed operating cost of battery system",
    )

    fixed_op_cost_expr += blk.battery_fixed_operating_cost

    # NOTE: we are assuming a fraction of the total battery replacement cost is incurred each year
    blk.battery_replacement_cost = pyo.Expression(
        expr=pyo.units.convert(
            (pv_batt_params.battery_replacement_cost_by_capacity * batt_size_kwh)
            / pv_batt_params.battery_replacement_frequency,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        ),
        doc="Replacement cost of battery system",
    )

    fixed_op_cost_expr += blk.battery_replacement_cost

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            fixed_op_cost_expr,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    var_op_cost_expr = 0

    blk.pv_variable_operating_cost = pyo.Expression(
        expr=pyo.units.convert(
            pv_batt_params.variable_operating_by_generation
            * blk.unit_model.electricity_annual,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        ),
        doc="Variable operating cost of PV system",
    )

    var_op_cost_expr += blk.pv_variable_operating_cost

    # TODO: need battery discharge flow as surrogate output to make this calculation
    # blk.battery_variable_operating_cost = pyo.Expression(
    #     expr=pyo.units.convert(
    #         pv_batt_params.battery_variable_operating_by_discharged
    #         * batt_size_kwh,
    #         to_units=blk.costing_package.base_currency
    #         / blk.costing_package.base_period,
    #     ),
    #     doc="Variable operating cost of battery system",
    # )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost == var_op_cost_expr
    )

    blk.costing_package.cost_flow(-1 * blk.unit_model.electricity, "electricity")
