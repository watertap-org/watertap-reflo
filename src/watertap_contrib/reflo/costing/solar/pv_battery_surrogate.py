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

import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)
from watertap_contrib.reflo.costing.solar.pv_surrogate import build_pv_surrogate_cost_param_block

# Costs are defaults from SAM 2024.12.12
# Energy Storage > Detailed PV-Battery > Single Owner

def build_pv_battery_surrogate_cost_param_block(blk):
    
    build_pv_surrogate_cost_param_block(blk)
    
    costing = blk.parent_block()
    
    blk.cost_per_kilowatt_battery = pyo.Var(
        initialize=233,
        units=costing.base_currency / pyo.units.kilowatt,
        bounds=(0, None),
        doc="Cost per kW for battery",
    )

    blk.cost_per_kilowatt_hour_battery = pyo.Var(
        initialize=252,
        units=costing.base_currency / pyo.units.kWh,
        bounds=(0, None),
        doc="Cost per kWh for battery",
    )

    blk.battery_fixed_operating_by_capacity = pyo.Var(
        initialize=31,
        units=costing.base_currency / (pyo.units.kWh * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of battery by capacity ($/kWhdc)",
    )

    blk.battery_variable_operating_by_discharged = pyo.Var(
        initialize=0,
        units=costing.base_currency / (pyo.units.MWh * costing.base_period),
        bounds=(0, None),
        doc="Variable operating cost of battery system per MWh dischraged",
    )

    blk.battery_replacement_cost_by_capacity = pyo.Var(
        initialize=0,
        units=costing.base_currency / (pyo.units.kWh * costing.base_period),
        bounds=(0, None),
        doc="Replacement cost of battery by capacity ($/kWhdc)",
    )

@register_costing_parameter_block(
    build_rule=build_pv_battery_surrogate_cost_param_block, parameter_block_name="pv_battery_surrogate"
)
def cost_pv_battery_surrogate(blk):
    
    blk.costing_package.has_electricity_generation = True
    global_params = blk.costing_package
    pv_batt_params = blk.costing_package.pv_battery_surrogate
    
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_capital_cost = pyo.Var(
        initialize=0,
        units=blk.costing_package.base_currency,
        bounds=(0, None),
        doc="Direct costs of PV system",
    )

    blk.indirect_capital_cost = pyo.Var(
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

    blk.direct_capital_cost_pv = pyo.Expression()

    blk.direct_capital_cost_constraint = pyo.Constraint(
        expr=blk.direct_capital_cost
        == pyo.units.convert(blk.unit_model.design_size, to_units=pyo.units.watt)
        * (
            pv_batt_params.cost_per_watt_module
            + pv_batt_params.cost_per_watt_inverter
            + pv_batt_params.cost_per_watt_other
        )
        + (
            pyo.units.convert(blk.unit_model.design_size, to_units=pyo.units.watt)
            * (
                pv_batt_params.cost_per_watt_module
                + pv_batt_params.cost_per_watt_inverter
                + pv_batt_params.cost_per_watt_other
            )
        )
        * (1 + pv_batt_params.contingency_frac_direct_capital_cost)  # BUG Check this
    )

    blk.indirect_capital_cost_constraint = pyo.Constraint(
        expr=blk.indirect_capital_cost
        == (
            pyo.units.convert(blk.unit_model.design_size, to_units=pyo.units.watt)
            * pv_batt_params.cost_per_watt_indirect
        )
        + (blk.unit_model.land_req * pv_batt_params.land_cost_per_acre)
    )

    blk.land_cost_constraint = pyo.Constraint(
        expr=blk.land_cost == (blk.unit_model.land_req * pv_batt_params.land_cost_per_acre)
    )

    blk.sales_tax_constraint = pyo.Constraint(
        expr=blk.sales_tax == blk.direct_capital_cost * global_params.sales_tax_frac
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.direct_capital_cost + blk.indirect_capital_cost + blk.sales_tax
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pv_batt_params.fixed_operating_by_capacity
        * pyo.units.convert(blk.unit_model.design_size, to_units=pyo.units.kW)
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == pv_batt_params.variable_operating_by_generation * blk.unit_model.annual_energy
    )

    blk.costing_package.cost_flow(-1 * blk.unit_model.electricity, "electricity")
