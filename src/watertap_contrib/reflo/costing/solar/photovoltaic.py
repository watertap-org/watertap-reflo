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

import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)


def build_photovoltaic_cost_param_block(blk):

    costing = blk.parent_block()

    blk.cost_per_watt_module = pyo.Var(
        initialize=0.41,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for solar module",
    )

    blk.cost_per_watt_inverter = pyo.Var(
        initialize=0.05,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for inverter",
    )

    blk.cost_per_watt_other = pyo.Var(
        initialize=0.1,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for other equipment, installation, and margin/overhead",
    )

    blk.cost_per_watt_indirect = pyo.Var(
        initialize=0.13,
        units=costing.base_currency / pyo.units.watt,
        bounds=(0, None),
        doc="Cost per watt for permitting, environmental studies, engineering, land prep, and grid interconnection",
    )

    blk.land_cost_per_acre = pyo.Var(
        initialize=4000,
        units=costing.base_currency / pyo.units.acre,
        bounds=(0, None),
        doc="Land cost per acre required",
    )

    blk.contingency_frac_direct_cost = pyo.Var(
        initialize=0.03,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs for contingency",
    )

    blk.tax_frac_direct_cost = pyo.Var(
        initialize=0.05,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs for sales tax",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=0,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of PV system per kW generated",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0,
        units=costing.base_currency / (pyo.units.MWh * costing.base_period),
        bounds=(0, None),
        doc="Annual operating cost of PV system per MWh generated",
    )


@register_costing_parameter_block(
    build_rule=build_photovoltaic_cost_param_block, parameter_block_name="photovoltaic"
)
def cost_pv(blk):

    global_params = blk.costing_package
    pv_params = blk.costing_package.photovoltaic
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=0,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Direct costs of PV system",
    )

    blk.indirect_cost = pyo.Var(
        initialize=0,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Indirect costs of PV system",
    )

    blk.sales_tax = pyo.Var(
        initialize=0,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Sales tax for PV system",
    )

    blk.system_capacity = pyo.Var(
        initialize=0,
        units=pyo.units.watt,
        bounds=(0, None),
        doc="DC system capacity for PV system",
    )

    blk.land_area = pyo.Var(
        initialize=0,
        units=pyo.units.acre,
        bounds=(0, None),
        doc="Land area required for PV system",
    )

    blk.annual_generation = pyo.Var(
        initialize=0,
        units=pyo.units.MWh,
        bounds=(0, None),
        doc="Annual electricity generation of PV system",
    )

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == blk.system_capacity
        * (
            pv_params.cost_per_watt_module
            + pv_params.cost_per_watt_inverter
            + pv_params.cost_per_watt_other
        )
        + (
            blk.system_capacity
            * (
                pv_params.cost_per_watt_module
                + pv_params.cost_per_watt_inverter
                + pv_params.cost_per_watt_other
            )
        )
        * pv_params.contingency_frac_direct_cost
    )

    blk.indirect_cost_constraint = pyo.Constraint(
        expr=blk.indirect_cost
        == (blk.system_capacity * pv_params.cost_per_watt_indirect)
        + (blk.land_area * pv_params.land_cost_per_acre)
    )

    blk.sales_tax_constraint = pyo.Constraint(
        expr=blk.sales_tax == blk.direct_cost * global_params.sales_tax_frac
    )

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.direct_cost + blk.indirect_cost + blk.sales_tax
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pv_params.fixed_operating_by_capacity
        * pyo.units.convert(blk.system_capacity, to_units=pyo.units.kW)
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == pv_params.variable_operating_by_generation * blk.annual_generation
    )

    blk.costing_package.cost_flow(blk.unit_model.electricity, "electricity")
