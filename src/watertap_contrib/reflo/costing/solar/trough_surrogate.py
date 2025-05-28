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


def build_trough_surrogate_cost_param_block(blk):

    costing = blk.parent_block()

    blk.cost_per_land_area = pyo.Var(
        initialize=4000,
        units=costing.base_currency / pyo.units.acre,
        bounds=(0, None),
        doc="Cost per acre of land",
    )

    blk.cost_per_total_aperture_area = pyo.Var(
        initialize=297,
        units=costing.base_currency / pyo.units.m**2,
        bounds=(0, None),
        doc="Cost per m2 of total aperture area (includes site improvement 16 $/m2, solar field 297 $/m2, HTF system 60 $/m2)",
    )

    blk.cost_per_heat_sink = pyo.Var(
        initialize=120,
        units=costing.base_currency / pyo.units.kW,
        bounds=(0, None),
        doc="Cost for expenses related to installation of the heat sink, including labor and equipment per kWh (thermal) heat load",
    )

    blk.cost_per_balance_of_plant = pyo.Var(
        initialize=90,
        units=costing.base_currency / pyo.units.kW,
        bounds=(0, None),
        doc="Cost per thermal kilowatt of heat sink capacity for expenses related to installation of the heat sink, including labor and equipment",
    )

    blk.cost_per_storage_capital = pyo.Var(
        initialize=62,
        units=costing.base_currency / pyo.units.kWh,
        bounds=(0, None),
        doc="Cost per kWh (thermal) for the trough plant thermal storage",
    )

    blk.contingency_frac_direct_cost = pyo.Var(
        initialize=0.07,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs for contingency",
    )

    blk.indirect_frac_direct_cost = pyo.Var(
        initialize=0.11,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs for indirect costs associated with engineer-procure-construction (EPC)",
    )

    blk.tax_frac_direct_cost = pyo.Var(
        initialize=0.05,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs applicable for sales tax",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=103758,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Fixed operating cost of trough plant in SAM. Not a function of electricity generated",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0.002,
        units=costing.base_currency / (pyo.units.kWh),
        bounds=(0, None),
        doc="Variable operating cost of trough plant per kWh generated",
    )


@register_costing_parameter_block(
    build_rule=build_trough_surrogate_cost_param_block,
    parameter_block_name="trough_surrogate",
)
def cost_trough_surrogate(blk):

    global_params = blk.costing_package
    trough_params = blk.costing_package.trough_surrogate
    trough = blk.unit_model
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=0,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Direct cost of trough plant",
    )

    blk.indirect_cost = pyo.Var(
        initialize=0,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Inirect cost of trough plant",
    )

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == (
            trough.total_aperture_area * trough_params.cost_per_total_aperture_area
            + pyo.units.convert(trough.heat_load, to_units=pyo.units.kW)
            * trough.hours_storage
            * trough_params.cost_per_storage_capital
            + pyo.units.convert(trough.heat_load, to_units=pyo.units.kW)
            * (
                trough_params.cost_per_heat_sink
                + trough_params.cost_per_balance_of_plant
            )
        )
        * (1 + trough_params.contingency_frac_direct_cost)
    )

    blk.indirect_cost_constraint = pyo.Constraint(
        expr=blk.indirect_cost
        == (
            blk.direct_cost * trough_params.indirect_frac_direct_cost
            + trough.land_area * trough_params.cost_per_land_area
        )
    )

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.indirect_cost + blk.direct_cost
        * (1 + trough_params.tax_frac_direct_cost * global_params.sales_tax_frac)
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == trough_params.fixed_operating_by_capacity
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == trough_params.variable_operating_by_generation * trough.heat_annual
    )

    blk.costing_package.cost_flow(
        trough.electricity,
        "electricity",
    )
    blk.costing_package.cost_flow(
        -1 * trough.heat,
        "heat",
    )
