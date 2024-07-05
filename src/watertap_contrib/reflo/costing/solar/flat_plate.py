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
)


def build_flat_plate_cost_param_block(blk):

    costing = blk.parent_block()

    blk.cost_per_area_collector = pyo.Var(
        initialize=600,
        units=costing.base_currency / pyo.units.m**2,
        bounds=(0, None),
        doc="Cost per area for solar collector",
    )

    blk.cost_per_volume_storage = pyo.Var(
        initialize=2000,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost per volume for thermal storage",
    )

    blk.land_cost_per_acre = pyo.Var(
        initialize=4000,
        units=costing.base_currency / pyo.units.acre,
        bounds=(0, None),
        doc="Land cost per acre required",
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
        initialize=16,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of flat plate plant per kW capacity",
    )


@register_costing_parameter_block(
    build_rule=build_flat_plate_cost_param_block,
    parameter_block_name="flat_plate",
)
def cost_flat_plate(blk):

    global_params = blk.costing_package
    flat_plate_params = blk.costing_package.flat_plate
    flat_plate = blk.unit_model
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_capital_cost = pyo.Var(
        initialize=1e4,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Direct cost of flat plate plant",
    )

    blk.indirect_capital_cost = pyo.Var(
        initialize=1e4,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Indirect costs of flat plate system",
    )

    blk.sales_tax = pyo.Var(
        initialize=1e2,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Sales tax for flat plate system",
    )

    if flat_plate.config.solar_model_type == "surrogate":
        blk.direct_cost_constraint = pyo.Constraint(
            expr=blk.direct_capital_cost
            == (
                (
                    flat_plate.collector_area_total
                    * flat_plate_params.cost_per_area_collector
                )
                + (
                    flat_plate.storage_volume
                    * flat_plate_params.cost_per_volume_storage
                )
            )
            * (1 + flat_plate_params.contingency_frac_direct_cost)
        )

    else:
        blk.direct_cost_constraint = pyo.Constraint(
            expr=blk.direct_capital_cost
            == (
                (
                    flat_plate.collector_area_total
                    * flat_plate_params.cost_per_area_collector
                )
            )
            * (1 + flat_plate_params.contingency_frac_direct_cost)
        )

    blk.land_area = pyo.Expression(
        expr=pyo.units.convert(flat_plate.collector_area_total, to_units=pyo.units.acre)
    )

    blk.indirect_cost_constraint = pyo.Constraint(
        expr=blk.indirect_capital_cost
        == blk.direct_capital_cost * flat_plate_params.indirect_frac_direct_cost
        + blk.land_area * flat_plate_params.land_cost_per_acre
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
        == flat_plate_params.fixed_operating_by_capacity
        * pyo.units.convert(flat_plate.heat_load, to_units=pyo.units.kW)
    )

    blk.costing_package.cost_flow(
        flat_plate.electricity,
        "electricity",
    )
    blk.costing_package.cost_flow(
        -1 * flat_plate.heat,
        "heat",
    )
