import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)


def build_flat_plate_surrogate_cost_param_block(blk):

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

    blk.tax_frac_direct_cost = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs applicable for sales tax",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=16,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of flat plate plant per kW capacity",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0,
        units=costing.base_currency / (pyo.units.MWh * costing.base_period),
        bounds=(0, None),
        doc="Variable operating cost of flat plate plant per MWh generated",
    )


@register_costing_parameter_block(
    build_rule=build_flat_plate_surrogate_cost_param_block,
    parameter_block_name="flat_plate_surrogate",
)
def cost_flat_plate(blk):

    global_params = blk.costing_package
    flat_plate_params = blk.costing_package.flat_plate_surrogate
    flat_plate = blk.unit_model
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.direct_cost = pyo.Var(
        initialize=1e4,
        units=blk.config.flowsheet_costing_block.base_currency,
        bounds=(0, None),
        doc="Direct cost of flat plate plant",
    )

    blk.indirect_cost = pyo.Var(
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

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == (
            (
                flat_plate.collector_area_total
                * flat_plate_params.cost_per_area_collector
            )
            + (flat_plate.storage_volume * flat_plate_params.cost_per_volume_storage)
        )
        * (1 + flat_plate_params.contingency_frac_direct_cost)
    )

    blk.indirect_cost_constraint = pyo.Constraint(
        expr=blk.indirect_cost
        == blk.direct_cost * flat_plate_params.indirect_frac_direct_cost
    )

    blk.sales_tax_constraint = pyo.Constraint(
        expr=blk.sales_tax == blk.direct_cost * global_params.sales_tax_frac
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.direct_cost + blk.indirect_cost + blk.sales_tax
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == flat_plate_params.fixed_operating_by_capacity
        * pyo.units.convert(flat_plate.heat_load, to_units=pyo.units.kW)
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == flat_plate_params.variable_operating_by_generation * flat_plate.heat_annual
    )

    # register the flows, e.g.,:
    blk.costing_package.cost_flow(
        flat_plate.electricity,
        "electricity",
    )
    blk.costing_package.cost_flow(
        flat_plate.heat,
        "heat",
    )
