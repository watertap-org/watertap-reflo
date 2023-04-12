import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)


def build_trough_surrogate_cost_param_block(blk):

    costing = blk.parent_block()

    blk.base_storage_hours = pyo.Param(
        mutable=True,
        initialize=6,
        units=pyo.units.hour,
        doc="Hours of storage that capital cost per capacity is based on",
    )

    blk.cost_per_capacity_capital = pyo.Var(
        initialize=560,
        units=costing.base_currency / pyo.units.kW,
        bounds=(0, None),
        doc="Cost per kW (thermal) for the trough plant, assuming six hours storage",
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

    blk.tax_frac_direct_cost = pyo.Var(
        initialize=0.05,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Fraction of direct costs applicable for sales tax",
    )

    blk.fixed_operating_by_capacity = pyo.Var(
        initialize=8,
        units=costing.base_currency / (pyo.units.kW * costing.base_period),
        bounds=(0, None),
        doc="Fixed operating cost of trough plant per kW capacity",
    )

    blk.variable_operating_by_generation = pyo.Var(
        initialize=0.001,
        units=costing.base_currency / (pyo.units.kWh * costing.base_period),
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

    blk.direct_cost_constraint = pyo.Constraint(
        expr=blk.direct_cost
        == (
            pyo.units.convert(trough.heat_load, to_units=pyo.units.kW)
            * (
                trough_params.cost_per_capacity_capital
                - trough_params.base_storage_hours
                * trough_params.cost_per_storage_capital
                + trough.hours_storage * trough_params.cost_per_storage_capital
            )
        )
        * (1 + trough_params.contingency_frac_direct_cost)
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.direct_cost
        * (1 + trough_params.tax_frac_direct_cost * global_params.sales_tax_frac)
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == trough_params.fixed_operating_by_capacity
        * pyo.units.convert(trough.heat_load, to_units=pyo.units.kW)
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
