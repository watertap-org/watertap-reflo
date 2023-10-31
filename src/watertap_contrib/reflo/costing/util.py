import pyomo.environ as pyo


def make_capital_cost_var(blk):
    blk.capital_cost = pyo.Var(
        initialize=1e5,
        # domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit capital cost",
    )


def make_variable_operating_cost_var(blk):
    blk.variable_operating_cost = pyo.Var(
        initialize=1e5,
        # domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit variable operating cost",
    )


def make_fixed_operating_cost_var(blk):
    blk.fixed_operating_cost = pyo.Var(
        initialize=1e5,
        # domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit fixed operating cost",
    )
