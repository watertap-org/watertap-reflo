import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing equations from:
# Kosmadakis G, Papapetrou M, Ortega-Delgado B, Cipollina A, Alarcón-Padilla D-C.
# "Correlations for estimating the specific capital cost of multi-effect distillation plants
#    considering the main design trends and operating conditions"
# doi: 10.1016/j.desal.2018.09.011


def build_air_stripping_cost_param_block(blk):

    costing = blk.parent_block()


@register_costing_parameter_block(
    build_rule=build_air_stripping_cost_param_block,
    parameter_block_name="air_stripping",
)
def cost_air_stripping(blk):

    ax = blk.costing_package.air_stripping
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.tower_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Aluminum tower cost",
    )

    blk.port_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Aluminum access port system cost",
    )

    blk.piping_liq_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Liquid inlet/outlet piping cost",
    )

    blk.piping_air_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Air inlet piping cost",
    )

    blk.tray_ring_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Tray ring cost",
    )

    blk.tower_internals_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Tray ring cost",
    )

    blk.packing_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Packing material cost",
    )

    blk.mist_eliminator_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Mist eliminator cost",
    )

    blk.pump_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Water pump cost",
    )

    blk.blower_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Air blower cost",
    )
    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
