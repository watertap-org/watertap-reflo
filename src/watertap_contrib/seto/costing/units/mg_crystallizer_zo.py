import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# TODO: add reference


def build_mg_crystallizer_zo_cost_param_block(blk):
    costing = blk.parent_block()

    blk.small_capital_cost_coefficient = pyo.Var(
        initialize=4000,
        units=costing.base_currency / (pyo.units.m**3 / pyo.units.day),
        bounds=(0, None),
        doc="Capital cost coefficient of small size crystallizer (<20 m3/day)",
    )

    blk.large_capital_cost_coefficient1 = pyo.Var(
        initialize=28531,
        units=costing.base_currency / (pyo.units.m**3 / pyo.units.day),
        bounds=(0, None),
        doc="Capital cost coefficient 1 of large size crystallizer (>20 m3/day)",
    )

    blk.large_capital_cost_coefficient2 = pyo.Var(
        initialize=18824,
        units=costing.base_currency / (pyo.units.m**3 / pyo.units.day),
        bounds=(0, None),
        doc="Capital cost coefficient 2 of small size crystallizer (>20 m3/day)",
    )

    blk.unit_operating_cost = pyo.Var(
        initialize=0.05,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Unit operating cost of 1 m3 of feed water",
    )

    blk.calcium_hydroxide_dosing_cost = pyo.Var(
        initialize=0.080,
        units=costing.base_currency / pyo.units.kg,
        bounds=(0, None),
        doc="Ca(OH)2 dosing cost (usually 80 - 300 $/ton)",
    )

    blk.hydromagnesite_recovery_revenue = pyo.Var(
        initialize=0.5,
        units=costing.base_currency / pyo.units.kg,
        bounds=(0, None),
        doc="Recovery revenue of hydromagnesite (usually 0.26 - 0.71 $/kg)",
    )

    blk.calcite_recovery_revenue = pyo.Var(
        initialize=0.4,
        units=costing.base_currency / pyo.units.kg,
        bounds=(0, None),
        doc="Recovery revenue of calcite (usually 0.30 - 0.42 $/kg)",
    )


@register_costing_parameter_block(
    build_rule=build_mg_crystallizer_zo_cost_param_block,
    parameter_block_name="mg_crystallizer",
)
def cost_mg_crystallizer(blk):

    mg_cryst_params = blk.costing_package.mg_crystallizer
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.include_metal_recovery = pyo.Var(
        initialize=0,
        units=pyo.units.dimensionless,
        doc="Include metal recovery (= 1), or not (=0)",
    )

    mg_cryst = blk.unit_model
    feed = mg_cryst.properties_in[0]

    feed_water_flow_rate = pyo.units.convert(
        feed.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.day
    )

    crystallizer_capital_cost = pyo.Expr_if(
        feed_water_flow_rate <= 20,
        mg_cryst_params.small_capital_cost_coefficient * feed_water_flow_rate,
        mg_cryst_params.large_capital_cost_coefficient1 * pyo.log(feed_water_flow_rate)
        + mg_cryst_params.large_capital_cost_coefficient2,
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == crystallizer_capital_cost
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == mg_cryst_params.unit_operating_cost
        * pyo.units.convert(
            feed.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.year
        )
        + mg_cryst_params.calcium_hydroxide_dosing_cost
        * pyo.units.convert(
            mg_cryst.calcium_hydroxide_mass, to_units=pyo.units.kg / pyo.units.year
        )
        - blk.include_metal_recovery
        * (
            mg_cryst_params.hydromagnesite_recovery_revenue
            * pyo.units.convert(
                mg_cryst.hydromagnesite_precipitated,
                to_units=pyo.units.kg / pyo.units.year,
            )
            + mg_cryst_params.calcite_recovery_revenue
            * pyo.units.convert(
                mg_cryst.calcite_precipitated_feed + mg_cryst.calcite_precipitated_dose,
                to_units=pyo.units.kg / pyo.units.year,
            )
        )
    )
