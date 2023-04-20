import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.seto.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing equations from:
# Kosmadakis G, Papapetrou M, Ortega-Delgado B, Cipollina A, Alarcón-Padilla D-C.
# "Correlations for estimating the specific capital cost of multi-effect distillation plants
#    considering the main design trends and operating conditions"
# doi: 10.1016/j.desal.2018.09.011


def build_lt_med_surrogate_cost_param_block(blk):

    costing = blk.parent_block()

    blk.cost_fraction_evaporator = pyo.Var(
        initialize=0.4,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Cost fraction of the evaporator",
    )

    blk.cost_fraction_maintenance = pyo.Var(
        initialize=0.02,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Fraction of capital cost for maintenance",
    )

    blk.cost_fraction_insurance = pyo.Var(
        initialize=0.005,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Fraction of capital cost for insurance",
    )

    blk.cost_chemicals_per_vol_dist = pyo.Var(
        initialize=0.04,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of chemicals per m3 distillate",
    )

    blk.cost_labor_per_vol_dist = pyo.Var(
        initialize=0.033,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of labor per m3 distillate",
    )

    blk.cost_misc_per_vol_dist = pyo.Var(
        initialize=0.033,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of labor per m3 distillate",
    )

    blk.cost_disposal_per_vol_brine = pyo.Var(
        initialize=0.02,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of disposal per m3 brine",
    )

    blk.specific_energy_consumption_electric = pyo.Var(
        initialize=1.5,
        units=pyo.units.kWh / pyo.units.m**3,
        bounds=(0, None),
        doc="Specific electric energy consumption",
    )

    blk.hours_thermal_storage = pyo.Var(
        initialize=0,
        units=pyo.units.hr,
        doc="Hours of thermal storage required",
    )

    blk.med_sys_A_coeff = pyo.Var(
        initialize=6291,
        units=pyo.units.dimensionless,
        doc="LT-MED system specific capital A coeff",
    )

    blk.med_sys_B_coeff = pyo.Var(
        initialize=-0.135,
        units=pyo.units.dimensionless,
        doc="LT-MED system specific capital B coeff",
    )

    blk.heat_exchanger_ref_area = pyo.Var(
        initialize=302.01,
        units=pyo.units.m**2 / (pyo.units.kg / pyo.units.s),
        doc="Specific heat exchanger area from reference plant",
    )

    blk.heat_exchanger_exp = pyo.Var(
        initialize=0.8,
        units=pyo.units.dimensionless,
        doc="Exponent for specific heat exchanger cost equation",
    )

    blk.fix_all_vars()


@register_costing_parameter_block(
    build_rule=build_lt_med_surrogate_cost_param_block,
    parameter_block_name="lt_med_surrogate",
)
def cost_lt_med_surrogate(blk):

    lt_med_params = blk.costing_package.lt_med_surrogate
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    lt_med = blk.unit_model
    feed = lt_med.feed_props[0]
    dist = lt_med.distillate_props[0]
    brine = lt_med.brine_props[0]
    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.membrane_system_cost = pyo.Var(
        initialize=100,
        units=base_currency,
        doc="Membrane system cost",
    )

    blk.evaporator_system_cost = pyo.Var(
        initialize=100,
        units=base_currency,
        doc="Evaporator system cost",
    )

    blk.med_specific_cost = pyo.Var(
        initialize=100,
        units=pyo.units.USD_2018 / (pyo.units.m**3 / pyo.units.day),
        doc="MED system cost per m3/day distillate",
    )

    blk.capacity = pyo.units.convert(
        dist.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.day
    )

    blk.annual_dist_production = pyo.units.convert(
        dist.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.year
    )
    blk.med_specific_cost_constraint = pyo.Constraint(
        expr=blk.med_specific_cost
        == (
            lt_med_params.med_sys_A_coeff
            * blk.capacity**lt_med_params.med_sys_B_coeff
        )
    )
    blk.membrane_system_cost_constraint = pyo.Constraint(
        expr=blk.membrane_system_cost
        == blk.capacity
        * (blk.med_specific_cost * (1 - lt_med_params.cost_fraction_evaporator))
    )

    blk.evaporator_system_cost_constraint = pyo.Constraint(
        expr=blk.evaporator_system_cost
        == blk.capacity
        * (
            blk.med_specific_cost
            * (
                lt_med_params.cost_fraction_evaporator
                * (
                    (
                        lt_med.specific_area_per_kg_s
                        / lt_med_params.heat_exchanger_ref_area
                    )
                    ** lt_med_params.heat_exchanger_exp
                )
            )
        )
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.membrane_system_cost + blk.evaporator_system_cost
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.annual_dist_production
        * (
            lt_med_params.cost_chemicals_per_vol_dist
            + lt_med_params.cost_labor_per_vol_dist
            + lt_med_params.cost_misc_per_vol_dist
        )
        + blk.capital_cost
        * (
            lt_med_params.cost_fraction_maintenance
            + lt_med_params.cost_fraction_insurance
        )
        + pyo.units.convert(
            brine.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.year
        )
        * lt_med_params.cost_disposal_per_vol_brine
    )

    blk.heat_flow = pyo.Expression(
        expr=lt_med.specific_energy_consumption_thermal
        * pyo.units.convert(blk.capacity, to_units=pyo.units.m**3 / pyo.units.hr)
    )
    blk.electricity_flow = pyo.Expression(
        expr=lt_med_params.specific_energy_consumption_electric
        * pyo.units.convert(blk.capacity, to_units=pyo.units.m**3 / pyo.units.hr)
    )

    blk.costing_package.cost_flow(blk.heat_flow, "heat")
    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
