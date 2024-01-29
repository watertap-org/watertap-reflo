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

# Costing equations from:
# Kosmadakis G, Papapetrou M, Ortega-Delgado B, Cipollina A, Alarcón-Padilla D-C.
# "Correlations for estimating the specific capital cost of multi-effect distillation plants
#    considering the main design trends and operating conditions"
# doi: 10.1016/j.desal.2018.09.011


def build_med_tvc_surrogate_cost_param_block(blk):

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

    blk.med_sys_A_coeff = pyo.Var(
        initialize=6291,
        units=pyo.units.dimensionless,
        doc="MED system specific capital A coeff",
    )

    blk.med_sys_B_coeff = pyo.Var(
        initialize=-0.135,
        units=pyo.units.dimensionless,
        doc="MED system specific capital B coeff",
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


@register_costing_parameter_block(
    build_rule=build_med_tvc_surrogate_cost_param_block,
    parameter_block_name="med_tvc_surrogate",
)
def cost_med_tvc_surrogate(blk):

    med_tvc_params = blk.costing_package.med_tvc_surrogate
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    med_tvc = blk.unit_model
    feed = med_tvc.feed_props[0]
    dist = med_tvc.distillate_props[0]
    brine = med_tvc.brine_props[0]
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
            med_tvc_params.med_sys_A_coeff
            * blk.capacity**med_tvc_params.med_sys_B_coeff
        )
    )
    blk.membrane_system_cost_constraint = pyo.Constraint(
        expr=blk.membrane_system_cost
        == blk.capacity
        * (blk.med_specific_cost * (1 - med_tvc_params.cost_fraction_evaporator))
    )

    blk.evaporator_system_cost_constraint = pyo.Constraint(
        expr=blk.evaporator_system_cost
        == blk.capacity
        * (
            blk.med_specific_cost
            * (
                med_tvc_params.cost_fraction_evaporator
                * (
                    (
                        med_tvc.specific_area_per_kg_s
                        / med_tvc_params.heat_exchanger_ref_area
                    )
                    ** med_tvc_params.heat_exchanger_exp
                )
            )
        )
    )
    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.membrane_system_cost + blk.evaporator_system_cost
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.annual_dist_production
        * (
            med_tvc_params.cost_chemicals_per_vol_dist
            + med_tvc_params.cost_labor_per_vol_dist
            + med_tvc_params.cost_misc_per_vol_dist
        )
        + blk.capital_cost
        * (
            med_tvc_params.cost_fraction_maintenance
            + med_tvc_params.cost_fraction_insurance
        )
        + pyo.units.convert(
            brine.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.year
        )
        * med_tvc_params.cost_disposal_per_vol_brine
    )

    blk.heat_flow = pyo.Expression(
        expr=med_tvc.specific_energy_consumption_thermal
        * pyo.units.convert(blk.capacity, to_units=pyo.units.m**3 / pyo.units.hr)
    )
    blk.electricity_flow = pyo.Expression(
        expr=med_tvc_params.specific_energy_consumption_electric
        * pyo.units.convert(blk.capacity, to_units=pyo.units.m**3 / pyo.units.hr)
    )

    blk.costing_package.cost_flow(blk.heat_flow, "heat")
    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
