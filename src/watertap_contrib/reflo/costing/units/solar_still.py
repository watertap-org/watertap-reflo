#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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


def build_solar_still_cost_param_block(blk):

    # Capital cost of solar still is function of number of stills
    # C = A * number_stills **b

    blk.cost_per_still_A_param = pyo.Var(
        initialize=300.65,
        units=pyo.units.USD_2020,
        bounds=(0, None),
        doc="Cost per still equation parameter A",
    )

    blk.cost_per_still_b_param = pyo.Var(
        initialize=-0.199,
        units=pyo.units.dimensionless,
        bounds=(None, 0),
        doc="Cost per still equation parameter b",
    )

    blk.cost_saltwater_pump = pyo.Var(
        initialize=750,
        units=pyo.units.USD_2020,
        bounds=(0, None),
        doc="Cost of single saltwater pump",
    )

    blk.cost_freshwater_pump = pyo.Var(
        initialize=150,
        units=pyo.units.USD_2020,
        bounds=(0, None),
        doc="Cost of single freshwater pump",
    )

    blk.cost_piping = pyo.Var(
        initialize=1.5,
        bounds=(0, None),
        units=pyo.units.USD_2020 / pyo.units.ft,
        doc="CPVC pipe per linear foot",
    )

    blk.pipe_length_param = pyo.Var(
        initialize=1.5,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Pipe length equation parameter",
    )
    blk.cost_feed_tank_base = pyo.Var(
        initialize=1601.1,
        bounds=(0, None),
        units=pyo.units.USD_2020,
        doc="Feed tank capital equation base",
    )

    blk.cost_feed_tank_exp = pyo.Var(
        initialize=0.6373,
        units=pyo.units.dimensionless,
        doc="Feed tank capital equation exponent",
    )

    blk.cost_dist_tank_base = pyo.Var(
        initialize=1680.8,
        units=pyo.units.USD_2020,
        doc="Distillate tank capital equation base",
    )

    blk.cost_dist_tank_exp = pyo.Var(
        initialize=0.6085,
        units=pyo.units.dimensionless,
        doc="Distillate tank capital equation exponent",
    )

    blk.cost_underground_tank_base = pyo.Var(
        initialize=411.78,
        units=pyo.units.USD_2020,
        doc="Underground tank capital equation base",
    )

    blk.cost_underground_tank_exp = pyo.Var(
        initialize=0.8693,
        units=pyo.units.dimensionless,
        doc="Underground tank capital equation exponent",
    )

    blk.cost_excavation_base = pyo.Var(
        initialize=646.33,
        bounds=(0, None),
        units=pyo.units.USD_2020,
        doc="Cost excavation base",
    )

    blk.cost_excavation_exp = pyo.Var(
        initialize=0.624,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Cost excavation exponent",
    )

    blk.fixed_opex_factor = pyo.Var(
        initialize=0.035,
        bounds=(0, None),
        units=pyo.units.year**-1,
        doc="Factor for calculating fixed operating costs as a fraction of CAPEX",
    )
    blk.fix_all_vars()


@register_costing_parameter_block(
    build_rule=build_solar_still_cost_param_block,
    parameter_block_name="solar_still",
)
def cost_solar_still(blk):
    ss_params = blk.costing_package.solar_still
    ss = blk.unit_model
    base_currency = blk.config.flowsheet_costing_block.base_currency
    base_period = blk.config.flowsheet_costing_block.base_period
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    dimensionless_water_yield = pyo.units.convert(
        ss.annual_water_yield * pyo.units.m**-1,
        to_units=pyo.units.dimensionless,
    )

    blk.pressure_drop = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.bar,
        doc="Pressure drop for pumps",
    )

    blk.max_sw_flow = pyo.Param(
        initialize=1728,
        mutable=True,
        units=pyo.units.m**3 / pyo.units.day,
        doc="Saltwater pump flow rate",
    )

    blk.max_fw_flow = pyo.Param(
        initialize=86.4,
        mutable=True,
        units=pyo.units.m**3 / pyo.units.day,
        doc="Freshwater pump flow rate",
    )

    blk.sw_pump_efficiency = pyo.Param(
        initialize=0.8,
        mutable=True,
        units=pyo.units.dimensionless,
        doc="Saltwater pump efficiency",
    )

    blk.fw_pump_efficiency = pyo.Param(
        initialize=0.8,
        mutable=True,
        units=pyo.units.dimensionless,
        doc="Freshwater pump efficiency",
    )

    blk.number_sw_pumps = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.dimensionless,
        doc="Number saltwater pumps",
    )

    blk.number_fw_pumps = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.dimensionless,
        doc="Number freshwater pumps",
    )

    blk.capital_cost_per_still = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Cost per solar still",
    )

    blk.sw_pump_power = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.kilowatt,
        doc="Saltwater pumping power per pump",
    )

    blk.fw_pump_power = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.kilowatt,
        doc="Freshwater pumping power per pump",
    )

    blk.pumping_power_required = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.kilowatt,
        doc="Total pumping power required",
    )

    blk.length_piping = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.ft,
        doc="Length of CPVC piping",
    )

    blk.capital_cost_solar_still = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of solar stills",
    )

    blk.capital_cost_sw_pumps = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of saltwater pumps",
    )

    blk.capital_cost_fw_pumps = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of freshwater pumps",
    )

    blk.capital_cost_feed_tank = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of feed tank",
    )

    blk.capital_cost_distillate_tank = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of distillate tank",
    )

    blk.capital_cost_underground_tank = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of underground tank",
    )

    blk.capital_cost_excavation = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of excavation",
    )

    blk.capital_cost_piping = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=base_currency,
        doc="Capital cost of piping",
    )

    blk.capital_cost_per_still_constraint = pyo.Constraint(
        expr=blk.capital_cost_per_still
        == pyo.units.convert(
            ss_params.cost_per_still_A_param
            * ss.number_stills**ss_params.cost_per_still_b_param,
            to_units=base_currency,
        )
    )

    blk.sw_pump_power_constraint = pyo.Constraint(
        expr=blk.sw_pump_power
        == pyo.units.convert(
            (blk.pressure_drop * blk.max_sw_flow) / blk.sw_pump_efficiency,
            to_units=pyo.units.kilowatt,
        )
    )

    blk.fw_pump_power_constraint = pyo.Constraint(
        expr=blk.fw_pump_power
        == pyo.units.convert(
            (blk.pressure_drop * blk.max_fw_flow) / blk.fw_pump_efficiency,
            to_units=pyo.units.kilowatt,
        )
    )

    blk.pumping_power_required_constraint = pyo.Constraint(
        expr=blk.pumping_power_required
        == blk.number_sw_pumps
        * pyo.units.convert(blk.sw_pump_power, to_units=pyo.units.kilowatt)
        + blk.number_fw_pumps
        * pyo.units.convert(blk.fw_pump_power, to_units=pyo.units.kilowatt)
    )

    blk.length_piping_constraint = pyo.Constraint(
        expr=blk.length_piping
        == pyo.units.convert(
            ss_params.pipe_length_param * ss.number_stills * ss.length_basin,
            to_units=pyo.units.ft,
        )
    )

    capital_cost_expr = 0

    blk.capital_cost_solar_still_constraint = pyo.Constraint(
        expr=blk.capital_cost_solar_still
        == pyo.units.convert(
            ss.number_stills * blk.capital_cost_per_still, to_units=base_currency
        )
    )

    capital_cost_expr += blk.capital_cost_solar_still

    blk.capital_cost_sw_pumps_constraint = pyo.Constraint(
        expr=blk.capital_cost_sw_pumps
        == pyo.units.convert(
            blk.number_sw_pumps * ss_params.cost_saltwater_pump, to_units=base_currency
        )
    )

    capital_cost_expr += blk.capital_cost_sw_pumps

    blk.capital_cost_fw_pumps_constraint = pyo.Constraint(
        expr=blk.capital_cost_fw_pumps
        == pyo.units.convert(
            blk.number_fw_pumps * ss_params.cost_freshwater_pump, to_units=base_currency
        )
    )

    capital_cost_expr += blk.capital_cost_fw_pumps

    blk.capital_cost_feed_tank_constraint = pyo.Constraint(
        expr=blk.capital_cost_feed_tank
        == pyo.units.convert(
            ss_params.cost_feed_tank_base
            * dimensionless_water_yield**ss_params.cost_feed_tank_exp,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_feed_tank

    blk.capital_cost_distillate_tank_constraint = pyo.Constraint(
        expr=blk.capital_cost_distillate_tank
        == pyo.units.convert(
            ss_params.cost_dist_tank_base
            * dimensionless_water_yield**ss_params.cost_dist_tank_exp,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_distillate_tank

    blk.capital_cost_underground_tank_constraint = pyo.Constraint(
        expr=blk.capital_cost_underground_tank
        == pyo.units.convert(
            ss_params.cost_underground_tank_base
            * dimensionless_water_yield**ss_params.cost_underground_tank_exp,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_underground_tank

    blk.capital_cost_excavation_constraint = pyo.Constraint(
        expr=blk.capital_cost_excavation
        == pyo.units.convert(
            ss_params.cost_excavation_base
            * dimensionless_water_yield**ss_params.cost_excavation_exp,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_excavation

    blk.capital_cost_piping_constraint = pyo.Constraint(
        expr=blk.capital_cost_piping
        == pyo.units.convert(
            blk.length_piping * ss_params.cost_piping,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_piping

    blk.costing_package.add_cost_factor(blk, None)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(capital_cost_expr, to_units=base_currency)
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.capital_cost * ss_params.fixed_opex_factor,
            to_units=base_currency / base_period,
        )
    )

    blk.costing_package.cost_flow(blk.pumping_power_required, "electricity")
