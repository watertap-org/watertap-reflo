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


def build_solar_still_cost_param_block(blk):

    blk.cost_solar_still = pyo.Var(
        initialize=50,
        units=pyo.units.USD_2020 / pyo.units.m**2,
        bounds=(0, None),
        doc="Cost of solar still per unit area",
    )

    blk.cost_saltwater_pump = pyo.Var(
        initialize=750, units=pyo.units.USD_2020, doc="Cost of single saltwater pump"
    )
    blk.cost_freshwater_pump = pyo.Var(
        initialize=150, units=pyo.units.USD_2020, doc="Cost of single freshwater pump"
    )
    blk.cost_piping = pyo.Var(
        initialize=1.5,
        units=pyo.units.USD_2020 / pyo.units.ft,
        doc="CPVC pipe per linear foot",
    )
    blk.cost_feed_tank_base = pyo.Var(
        initialize=1601, units=pyo.units.USD_2020, doc="Feed tank capital equation base"
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
        initialize=646.33, units=pyo.units.USD_2020, doc="Cost excavation base"
    )
    blk.cost_excavation_exp = pyo.Var(
        initialize=0.624, units=pyo.units.dimensionless, doc="Cost excavation exponent"
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
    make_capital_cost_var(blk)
    # make_fixed_operating_cost_var(blk)
    blk.water_yield_total = (ss.water_yield * ss.total_area) / ss.properties_in[
        0
    ].dens_mass_phase["Liq"]

    dimensionless_water_yield_units = pyo.units.year * pyo.units.m**-3
    dimensionless_water_yield = pyo.units.convert(
        blk.water_yield_total * dimensionless_water_yield_units,
        to_units=pyo.units.dimensionless,
    )
    # dimensionless_water_yield = pyo.units.convert(
    #     ss.water_yield * dimensionless_water_yield_units,
    #     to_units=pyo.units.dimensionless,
    # )

    blk.pressure_drop = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.bar,
        doc="Pressure drop for pumps",
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

    blk.max_sw_flow = pyo.Param(
        initialize=72,
        mutable=True,
        units=pyo.units.m**3 / pyo.units.hr,
        doc="Max volumetric flow rate for saltwater pump",
    )

    blk.max_fw_flow = pyo.Param(
        initialize=3.6,
        mutable=True,
        units=pyo.units.m**3 / pyo.units.hr,
        doc="Max volumetric flow rate for freshwater pump",
    )

    blk.sw_pump_power = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.hp,
        doc="Saltwater pump power reqiured",
    )

    blk.fw_pump_power = pyo.Param(
        initialize=1,
        mutable=True,
        units=pyo.units.hp,
        doc="Freshwater pump power required",
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

    blk.number_sw_pumps = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Number saltwater pumps",
    )

    blk.number_fw_pumps = pyo.Var(
        initialize=1,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Number freshwater pumps",
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

    # blk.sw_pump_power_constraint = pyo.Constraint(
    #     expr=blk.sw_pump_power
    #     == pyo.units.convert(
    #         (blk.pressure_drop * blk.max_sw_flow) / blk.sw_pump_efficiency,
    #         to_units=pyo.units.kilowatt,
    #     )
    # )

    # blk.fw_pump_power_constraint = pyo.Constraint(
    #     expr=blk.fw_pump_power
    #     == pyo.units.convert(
    #         (blk.pressure_drop * blk.max_fw_flow) / blk.fw_pump_efficiency,
    #         to_units=pyo.units.kilowatt,
    #     )
    # )

    blk.pumping_power_required_constraint = pyo.Constraint(
        expr=blk.pumping_power_required
        == blk.number_sw_pumps * pyo.units.convert(blk.sw_pump_power, to_units=pyo.units.kilowatt)
        + blk.number_fw_pumps * pyo.units.convert(blk.fw_pump_power, to_units=pyo.units.kilowatt)
    )

    blk.length_piping_constraint = pyo.Constraint(
        expr=blk.length_piping == pyo.units.convert(ss.number_stills * ss.still_length, to_units=pyo.units.ft)
    )

    blk.number_sw_pumps_constraint = pyo.Constraint(
        expr=blk.number_sw_pumps
        == pyo.units.convert(
            ss.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyo.units.m**3 / pyo.units.hr,
        )
        / blk.max_sw_flow
    )

    blk.number_fw_pumps_constraint = pyo.Constraint(
        expr=blk.number_fw_pumps
        == pyo.units.convert(
            ss.properties_out[0].flow_vol_phase["Liq"],
            to_units=pyo.units.m**3 / pyo.units.hr,
        )
        / blk.max_fw_flow
    )

    capital_cost_expr = 0

    blk.capital_cost_solar_still_constraint = pyo.Constraint(
        expr=blk.capital_cost_solar_still
        == pyo.units.convert(
            ss.total_area * ss_params.cost_solar_still, to_units=base_currency
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
        == pyo.units.convert(blk.length_piping * ss_params.cost_piping,
            to_units=base_currency,
        )
    )

    capital_cost_expr += blk.capital_cost_piping

    blk.costing_package.add_cost_factor(blk, None)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(capital_cost_expr, to_units=base_currency)
    )

    blk.costing_package.cost_flow(blk.pumping_power_required, "electricity")
