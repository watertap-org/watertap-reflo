#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)
from watertap.costing.unit_models.pump import build_low_pressure_pump_cost_param_block
from watertap.costing.util import register_costing_parameter_block


"""
Developed from 

U.S. Dept. of Interior & Michael C. Mickley (2006)
"Membrane Concentrate Disposal: Practices and Regulation"
Desalination and Water Purification Research and Development Program Report No. 123 (Second Edition)
Chapter 10: Evaporation Pond Disposal

Loh, H. P., Lyons, Jennifer, White, Charles W. (2002)
"Process Equipment Cost Estimation, Final Report"
https://doi.org/10.2172/797810

J. Gilron, E. Ramon, N. Assaf and O. Kedem (2019)
Current Trends and Future Developments on (Bio-) Membranes
"Chapter 9: Wind-Aided Intensified Evaporation (WAIV)"
https://doi.org/10.1016/B978-0-12-813551-8.00009-7

"""


def build_recovered_solids_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0,
        doc="Revenue from recovered salt",
        units=pyo.units.USD_2023 / pyo.units.kg,
    )

    costing = blk.parent_block()
    costing.register_flow_type("recovered_solids", blk.cost)


def build_waiv_cost_param_block(blk):

    blk.geotextile_membrane_cost = pyo.Var(
        initialize=5.77,
        bounds=(0, None),
        units=pyo.units.USD_2023 / pyo.units.m**2,
        doc="Cost of geotextile membrane fabric for WAIV module",
    )
    blk.geotextile_membrane_lifetime = pyo.Var(
        initialize=6,
        bounds=(0, None),
        units=pyo.units.year,
        doc="Lifetime of geotextile membranes (replacement frequency)",
    )
    # TODO: make land cost a global costing parameter
    blk.land_cost = pyo.Var(
        initialize=5000,  # From Mickley (2006), Chap. 10
        bounds=(0, None),
        units=pyo.units.USD_2001 / pyo.units.acre,
        doc="Land purchase cost per acre",
    )
    blk.land_clearing_cost = pyo.Var(
        initialize=2000,  # for sparsley wooded areas; from Mickley (2006), Chap. 10
        bounds=(0, None),
        units=pyo.units.USD_2001 / pyo.units.acre,
        doc="Land clearing cost per acre",
    )
    blk.foundation_cost = pyo.Var(
        initialize=217800,
        bounds=(0, None),
        units=pyo.units.USD_2023 / pyo.units.acre,
        doc="Foundation cost per acre",
    )
    blk.process_equipment_frac_cost = pyo.Var(
        initialize=0.525,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Cost of process equipment for WAIV system as fraction of geotextile membrane cost",
    )
    blk.labor_frac_cost = pyo.Var(
        initialize=0.155,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Cost of labor for WAIV system construction as fraction of process equipment cost",
    )
    blk.other_operating_cost_param = pyo.Var(
        initialize=85000 / 144000,
        bounds=(0, None),
        units=(pyo.units.USD_2015 / pyo.units.year)
        / (pyo.units.gallons / pyo.units.day),
        doc="Other operating costs (chemical, membranes, filter media) per gal/day",
    )
    blk.feed_tank_sizing_factor = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=(pyo.units.gallons) / (pyo.units.Mgallons / pyo.units.day),
        doc="Sizing factor to scale feed tank size required based on feed flow rate",
    )
    blk.feed_tank_cost_eq_base = pyo.Var(
        initialize=249.48,
        bounds=(0, None),
        units=pyo.units.USD_1998,
        doc="Base parameter for feed tank cost curve",
    )
    blk.feed_tank_cost_eq_exp = pyo.Var(
        initialize=0.5737,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Exponent for feed tank cost curve",
    )


@register_costing_parameter_block(
    build_rule=build_low_pressure_pump_cost_param_block,
    parameter_block_name="low_pressure_pump",
)
@register_costing_parameter_block(
    build_rule=build_recovered_solids_cost_param_block,
    parameter_block_name="recovered_solids",
)
@register_costing_parameter_block(
    build_rule=build_waiv_cost_param_block,
    parameter_block_name="waiv",
)
def cost_waiv(blk):
    """
    WAIV system cost
    """

    waiv_params = blk.costing_package.waiv

    blk.geotextile_membrane_capital_cost = pyo.Var(
        initialize=5.77,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Geotextile membrane fabric capital cost for WAIV module",
    )

    blk.land_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Land capital cost",
    )

    blk.land_clearing_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Land clearing capital cost",
    )

    blk.foundation_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Foundation capital cost",
    )

    blk.process_equipment_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Process equipment (structural steel, piping, instruments, electrical, painting, misc.) capital cost",
    )

    blk.feed_tank_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Feed tank capital cost",
    )

    blk.pump_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Pump capital cost",
    )

    blk.recovered_solids_handling_operating_cost = pyo.Var(
        initialize=0.1,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.year,
        doc="Operating cost to handle recovered solids",
    )

    blk.geotextile_membrane_replacement_operating_cost = pyo.Var(
        initialize=10000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.year,
        doc="Operating cost for replacement of geotextile membranes",
    )

    blk.other_operating_cost = pyo.Var(
        initialize=10000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.year,
        doc="Other operating cost (chemical, membrane, filter media, maintenance)",
    )

    blk.feed_tank_volume_required = pyo.Expression(
        expr=pyo.units.convert(
            waiv_params.feed_tank_sizing_factor
            * blk.unit_model.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyo.units.gallons,
        )
    )
    feed_tank_volume_required_dimensionless = pyo.units.convert(
        blk.feed_tank_volume_required * pyo.units.gallons**-1,
        to_units=pyo.units.dimensionless,
    )

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_fixed_operating_cost_var(blk)

    capital_cost_expr = 0

    blk.geotextile_membrane_capital_cost_constraint = pyo.Constraint(
        expr=blk.geotextile_membrane_capital_cost
        == pyo.units.convert(
            waiv_params.geotextile_membrane_cost
            * blk.unit_model.geotextile_area_required,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.geotextile_membrane_capital_cost

    blk.land_capital_cost_constraint = pyo.Constraint(
        expr=blk.land_capital_cost
        == pyo.units.convert(
            waiv_params.land_cost * blk.unit_model.total_land_area_required,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.land_capital_cost

    blk.land_clearing_capital_cost_constraint = pyo.Constraint(
        expr=blk.land_clearing_capital_cost
        == pyo.units.convert(
            waiv_params.land_clearing_cost * blk.unit_model.total_land_area_required,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.land_clearing_capital_cost

    blk.feed_tank_capital_cost_constraint = pyo.Constraint(
        expr=blk.feed_tank_capital_cost
        == pyo.units.convert(
            waiv_params.feed_tank_cost_eq_base
            * feed_tank_volume_required_dimensionless
            ** waiv_params.feed_tank_cost_eq_exp,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.feed_tank_capital_cost

    blk.pump_capital_cost_constraint = pyo.Constraint(
        expr=blk.pump_capital_cost
        == pyo.units.convert(
            blk.costing_package.low_pressure_pump.cost
            * blk.unit_model.properties_in[0].flow_vol_phase["Liq"],
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.pump_capital_cost

    blk.foundation_capital_cost_constraint = pyo.Constraint(
        expr=blk.foundation_capital_cost
        == pyo.units.convert(
            waiv_params.foundation_cost * blk.unit_model.total_land_area_required,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.foundation_capital_cost

    blk.process_equipment_capital_cost_constraint = pyo.Constraint(
        expr=blk.process_equipment_capital_cost
        == pyo.units.convert(
            blk.geotextile_membrane_capital_cost
            * waiv_params.process_equipment_frac_cost
            * (1 + waiv_params.labor_frac_cost),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.process_equipment_capital_cost

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == capital_cost_expr
    )

    fixed_operating_cost_expr = 0

    blk.other_operating_cost_constraint = pyo.Constraint(
        expr=blk.other_operating_cost
        == pyo.units.convert(
            waiv_params.other_operating_cost_param
            * blk.unit_model.properties_in[0].flow_vol_phase["Liq"],
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    fixed_operating_cost_expr += blk.other_operating_cost

    blk.geotextile_membrane_replacement_operating_cost_constraint = pyo.Constraint(
        expr=blk.geotextile_membrane_replacement_operating_cost
        == pyo.units.convert(
            blk.geotextile_membrane_capital_cost
            / waiv_params.geotextile_membrane_lifetime,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    fixed_operating_cost_expr += blk.geotextile_membrane_replacement_operating_cost

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == fixed_operating_cost_expr
    )

    if blk.unit_model.config.terminal_process:
        # assume all TDS is recovered as solids
        blk.costing_package.cost_flow(
            blk.unit_model.properties_in[0].flow_mass_phase_comp["Liq", "TDS"],
            "recovered_solids",
        )

    if not blk.unit_model.config.terminal_process:
        # assume difference between inlet and outlet is recovered as solids
        blk.costing_package.cost_flow(
            (
                blk.unit_model.properties_in[0].flow_mass_phase_comp["Liq", "TDS"]
                - blk.unit_model.properties_out[0].flow_mass_phase_comp["Liq", "TDS"]
            ),
            "recovered_solids",
        )
