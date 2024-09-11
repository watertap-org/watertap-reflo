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

# Costing numbers are from data provided by Trevi (2021),
# and should be applied to a Trevi FO flowsheet.
# Mass and energy flows should be calculated from a Trevi FO flowsheet.


def build_forward_osmosis_cost_param_block(blk):

    costing = blk.parent_block()

    blk.base_unit_cost = pyo.Var(
        initialize=26784,
        units=costing.base_currency / (pyo.units.m**3 / pyo.units.day),
        bounds=(0, None),
        doc="Base price of Trevi FO system ($/m3/day), including membranes, heat exchangers, \
            construction, draw solution, coalescers, structural, polishing, pipes, plumbing, \
            pre-filtration, control, pumps, instrumentation, valves, CIP and tanks",
    )

    blk.unit_cost_index = pyo.Var(
        initialize=-0.428,
        units=pyo.units.dimensionless,
        bounds=(None, 0),
        doc="Scaling factor of Trevi FO system capital cost",
    )

    blk.base_labor_cost = pyo.Var(
        initialize=0.4757,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Base price of the labor cost ($/m3)",
    )

    blk.labor_cost_index = pyo.Var(
        initialize=-0.178,
        units=pyo.units.dimensionless,
        bounds=(None, 0),
        doc="Scaling factor of Trevi FO system capital cost",
    )

    blk.specific_energy_consumption_electric = pyo.Var(
        initialize=1,
        units=pyo.units.kWh / pyo.units.m**3,
        bounds=(0, None),
        doc="Specific electric energy consumption (kWh/m3)",
    )

    blk.cost_chemicals = pyo.Var(
        initialize=0.07,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of chemicals per m3 product",
    )

    blk.cost_durable_goods = pyo.Var(
        initialize=0.05,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of durable goods per m3 product",
    )

    blk.cost_disposal = pyo.Var(
        initialize=0.02,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Cost of disposal per m3 brine",
    )

    blk.cost_fraction_insurance = pyo.Var(
        initialize=0.005,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Fraction of capital cost for insurance",
    )


@register_costing_parameter_block(
    build_rule=build_forward_osmosis_cost_param_block,
    parameter_block_name="forward_osmosis",
)
def cost_forward_osmosis(blk):

    fo_params = blk.costing_package.forward_osmosis
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    fo = blk.unit_model
    test_fs = fo.flowsheet()

    brine = fo.brine_props[0]

    blk.costing_package.add_cost_factor(blk, None)
    blk.annual_dist_production = pyo.units.convert(
        test_fs.system_capacity, to_units=pyo.units.m**3 / pyo.units.year
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == fo_params.base_unit_cost * test_fs.system_capacity**fo_params.unit_cost_index
    )

    blk.labor_cost = pyo.Expression(
        expr=fo_params.base_labor_cost
        * test_fs.system_capacity**fo_params.labor_cost_index
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=(
            blk.fixed_operating_cost
            == blk.annual_dist_production
            * (blk.labor_cost + fo_params.cost_chemicals + fo_params.cost_durable_goods)
            + blk.capital_cost * fo_params.cost_fraction_insurance
            + pyo.units.convert(
                brine.flow_vol_phase["Liq"], to_units=pyo.units.m**3 / pyo.units.year
            )
            * fo_params.cost_disposal
        )
    )

    blk.thermal_energy_flow = pyo.Expression(
        expr=test_fs.specific_energy_consumption_thermal
        * pyo.units.convert(
            test_fs.system_capacity, to_units=pyo.units.m**3 / pyo.units.hr
        )
    )

    blk.electricity_flow = pyo.Expression(
        expr=fo_params.specific_energy_consumption_electric
        * pyo.units.convert(
            test_fs.system_capacity, to_units=pyo.units.m**3 / pyo.units.hr
        )
    )

    blk.costing_package.cost_flow(blk.thermal_energy_flow, "heat")
    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
