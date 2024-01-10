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
import idaes.core.util.scaling as iscale
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)

# Costing equations from:
# Hitsov, I., De Sitter, K., Dotremont, C., & Nopens, I. (2018).
# Economic modelling and model-based process optimization of
# membrane distillation. Desalination, 436, 125-143.
# doi: 10.1016/j.desal.2018.01.038


def build_vagmd_surrogate_cost_param_block(blk):

    costing = blk.parent_block()

    blk.base_module_cost = pyo.Var(
        initialize=2.34e3,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Base price of AGMD module assembly ($/base capacity)",
    )

    blk.base_module_capacity = pyo.Var(
        initialize=3,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Base capacity of module assembly (modules)",
    )

    blk.module_cost_index = pyo.Var(
        initialize=0.8,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Module assembly cost scaling index",
    )

    blk.membrane_cost = pyo.Var(
        initialize=90,
        units=costing.base_currency / pyo.units.m**2,
        bounds=(0, None),
        doc="Base price of AGMD membrane",
    )

    blk.heat_exchanger_eff = pyo.Var(
        initialize=0.85,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Heat exchanger efficiency",
    )

    blk.heat_transfer_coefficient = pyo.Var(
        initialize=2.5,
        units=pyo.units.kW / pyo.units.K / pyo.units.m**2,
        bounds=(0, None),
        doc="Heat transfer coeffcient of the heat exchanger",
    )

    blk.base_HX_cost = pyo.Var(
        initialize=420,
        units=costing.base_currency / pyo.units.m**2,
        bounds=(0, None),
        doc="Base price of heat excahnger ($/m2)",
    )

    blk.base_endplates_cost = pyo.Var(
        initialize=1.02e3,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Base price of heat excahnger endplates ($/m2 HX)",
    )

    blk.base_endplates_capacity = pyo.Var(
        initialize=10,
        units=pyo.units.m**2,
        bounds=(0, None),
        doc="Base capacity of heat excahnger endplates (m2 HX)",
    )

    blk.endplates_cost_index = pyo.Var(
        initialize=0.6,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Endplates cost scaling index",
    )

    blk.base_housing_rack_cost = pyo.Var(
        initialize=6e3,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Base price of housing rack ($/base capacity)",
    )

    blk.base_housing_rack_capacity = pyo.Var(
        initialize=3,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Base capacity of housing rack(modules)",
    )

    blk.housing_rack_cost_index = pyo.Var(
        initialize=0.6,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Housing rack cost scaling index",
    )

    blk.base_tank_cost = pyo.Var(
        initialize=6e3,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Base price of tank ($/base capacity)",
    )

    blk.base_tank_capacity = pyo.Var(
        initialize=3,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Base capacity of tank(modules)",
    )

    blk.tank_cost_index = pyo.Var(
        initialize=0.5,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Tank cost scaling index",
    )

    blk.base_pump_cost = pyo.Var(
        initialize=3.6e3,
        units=costing.base_currency / (pyo.units.m**3 / pyo.units.h),
        bounds=(0, None),
        doc="Base price of pump ($/base capacity)",
    )

    blk.base_pump_capacity = pyo.Var(
        initialize=5,
        units=pyo.units.m**3 / pyo.units.h,
        bounds=(0, None),
        doc="Base capacity of pump(m3/h)",
    )

    blk.pump_cost_index = pyo.Var(
        initialize=0.5,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Pump cost scaling index",
    )

    blk.base_other_cost = pyo.Var(
        initialize=18e3,
        units=costing.base_currency,
        bounds=(0, None),
        doc="Base price of other capital ($/base capacity)",
    )

    blk.base_other_capacity = pyo.Var(
        initialize=3,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Base capacity of other capital(modules)",
    )

    blk.other_cost_index = pyo.Var(
        initialize=0.3,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Tank cost scaling index",
    )

    blk.cost_fraction_maintenance = pyo.Var(
        initialize=0.013,
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

    blk.membrane_replacement_cost = pyo.Var(
        initialize=0.22,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Membrane replacement cost ($/m3)",
    )

    blk.specific_operational_cost = pyo.Var(
        initialize=0.1,
        units=costing.base_currency / pyo.units.m**3,
        bounds=(0, None),
        doc="Specific operational cost ($/m3)",
    )


@register_costing_parameter_block(
    build_rule=build_vagmd_surrogate_cost_param_block,
    parameter_block_name="vagmd_surrogate",
)
def cost_vagmd_surrogate(blk):

    vagmd_params = blk.costing_package.vagmd_surrogate
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    vagmd = blk.unit_model
    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.module_cost = pyo.Var(
        initialize=100000,
        bounds=(0, None),
        units=base_currency,
        doc="MD module cost",
    )

    blk.other_capital_cost = pyo.Var(
        initialize=10000,
        bounds=(0, None),
        units=base_currency,
        doc="""Other capital cost (includes housing rack,
                tank, controller, cabling and programming)""",
    )

    blk.annual_dist_production = pyo.units.convert(
        vagmd.system_capacity, to_units=pyo.units.m**3 / pyo.units.year
    )

    blk.module_cost_constraint = pyo.Constraint(
        expr=blk.module_cost
        == vagmd_params.base_module_cost
        * vagmd_params.base_module_capacity
        * (vagmd.num_modules / vagmd_params.base_module_capacity)
        ** vagmd_params.module_cost_index
        + vagmd_params.membrane_cost * vagmd.module_area * vagmd.num_modules
    )

    blk.other_capital_cost_constraint = pyo.Constraint(
        expr=blk.other_capital_cost
        == vagmd_params.base_housing_rack_cost
        * (vagmd.num_modules / vagmd_params.base_housing_rack_capacity)
        ** vagmd_params.housing_rack_cost_index
        + vagmd_params.base_tank_cost
        * (vagmd.num_modules / vagmd_params.base_tank_capacity)
        ** vagmd_params.tank_cost_index
        + vagmd_params.base_other_cost
        * (vagmd.num_modules / vagmd_params.base_other_capacity)
        ** vagmd_params.other_cost_index
    )

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.module_cost + blk.other_capital_cost
    )

    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            vagmd.system_capacity, to_units=pyo.units.m**3 / pyo.units.year
        )
        * (
            vagmd_params.membrane_replacement_cost
            + vagmd_params.specific_operational_cost
        )
        + blk.capital_cost
        * (
            vagmd_params.cost_fraction_maintenance
            + vagmd_params.cost_fraction_insurance
        )
    )

    blk.costing_package.cost_flow(vagmd.thermal_power_requirement, "heat")
    blk.costing_package.cost_flow(vagmd.elec_power_requirement, "electricity")

    if iscale.get_scaling_factor(blk.module_cost) is None:
        iscale.set_scaling_factor(blk.module_cost, 1e-5)

    if iscale.get_scaling_factor(blk.other_capital_cost) is None:
        iscale.set_scaling_factor(blk.other_capital_cost, 1e-5)

    if iscale.get_scaling_factor(blk.capital_cost) is None:
        iscale.set_scaling_factor(blk.capital_cost, 1e-5)
