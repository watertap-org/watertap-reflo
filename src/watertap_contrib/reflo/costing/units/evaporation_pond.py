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
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)
from watertap.costing.util import register_costing_parameter_block


costing_params_dict = {
    "dike_capital_cost": {12: (29471, -0.317), 8: (19000, -0.369), 4: (84624, -0.431)},
    "nominal_liner_capital_cost": {
        12: (1176, 0.2021),
        8: (15151, 0.1516),
        4: (20255, 0.0912),
    },
    "fence_capital_cost": {12: (13505, -0.455), 8: (11613, -0.423), 4: (1024, -0.398)},
    "road_capital_cost": {12: (2244, -0.454), 8: (1943, -0.424), 4: (1712.6, -0.399)},
}


def build_evaporation_pond_cost_param_block(blk):

    blk.basis_year = 2001
    blk.basis_currency = getattr(pyo.units, f"USD_{blk.basis_year}")

    blk.liner_thickness_base = pyo.Var(
        initialize=60,
        units=pyo.units.mil,
        doc="Basis for adjusting liner cost based on thickness",
    )

    blk.dike_capital_cost_base = pyo.Var(
        initialize=1e3,
        bounds=(0, None),
        units=blk.basis_currency / pyo.units.acre,
        doc="Dike capital cost per acre base",
    )

    blk.dike_capital_cost_exp = pyo.Var(
        initialize=-0.7,
        bounds=(None, 0),
        units=pyo.units.dimensionless,
        doc="Dike capital cost per acre exponent",
    )

    blk.nominal_liner_capital_cost_base = pyo.Var(
        initialize=1e3,
        bounds=(0, None),
        units=blk.basis_currency / pyo.units.acre,
        doc="Nominal liner capital cost per acre base",
    )

    blk.nominal_liner_capital_cost_exp = pyo.Var(
        initialize=0.7,
        bounds=(0, None),
        units=pyo.units.dimensionless,
        doc="Nominal liner capital cost per acre exponent",
    )

    blk.fence_capital_cost_base = pyo.Var(
        initialize=1e3,
        bounds=(0, None),
        units=blk.basis_currency / pyo.units.acre,
        doc="Fence capital cost per acre base",
    )

    blk.fence_capital_cost_exp = pyo.Var(
        initialize=-0.7,
        bounds=(None, 0),
        units=pyo.units.dimensionless,
        doc="Fence capital cost per acre exponent",
    )

    blk.road_capital_cost_base = pyo.Var(
        initialize=1e3,
        bounds=(0, None),
        units=blk.basis_currency / pyo.units.acre,
        doc="Road capital cost per acre base",
    )

    blk.road_capital_cost_exp = pyo.Var(
        initialize=-0.7,
        bounds=(None, 0),
        units=pyo.units.dimensionless,
        doc="Road capital cost per acre exponent",
    )

    blk.land_cost = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=pyo.units.USD_2001 / pyo.units.acre,
        doc="Land cost per acre",
    )

    blk.land_clearing_cost = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=pyo.units.USD_2001 / pyo.units.acre,
        doc="Land clearing cost per acre",
    )

    blk.liner_thickness = pyo.Var(
        initialize=60,
        bounds=(20, 120),
        units=pyo.units.mil,
        doc="Liner thickness",
    )

    blk.liner_replacement_frequency = pyo.Var(
        initialize=20,
        bounds=(0, None),
        units=pyo.units.year,
        doc="Liner replacement frequency",
    )

    blk.precipitate_handling_cost = pyo.Var(
        initialize=0,
        bounds=(None, None),
        units=blk.parent_block().base_currency / pyo.units.kg,
        doc="Cost to handle precipitated solids",
    )


@register_costing_parameter_block(
    build_rule=build_evaporation_pond_cost_param_block,
    parameter_block_name="evaporation_pond",
)
def cost_evaporation_pond(blk):
    """
    Evaporation pond cost
    """

    dike_height = blk.unit_model.config.dike_height
    pond_params = blk.costing_package.evaporation_pond

    for cv, d in costing_params_dict.items():
        cvps = d[dike_height]  # costing var params
        for i, n in zip([0, 1], ["base", "exp"]):
            v = getattr(pond_params, f"{cv}_{n}")
            v.fix(cvps[i])

    blk.basis_year = 2001
    blk.basis_currency = getattr(pyo.units, f"USD_{blk.basis_year}")
    evap_area_acre_dim = pyo.units.convert(
        blk.unit_model.evaporative_area_acre * pyo.units.acre**-1,
        to_units=pyo.units.dimensionless,
    )

    total_area_required_acre = pyo.units.convert(
        blk.unit_model.evaporation_pond_area * blk.unit_model.number_evaporation_ponds,
        to_units=pyo.units.acre,
    )

    blk.dike_cost_per_acre = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.acre,
        doc="Dike capital cost per acre",
    )

    blk.nominal_liner_cost_per_acre = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.acre,
        doc="Nominal liner cost per acre",
    )

    blk.liner_cost_per_acre = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.acre,
        doc="Liner cost per acre",
    )

    blk.fence_cost_per_acre = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.acre,
        doc="Fence cost per acre",
    )

    blk.road_cost_per_acre = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.acre,
        doc="Road cost per acre",
    )

    ########

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

    blk.dike_capital_cost = pyo.Var(
        initialize=500000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Dike capital cost",
    )

    blk.liner_capital_cost = pyo.Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Liner capital cost",
    )

    blk.fence_capital_cost = pyo.Var(
        initialize=1e5,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Fence capital cost",
    )

    blk.road_capital_cost = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Road capital cost",
    )

    blk.evaporation_enhancement_capital_cost = pyo.Var(
        initialize=0,
        bounds=(0, None),
        units=blk.costing_package.base_currency,
        doc="Capital cost for evaporation enhancement technology",
    )

    blk.precipitate_handling_operating_cost = pyo.Var(
        initialize=0.1,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.year,
        doc="Operating cost to handle precipitated solids",
    )

    blk.liner_replacement_operating_cost = pyo.Var(
        initialize=5000,
        bounds=(0, None),
        units=blk.costing_package.base_currency / pyo.units.year,
        doc="Operating cost to replace liner",
    )

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)
    make_fixed_operating_cost_var(blk)

    @blk.Expression(doc="Total cost per acre")
    def total_cost_per_acre(b):
        return (
            pyo.units.convert(
                pond_params.land_cost,
                to_units=b.costing_package.base_currency / pyo.units.acre,
            )
            + pyo.units.convert(
                pond_params.land_clearing_cost,
                to_units=b.costing_package.base_currency / pyo.units.acre,
            )
            + b.dike_cost_per_acre
            + b.liner_cost_per_acre
            + b.fence_cost_per_acre
            + b.road_cost_per_acre
        )

    @blk.Constraint(doc="Dike cost per acre")
    def dike_cost_per_acre_constraint(b):
        return b.dike_cost_per_acre == pyo.units.convert(
            pond_params.dike_capital_cost_base
            * evap_area_acre_dim**pond_params.dike_capital_cost_exp,
            to_units=b.costing_package.base_currency / pyo.units.acre,
        )

    @blk.Constraint(doc="Nominal liner cost per acre")
    def nominal_liner_cost_per_acre_constraint(b):
        return b.nominal_liner_cost_per_acre == pyo.units.convert(
            pond_params.nominal_liner_capital_cost_base
            * evap_area_acre_dim**pond_params.nominal_liner_capital_cost_exp,
            to_units=b.costing_package.base_currency / pyo.units.acre,
        )

    @blk.Constraint(doc="Liner cost per acre")
    def liner_cost_per_acre_constraint(b):
        return b.liner_cost_per_acre == pyo.units.convert(
            b.nominal_liner_cost_per_acre
            * (pond_params.liner_thickness / pond_params.liner_thickness_base),
            to_units=b.costing_package.base_currency / pyo.units.acre,
        )

    @blk.Constraint(doc="Fence cost per acre")
    def fence_cost_per_acre_constraint(b):
        return b.fence_cost_per_acre == pyo.units.convert(
            pond_params.fence_capital_cost_base
            * evap_area_acre_dim**pond_params.fence_capital_cost_exp,
            to_units=b.costing_package.base_currency / pyo.units.acre,
        )

    @blk.Constraint(doc="Road cost per acre")
    def road_cost_per_acre_constraint(b):
        return b.road_cost_per_acre == pyo.units.convert(
            pond_params.road_capital_cost_base
            * evap_area_acre_dim**pond_params.road_capital_cost_exp,
            to_units=b.costing_package.base_currency / pyo.units.acre,
        )

    capital_cost_expr = 0

    @blk.Constraint(doc="Capital cost for land")
    def land_capital_cost_constraint(b):
        return b.land_capital_cost == pyo.units.convert(
            total_area_required_acre * pond_params.land_cost,
            to_units=b.costing_package.base_currency,
        )

    capital_cost_expr += blk.land_capital_cost

    @blk.Constraint(doc="Capital cost for land clearing")
    def land_clearing_capital_cost_constraint(b):
        return b.land_clearing_capital_cost == pyo.units.convert(
            total_area_required_acre * pond_params.land_clearing_cost,
            to_units=b.costing_package.base_currency,
        )

    capital_cost_expr += blk.land_clearing_capital_cost

    @blk.Constraint(doc="Dike capital cost")
    def dike_capital_cost_constraint(b):
        return b.dike_capital_cost == b.dike_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.dike_capital_cost

    @blk.Constraint(doc="Liner capital cost")
    def liner_capital_cost_constraint(b):
        return b.liner_capital_cost == b.liner_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.liner_capital_cost

    @blk.Constraint(doc="Fence capital cost")
    def fence_capital_cost_constraint(b):
        return b.fence_capital_cost == b.fence_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.fence_capital_cost

    @blk.Constraint(doc="Road capital cost")
    def road_capital_cost_constraint(b):
        return b.road_capital_cost == b.road_cost_per_acre * total_area_required_acre

    capital_cost_expr += blk.road_capital_cost

    @blk.Constraint(doc="Evaporation enhancement technology capital cost")
    def evaporation_enhancement_capital_cost_constraint(b):
        return (
            b.evaporation_enhancement_capital_cost
            == 0 * b.costing_package.base_currency
        )  # placeholder

    capital_cost_expr += blk.evaporation_enhancement_capital_cost

    @blk.Constraint(doc="Capital cost for evaporation pond")
    def capital_cost_constraint(b):
        return b.capital_cost == capital_cost_expr

    fixed_operating_cost_expr = 0

    @blk.Constraint(doc="Liner replacement cost")
    def liner_replacement_operating_cost_constraint(b):
        return b.liner_replacement_operating_cost == pyo.units.convert(
            b.liner_capital_cost / pond_params.liner_replacement_frequency,
            to_units=b.costing_package.base_currency / b.costing_package.base_period,
        )

    fixed_operating_cost_expr += blk.liner_replacement_operating_cost

    @blk.Constraint(doc="Solids handling operating cost")
    def precipitate_handling_operating_cost_constraint(b):
        return b.precipitate_handling_operating_cost == pyo.units.convert(
            b.unit_model.mass_flow_precipitate * pond_params.precipitate_handling_cost,
            to_units=b.costing_package.base_currency / b.costing_package.base_period,
        )

    fixed_operating_cost_expr += blk.precipitate_handling_operating_cost

    @blk.Constraint(doc="Fixed operating cost for evaporation pond")
    def fixed_operating_cost_constraint(b):
        return b.fixed_operating_cost == fixed_operating_cost_expr
