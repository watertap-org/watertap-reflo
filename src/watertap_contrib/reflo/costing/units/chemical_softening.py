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


from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Var,
    Param,
    units as pyunits,
)

from watertap.costing.util import register_costing_parameter_block
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_lime_cost_param_block(blk):

    blk.cost = Param(
        initialize=0.10,
        units=pyunits.USD_2021 / pyunits.kg,
        doc="Cost of CaO $/kg",  # taken from CatCost, 9/2024
    )

    blk.purity = Param(
        mutable=True,
        initialize=1,
        doc="CaO purity",
        units=pyunits.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("lime", blk.cost / blk.purity)


def build_soda_ash_cost_param_block(blk):

    blk.cost = Param(
        initialize=0.24,
        units=pyunits.USD_2021 / pyunits.kg,
        doc="Cost of Na2CO3 $/kg",  # taken from CatCost, 9/2024
    )

    blk.purity = Param(
        mutable=True,
        initialize=1,
        doc="Na2CO3 purity",
        units=pyunits.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("soda_ash", blk.cost / blk.purity)


def build_mgcl2_cost_param_block(blk):

    blk.cost = Param(
        initialize=0.49,
        units=pyunits.USD_2020 / pyunits.kg,
        doc="Cost of MgCl2 $/kg",  # taken from CatCost, 9/2024
    )

    blk.purity = Param(
        mutable=True,
        initialize=1,
        doc="MgCl2 purity",
        units=pyunits.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("mgcl2", blk.cost / blk.purity)


def build_co2_cost_param_block(blk):

    blk.cost = Param(
        initialize=0.38,
        units=pyunits.USD_2020 / pyunits.kg,
        doc="Cost of CO2 $/kg",
    )

    blk.purity = Param(
        mutable=True,
        initialize=1,
        doc="CO2 purity",
        units=pyunits.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("co2", blk.cost / blk.purity)


def build_chemical_softening_cost_param_block(blk):
    """
    Parameters and variables to be used in the costing model
    References :
    Cost functions- http://hdl.handle.net/10106/4924
    Material costs- https://www.chemanalyst.com/
    """
    costing = blk.parent_block()

    blk.sed_basin_depth = Var(
        initialize=4.5,
        bounds=(3, 6),
        units=pyunits.m,
        doc="Depth of sedimentation basin",
    )

    blk.sludge_disposal_cost = Var(
        initialize=35,
        units=costing.base_currency / pyunits.tons,
        doc="Cost of sludge disposal $/tonne",
    )

    # Costing equation coefficients - Capital cost

    # Rapid Mix, G=300 s-1
    # Adapted from Equation 4.22 in Sharma
    # C = A*x**b
    # Original bounds = (100, 20000)

    blk.mixer_capital_coeff = Param(
        initialize=1219.22948,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Capital coefficient for mixer",
    )
    blk.mixer_capital_exponent = Param(
        initialize=0.730785,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Capital exponent for mixer",
    )

    blk.floc_tank_capital_coeff = Param(
        initialize=673894,
        mutable=True,
        units=pyunits.USD_2007 * pyunits.Mgallons**-1,
        doc="Coefficient of flocculation tank capital cost equation",
    )

    blk.floc_tank_capital_intercept = Param(
        initialize=217222,
        mutable=True,
        units=pyunits.USD_2007,
        doc="Intercept in the flocculation tank capital cost equation",
    )

    # Sedimentation, circular clarifiers
    # Adapted from Equation 4.32 in Sharma
    # C = A*x**b
    # Original bounds = (707, 31416)

    blk.sed_capital_coeff = Param(
        initialize=969.38,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Capital coefficient for sedimentation basin",
    )
    blk.sed_capital_exponent = Param(
        initialize=0.75539,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Capital exponent for sedimentation basin",
    )

    # Recarbonation basin
    # Adapted from Equation 4.60 in Sharma to fit
    # C = A*x**b
    # Original bounds = (770, 35200)

    blk.recarb_basin_capital_coeff = Param(
        initialize=56.60714729,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Capital coefficient for recarbonation basin",
    )
    blk.recarb_basin_capital_exponent = Param(
        initialize=0.81337136,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Capital coefficient for recarbonation basin",
    )

    # Recarbonation basin, CO2 source (Liquid CO2 as CO2 source)
    # Adapted from Equation 4.61 in Sharma to fit
    # C = A*x**b
    # Original bounds = (380, 15000)

    blk.recarb_source_capital_coeff = Param(
        initialize=549.13770574,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Capital coefficient for recarbonation basin",
    )
    blk.recarb_source_capital_exponent = Param(
        initialize=0.75494096,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Capital exponent for recarbonation basin",
    )

    # Lime feed system capital cost
    # Adapted from Equation 4.59 in Sharma to fit
    # C = A*x**b
    # Original bounds = (1000, 10000)

    blk.lime_feed_system_capital_coeff = Param(
        initialize=21420.1311,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Capital coefficient for lime feed system",
    )

    blk.lime_feed_system_capital_exponent = Param(
        initialize=0.311419,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Capital exponent for lime feed system",
    )

    # Admin capital costs
    # Equation 4.94 in Sharma

    blk.admin_capital_coeff = Param(
        initialize=69195,
        mutable=True,
        units=pyunits.USD_2009,
        doc="Coefficient in lime feed system capital cost",
    )

    blk.admin_capital_exp = Param(
        initialize=0.5523,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Exponent in lime feed system capital cost",
    )

    # Costing equation coefficients - Operational cost

    # Mixer operational costs
    # Adapted from Equation 4.144 in Sharma to fit
    # C = A*x**b
    # Original bounds = (100, 20000)

    blk.mixer_op_coeff = Param(
        initialize=78.98,
        mutable=True,
        units=pyunits.USD_2009 / pyunits.year,
        doc="Opex coefficient for mixers",
    )

    blk.mixer_op_exponent = Param(
        initialize=0.776291,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Opex exponent for mixers",
    )

    # Floc operational costs
    # Adapted from Equation 4.148 in Sharma to fit
    # C = A*x**b
    # Original bounds = (1800, 1000000)

    blk.floc_op_coeff = Param(
        initialize=0.87966855,
        mutable=True,
        units=pyunits.USD_2009 / pyunits.year,
        doc="Opex coefficient for flocculators",
    )

    blk.floc_op_exponent = Param(
        initialize=0.89292125,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Opex exponent for flocculators",
    )

    # Sedimentation operational costs
    # Adapted from Equation 4.151 in Sharma to fit
    # C = m*x + b
    # Original bounds = (30, 200)

    blk.sed_op_slope = Param(
        initialize=1.54138506,
        mutable=True,
        units=pyunits.USD_2009 / pyunits.year,
        doc="Opex slope for sedimenation basin",
    )

    blk.sed_op_intercept = Param(
        initialize=6880,
        mutable=True,
        units=pyunits.USD_2009 / pyunits.year,
        doc="Opex intercept for sedimentation basin",
    )

    # Recarbonation basin operational costs
    # Adapted from Equation 4.171 in Sharma to fit
    # C = A*x**b
    # Original bounds = (30, 200)

    blk.recarb_basin_op_coeff = Param(
        initialize=1245.55473,
        mutable=True,
        units=pyunits.USD_2009 / pyunits.year,
        doc="Opex coefficient for recarbonation basin",
    )

    blk.recarb_basin_op_exponent = Param(
        initialize=0.3806758,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Opex exponent for recarbonation basin",
    )

    # Lime feed system operational costs
    # Equation 4.170 in Sharma
    # DOES NOT INCLUDE CHEMICAL COSTS (see Table C.1 in Sharma)

    blk.lime_feed_system_op_coeff = Param(
        initialize=4616.7,
        mutable=True,
        doc="Coefficient in lime feed system operational cost",
    )

    blk.lime_feed_system_op_exp = Param(
        initialize=0.4589,
        mutable=True,
        doc="Constant in lime feed system operational cost",
    )

    # Admin operational costs
    # Equation 4.196 in Sharma

    blk.admin_op_coeff = Param(
        initialize=88589,
        mutable=True,
        doc="Coefficient in lime feed system operational cost",
    )

    blk.admin_op_exp = Param(
        initialize=0.4589,
        mutable=True,
        doc="Exponent in lime feed system operational cost",
    )


@register_costing_parameter_block(
    build_rule=build_lime_cost_param_block,
    parameter_block_name="lime",
)
@register_costing_parameter_block(
    build_rule=build_soda_ash_cost_param_block,
    parameter_block_name="soda_ash",
)
@register_costing_parameter_block(
    build_rule=build_mgcl2_cost_param_block,
    parameter_block_name="mgcl2",
)
@register_costing_parameter_block(
    build_rule=build_co2_cost_param_block,
    parameter_block_name="co2",
)
@register_costing_parameter_block(
    build_rule=build_chemical_softening_cost_param_block,
    parameter_block_name="chemical_softening",
)
def cost_chemical_softening(blk):
    """
    Capital and operating costs for chemical softening
    """
    chem_soft = blk.costing_package.chemical_softening
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.mixer_power = Var(
        initialize=1000,
        units=pyunits.W,
        domain=NonNegativeReals,
        doc="Power consumption for rapid mixer",
    )

    blk.floc_power = Var(
        initialize=100,
        units=pyunits.W,
        domain=NonNegativeReals,
        doc="Power consumption for flocculation tank",
    )

    blk.electricity_flow = Var(
        initialize=100,
        domain=NonNegativeReals,
        units=pyunits.kW,
        doc="Total electricity consumption",
    )

    # Capital cost components

    blk.mix_tank_capital_cost = Var(
        initialize=1e6,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital Cost of mixing tank",
    )

    blk.floc_tank_capital_cost = Var(
        initialize=1e6,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of flocculation tank",
    )

    blk.sed_basin_capital_cost = Var(
        initialize=1e6,
        doc="Capital cost of sedimentation basin",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.recarb_basin_capital_cost = Var(
        initialize=1e6,
        doc="Capital cost of recarbonation basin",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.recarb_basin_source_capital_cost = Var(
        initialize=1e6,
        doc="Capital cost of recarbonation tank - liquid CO2 as CO2 source",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.lime_feed_system_capital_cost = Var(
        initialize=1e6,
        doc="Capital cost of lime feed system",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.admin_capital_cost = Var(
        initialize=1e6,
        doc="Capital cost of administrative, laboratory, and maintenance building",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    # Operating cost components

    blk.mixer_op_cost = Var(
        initialize=1e4,
        doc="Operational cost of mixing tank - energy consumption and labor",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.floc_tank_op_cost = Var(
        initialize=1e4,
        doc="Operational cost of flocculation tank - energy consumption and labor",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.sed_basin_op_cost = Var(
        initialize=1e4,
        doc="Operational cost of sedimentation basin - energy consumption, labor, and maintenance",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.recarb_basin_op_cost = Var(
        initialize=1e4,
        doc="Operationg cost of recarbonation basin - energy consumption, labor, and maintenance",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.lime_feed_op_cost = Var(
        initialize=1e4,
        doc="Operation cost of lime feed system - energy consumption, labor, and maintenance",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.lime_sludge_mngt_op_cost = Var(
        initialize=1e5,
        doc="Operational cost of lime sludge management",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    blk.admin_op_cost = Var(
        initialize=1e5,
        doc="Operational cost of administrative, laboratory and maintenance building",
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
    )

    # Capital cost component constraints

    capital_cost_expr = 0

    # Mixing tank
    blk.volume_mixer_ft3_dimensionless = pyunits.convert(
        blk.unit_model.volume_mixer * pyunits.ft**-3, to_units=pyunits.dimensionless
    )
    blk.mix_tank_capital_cost_constraint = Constraint(
        expr=blk.mix_tank_capital_cost
        == pyunits.convert(
            chem_soft.mixer_capital_coeff
            * blk.volume_mixer_ft3_dimensionless**chem_soft.mixer_capital_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.mix_tank_capital_cost * blk.unit_model.number_mixers

    # Flocculation tank
    blk.floc_tank_capital_cost_constraint = Constraint(
        expr=blk.floc_tank_capital_cost
        == pyunits.convert(
            (
                chem_soft.floc_tank_capital_coeff
                * pyunits.convert(blk.unit_model.volume_floc, to_units=pyunits.Mgallon)
                + chem_soft.floc_tank_capital_intercept
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.floc_tank_capital_cost * blk.unit_model.number_floc

    # Sedimentation basin
    blk.sed_basin_effective_settling_area_ft2_dimensionless = pyunits.convert(
        blk.unit_model.volume_sed / chem_soft.sed_basin_depth * pyunits.ft**-2,
        to_units=pyunits.dimensionless,
    )
    blk.sed_basin_capital_cost_constraint = Constraint(
        expr=blk.sed_basin_capital_cost
        == pyunits.convert(
            (
                chem_soft.sed_capital_coeff
                * blk.sed_basin_effective_settling_area_ft2_dimensionless
                ** chem_soft.sed_capital_exponent
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.sed_basin_capital_cost

    # Recarbonation basin
    blk.recarb_basin_vol_ft3_dimensionless = pyunits.convert(
        blk.unit_model.volume_recarb * pyunits.ft**-3, to_units=pyunits.dimensionless
    )
    blk.recarb_basin_capital_cost_constraint = Constraint(
        expr=blk.recarb_basin_capital_cost
        == pyunits.convert(
            chem_soft.recarb_basin_capital_coeff
            * blk.recarb_basin_vol_ft3_dimensionless
            ** chem_soft.recarb_basin_capital_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.recarb_basin_capital_cost

    # Recarbonation source cost
    blk.recarb_source_first_basin_lb_day_dimensionless = pyunits.convert(
        blk.unit_model.CO2_first_basin * pyunits.day * pyunits.lb**-1,
        to_units=pyunits.dimensionless,
    )
    blk.recarb_source_second_basin_lb_day_dimensionless = pyunits.convert(
        blk.unit_model.CO2_second_basin * pyunits.day * pyunits.lb**-1,
        to_units=pyunits.dimensionless,
    )
    blk.recarb_basin_source_capital_cost_constraint = Constraint(
        expr=blk.recarb_basin_source_capital_cost
        == pyunits.convert(
            chem_soft.recarb_source_capital_coeff
            * (
                blk.recarb_source_first_basin_lb_day_dimensionless
                + blk.recarb_source_second_basin_lb_day_dimensionless
            )
            ** chem_soft.recarb_source_capital_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.recarb_basin_source_capital_cost

    # Lime feed system cost
    blk.lime_dosing_lb_hr_dimensionless = pyunits.convert(
        blk.unit_model.CaO_dosing * pyunits.hr * pyunits.lb**-1,
        to_units=pyunits.dimensionless,
    )
    blk.lime_feed_system_capital_cost_constraint = Constraint(
        expr=blk.lime_feed_system_capital_cost
        == pyunits.convert(
            chem_soft.lime_feed_system_capital_coeff
            * blk.lime_dosing_lb_hr_dimensionless
            ** chem_soft.lime_feed_system_capital_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.lime_feed_system_capital_cost

    # Admin cost
    flow_vol_mgd_dimensionless = pyunits.convert(
        blk.unit_model.properties_in[0].flow_vol_phase["Liq"]
        * pyunits.day
        * pyunits.Mgallons**-1,
        to_units=pyunits.dimensionless,
    )
    blk.admin_capital_cost_constraint = Constraint(
        expr=blk.admin_capital_cost
        == pyunits.convert(
            chem_soft.admin_capital_coeff
            * flow_vol_mgd_dimensionless**chem_soft.admin_capital_exp,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.admin_capital_cost

    # Operational cost component constraints

    op_cost_expr = 0

    # Mixing tank
    blk.mixer_volume_ft3_dimensionless = pyunits.convert(
        blk.unit_model.volume_mixer * pyunits.ft**-3, to_units=pyunits.dimensionless
    )
    blk.mix_tank_op_cost_constraint = Constraint(
        expr=blk.mixer_op_cost
        == pyunits.convert(
            chem_soft.mixer_op_coeff
            * blk.mixer_volume_ft3_dimensionless**chem_soft.mixer_op_exponent,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    op_cost_expr += blk.mixer_op_cost

    # Flocculation tank
    blk.floc_vol_ft3_dimensionless = pyunits.convert(
        blk.unit_model.volume_floc * pyunits.ft**-3, to_units=pyunits.dimensionless
    )
    blk.floc_tank_op_cost_constraint = Constraint(
        expr=blk.floc_tank_op_cost
        == pyunits.convert(
            chem_soft.floc_op_coeff
            * blk.floc_vol_ft3_dimensionless**chem_soft.floc_op_exponent,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    op_cost_expr += blk.floc_tank_op_cost

    # Sedimentation basin
    blk.sed_surf_area_ft2_dimensionless = pyunits.convert(
        (blk.unit_model.volume_sed / chem_soft.sed_basin_depth) * pyunits.ft**-2,
        to_units=pyunits.dimensionless,
    )
    blk.sed_basin_op_cost_constraint = Constraint(
        expr=blk.sed_basin_op_cost
        == pyunits.convert(
            chem_soft.sed_op_intercept
            + chem_soft.sed_op_slope * blk.sed_surf_area_ft2_dimensionless,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    op_cost_expr += blk.sed_basin_op_cost

    # Recarbonation basin
    blk.recarb_basin_op_cost_constraint = Constraint(
        expr=blk.recarb_basin_op_cost
        == pyunits.convert(
            chem_soft.recarb_basin_op_coeff
            * (
                blk.recarb_source_first_basin_lb_day_dimensionless
                + blk.recarb_source_second_basin_lb_day_dimensionless
            )
            ** chem_soft.recarb_basin_op_exponent,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    op_cost_expr += blk.recarb_basin_op_cost

    # Lime feed
    blk.lime_feed_op_cost_constraint = Constraint(
        expr=blk.lime_feed_op_cost
        == (
            chem_soft.lime_feed_system_op_coeff
            * (
                pyunits.convert(
                    blk.unit_model.CaO_dosing, to_units=pyunits.lb / pyunits.day
                )
            )
            ** chem_soft.lime_feed_system_op_exp
        )
    )

    op_cost_expr += blk.lime_feed_op_cost

    # Lime sludge management
    blk.lime_sludge_mngt_op_cost_constraint = Constraint(
        expr=blk.lime_sludge_mngt_op_cost
        == (
            pyunits.convert(
                blk.unit_model.sludge_prod, to_units=pyunits.tons / pyunits.year
            )
            * chem_soft.sludge_disposal_cost
        )
    )

    op_cost_expr += blk.lime_sludge_mngt_op_cost

    # Admin cost
    blk.admin_op_cost_constraint = Constraint(
        expr=blk.admin_op_cost
        == (
            chem_soft.admin_op_coeff
            * (
                (
                    pyunits.convert(
                        blk.unit_model.properties_in[0].flow_vol_phase["Liq"],
                        to_units=pyunits.megagallon / pyunits.day,
                    )
                )
                ** chem_soft.admin_op_exp
            )
        )
    )

    op_cost_expr += blk.admin_op_cost

    # Sum of all capital costs
    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = Constraint(expr=blk.capital_cost == capital_cost_expr)

    # Sum of all operational costs

    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost == op_cost_expr
    )

    blk.mixer_power_constraint = Constraint(
        expr=blk.mixer_power
        == (blk.unit_model.vel_gradient_mix**2)
        * blk.unit_model.volume_mixer
        * blk.unit_model.properties_in[0].params.visc_d_phase["Liq"]
    )

    blk.floc_power_constraint = Constraint(
        expr=blk.floc_power
        == (blk.unit_model.vel_gradient_floc**2)
        * blk.unit_model.volume_floc
        * blk.unit_model.properties_in[0].params.visc_d_phase["Liq"]
    )

    blk.electricity_flow_constraint = Constraint(
        expr=blk.electricity_flow
        == pyunits.convert((blk.mixer_power + blk.floc_power), to_units=pyunits.kW)
    )

    @blk.Expression(doc="Costing block lime dosing")
    def cao_dosing(b):
        return pyunits.convert(
            b.unit_model.CaO_dosing
            + pyunits.convert(
                b.unit_model.excess_CaO
                * b.unit_model.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.day,
            ),
            to_units=pyunits.kg / pyunits.year,
        )

    @blk.Expression(doc="Costing block MgCl2 dosing")
    def mgcl2_dosing(b):
        return pyunits.convert(
            b.unit_model.MgCl2_dosing
            * b.unit_model.MgCl2_mw
            / b.unit_model.Mg_mw  # to convert to MgCl2 solid
            * b.unit_model.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyunits.kg / pyunits.year,
        )

    @blk.Expression(doc="Costing block CO2 dosing")
    def co2_dosing(b):
        return pyunits.convert(
            b.unit_model.CO2_first_basin + b.unit_model.CO2_second_basin,
            to_units=pyunits.kg / pyunits.year,
        )

    blk.costing_package.cost_flow(blk.electricity_flow, "electricity")

    blk.costing_package.cost_flow(blk.cao_dosing, "lime")
    blk.costing_package.cost_flow(blk.unit_model.Na2CO3_dosing, "soda_ash")
    blk.costing_package.cost_flow(blk.mgcl2_dosing, "mgcl2")
    blk.costing_package.cost_flow(blk.co2_dosing, "co2")
