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
        initialize=0.171,
        units=pyunits.USD_2020 / pyunits.kg,
        doc="Cost of CaO $/kg",
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
        initialize=0.65,
        units=pyunits.USD_2020 / pyunits.kg,
        doc="Cost of Na2CO3 $/kg",
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
        initialize=0.55,
        units=pyunits.USD_2020 / pyunits.kg,
        doc="Cost of MgCl2 $/kg",
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
    blk.mix_tank_capital_coeff_1 = Param(
        initialize=0.0002,
        mutable=True,
        doc="Coefficient of first term in mixer tank capital cost equation",
    )

    blk.mix_tank_capital_coeff_2 = Param(
        initialize=22.776,
        mutable=True,
        doc="Coefficient of second term in mixer tank capital cost equation",
    )

    blk.mix_tank_capital_exp_1 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of the first term in the mixer tank capital cost equation",
    )

    blk.mix_tank_capital_constant = Param(
        initialize=28584,
        mutable=True,
        doc="Constant in the mixer tank capital cost equation",
    )

    blk.floc_tank_capital_coeff = Param(
        initialize=673894,
        mutable=True,
        doc="Coefficient of flocculation tank capital cost equation",
    )

    blk.floc_tank_capital_constant = Param(
        initialize=217222,
        mutable=True,
        doc="Constant in the flocculation tank capital cost equation",
    )

    blk.sed_basin_capital_coeff_1 = Param(
        initialize=-0.0005,
        mutable=True,
        doc="Coefficient of first term in sedimentation basin capital cost equation",
    )

    blk.sed_basin_capital_coeff_2 = Param(
        initialize=86.89,
        mutable=True,
        doc="Coefficient of second term in sedimentation basin capital cost equation",
    )

    blk.sed_basin_capital_exp_1 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of the first term in the sedimentation basin capital cost equation",
    )

    blk.sed_basin_capital_constant = Param(
        initialize=182801,
        mutable=True,
        doc="Constant in the sedimentation basin capital cost equation",
    )

    blk.recarb_basin_capital_coeff_1 = Param(
        initialize=4e-9,
        mutable=True,
        doc="Coefficient of first term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_capital_coeff_2 = Param(
        initialize=-0.0002,
        mutable=True,
        doc="Coefficient of second term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_capital_coeff_3 = Param(
        initialize=10.027,
        mutable=True,
        doc="Coefficient of third term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_capital_exp_1 = Param(
        initialize=3,
        mutable=True,
        doc="Exponent of first term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_capital_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of second term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_capital_constant = Param(
        initialize=19287,
        mutable=True,
        doc="Constant in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_source_capital_coeff_1 = Param(
        initialize=9e-8,
        mutable=True,
        doc="Coefficient of first term in recarbonation basin capital cost equation",
    )

    blk.recarb_basin_source_capital_coeff_2 = Param(
        initialize=-0.001,
        mutable=True,
        doc="Coefficient of second term in recarbonation basin source capital cost equation",
    )

    blk.recarb_basin_source_capital_coeff_3 = Param(
        initialize=42.578,
        mutable=True,
        doc="Coefficient of third term in recarbonation basin source capital cost equation",
    )

    blk.recarb_basin_source_capital_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of second term in recarbonation basin source capital cost equation",
    )

    blk.recarb_basin_source_capital_constant = Param(
        initialize=130812,
        mutable=True,
        doc="Constant in recarbonation basin source capital cost equation",
    )

    blk.lime_feed_system_capital_coeff = Param(
        initialize=20.065,
        mutable=True,
        doc="Coefficient in lime feed system capital cost",
    )

    blk.lime_feed_system_capital_constant = Param(
        initialize=193268, mutable=True, doc="Constant in lime feed system capital cost"
    )

    blk.admin_capital_coeff = Param(
        initialize=69195,
        mutable=True,
        doc="Coefficient in lime feed system capital cost",
    )

    blk.admin_capital_exp = Param(
        initialize=0.5523, mutable=True, doc="Exponent in lime feed system capital cost"
    )

    # Costing equation coefficients - Operational cost
    blk.mix_tank_op_coeff_1 = Param(
        initialize=-3e-8,
        mutable=True,
        doc="Coefficient of first term in mixer tank operational cost equation",
    )

    blk.mix_tank_op_coeff_2 = Param(
        initialize=0.0008,
        mutable=True,
        doc="Coefficient of second term in mixer tank operational cost equation",
    )

    blk.mix_tank_op_coeff_3 = Param(
        initialize=2.8375,
        mutable=True,
        doc="Coefficient of third term in mixer tank operational cost equation",
    )

    blk.mix_tank_op_exp_1 = Param(
        initialize=3,
        mutable=True,
        doc="Exponent of the first term in the mixer tank operational cost equation",
    )

    blk.mix_tank_op_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of the second term in the mixer tank operational cost equation",
    )

    blk.mix_tank_op_constant = Param(
        initialize=22588,
        mutable=True,
        doc="Constant in the mixer tank operational cost equation",
    )

    blk.floc_tank_op_coeff_1 = Param(
        initialize=3e-13,
        mutable=True,
        doc="Coefficient of the first term in the flocculation tank capital cost equation",
    )

    blk.floc_tank_op_coeff_2 = Param(
        initialize=-4e-7,
        mutable=True,
        doc="Coefficient of the second term in the flocculation tank capital cost equation",
    )

    blk.floc_tank_op_coeff_3 = Param(
        initialize=0.318,
        mutable=True,
        doc="Coefficient of the third term in the flocculation tank capital cost equation",
    )

    blk.floc_tank_op_exp_1 = Param(
        initialize=3,
        mutable=True,
        doc="Exponent of the first term in the flocculation tank capital cost equation",
    )

    blk.floc_tank_op_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of the second term in the flocculation tank capital cost equation",
    )

    blk.floc_tank_op_constant = Param(
        initialize=6040,
        mutable=True,
        doc="Constant in the flocculation tank capital cost equation",
    )

    blk.sed_basin_op_coeff_1 = Param(
        initialize=7e-10,
        mutable=True,
        doc="Coefficient of first term in sedimentation basin operational cost equation",
    )

    blk.sed_basin_op_coeff_2 = Param(
        initialize=-0.00005,
        mutable=True,
        doc="Coefficient of second term in sedimentation basin operational cost equation",
    )

    blk.sed_basin_op_coeff_3 = Param(
        initialize=1.5908,
        mutable=True,
        doc="Coefficient of second term in sedimentation basin operational cost equation",
    )

    blk.sed_basin_op_exp_1 = Param(
        initialize=3,
        mutable=True,
        doc="Exponent of the first term in the sedimentation basin operational cost equation",
    )

    blk.sed_basin_op_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of the first term in the sedimentation basin operational cost equation",
    )

    blk.sed_basin_op_constant = Param(
        initialize=6872,
        mutable=True,
        doc="Constant in the sedimentation basin operational cost equation",
    )

    blk.recarb_basin_op_coeff_1 = Param(
        initialize=1e-8,
        mutable=True,
        doc="Coefficient of first term in recarbonation basin operational cost equation",
    )

    blk.recarb_basin_op_coeff_2 = Param(
        initialize=-0.0004,
        mutable=True,
        doc="Coefficient of second term in recarbonation basin operational cost equation",
    )

    blk.recarb_basin_op_coeff_3 = Param(
        initialize=6.19,
        mutable=True,
        doc="Coefficient of third term in recarbonation basin operational cost equation",
    )

    blk.recarb_basin_op_exp_1 = Param(
        initialize=3,
        mutable=True,
        doc="Exponent of first term in recarbonation basin operational cost equation",
    )

    blk.recarb_basin_op_exp_2 = Param(
        initialize=2,
        mutable=True,
        doc="Exponent of second term in recarbonation basin operational cost equation",
    )

    blk.recarb_basin_op_constant = Param(
        initialize=10265,
        mutable=True,
        doc="Constant in recarbonation basin operational cost equation",
    )

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
    parameter_block_name="soda ash",
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

    blk.mix_tank_op_cost = Var(
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
    blk.mix_tank_capital_cost_constraint = Constraint(
        expr=blk.mix_tank_capital_cost
        == (
            chem_soft.mix_tank_capital_coeff_1
            * (pyunits.convert(blk.unit_model.volume_mixer, to_units=pyunits.ft**3))
            ** chem_soft.mix_tank_capital_exp_1
        )
        + (
            chem_soft.mix_tank_capital_coeff_2
            * (pyunits.convert(blk.unit_model.volume_mixer, to_units=pyunits.ft**3))
            + chem_soft.mix_tank_capital_constant
        )
    )

    capital_cost_expr += blk.mix_tank_capital_cost

    # Flocculation tank
    blk.floc_tank_capital_cost_constraint = Constraint(
        expr=blk.floc_tank_capital_cost
        == (
            chem_soft.floc_tank_capital_coeff
            * pyunits.convert(blk.unit_model.volume_floc, to_units=pyunits.megagallon)
            + chem_soft.floc_tank_capital_constant
        )
    )

    capital_cost_expr += blk.floc_tank_capital_cost

    # Sedimentation basin
    blk.sed_basin_capital_cost_constraint = Constraint(
        expr=blk.sed_basin_capital_cost
        == (
            chem_soft.sed_basin_capital_coeff_1
            * (
                pyunits.convert(blk.unit_model.volume_sed, to_units=pyunits.ft**3)
                / pyunits.convert(chem_soft.sed_basin_depth, to_units=pyunits.ft)
            )
            ** chem_soft.sed_basin_capital_exp_1
            + chem_soft.sed_basin_capital_coeff_2
            * pyunits.convert(blk.unit_model.volume_sed, to_units=pyunits.ft**3)
            / pyunits.convert(chem_soft.sed_basin_depth, to_units=pyunits.ft)
            + chem_soft.sed_basin_capital_constant
        )
    )

    capital_cost_expr += blk.sed_basin_capital_cost

    # Recarbonation basin
    blk.recarb_basin_capital_cost_constraint = Constraint(
        expr=blk.recarb_basin_capital_cost
        == (
            chem_soft.recarb_basin_capital_coeff_1
            * (pyunits.convert(blk.unit_model.volume_recarb, to_units=pyunits.ft**3))
            ** chem_soft.recarb_basin_capital_exp_1
            + chem_soft.recarb_basin_capital_coeff_2
            * (pyunits.convert(blk.unit_model.volume_recarb, to_units=pyunits.ft**3))
            ** chem_soft.recarb_basin_capital_exp_2
            + chem_soft.recarb_basin_capital_coeff_3
            * (pyunits.convert(blk.unit_model.volume_recarb, to_units=pyunits.ft**3))
            + chem_soft.recarb_basin_capital_constant
        )
    )

    capital_cost_expr += blk.recarb_basin_capital_cost

    # Recarbonation source cost
    blk.recarb_basin_source_capital_cost_constraint = Constraint(
        expr=blk.recarb_basin_source_capital_cost
        == (
            chem_soft.recarb_basin_source_capital_coeff_1
            * (
                pyunits.convert(
                    blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,
                    to_units=pyunits.lb / pyunits.day,
                )
            )
            + chem_soft.recarb_basin_source_capital_coeff_2
            * (
                pyunits.convert(
                    blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,
                    to_units=pyunits.lb / pyunits.day,
                )
            )
            ** chem_soft.recarb_basin_source_capital_exp_2
            + chem_soft.recarb_basin_source_capital_coeff_3
            * (
                pyunits.convert(
                    blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,
                    to_units=pyunits.lb / pyunits.day,
                )
            )
            + chem_soft.recarb_basin_source_capital_constant
        )
    )

    capital_cost_expr += blk.recarb_basin_source_capital_cost

    # Lime feed system cost
    blk.lime_feed_system_capital_cost_constraint = Constraint(
        expr=blk.lime_feed_system_capital_cost
        == (
            chem_soft.lime_feed_system_capital_coeff
            * (
                pyunits.convert(
                    blk.unit_model.CaO_dosing, to_units=pyunits.lb / pyunits.hour
                )
            )
            + chem_soft.lime_feed_system_capital_constant
        )
    )

    capital_cost_expr += blk.lime_feed_system_capital_cost

    # Admin cost
    blk.admin_capital_cost_constraint = Constraint(
        expr=blk.admin_capital_cost
        == chem_soft.admin_capital_coeff
        * (
            (
                pyunits.convert(
                    blk.unit_model.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.megagallon / pyunits.day,
                )
            )
            ** chem_soft.admin_capital_exp
        )
    )

    capital_cost_expr += blk.admin_capital_cost

    # Operational cost component constraints

    op_cost_expr = 0

    # Mixing tank
    blk.mix_tank_op_cost_constraint = Constraint(
        expr=blk.mix_tank_op_cost
        == (
            chem_soft.mix_tank_op_coeff_1
            * (pyunits.convert(blk.unit_model.volume_mixer, to_units=pyunits.ft**3))
            ** chem_soft.mix_tank_op_exp_1
        )
        + (
            chem_soft.mix_tank_op_coeff_2
            * (pyunits.convert(blk.unit_model.volume_mixer, to_units=pyunits.ft**3))
            ** chem_soft.mix_tank_op_exp_2
        )
        + (
            chem_soft.mix_tank_op_coeff_3
            * (pyunits.convert(blk.unit_model.volume_mixer, to_units=pyunits.ft**3))
        )
        + chem_soft.mix_tank_op_constant
    )

    op_cost_expr += blk.mix_tank_op_cost

    # Flocculation tank
    blk.floc_tank_op_cost_constraint = Constraint(
        expr=blk.floc_tank_op_cost
        == (
            (
                chem_soft.floc_tank_op_coeff_1
                * (
                    pyunits.convert(
                        blk.unit_model.volume_floc, to_units=pyunits.ft**3
                    )
                )
                ** chem_soft.floc_tank_op_exp_1
            )
            + (
                chem_soft.floc_tank_op_coeff_2
                * (
                    pyunits.convert(
                        blk.unit_model.volume_floc, to_units=pyunits.ft**3
                    )
                )
                ** chem_soft.floc_tank_op_exp_2
            )
            + (
                chem_soft.floc_tank_op_coeff_3
                * (
                    pyunits.convert(
                        blk.unit_model.volume_floc, to_units=pyunits.ft**3
                    )
                )
                + chem_soft.floc_tank_op_constant
            )
        )
    )

    op_cost_expr += blk.floc_tank_op_cost

    # Sedimentation basin
    blk.sed_basin_op_cost_constraint = Constraint(
        expr=blk.sed_basin_op_cost
        == (
            chem_soft.sed_basin_op_coeff_1
            * (
                (pyunits.convert(blk.unit_model.volume_sed, to_units=pyunits.ft**3))
                / pyunits.convert(chem_soft.sed_basin_depth, to_units=pyunits.ft)
            )
            ** chem_soft.sed_basin_op_exp_1
        )
        + (
            chem_soft.sed_basin_op_coeff_2
            * (
                (pyunits.convert(blk.unit_model.volume_sed, to_units=pyunits.ft**3))
                / pyunits.convert(chem_soft.sed_basin_depth, to_units=pyunits.ft)
            )
            ** chem_soft.sed_basin_op_exp_2
        )
        + (
            chem_soft.sed_basin_op_coeff_3
            * (pyunits.convert(blk.unit_model.volume_sed, to_units=pyunits.ft**3))
            / pyunits.convert(chem_soft.sed_basin_depth, to_units=pyunits.ft)
            + chem_soft.sed_basin_op_constant
        )
    )

    op_cost_expr += blk.sed_basin_op_cost

    # Recarbonation basin
    blk.recarb_basin_op_cost_constraint = Constraint(
        expr=blk.recarb_basin_op_cost
        == (
            (
                chem_soft.recarb_basin_op_coeff_1
                * (
                    pyunits.convert(
                        blk.unit_model.CO2_first_basin
                        + blk.unit_model.CO2_second_basin,
                        to_units=pyunits.lb / pyunits.day,
                    )
                )
                ** chem_soft.recarb_basin_op_exp_1
            )
            + (
                chem_soft.recarb_basin_op_coeff_2
                * (
                    pyunits.convert(
                        blk.unit_model.CO2_first_basin
                        + blk.unit_model.CO2_second_basin,
                        to_units=pyunits.lb / pyunits.day,
                    )
                )
                ** chem_soft.recarb_basin_op_exp_2
            )
            + (
                chem_soft.recarb_basin_op_coeff_3
                * (
                    pyunits.convert(
                        blk.unit_model.CO2_first_basin
                        + blk.unit_model.CO2_second_basin,
                        to_units=pyunits.lb / pyunits.day,
                    )
                )
                + chem_soft.recarb_basin_op_constant
            )
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
            * b.unit_model.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyunits.kg / pyunits.year,
        )

    @blk.Expression(doc="Costing block co2 dosing")
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
