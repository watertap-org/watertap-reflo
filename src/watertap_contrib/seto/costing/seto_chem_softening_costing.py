from watertap.core.zero_order_costing import ZeroOrderCostingData
from idaes.core import declare_process_block_class

from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    NonNegativeReals,
    log,
    Var,
    Param,
    Set,
    Suffix,
    value,
    check_optimal_termination,
    units as pyunits,
)

from idaes.core.util.misc import StrEnum
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


# TODO: Remove the hard coded constants in the costing equations

class SofteningProcedureType(StrEnum): 
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


def build_chem_softening_cost_param_block(blk):
    '''
        Parameters and variables to be used in the costing model
    '''
    costing = blk.parent_block()
    #
    # blk.unit_cost = Var(
    #     initialize = 10,
    #     doc = 'sample',
    #     units = dimensionless,
    # )


    blk.sed_basin_depth = Var(
        initialize=4.5,
        bounds=(3, 6),
        units= pyunits.m,
        doc='Depth of sedimentation basin'
        )

    blk.sludge_disposal_cost = Var (
        intialize = 35,
        units = costing.base_currency/pyunits.tons,
        doc = "Cost of sludge disposal $/tonne"
    )

    blk.MgCl2_cost = Var (
        intialize = 1.5,
        units = costing.base_currency/pyunits.kg,
        doc = "Cost of MgCl2 $/kg"
    )

    blk.Na2CO3_cost = Var (
        intialize = 0.65,
        units = costing.base_currency/pyunits.kg,
        doc = "Cost of soda ash $/kg"
    )


    blk.fix_all_vars()

@register_costing_parameter_block(
    build_rule = build_chem_softening_cost_param_block,
    parameter_block_name = "chemical_softening"
    )

def cost_chem_softening(blk, softeningType = SofteningProcedureType.single_stage_lime):
    '''
        Function for costing chemical softening unit for 
        1. Lime softening
    '''
    if (softeningType == blk.config.softening_train_type):
        cost_chem_softening_single_stage_lime(blk)




def cost_chem_softening_single_stage_lime(blk):
    '''
        Capital cost for single stage lime softening
    '''

    chem_soft = blk.costing_package.chem_softening
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    # Capital cost components 
    blk.mix_tank_capital_cost = Var(
        initialize = 1e6,
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc = 'Capital Cost of mixing tank'
    )

    blk.floc_tank_capital_cost = Var(
        initialize = 1e6,
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc = 'Capital cost of flocculation tank'
    )

    blk.sed_basin_capital_cost = Var(
        initialize = 1e6,
        doc = 'Capital cost of sedimentation basin',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.recarb_basin_capital_cost = Var(
        initialize = 1e6,
        doc = 'Capital cost of recarbonation basin',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.recarb_basin_source_capital_cost = Var(
        initialize = 1e6,
        doc = 'Capital cost of recarbonation tank - liquid CO2 as CO2 source',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.lime_feed_system_capital_cost = Var(
        initialize = 1e6,
        doc = 'Capital cost of lime feed system',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.admin_capital_cost = Var(
        initialize = 1e6,
        doc = 'Capital cost of administrative, laboratory, and maintenance building',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    # Operating cost components
   
    blk.mix_tank_op_cost = Var(
        initialize = 1e4,
        doc = 'Operational cost of mixing tank - energy consumption and labor',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.floc_tank_op_cost = Var(
        initialize = 1e4,
        doc = 'Operational cost of flocculation tank - energy consumption and labor',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.sed_basin_op_cost = Var(
        initialize = 1e4,
        doc = 'Operational cost of sedimentation basin - energy consumption, labor, and maintenance',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.recarb_basin_op_cost = Var(
        initialize = 1e4,
        doc = 'Operationg cost of recarbonation basin - energy consumption, labor, and maintenance',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.lime_feed_op_cost= Var(
        initialize = 1e4,
        doc = 'Operation cost of lime feed system - energy consumption, labor, and maintenance',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    blk.soda_ash_add_op_cost= Var(
        initialize = 1e4,
        doc = 'Operation cost of adding soda ash - chemical addition cost',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    blk.mg_add_op_cost= Var(
        initialize = 1e4,
        doc = 'Operational cost of adding Mg - chemical addition cost',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    blk.lime_sludge_mngt_op_cost= Var(
        initialize = 1e5,
        doc = 'Operational cost of lime sludge management',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    blk.admin_op_cost= Var(
        initialize = 1e5,
        doc = 'Operational cost of administrative, laboratory and maintenance building',
        domain= NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    #region Capital cost component constraints

    # Mixing tank
    if blk.unit_model.vel_gradient_mix <= 300:
        blk.mix_tank_capital_cost_constraint = Constraint(
            expr = blk.mix_tank_capital_cost 
            ==  pyunits.convert( 
                        (0.0002 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) ** 2) + 
                        (22.776 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) + 28584 )
            , to_units = blk.costing_package.base_currency
            )
        )
    elif blk.unit_model.vel_gradient_mix > 300 and blk.unit_model.vel_gradient_mix <= 600:
        blk.mix_tank_capital_cost_constraint = Constraint(
            expr = blk.mix_tank_capital_cost 
            ==  pyunits.convert( 
                        (0.0002 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) ** 2) + 
                        (29.209 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) + 30388)
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_mix > 600 and blk.unit_model.vel_gradient_mix <= 1000:
        blk.mix_tank_capital_cost_constraint = Constraint(
            expr = blk.mix_tank_capital_cost 
            ==  pyunits.convert(  
                        (0.0002 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) ** 2) + 
                        (55.443* (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) + 29756)
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('Error - Out of bounds (300 - 1000 s-1)')


    # Flocculation tank
    if blk.unit_model.vel_gradient_floc <= 20:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( (566045 * (blk.unit_model.volume_floc / 3785) + 224745)
            , to_units = blk.costing_package.base_currency
            )
        )
    elif blk.unit_model.vel_gradient_floc > 20 and blk.unit_model.vel_gradient_floc <= 50:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( (673894 * (blk.unit_model.volume_floc / 3785) + 217222)
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_floc > 50 and blk.unit_model.vel_gradient_floc <= 80:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( (952902 * (blk.unit_model.volume_floc / 3785) + 177335)
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('Error - Out of bounds (20-80 s-1')


    # Sedimentation basin 
    blk.sed_basin_capital_cost_constraint = Constraint(
        expr = blk.sed_basin_capital_cost 
        == pyunits.convert(
                (- 0.0005 * ( pyunits.convert(blk.unit_model.volume_sed,to_units=pyunits.ft**3) / pyunits.convert(chem_soft.sed_basin_depth,to_units=pyunits.ft) ) ** 2
                + 86.89 * (blk.unit_model.volume_sed * 35.315) / pyunits.convert(chem_soft.sed_basin_depth,to_units=pyunits.ft) + 182801 )
            , to_units = blk.costing_package.base_currency 
        )
    )

    # Recarbonation basin
    blk.recarb_basin_capital_cost_constraint =  Constraint(
        expr = blk.recarb_basin_capital_cost 
        == pyunits.convert(
                    ( 4e-9 * (pyunits.convert(blk.unit_model.volume_recarb,to_units=pyunits.ft**3))**3 -
                    0.0002 * (pyunits.convert(blk.unit_model.volume_recarb,to_units=pyunits.ft**3)) ** 2 +
                    10.027 * (pyunits.convert(blk.unit_model.volume_recarb,to_units=pyunits.ft**3)) + 19287 )
            ,to_units = blk.costing_package.base_currency
        )
    )

    # Recarbonation source cost
    blk.recarb_basin_source_capital_cost_constraint = Constraint(
        expr = blk.recarb_basin_source_capital_cost 
        == pyunits.convert(
                    (9e-8 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin ,to_units = pyunits.lb / pyunits.day)) - 
                    0.001 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin ,to_units = pyunits.lb / pyunits.day)) ** 2 + 
                    42.578 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,to_units = pyunits.lb / pyunits.day)) + 130812 )
            , to_units = blk.costing_package.base_currency
        )
    )

    # Lime feed system cost
    blk.lime_feed_system_capital_cost_constraint = Constraint(
        expr = blk.lime_feed_system_capital_cost
        == pyunits.convert(
             (20.065 * (pyunits.convert(blk.unit_model.CaOH,to_units = pyunits.lb / pyunits.hour)) + 193268)
            , to_units = blk.costing_package.base_currency
        )
    )

    # Admin cost 
    blk.admin_capital_cost_constraint = Constraint(
        expr = blk.admin_capital_cost
        == pyunits.convert(
             (69195 * ((blk.unit_model.properties_in[0].flow_vol / 3785)**0.5523))
            , to_units = blk.costing_package.base_currency
        )
    )

    #endregion

    #region Operational cost component constraints

    # Mixing tank
    if blk.unit_model.vel_gradient_mix <= 300:
        blk.mix_tank_op_cost_constraint = Constraint(
            expr = blk.mix_tank_op_cost 
            ==  pyunits.convert( 
                        (-3e-8  * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3))** 3) + 
                        (0.0008 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3))** 2) + 
                        (2.8375 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3))) + 22588
            , to_units = blk.costing_package.base_currency
            )
        )
    elif blk.unit_model.vel_gradient_mix > 300 and blk.unit_model.vel_gradient_mix <= 600:
        blk.mix_tank_op_cost_constraint = Constraint(
            expr = blk.mix_tank_op_cost 
            ==  pyunits.convert( 
                        (-3e-8  * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) ** 3) + 
                        (0.0008 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**2))** 2) + 
                        (7.8308 * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) + 22588)
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_mix > 600 and blk.unit_model.vel_gradient_mix <= 1000:
        blk.mix_tank_op_cost_constraint = Constraint(
            expr = blk.mix_tank_op_cost 
            ==  pyunits.convert(  
                        (36.096  * (pyunits.convert(blk.unit_model.volume_mixer,to_units=pyunits.ft**3)) ** 2 + 18928)
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('Error - Out of bounds (300 - 1000 s-1)')


    # Flocculation tank
    if blk.unit_model.vel_gradient_floc <= 20:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( 
                        ( (3e-13 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) ** 3) +
                        (-5e-7  * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) ** 2) + 
                        (0.2757 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) + 6594) )
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_floc > 20 and blk.unit_model.vel_gradient_floc <= 50:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( 
                        ( (3e-13 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) ** 3) +
                        (-4e-7  * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) ** 2) + 
                        (0.318 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) + 6040) )
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_floc > 50 and blk.unit_model.vel_gradient_floc <= 80:
        blk.floc_tank_capital_cost_constraint = Constraint(
            expr = blk.floc_tank_capital_cost 
            ==  pyunits.convert( 
                        ( (-3e-7 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) ** 2) + 
                        (0.5692 * (pyunits.convert(blk.unit_model.volume_floc,to_units=pyunits.ft**3)) + 6748) )
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('Error - Out of bounds (20-80 s-1')

    
    blk.sed_basin_op_cost_constraint = Constraint(
        expr = blk.sed_basin_op_cost
        == pyunits.convert(
            (((7e-10 * pyunits.convert(blk.unit_model.volume_sed,to_units=pyunits.ft**3))/pyunits.convert(chem_soft.sed_basin_depth,to_units=pyunits.ft)) ** 3) - 
            (((0.00005 * pyunits.convert(blk.unit_model.volume_sed,to_units=pyunits.ft**3))/pyunits.convert(chem_soft.sed_basin_depth,to_units=pyunits.ft)) ** 2) +
            ((1.5908 * pyunits.convert(blk.unit_model.volume_sed,to_units=pyunits.ft**3))/pyunits.convert(chem_soft.sed_basin_depth,to_units=pyunits.ft) + 6872)
            , to_units = blk.costing_package.base_currency
        )
    ) 


    blk.recarb_basin_op_cost_constraint = Constraint(
        expr = blk.recarb_basin_op_cost
        == pyunits.convert(
            ((1e-8 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,to_units = pyunits.lb / pyunits.day))**3) - 
            (0.0004 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,to_units = pyunits.lb / pyunits.day)) ** 2) + 
            (6.19 * (pyunits.convert(blk.unit_model.CO2_first_basin + blk.unit_model.CO2_second_basin,to_units = pyunits.lb / pyunits.day)) + 10265))
            , to_units = blk.costing_package.base_currency
        )
    ) 

    blk.lime_feed_op_cost_constraint = Constraint(
        expr = blk.lime_feed_op_cost
        == pyunits.convert(
             (4616.7 * (pyunits.convert(blk.unit_model.CaOH, to_units = pyunits.lb/pyunits.day))**0.4589)
            , to_units = blk.costing_package.base_currency
        )
    )  


    blk.soda_ash_add_op_cost_constraint = Constraint(
        expr = blk.soda_ash_add_op_cost
        == pyunits.convert(
            (pyunits.convert(blk.unit_model.Na2CO3, to_units = pyunits.kg/pyunits.year)*chem_soft.NaCO3_cost)
            , to_units = blk.costing_package.base_currency
        )
    ) 

 
    blk.mg_add_op_cost_constraint = Constraint(
        expr = blk.mg_add_op_cost
        == pyunits.convert(
            (pyunits.convert(blk.unit_model.MgCl2, to_units = pyunits.kg/pyunits.year)*chem_soft.MgCl2_cost)
            , to_units = blk.costing_package.base_currency
        )
    ) 


    blk.lime_sludge_mngt_op_cost_constraint = Constraint(
        expr = blk.lime_sludge_mngt_op_cost
        == pyunits.convert(
            (pyunits.convert(blk.unit_model.sludge_prod, to_units = pyunits.tons/pyunits.year)* chem_soft.sludge_disposal_cost)
            , to_units = blk.costing_package.base_currency
        )
    ) 


    blk.admin_op_cost_constraint = Constraint(
        expr = blk.admin_op_cost
        == pyunits.convert(
             (88589 * ((blk.unit_model.properties_in[0].flow_vol / 3785)**0.4589))
            , to_units = blk.costing_package.base_currency
        )
    ) 

    #endregion


    # Sum of all capital costs
    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == blk.mix_tank_capital_cost + blk.floc_tank_capital_cost + blk.sed_basin_capital_cost + blk.recarb_basin_capital_cost + 
        blk.recarb_basin_source_capital_cost + blk.lime_feed_system_capital_cost + blk.admin_capital_cost 
    )


    # Sum of all operational costs
    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        ==  blk.mix_tank_op_cost + blk.floc_tank_op_cost + blk.sed_basin_op_cost + blk.recarb_basin_op_cost +
        blk.lime_feed_op_cost + blk.soda_ash_add_op_cost + blk.mg_add_op_cost + blk.lime_sludge_mngt_op_cost + blk.admin_op_cost
    )