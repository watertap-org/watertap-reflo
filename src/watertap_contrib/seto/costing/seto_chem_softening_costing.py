from watertap.core.zero_order_costing import ZeroOrderCostingData
from idaes.core import declare_process_block_class
import pyomo.environ as pyo
from idaes.core.util.misc import StrEnum
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


def build_chem_softening_cost_param_block(blk):
    '''
        Parameters and variables to be used in the costing model
    '''
    #
    # blk.unit_cost = pyo.Var(
    #     initialize = 10,
    #     doc = 'sample',
    #     units = pyo.dimensionless,
    # )


    blk.sed_basin_depth = pyo.Var(
        initialize=4.5,
        bounds=(3, 6),
        units=pyo.units.m,
        doc='Depth of sedimentation basin')





@register_costing_parameter_block(
    build_rule = build_chem_softening_cost_param_block,
    parameter_block_name = "chemical_softening"
    )

def cost_chem_softening(blk, softeningType = SofteningProcedureType.single_stage_lime):
    '''
        Function for costing chemical softening unit for 
        1. Lime softening
    '''
    if (softeningType == SofteningProcedureType.single_stage_lime):
        cost_chem_softening_single_stage_lime(blk)






def cost_chem_softening_single_stage_lime(blk):
    '''
        Capital cost for single stage lime softening
    '''

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    # Capital cost components 
    blk.mix_tank_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of mixing tank',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.floc_tank_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of flocculation tank',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.sed_basin_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of sedimentation basin',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.recarb_basin_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of recarbonation basin',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.recarb_basin_source_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of recarbonation tank - liquid CO2 as CO2 source',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.lime_feed_system_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of lime feed system',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )

    blk.admin_cost = pyo.Var(
        initialize = 1e5,
        doc = 'Cost of administrative, laboratory, and maintenance building',
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency
    )


    # Capital cost component constraints

    # Mixing tank
    if blk.unit_model.vel_gradient_mix <= 300:
        blk.mix_tank_cost_constraint = pyo.Constraint(
            expr = blk.mix_tank_cost 
            ==  pyo.units.convert( (0.0002 * (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) ** 2) + 22.776 * (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) + 28584
            , to_units = blk.costing_package.base_currency
            )
        )
    elif blk.unit_model.vel_gradient_mix > 300 and blk.unit_model.vel_gradient_mix <= 600:
        blk.mix_tank_cost_constraint = pyo.Constraint(
            expr = blk.mix_tank_cost 
            ==  pyo.units.convert( (0.0002 * (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) ** 2) + 29.209 * (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) + 30388
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_mix > 600 and blk.unit_model.vel_gradient_mix <= 1000:
        blk.mix_tank_cost_constraint = pyo.Constraint(
            expr = blk.mix_tank_cost 
            ==  pyo.units.convert(  (0.0002 * (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) ** 2) + 55.443* (pyo.units.convert(blk.unit_model.volume_mixer,to_units=pyo.units.ft**3)) + 29756
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('ERROR - OUT OF BOUNDS (300 - 1000 s-1)')


    # Flocculation tank
    if blk.unit_model.vel_gradient_floc <= 20:
        blk.floc_tank_cost_constraint = pyo.Constraint(
            expr = blk.floc_tank_cost 
            ==  pyo.units.convert( (566045 * (blk.unit_model.volume_floc / 3785) + 224745)
            , to_units = blk.costing_package.base_currency
            )
        )
    elif blk.unit_model.vel_gradient_floc > 20 and blk.unit_model.vel_gradient_floc <= 50:
        blk.floc_tank_cost_constraint = pyo.Constraint(
            expr = blk.floc_tank_cost 
            ==  pyo.units.convert( (673894 * (blk.unit_model.volume_floc / 3785) + 217222)
            , to_units = blk.costing_package.base_currency
            )
        )

    elif blk.unit_model.vel_gradient_floc > 50 and blk.unit_model.vel_gradient_floc <= 80:
        blk.floc_tank_cost_constraint = pyo.Constraint(
            expr = blk.floc_tank_cost 
            ==  pyo.units.convert( (952902 * (blk.unit_model.volume_floc / 3785) + 177335)
            , to_units = blk.costing_package.base_currency
            )
        )
    else:
        print('Error - Out of bounds (20-80 s-1')


    # Sedimentation basin 
    blk.sed_basin_cost_constraint = pyo.Constraint(
        expr = blk.sed_basin_cost 
        == pyo.units.convert(
            (- 0.0005 * ( pyo.units.convert(blk.unit_model.volume_sed,to_units=pyo.units.ft**3) / pyo.units.convert(blk.chemical_softening.sed_basin_depth,to_units=pyo.units.ft) ) ** 2
            + 86.89 * (blk.unit_model.volume_sed * 35.315) / pyo.units.convert(blk.chemical_softening.sed_basin_depth,to_units=pyo.units.ft) + 182801 )
            , to_units = blk.costing_package.base_currency 
        )
    )

    # Recarbonation basin
    blk.recarb_basin_cost_constraint =  pyo.Constraint(
        expr = blk.sed_basin_cost 
        == pyo.units.convert(
            ( 4e-9 * (pyo.units.convert(blk.unit_model.volume_recarb,to_units=pyo.units.ft**3))**3 
            - 0.0002 * (pyo.units.convert(blk.unit_model.volume_recarb,to_units=pyo.units.ft**3)) ** 2 
            + 10.027 * (pyo.units.convert(blk.unit_model.volume_recarb,to_units=pyo.units.ft**3)) + 19287 )
            , to_units = blk.costing_package.base_currency
        )
    )

    # Recarbonation source cost


    # Lime feed system cost


    # Admin cost 



    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.mix_tank_cost + blk.floc_tank_cost + blk.sed_basin_cost + blk.recarb_basin_cost + blk.recarb_basin_source_cost + blk.lime_feed_system_cost + blk.admin_cost 
    )


    
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        ==
    )