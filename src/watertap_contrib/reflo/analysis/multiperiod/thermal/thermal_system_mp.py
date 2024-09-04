from pyomo.environ import (
    ConcreteModel, 
    Var, 
    value, 
    Objective, 
    assert_optimal_termination,
    units as pyunits
)
from idaes.core.util.model_statistics import *
from idaes.core import FlowsheetBlock
import logging
import pandas as pd
import datetime
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from watertap_contrib.reflo.analysis.multiperiod.thermal.thermal_process_steady_state_flowsheet import (
    create_mp_steady_state, 
    fix_dof_and_initialize,
    print_results
)

from watertap.core.solvers import get_solver


# Can optimize storage tank volume
def unfix_dof(m):
    m.fs.previous_hx_solar_hot_outlet_temperature.unfix()
    m.fs.previous_fpc_outlet_temperature.unfix()
    m.fs.previous_tes_tank_temp.unfix()
    m.fs.previous_hx_solar_cold_outlet_temperature.unfix()
    m.fs.previous_process_outlet_temperature.unfix()
    m.fs.previous_acc_grid_duty.unfix()
    m.fs.previous_grid_duty.unfix()

    return None

def get_variable_pairs(t1,t2):
    """
    This function returns pairs of variables that need to be connected across two time periods

    Args:
        t1: current time block
        t2: next time block

    Returns:
        None
    """

    return [

        (t1.fs.hx_solar.hot_side_outlet.temperature[0], t2.fs.previous_hx_solar_hot_outlet_temperature),
        (t1.fs.fpc.outlet.temperature[0], t2.fs.previous_fpc_outlet_temperature),
        (t1.fs.tes.tes_temperature[0], t2.fs.previous_tes_tank_temp),
        (t1.fs.hx_solar.cold_side_outlet.temperature[0], t2.fs.previous_hx_solar_cold_outlet_temperature),
        (t1.fs.tes.tes_process_outlet.temperature[0], t2.fs.previous_process_outlet_temperature),
        (t1.fs.acc_grid_duty, t2.fs.previous_acc_grid_duty),
        (t1.fs.grid_heater.heat_duty[0], t2.fs.previous_grid_duty),
        
        ]

def create_multiperiod_thermal_model(
        n_time_points = 3,
        # 24-hr GHI in Phoenix, AZ on June 18th (W/m2)
        GHI = [0, 0, 0, 0, 0, 23, 170, 386, 596, 784, 939, 1031, 1062, 1031, 938, 790, 599, 383, 166, 31, 0, 0, 0, 0],
):

    """
    This function creates a multi-period pv battery flowsheet object. This object contains
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """ 

    mp = MultiPeriodModel(
        n_time_points = n_time_points,
        process_model_func = create_mp_steady_state,
        linking_variable_func = get_variable_pairs,
        initialization_func = fix_dof_and_initialize,
        unfix_dof_func = unfix_dof,
        outlvl = logging.WARNING,
    )


    model_data = {
        t: {
            "GHI" : GHI[t],
            "elec_price" : 0.07 ,
            "mass_fr_fpc" : 0.05,
            "mass_fr_tes_hx_solar" : 0.1,
            "mass_fr_tes_process" : 0.05,
        }
        for t in range(n_time_points)
    }

    flowsheet_options = {
            "dt": 3600,
            "GHI": 0 ,
            "mass_fr_fpc" : 0.05,
            "mass_fr_tes_hx_solar" : 0.1,
            "mass_fr_tes_process" : 0.05, 
            "process_inlet_temp": 60, 
        }


    # create the multiperiod object
    mp.build_multi_period_model(
        model_data_kwargs = model_data,
        flowsheet_options= flowsheet_options,
        initialization_options= None, 
        # {
        #     "GHI" : 0,
        #     "mass_fr_fpc" : 0.05,
        #     "mass_fr_tes_hx_solar" : 0.1,
        #     "mass_fr_tes_process" : 0.05,
        # },
        unfix_dof_options=None,
    )


    active_blocks = mp.get_active_process_blocks()


    # Initialize and unfix dof for each period
    solver = get_solver()
    for blk in active_blocks:
        fix_dof_and_initialize(
            blk=blk,
            dt= 3600,
            GHI= 0 ,
            mass_fr_tes_hx_solar=0.1,
            mass_fr_tes_process= 0.05, 
            process_inlet_temp=40, 
        )
        print('dof before init: ', degrees_of_freedom(blk))
        result = solver.solve(blk)
        unfix_dof(m=blk)

    
    #Initialize the first step
    active_blocks[0].fs.previous_hx_solar_hot_outlet_temperature.fix(41+273.15)
    active_blocks[0].fs.previous_fpc_outlet_temperature.fix(41+273.15)
    active_blocks[0].fs.previous_tes_tank_temp.fix(39+273.15)
    active_blocks[0].fs.previous_hx_solar_cold_outlet_temperature.fix(41+273.15)
    active_blocks[0].fs.previous_process_outlet_temperature.fix(36+273.15)
    active_blocks[0].fs.previous_grid_duty.fix(0)
    active_blocks[0].fs.previous_acc_grid_duty.fix(0)

    return mp



if __name__ == "__main__":

    mp = create_multiperiod_thermal_model()
    solver = get_solver()
    results = solver.solve(mp)

    # assert_optimal_termination(results)

    print('\nStep 1')
    print_results(mp.blocks[0].process)

    print('\nStep 2')
    print_results(mp.blocks[1].process)

    print('\nStep 3')
    print_results(mp.blocks[2].process)

    # print(mp.blocks[0].process.fs.previous_hx_solar_hot_outlet_temperature.value-273.15)