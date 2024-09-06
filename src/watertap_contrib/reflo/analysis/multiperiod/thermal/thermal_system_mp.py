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

from idaes.core.solvers.get_solver import get_solver
import matplotlib.pyplot as plt
import numpy as np


# Can optimize storage tank volume
def unfix_dof(blk, process_outlet_temperature = None):
    blk.fs.previous_hx_solar_hot_outlet_temperature.unfix()
    blk.fs.previous_fpc_outlet_temperature.unfix()
    blk.fs.previous_tes_tank_temp.unfix()
    blk.fs.previous_hx_solar_cold_outlet_temperature.unfix()
    blk.fs.previous_acc_grid_duty.unfix()
    blk.fs.previous_grid_duty.unfix()

    # Checks if the process outlet temperature is fixed and provided. This is the steady-state operation of the water treatment process
    if process_outlet_temperature!= None:
        blk.fs.previous_process_outlet_temperature.unfix()

    return None

def get_variable_pairs_steady_treatment(t1,t2):
    """
    This function returns pairs of variables that need to be connected across two time periods
    This function assumes that water treatment process is operating at steady state and the stream return to the TES has a constant temperature

    Args:
        t1: current time block
        t2: next time block

    Returns:
        List of tuples with variable pairing
    """
    return [

        (t1.fs.hx_solar.hot_side_outlet.temperature[0], t2.fs.previous_hx_solar_hot_outlet_temperature),
        (t1.fs.fpc.outlet.temperature[0], t2.fs.previous_fpc_outlet_temperature),
        (t1.fs.tes.tes_temperature[0], t2.fs.previous_tes_tank_temp),
        (t1.fs.hx_solar.cold_side_outlet.temperature[0], t2.fs.previous_hx_solar_cold_outlet_temperature),
        
        # In steady-state treatment case- the process outlet temperature is fixed
        # (t1.fs.tes.tes_process_outlet.temperature[0], t2.fs.previous_process_outlet_temperature),
        
        (t1.fs.acc_grid_duty, t2.fs.previous_acc_grid_duty),
        (t1.fs.grid_heater.heat_duty[0], t2.fs.previous_grid_duty),
        
        ]

def get_variable_pairs(t1,t2):
    """
    This function returns pairs of variables that need to be connected across two time periods

    Args:
        t1: current time block
        t2: next time block

    Returns:
        List of tuples with variable pairing
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
        GHI = [10, 20, 30, 0, 0, 23, 170, 386, 596, 784, 939, 1031, 
               1062, 1031, 938, 790, 599, 383, 166, 31, 0, 0, 0, 0],
        initial_temperature = {
        'solar_hx_hot_outlet': 49,
        'fpc_outlet': 49,
        'tes_tank': 49,
        'solar_hx_cold_outlet':49,
        'process_outlet': 49
    },
        process_inlet_temperature = 70,
        process_outlet_temperature = None,

        tank_vol = 2,
        fpc_collector_area = 3 
):

    """
    This function creates a multi-period thermal flowsheet object. This object contains
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period fpc-tes flowsheet model
    """ 

    if process_outlet_temperature!= None:

        mp = MultiPeriodModel(
            n_time_points = n_time_points,
            process_model_func = create_mp_steady_state,
            linking_variable_func = get_variable_pairs_steady_treatment,
            initialization_func = fix_dof_and_initialize,
            unfix_dof_func = unfix_dof,
            outlvl = logging.WARNING,
        )
    
    else:
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
            "dt" : 3600,
            "GHI" : GHI[t],
            "mass_fr_fpc" : 0.05,
            "mass_fr_tes_hx_solar" : 0.1,
            "mass_fr_tes_process" : 0.05,
        }
        for t in range(n_time_points)
    }

    flowsheet_options = {
            "mass_fr_fpc" : 0.05,
            "mass_fr_tes_hx_solar" : 0.1,
            "mass_fr_tes_process" : 0.05, 
            "tank_vol" : tank_vol, 
            "fpc_collector_area": fpc_collector_area,
            "process_inlet_temp": process_inlet_temperature, # This the set point the grid heater needs to achieve
            "process_outlet_temp": process_outlet_temperature # This is the fixed temperature exiting the process
        }


    # create the multiperiod object
    # model_data_kwargs-Passes the inputs required for the process_model_func
    # flowsheet_options-Passes the variables required for the initializtion function

    mp.build_multi_period_model(
        model_data_kwargs = model_data,
        flowsheet_options= flowsheet_options,
        initialization_options= None,
        # {
            # "GHI" : 0,
        #     "mass_fr_fpc" : 0.05,
        #     "mass_fr_tes_hx_solar" : 0.1,
        #     "mass_fr_tes_process" : 0.05,
        # },
        unfix_dof_options=None,
    )

    active_blocks = mp.get_active_process_blocks()

    # Initialize the first step
    active_blocks[0].fs.previous_hx_solar_hot_outlet_temperature.fix(initial_temperature['solar_hx_hot_outlet'] + 273.15)
    active_blocks[0].fs.previous_fpc_outlet_temperature.fix(initial_temperature['fpc_outlet'] + 273.15)
    active_blocks[0].fs.previous_tes_tank_temp.fix(initial_temperature['tes_tank'] + 273.15)
    active_blocks[0].fs.previous_hx_solar_cold_outlet_temperature.fix(initial_temperature['solar_hx_cold_outlet'] + 273.15)

    if process_outlet_temperature == None:
        active_blocks[0].fs.previous_process_outlet_temperature.fix(initial_temperature['process_outlet'] + 273.15)

    active_blocks[0].fs.previous_grid_duty.fix(0)
    active_blocks[0].fs.previous_acc_grid_duty.fix(0)

    # Initialize and unfix dof for each period
    solver = get_solver()
    for blk in active_blocks:
        print('\nDegrees of freedom before initialization function \n', degrees_of_freedom(blk))
        fix_dof_and_initialize(
            blk=blk,
            mass_fr_tes_hx_solar=0.1,
            mass_fr_tes_process= 0.05, 
            tank_vol = tank_vol,
            fpc_collector_area = fpc_collector_area,
            process_inlet_temp=process_inlet_temperature, 
            process_outlet_temp=process_outlet_temperature
        )
        result = solver.solve(blk)
        print('\nDegrees of freedom after initialization function \n', degrees_of_freedom(blk))
        unfix_dof(blk)
        print('\nDegrees of freedom after unfix_dof function \n', degrees_of_freedom(blk))

    # Initialize the first step
    active_blocks[0].fs.previous_hx_solar_hot_outlet_temperature.fix(initial_temperature['solar_hx_hot_outlet'] + 273.15)
    active_blocks[0].fs.previous_fpc_outlet_temperature.fix(initial_temperature['fpc_outlet'] + 273.15)
    active_blocks[0].fs.previous_tes_tank_temp.fix(initial_temperature['tes_tank'] + 273.15)
    active_blocks[0].fs.previous_hx_solar_cold_outlet_temperature.fix(initial_temperature['solar_hx_cold_outlet'] + 273.15)

    if process_outlet_temperature == None:
        active_blocks[0].fs.previous_process_outlet_temperature.fix(initial_temperature['process_outlet'] + 273.15)

    active_blocks[0].fs.previous_grid_duty.fix(0)
    active_blocks[0].fs.previous_acc_grid_duty.fix(0)

    return mp


if __name__ == "__main__":

    n_time_points = 72
    GHI = [0, 0, 0, 0, 0, 23, 170, 386, 596, 784, 939, 1031, 
               1062, 1031, 938, 790, 599, 383, 166, 31, 0, 0, 0, 0,]
    
    GHI = np.tile(GHI,3)

    process_inlet_temperature = 50
    process_outlet_temperature = 45

    fpc_collector_area = 3  #m2
    tank_vol = 3 # m3
    
    # Initial temperatures in C
    initial_temperature = {
        'solar_hx_hot_outlet': 49,
        'fpc_outlet': 49,
        'tes_tank': 49,
        'solar_hx_cold_outlet':49,
        'process_outlet': 49
    }


    mp = create_multiperiod_thermal_model(        
        n_time_points = n_time_points,
        # 24-hr GHI in Phoenix, AZ on June 18th (W/m2)
        GHI = GHI,
        initial_temperature = initial_temperature,

        process_inlet_temperature = process_inlet_temperature,
        process_outlet_temperature = process_outlet_temperature,

        tank_vol = tank_vol,
        fpc_collector_area =fpc_collector_area
        
        )
    
    # print('\nDegrees of freedom after initialization', degrees_of_freedom(mp))

    solver = get_solver()
    results = solver.solve(mp)

    assert_optimal_termination(results)

    print('\nStep 1')
    print_results(mp.blocks[0].process)

    print('\nStep 2')
    print_results(mp.blocks[1].process)

    print('\nStep 3')
    print_results(mp.blocks[2].process)


    tes_temp = []
    grid_duty = []

    fpc_inlet = []
    fpc_outlet = []

    hx_hot_inlet = []
    hx_hot_outlet = []

    hx_cold_inlet = []
    hx_cold_outlet = []

    process_outlet = []
    grid_heater_outlet = []

    # fig, ((ax00,ax0,ax1), (ax2, ax3)) = plt.subplots(3,3)
    fig = plt.figure()
    gs = fig.add_gridspec(3,2)
    ax00 = fig.add_subplot(gs[0, :])
    ax0 = fig.add_subplot(gs[1, 0])
    ax1 = fig.add_subplot(gs[1, 1])
    ax2 = fig.add_subplot(gs[2, 0])
    ax3 = fig.add_subplot(gs[2, 1])

    for n in range(0,n_time_points):
        tes_temp.append(mp.blocks[n].process.fs.tes.tes_temperature[0].value-273.15)
        process_outlet.append(mp.blocks[n].process.fs.tes.tes_process_inlet.temperature[0]()-273.15)
        grid_heater_outlet.append(mp.blocks[n].process.fs.grid_heater.outlet.temperature[0]() - 273.15)

        fpc_inlet.append(mp.blocks[n].process.fs.fpc.inlet.temperature[0].value-273.15)
        fpc_outlet.append(mp.blocks[n].process.fs.fpc.outlet.temperature[0]()-273.15)

        hx_hot_inlet.append(mp.blocks[n].process.fs.hx_solar.hot_side_inlet.temperature[0].value - 273.15)
        hx_hot_outlet.append(mp.blocks[n].process.fs.hx_solar.hot_side_outlet.temperature[0].value - 273.15)

        hx_cold_inlet.append(mp.blocks[n].process.fs.hx_solar.cold_side_inlet.temperature[0].value - 273.15)
        hx_cold_outlet.append(mp.blocks[n].process.fs.hx_solar.cold_side_outlet.temperature[0].value - 273.15)

        grid_duty.append(mp.blocks[n].process.fs.grid_heater.heat_duty[0].value)

    ax00.plot(range(0,n_time_points),GHI[0:n_time_points], label = 'GHI', marker = 'o')

    ax0.plot(range(0,n_time_points),fpc_inlet, label = 'FPC inlet', marker = 'o')
    ax0.plot(range(0,n_time_points),fpc_outlet, label = 'FPC outlet', marker = 'o')

    ax1.plot(range(0,n_time_points),hx_hot_inlet, label = 'Solar HX hot inlet', marker = 'o')
    ax1.plot(range(0,n_time_points),hx_hot_outlet, label = 'Solar HX hot outlet', marker = 'o')

    ax1.plot(range(0,n_time_points),hx_cold_inlet, label = 'Solar HX cold inlet', marker = 'o')
    ax1.plot(range(0,n_time_points),hx_cold_outlet, label = 'Solar HX cold outlet', marker = 'o')

    ax2.plot(range(0,n_time_points),tes_temp, label = 'TES temp', marker = 'o')
    ax2.plot(range(0,n_time_points),process_outlet, label = 'Process to TES', marker = 'o')
    ax2.plot(range(0,n_time_points),grid_heater_outlet, label = 'Grid heater outlet', marker = 'o')

    ax3.plot(range(0,n_time_points),grid_duty, marker = 'o')

    ax00.set_title('GHI',fontsize = 20)
    ax0.set_title('FPC',fontsize = 20)
    ax1.set_title('Solar HX',fontsize = 20)
    ax2.set_title('TES',fontsize = 20)
    ax3.set_title('Additional grid duty',fontsize = 20)


    ax00.set_ylabel(r'W/m${^2}$',fontsize = 20)
    ax0.set_ylabel(r'Temperature (${^o}$C)',fontsize = 20)
    ax1.set_ylabel(r'Temperature (${^o}$C)',fontsize = 20)
    ax2.set_ylabel(r'Temperature (${^o}$C)',fontsize = 20)
    ax3.set_ylabel('Heat Duty (W)',fontsize = 20)


    ax0.set_ylim([40,70])
    ax1.set_ylim([40,70])
    ax3.set_ylim([-200,800])

    # Function to find xmin and xmax for the axvspan function
    # Get first non-zero point and first zero point
    xmin_list = []
    xmax_list = []
    for n in range(0,len(GHI)-1):
        if GHI[n]<=0 and GHI[n+1]>0:
            xmin_list.append(n)
        if GHI[n]>0 and GHI[n+1]<=0:
            xmax_list.append(n)

    print(xmin_list)
    print(xmax_list)

    for ax in (ax00,ax0,ax1, ax2, ax3):
        ax.legend()
        ax.set_xlabel('Time (h)',fontsize = 20)
        ax.tick_params(axis='both', which='major', labelsize=20)

        for xmin, xmax in zip(xmin_list,xmax_list):
            ax.axvspan(xmin,xmax, facecolor ='gold',alpha=0.3)

    # gs.tight_layout(fig, pad = 0.1)
    plt.subplots_adjust(left=0.1,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.2, 
                    hspace=0.5)
    plt.show()

    # Calculating the average grid heat duty
    grid_duty = np.array(grid_duty)
    avg_grid_duty = grid_duty[grid_duty>0].mean()
                              
    print('Average grid:',avg_grid_duty)
                              