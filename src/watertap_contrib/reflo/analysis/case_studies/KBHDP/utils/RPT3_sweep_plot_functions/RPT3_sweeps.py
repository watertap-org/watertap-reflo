from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_3 import (
    main,
)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import matplotlib.patches as mpatches


def sweep(
    folder_path,
    folder_name,
    sweep_type,
    water_recovery=0.5,
    heat_price=0.07,
    electricity_price=0.07,
    frac_heat_from_grid=0.01,
    hours_storage=6,
    cost_per_area_collector=700,
    cost_per_volume_storage=1000,
    dwi_lcow = None
):

    m = main(
        water_recovery=water_recovery,
        heat_price=heat_price,
        electricity_price=electricity_price,
        grid_frac_heat=frac_heat_from_grid,
        hours_storage=hours_storage,
        cost_per_area_collector=cost_per_area_collector,
        cost_per_volume_storage=cost_per_volume_storage,
        dwi_lcow = dwi_lcow
    )
    save_kbhdp_rpt3_results(m, folder_path, folder_name=folder_name, sweep_type=sweep_type)


def save_kbhdp_rpt3_results(m, folder_path, folder_name=None, sweep_type=None):

    results_df = pd.DataFrame(
        columns=[
            "water_recovery",
            "heat_price",
            "LCOH",
            "hours_storage",
            "fpc_heat_load",
            "cost_per_area_collector",
            "cost_per_volume_storage",
            "frac_heat_from_grid",
            "product_annual_production",
            "utilization_factor",
            "capital_recovery_factor",
            "unit",
            "cost_component",
            "cost",
            "norm_cost_component",
        ]
    )

    capex_output = {
        "FPC": value(m.fs.energy.FPC.costing.capital_cost)
        * value(m.fs.costing.capital_recovery_factor),
        "MD": value(
            m.fs.treatment.md.unit.mp.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.capital_cost
        )
        * value(m.fs.costing.capital_recovery_factor),
        "DWI": 0,
        "Heat": 0,
        "Electricity": 0,
    }

    fixed_opex_output = {
        "FPC": value(m.fs.energy.FPC.costing.fixed_operating_cost)
        + value(m.fs.energy.FPC.costing.capital_cost)
        * value(m.fs.energy.costing.maintenance_labor_chemical_factor),
        "MD": value(
            m.fs.treatment.md.unit.mp.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.fixed_operating_cost
        )
        + value(
            m.fs.treatment.md.unit.mp.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.capital_cost
        )
        * value(m.fs.treatment.costing.maintenance_labor_chemical_factor),
        "DWI": 0,
        "Heat": 0,
        "Electricity": 0,
    }
    variable_opex_output = {
        "FPC": 0,
        "MD": 0,
        "DWI": value(m.fs.treatment.dwi.unit.costing.variable_operating_cost),
        "Heat": value(m.fs.costing.total_heat_operating_cost),
        "Electricity": value(m.fs.costing.total_electric_operating_cost),
    }

    for unit in ["FPC", "MD", "DWI", "Heat", "Electricity"]:
        # Add fixed_opex
        temp = {
            "water_recovery": value(m.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.FPC.hours_storage),
            "fpc_heat_load": value(m.fs.energy.FPC.heat_load),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "cost_per_area_collector":value(m.fs.energy.costing.flat_plate.cost_per_area_collector),
            "cost_per_volume_storage": value(m.fs.energy.costing.flat_plate.cost_per_volume_storage),
            "product_annual_production": value(
                m.fs.costing.annual_water_production
            ),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "fixed_opex",
            "cost": fixed_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

        # Add variable opex
        temp = {
            "water_recovery": value(m.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.FPC.hours_storage),
            "fpc_heat_load": value(m.fs.energy.FPC.heat_load),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "cost_per_area_collector":value(m.fs.energy.costing.flat_plate.cost_per_area_collector),
            "cost_per_volume_storage": value(m.fs.energy.costing.flat_plate.cost_per_volume_storage),
            "product_annual_production": value(
                m.fs.costing.annual_water_production
            ),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "variable_opex",
            "cost": variable_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

        # Add opex
        temp = {
            "water_recovery": value(m.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.FPC.hours_storage),
            "fpc_heat_load": value(m.fs.energy.FPC.heat_load),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "cost_per_area_collector":value(m.fs.energy.costing.flat_plate.cost_per_area_collector),
            "cost_per_volume_storage": value(m.fs.energy.costing.flat_plate.cost_per_volume_storage),
            "product_annual_production": value(
                m.fs.costing.annual_water_production
            ),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "opex",
            "cost": variable_opex_output[unit] + fixed_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

        # Add capex
        temp = {
            "water_recovery": value(m.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.FPC.hours_storage),
            "fpc_heat_load": value(m.fs.energy.FPC.heat_load),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "cost_per_area_collector":value(m.fs.energy.costing.flat_plate.cost_per_area_collector),
            "cost_per_volume_storage": value(m.fs.energy.costing.flat_plate.cost_per_volume_storage),
            "product_annual_production": value(
                m.fs.costing.annual_water_production
            ),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "capex",
            "cost": capex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

    results_df["norm_cost_component"] = (
        results_df["cost"]
        / results_df["product_annual_production"]
        / results_df["utilization_factor"]
    )

    if folder_name != None:

        file_name = "RPT3_" + sweep_type + "_" + str(results_df[sweep_type][0])
        results_df.to_csv(
            folder_path
            + folder_name
            + "\\"
            + file_name
            + ".csv"
        )


def plot_results(folder_name, sweep_type, xaxis_label=None):
    # results_folder_path = r'C:\Users\mhardika\Documents\SETO\Case Studies\RPT3\RPT3_results\\' + folder_name
    # Read output files
    results_dict = {}
    for filename in os.listdir(results_folder_path):
            filepath = os.path.join(results_folder_path, filename)
            file = os.path.splitext(filepath)[0]
            sweep_value = file.split('_')[-1]
            # print(sweep_value)
            results_dict[sweep_value] = filepath

    results = []
    for key in results_dict.keys():
            temp = pd.read_csv(results_dict[key],index_col=[0])
            try:
                results = pd.concat([results,temp])
            except TypeError:
                results = temp

    # Plot the results
    pivot_results = results.pivot(index=['unit','cost_component'],columns=sweep_type,values='norm_cost_component').reset_index()

    # Create a list of y
    hatch_list = []  # Hatched if opex
    label_list = []
    color_list = []
    
    sweep_var = sorted(results[sweep_type].unique())
    col = sns.color_palette("Paired", len(pivot_results['unit'].unique()))   # color for every unit
    y = np.array([]).reshape(0,len(sweep_var))

    i=0
    for unit in ['MD','DWI','FPC','Heat','Electricity']:
        print(unit)
        row = pivot_results[pivot_results['unit']==unit]
        # Add the capex
        y = np.vstack((y,row[row['cost_component']=='capex'][sweep_var]))
        # Add the opex
        y = np.vstack((y,row[row['cost_component']=='opex'][sweep_var]))
        label_list.append(unit)
        label_list.append(unit)
        hatch_list.append('')
        hatch_list.append('//')
        if unit in ['Heat','Electricity']:
             hatch_list[-1] = ''
        color_list.append(col[i])
        color_list.append(col[i])
        i=i+1

    print(hatch_list)
    # Plot results
    sns.set_theme('talk')
    sns.set_style("ticks")

    fig, ax = plt.subplots(figsize=(6,6))
    a = ax.stackplot(sweep_var, y, labels=label_list, colors=color_list, linestyle ='-',edgecolor='black')
    
    ax.set_ylim(0,10)
    for stack, hatch in zip(a, hatch_list):
        stack.set_hatch(hatch)

    ax2=ax.twinx()
    a2 = ax2.scatter(results[sweep_type],results['fpc_heat_load'],c='black',label='FPC Heat Load')
    ax2.set_ylim(round(results['fpc_heat_load'].min()/10)*10-10,round(results['fpc_heat_load'].max()/10)*10+20)

    capex_handle = mpatches.Patch( facecolor='white', edgecolor='black',hatch='',label='Capex')
    opex_handle = mpatches.Patch( facecolor='white', edgecolor='black',hatch='//',label='Opex')
    
    a2_handle, a2_label = ax2.get_legend_handles_labels()
    
    handles = [capex_handle,opex_handle] +  a[0::2] + a2_handle
    labels = ['Capex','Opex'] + label_list[0::2] + a2_label
    # ax.legend(a[1::2],label_list[1::2],ncols=2)
    ax.legend(handles,labels,ncols=2,frameon=False)
    ax.set_xlabel(xaxis_label)
    ax.set_ylabel(r'LCOW (\$/m${^3}$)')
    ax2.set_ylabel(r'FPC Heat Load (MW)')
    plt.show()

if __name__ == "__main__":
    # Recovery
    water_recovery_sweep = [0.5, 0.6, 0.7, 0.8, 0.9]
    # $/kwh
    heat_price_sweep = [0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04]

    hours_storage_sweep = [4, 8, 12, 16, 20, 24]

    grid_frac_sweep = [0, 0.2, 0.4, 0.6, 0.8]

    cost_per_area_collector_sweep = [300, 400, 500, 600, 700]

    cost_per_volume_storage_sweep = [1000,1500,2000,2500]
    
    dwi_lcow_sweep = [0, 0.02, 0.04, 0.06, 0.08, 0.1]


    #############################################################################################
    # Select sweep type
    #############################################################################################
    folder_name = "cost_per_volume_storage_rec_0.9_grid_frac_0.5_lcot"
    sweep_type = "cost_per_volume_storage"

    folder_path = r"C:\Users\mhardika\Documents\SETO\Case Studies\RPT3\RPT3_results\\"

    # Create folder if it doesn't exist
    results_folder_path = (
        folder_path
        + folder_name
    )
    os.makedirs(results_folder_path, exist_ok=True)

    if sweep_type == "water_recovery":
        for recovery in water_recovery_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=recovery,
                heat_price=0.01,
                frac_heat_from_grid=0.99,
                hours_storage=24,
            )

        plot_results(
            folder_name=folder_name, sweep_type=sweep_type, xaxis_label="Water Recovery"
        )

    elif sweep_type == "heat_price":
        for heat_price in heat_price_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.5,
                heat_price=heat_price,
                frac_heat_from_grid=0.5,
                hours_storage=6,
            )

        plot_results(
            folder_name=folder_name,
            sweep_type=sweep_type,
            xaxis_label="Heat Price ($/kWh)",
        )

    elif sweep_type == "hours_storage":
        for h in hours_storage_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.9,
                heat_price=0.01,
                frac_heat_from_grid=0.5,
                hours_storage=h,
            )

        plot_results(
            folder_name=folder_name,
            sweep_type=sweep_type,
            xaxis_label="Hours Storage (h)",
        )

    elif sweep_type == "frac_heat_from_grid":
        for grid_frac in grid_frac_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.9,
                heat_price=0.08,
                frac_heat_from_grid=grid_frac,
                hours_storage=24,
                run_optimization=False,
            )

        plot_results(
            folder_name=folder_name, sweep_type=sweep_type, xaxis_label="Grid Fraction"
        )
    
    elif sweep_type == "cost_per_area_collector":
        for cost_per_area_collector in cost_per_area_collector_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.9,
                heat_price=0.01,
                frac_heat_from_grid=0.5,
                hours_storage=24,
                cost_per_area_collector=cost_per_area_collector,
                cost_per_volume_storage=2000
                
            )

        plot_results(
            folder_name=folder_name, sweep_type=sweep_type, xaxis_label="Cost per Area Collector"
        )
    
    elif sweep_type == "cost_per_volume_storage":
        for cost_per_volume_storage in cost_per_volume_storage_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.9,
                heat_price=0.01,
                frac_heat_from_grid=0.5,
                hours_storage=24,
                cost_per_area_collector=600,
                cost_per_volume_storage=cost_per_volume_storage
                
            )

        plot_results(
            folder_name=folder_name, sweep_type=sweep_type, xaxis_label="Cost per Storage Volume"
        )
    
    elif sweep_type == "dwi_lcow":
        for dwi_lcow in dwi_lcow_sweep:
            sweep(
                folder_path = folder_path,
                folder_name=folder_name,
                sweep_type=sweep_type,
                water_recovery=0.9,
                heat_price=0.01,
                frac_heat_from_grid=0.5,
                hours_storage=24,
                cost_per_area_collector=600,
                cost_per_volume_storage=2000,
                dwi_lcow = dwi_lcow
                
            )

        plot_results(
            folder_name=folder_name, sweep_type=sweep_type, xaxis_label="DWI LCOW"
        )
    