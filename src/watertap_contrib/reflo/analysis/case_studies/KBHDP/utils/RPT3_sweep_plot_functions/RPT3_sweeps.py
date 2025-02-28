from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_3 import (
    main,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.RPT3_sweep_plot_functions.case_study_plotting import *


def plot_case_study(df,xcol,ax_dict):

    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "MD": "fs.treatment.md.unit",
        "DWI": "fs.treatment.dwi.unit.costing",
        "FPC": "fs.energy.FPC.costing",
    }
    agg_flows = {
        "Electricity":"electric",
        "Heat":"heat",
    }

    ax_dict = dict(xlabel=ax_dict, ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        unit_dict=unit_dict,
        global_costing_blk = "fs.treatment.costing",
        agg_flows=agg_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
    )

    plt.show()


if __name__ == "__main__":

    sweep_dict = {
    'water_recovery':np.linspace(0.4,0.7,4),
    'heat_price': np.linspace(0.0083,0.02075,4),     # $/kwh
    'grid_frac_heat':np.linspace(0.5,0.9,5),
    'dwi_lcow':np.linspace(0.0294,0.073375,4),     # $/m3 treated water
    'cost_per_area_collector':np.linspace(300,750,4),  # $/m2
    'cost_per_volume_storage': np.linspace(1000,2500,4),    # $/m3
    
    }   
    
    input_dict = {
        'water_recovery':0.7,
        'heat_price':0.0166,
        "electricity_price":0.04989,
        'grid_frac_heat':0.5,
        'hours_storage':24,
        'cost_per_area_collector':600,
        'cost_per_volume_storage':2000,
        'dwi_lcow': 0.0587
    }


    #############################################################################################
    # Select sweep type
    #############################################################################################
    
    sweep_type = "heat_price"
    only_plot = False
    # only_plot = True

    xcol_dict = {
        "water_recovery":"fs.water_recovery",
        "heat_price": "fs.costing.heat_cost_buy",
        "hours_storage":"fs.energy.FPC.hours_storage",
        "grid_frac_heat":"fs.costing.frac_heat_from_grid",
        "dwi_lcow": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        "cost_per_area_collector":"fs.energy.costing.flat_plate.cost_per_area_collector",
        "cost_per_volume_storage":"fs.energy.costing.flat_plate.cost_per_volume_storage",
    }

    ax_dict = {
        "water_recovery": "Water Recovery (%)",
        "heat_price": "Heat Price ($/kWh)",
        "hours_storage": "Hours Storage (h)",
        "grid_frac_heat": "Grid Fraction (Heat)",
        "dwi_lcow": "DWI LCOW ($/m3)",
        "cost_per_area_collector":"Collector Cost (\$/m$^2$)",
        "cost_per_volume_storage":"Thermal Storage Cost (\$/m$^3$)",
    }


    skips = [
        "bpe_",
        "dh_vap_w_param",
        "cp_phase_param",
        "pressure_sat_param",
        "enth_mass_param",
        "osm_coeff_param",
        "diffus_param",
        "visc_d_param",
        "diffus_phase",
        "dens_mass_param",
        "therm_cond_phase_param",
        "TIC",
        "TPEC",
        "blocks[",
        "yearly_heat_production",
        "yearly_electricity_production",
        "cp_vap_param",
        "cp_mass_phase",
    ]
 
    if only_plot==False:
        xcol = xcol_dict[sweep_type]
        m = main(
                water_recovery= input_dict['water_recovery'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                grid_frac_heat=input_dict['grid_frac_heat'],
                hours_storage=input_dict['hours_storage'],
                cost_per_area_collector= input_dict['cost_per_area_collector'],
                cost_per_volume_storage= input_dict['cost_per_volume_storage'],
                dwi_lcow=input_dict['dwi_lcow']
                )
        
        results_dict_test = build_results_dict(m, skips=skips)

        
        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = main(
                water_recovery= input_dict['water_recovery'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                grid_frac_heat=input_dict['grid_frac_heat'],
                hours_storage=input_dict['hours_storage'],
                cost_per_area_collector= input_dict['cost_per_area_collector'],
                cost_per_volume_storage= input_dict['cost_per_volume_storage'],
                dwi_lcow=input_dict['dwi_lcow']
                )
            
            results_dict_test = results_dict_append(m, results_dict_test)

        if sweep_type!="water_recovery":
            rec = str(input_dict['water_recovery'])
        else:
            rec = 'var'

        if sweep_type!="grid_frac_heat":
            grid_frac = str(input_dict['grid_frac_heat'])
        else:
            grid_frac = 'var'

        df = pd.DataFrame.from_dict(results_dict_test)
        filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//" + sweep_type + "_grid_frac_"+grid_frac+"_recovery_"+rec+".csv"
        df.to_csv(filename)
        # df_T= pd.DataFrame.from_dict(results_dict_test, orient='index')
        # df_T.to_csv("/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//"+sweep_type+ "_T.csv")

    if sweep_type!="water_recovery":
        rec = str(input_dict['water_recovery'])
    else:
        rec = 'var'

    if sweep_type!="grid_frac_heat":
        grid_frac = str(input_dict['grid_frac_heat'])
    else:
        grid_frac = 'var'

    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//"+ sweep_type + "_grid_frac_" + grid_frac + "_recovery_" + rec + ".csv" 
    # filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//"+ sweep_type + ".csv" 
    df = pd.read_csv(filename).drop(columns="Unnamed: 0")
    plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type])