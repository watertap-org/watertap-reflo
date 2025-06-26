from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.code_demos.REFLO_demo.flowsheets.demo_kbhdp_mld import (
    demo2_mld_grid_only,
    demo2_mld_solar,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from watertap_contrib.reflo.code_demos.REFLO_demo.flowsheets.sweep_functions.case_study_plotting import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.results_dict import *


def plot_case_study(df,xcol,ax_dict, unit_dict, agg_flows):
    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    ax_dict = dict(xlabel=ax_dict, ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        unit_dict=unit_dict,
        costing_blk="fs.costing",
        agg_flows=agg_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
    )
    
    fig1 = plot_elec(df,flow_col,ax_dict)
    fig2 = plot_heat(df,flow_col,ax_dict)

    plt.show()


def sweep_kbhdp_mld_grid_only(sweep_type = "ro_water_recovery", only_plot = False):
    sweep_dict = {
    'md_water_recovery':[0.3,0.4,0.5,0.6,0.7],
    }   
    
    input_dict = {
        'ro_water_recovery':0.8,
        'md_water_recovery':0.5,
        'heat_price':0.00894,
        'electricity_price':0.04989,
        'dwi_lcow':0.58,

    }

    #############################################################################################
    # Select sweep type
    #############################################################################################

    xcol_dict = {
       "md_water_recovery":"md_water_recovery",
    }
    ax_dict = {
        "md_water_recovery": "MD Water Recovery (%)",
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
        m = demo2_mld_grid_only(
                ro_recovery=input_dict["ro_water_recovery"],
                md_water_recovery = input_dict['md_water_recovery'],
                heat_price=input_dict["heat_price"],
                electricity_price=input_dict["electricity_price"],
                dwi_lcow=input_dict["dwi_lcow"],
                )
        
        results_dict_test = build_results_dict(m, skips=skips)

        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = demo2_mld_grid_only(    
                ro_recovery=input_dict["ro_water_recovery"],
                md_water_recovery = input_dict['md_water_recovery'],
                heat_price=input_dict["heat_price"],
                electricity_price=input_dict["electricity_price"],
                dwi_lcow=input_dict["dwi_lcow"],
                )
            
            results_dict_test = results_dict_append(m, results_dict_test)

        df = pd.DataFrame.from_dict(results_dict_test,orient='index')
        filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/kbhdp_mld_grid_only_" + sweep_type + ".csv"
        df.to_csv(filename)
    
    unit_dict = {
        "EC": "fs.treatment.EC.ec.costing",
        "UF": "fs.treatment.UF.unit.costing",
        "Pump": "fs.treatment.pump.costing",
        "RO": "fs.treatment.RO.stage[1].module.costing",
        "MD": "fs.treatment.md.unit",
        "DWI" : "fs.treatment.dwi.unit.costing",
    }
    agg_flows = {
        "Electricity":"electricity",
        "Heat":"heat",
        "Aluminum":"aluminum",
    }

    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/kbhdp_mld_grid_only_" + sweep_type + ".csv"
    df = pd.read_csv(filename, index_col="Unnamed: 0")
    df = df.transpose()
    df['md_water_recovery'] = df["fs.water_recovery"]*100
    plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type],unit_dict=unit_dict,agg_flows=agg_flows)

def sweep_kbhdp_mld_solar(sweep_type = "ro_water_recovery", only_plot = False):
    
    sweep_dict = {
    'md_water_recovery':[0.3,0.4,0.5,0.6,0.7],
    'frac_heat_from_grid': np.linspace(0.3,0.7,5)
    }   
    
    input_dict = {
        'ro_water_recovery':0.8,
        'md_water_recovery':0.5,
        'heat_price':0.00894,
        'electricity_price':0.04989,
        'frac_heat_from_grid':0.5,
        'frac_elec_from_grid':0.5,
        'cost_per_watt_installed':1.6,
        'cost_per_total_aperture_area':297,
        'cost_per_storage_capital':62,
        'dwi_lcow':0.58,
        
    }

    #############################################################################################
    # Select sweep type
    #############################################################################################

    xcol_dict = {
        "md_water_recovery":"md_water_recovery",
        "frac_heat_from_grid":"fs.costing.frac_heat_from_re"
    }
    ax_dict = {
        "md_water_recovery": "MD Water Recovery (%)",
        "frac_heat_from_grid":"Solar Energy (%)"
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

        m = demo2_mld_solar(
                ro_recovery=input_dict['ro_water_recovery'],
                md_water_recovery=input_dict['md_water_recovery'],
                frac_heat_from_grid=input_dict['frac_heat_from_grid'],
                frac_elec_from_grid=input_dict['frac_elec_from_grid'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                cost_per_total_aperture_area=input_dict['cost_per_total_aperture_area'],
                cost_per_storage_capital=input_dict['cost_per_storage_capital'],
                cost_per_watt_installed=input_dict['cost_per_watt_installed'],
                dwi_lcow=input_dict['dwi_lcow'],
                )
        
        results_dict_test = build_results_dict(m, skips=skips)

        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = demo2_mld_solar(    
                ro_recovery=input_dict['ro_water_recovery'],
                md_water_recovery=input_dict['md_water_recovery'],
                frac_heat_from_grid=input_dict['frac_heat_from_grid'],
                frac_elec_from_grid=input_dict['frac_elec_from_grid'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                cost_per_total_aperture_area=input_dict['cost_per_total_aperture_area'],
                cost_per_storage_capital=input_dict['cost_per_storage_capital'],
                cost_per_watt_installed=input_dict['cost_per_watt_installed'],
                dwi_lcow=input_dict['dwi_lcow'],
                )
            
            results_dict_test = results_dict_append(m, results_dict_test)

        df = pd.DataFrame.from_dict(results_dict_test,orient='index')
        filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/kbhdp_mld_solar_" + sweep_type + ".csv"
        df.to_csv(filename)
    
    unit_dict = {
        "EC": "fs.treatment.EC.ec.costing",
        "UF": "fs.treatment.UF.unit.costing",
        "Pump": "fs.treatment.pump.costing",
        "RO": "fs.treatment.RO.stage[1].module.costing",
        "DWI" : "fs.treatment.dwi.unit.costing",
        "MD": "fs.treatment.md.unit",
        "PV": "fs.energy.pv.costing",
        "CST": "fs.energy.cst.unit.costing",
    }
    agg_flows = {
        "Electricity":"electric",
        "Heat":"heat",
        "Aluminum":"aluminum",
    }

    

    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/kbhdp_mld_solar_" + sweep_type + ".csv"
    df = pd.read_csv(filename, index_col="Unnamed: 0")
    df = df.transpose()
    try:
        df["fs.costing.frac_elec_from_re"] = (1-df["fs.costing.frac_elec_from_grid"])*100
    except:
        pass
    try:
        df["fs.costing.frac_heat_from_re"] = (1-df["fs.costing.frac_heat_from_grid"])*100
        df['md_water_recovery'] = df["fs.water_recovery"]*100
    except:
        pass
    plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type],unit_dict=unit_dict,agg_flows=agg_flows)

if __name__ == "__main__":

    # sweep_kbhdp_mld_grid_only(sweep_type = "md_water_recovery", only_plot = True)
    # sweep_kbhdp_mld_solar(sweep_type = "md_water_recovery", only_plot = True)
    sweep_kbhdp_mld_solar(sweep_type = "frac_heat_from_grid", only_plot = True)