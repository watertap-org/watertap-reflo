from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_ZLD_MH import (
    zld_main,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.ZLD_sweep_plot_functions.case_study_plotting import *


def plot_case_study(df,xcol,ax_dict):

    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "EC": "fs.treatment.EC.ec.costing",
        "UF": "fs.treatment.UF.unit.costing",
        "Pump": "fs.treatment.pump.costing",
        "RO": "fs.treatment.RO.stage[1].module.costing",
        "MD": "fs.treatment.md.unit",
        "MEC": "fs.treatment.mec.unit.costing",
        "CST": "fs.energy.cst.unit.costing",
        "PV": "fs.energy.pv.costing",
    }
    agg_flows = {
        "Electricity":"electric",
        "Heat":"heat",
        "Aluminum":"aluminum",
    }

    revenue_flows = {
        "NaCl_recovered":"NaCl_recovered"
    }

    ax_dict = dict(xlabel=ax_dict, ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        unit_dict=unit_dict,
        global_costing_blk = "fs.treatment.costing",
        agg_flows=agg_flows,
        revenue_flows =revenue_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
    )

    plt.show()


if __name__ == "__main__":

    sweep_dict = {
    'ro_water_recovery':[0.6,0.8],
    'nacl_recovery_price':[0,-0.012,-0.024],
    'grid_frac_heat':np.linspace(0.2,0.9,8),
    'cst_cost_per_total_aperture_area':np.linspace(148.5,445.5,4),
    'cst_cost_per_storage_capital':np.linspace(31,93,4),
    'pv_cost_per_watt_installed':np.linspace(0.8,2.4,4),
    }   
    
    input_dict = {
        'ro_water_recovery':0.8,
        'md_water_recovery':0.78,
        'nacl_recovery_price':0,
        'heat_price':0.00894,
        'electricity_price':0.04989,
        'grid_frac_heat': 0.5,
        'cst_cost_per_total_aperture_area': 297,
        'cst_cost_per_storage_capital': 62,
        'pv_cost_per_watt_installed': 1.6
    }


    #############################################################################################
    # Select sweep type
    #############################################################################################
    
    sweep_type = "grid_frac_heat"
    only_plot = False
    # only_plot = True

    xcol_dict = {
        "ro_water_recovery":"fs.treatment.ro_water_recovery",
        "md_water_recovery":"water_recovery",
        "nacl_recovery_price":"fs.treatment.costing.nacl_recovered.cost",
        'cst_cost_per_total_aperture_area':'fs.energy.costing.trough_surrogate.cost_per_total_aperture_area',
        'cst_cost_per_storage_capital':'fs.energy.costing.trough_surrogate.cost_per_storage_capital',
        'pv_cost_per_watt_installed':'fs.energy.costing.pv_surrogate.cost_per_watt_installed',
        'grid_frac_heat':'fs.costing.frac_elec_from_grid'
    }

    ax_dict = {
        "ro_water_recovery": "RO Water Recovery (%)",
        "md_water_recovery": "MD Water Recovery (%)",
        "nacl_recovery_price": "NaCl Recovery Price ($/kg)",
        'cst_cost_per_total_aperture_area':"Cost per Total Aperture Area ($/m2)",
        'cst_cost_per_storage_capital':"Cost per Thermal Storage Capacity ($/kWh)",
        'pv_cost_per_watt_installed':"Cost per Watt Installed ($/W)",
        'grid_frac_heat':'Grid Fraction'
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
        m = zld_main(
                ro_recovery= input_dict['ro_water_recovery'],
                md_water_recovery= input_dict['md_water_recovery'],
                nacl_recovery_price = input_dict['nacl_recovery_price'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                grid_frac_heat = input_dict["grid_frac_heat"],
                cost_per_total_aperture_area=input_dict['cst_cost_per_total_aperture_area'],
                cost_per_storage_capital=input_dict['cst_cost_per_storage_capital'],
                cost_per_watt_installed = input_dict['pv_cost_per_watt_installed'],
                )
        
        results_dict_test = build_results_dict(m, skips=skips)

        
        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = zld_main(
                ro_recovery= input_dict['ro_water_recovery'],
                md_water_recovery = input_dict['md_water_recovery'],
                nacl_recovery_price = input_dict['nacl_recovery_price'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                grid_frac_heat = input_dict["grid_frac_heat"],
                cost_per_total_aperture_area=input_dict['cst_cost_per_total_aperture_area'],
                cost_per_storage_capital=input_dict['cst_cost_per_storage_capital'],
                cost_per_watt_installed = input_dict['pv_cost_per_watt_installed'],
                )
            
            results_dict_test = results_dict_append(m, results_dict_test)


        df = pd.DataFrame.from_dict(results_dict_test)
        filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/ZLD_sweep_results/kbhdp_ZLD_" + sweep_type + ".csv"
        df.to_csv(filename)
        # df_T= pd.DataFrame.from_dict(results_dict_test, orient='index')
        # df_T.to_csv("/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//"+sweep_type+ "_T.csv")


    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/ZLD_sweep_results/kbhdp_ZLD_"+ sweep_type + ".csv" 
    # filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/RPT3_sweep_results//"+ sweep_type + ".csv" 
    df = pd.read_csv(filename).drop(columns="Unnamed: 0")
    plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type])