from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.analysis.case_studies.permian.permian_ST2_MD import (
    run_permian_st2_md,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from watertap_contrib.reflo.analysis.case_studies.permian.utils.results_dict import *
from watertap_contrib.reflo.analysis.case_studies.permian.utils.permian_MD_sweep_plot_functions.case_study_plotting import *


def plot_case_study(df,xcol,ax_dict):

    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
        "EC": "fs.treatment.EC.unit.costing",
        "CF": "fs.treatment.cart_filt.unit.costing",
        "MD": "fs.treatment.md.unit",
        "MEC": "fs.treatment.mec.unit.costing",
        "CST": "fs.energy.cst.unit.costing"
    }
    agg_flows = {
        "Electricity":"electric",
        "Heat":"heat",
        "H2O2":"hydrogen_peroxide",
        "Aluminum":"aluminum",
        # "Steam":"steam"
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
    'water_recovery':[0.5], #np.linspace(0.2,0.5,4),
    'heat_price':np.linspace(0.00447,0.011175,4),     # $/kwh
    'grid_frac_heat':np.linspace(0.5,0.9,4),
    }   
    
    input_dict = {
        'Qin': 5, 
        'tds': 130,
        'water_recovery':0.5,
        'grid_frac_heat':1,
        'heat_price':0.00894,
        "electricity_price":0.04346,
    }

    permian_cryst_config = {
        "operating_pressures": [0.45, 0.25, 0.208, 0.095], # Operating pressure of each effect (bar)
        "nacl_yield": 0.9, # Yield
        "heat_transfer_coefficient": 0.13
        }


    #############################################################################################
    # Select sweep type
    #############################################################################################
    
    sweep_type = "water_recovery"
    only_plot = False
    # only_plot = True

    if input_dict['grid_frac_heat'] == 1:
        treatment_only = True
    else:
        treatment_only = False 
    

    xcol_dict = {
        "water_recovery":"fs.water_recovery",
        "heat_price": "fs.costing.heat_cost_buy",
        "hours_storage": "fs.energy.FPC.hours_storage",
        "grid_frac_heat": "fs.costing.frac_heat_from_grid"
    }

    ax_dict = {
        "water_recovery": "MD Water Recovery (%)",
        "heat_price": "Heat Price ($/kWh)",
        "hours_storage": "Hours Storage (h)",
        "grid_frac_heat": "Grid Fraction (Heat)"
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
        m = run_permian_st2_md(
                Qin=input_dict['Qin'], 
                tds=input_dict['tds'], 
                grid_frac_heat=input_dict['grid_frac_heat'],
                water_recovery= input_dict['water_recovery'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                permian_cryst_config=permian_cryst_config,
                treatment_only=treatment_only
                )
        
        results_dict_test = build_results_dict(m, skips=skips)

        
        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = run_permian_st2_md(
                Qin=input_dict['Qin'], 
                tds=input_dict['tds'], 
                grid_frac_heat=input_dict['grid_frac_heat'],
                water_recovery= input_dict['water_recovery'],
                heat_price=input_dict['heat_price'],
                electricity_price=input_dict['electricity_price'],
                permian_cryst_config=permian_cryst_config,
                treatment_only=treatment_only
                )
            
            results_dict_test = results_dict_append(m, results_dict_test)

        df = pd.DataFrame.from_dict(results_dict_test)

        if sweep_type != "grid_frac_heat":
            grid_frac_var = str(input_dict["grid_frac_heat"])
        else:
            grid_frac_var = "var"

        if sweep_type != "water_recovery":
            rec_var = str(input_dict["water_recovery"])
        else:
            rec_var = "var"

        filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/permian-case-study-md/ST2_MD_sweep_results//permian_ST2_MD_"+ sweep_type + "_grid_frac_" + grid_frac_var + "_recovery_" + rec_var + "_check3.csv"
        df.to_csv(filename)
        # df_T= pd.DataFrame.from_dict(results_dict_test, orient='index')
        # df_T.to_csv("/Users/mhardika/Documents/watertap-seto/Mukta-Work//permian-case-study-md/ST1_MD_sweep_results//"+"grid_frac_heat_0.5"+ "_T.csv")

    if sweep_type != "grid_frac_heat":
            grid_frac_var = str(input_dict["grid_frac_heat"])
    else:
        grid_frac_var = "var"

    if sweep_type != "water_recovery":
        rec_var = str(input_dict["water_recovery"])
    else:
        rec_var = "var"

    # filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work//permian-case-study-md/ST2_MD_sweep_results//permian_ST2_MD_"+sweep_type + "_grid_frac_" + grid_frac_var + "_recovery_" + rec_var + ".csv"
    # df = pd.read_csv(filename).drop(columns="Unnamed: 0")
    # plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type])