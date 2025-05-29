from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory, value
import pandas as pd
from watertap_contrib.reflo.analysis.case_studies.permian.permian_RPT1_MD import (
    run_permian_st1_md,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
from watertap_contrib.reflo.analysis.case_studies.permian.utils.results_dict import *
from watertap_contrib.reflo.analysis.case_studies.permian.utils.permian_MD_sweep_plot_functions.case_study_plotting import *


def plot_case_study(df, xcol, ax_dict):

    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
        "EC": "fs.treatment.EC.unit.costing",
        "CF": "fs.treatment.cart_filt.unit.costing",
        "MD": "fs.treatment.md.unit",
        "DWI": "fs.treatment.DWI.unit.costing",
        "CST": "fs.energy.cst.unit.costing",
    }
    agg_flows = {
        "Electricity": "electric",
        "Heat": "heat",
        "H2O2": "hydrogen_peroxide",
        "Aluminum": "aluminum",
    }

    ax_dict = dict(xlabel=ax_dict, ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        unit_dict=unit_dict,
        global_costing_blk="fs.treatment.costing",
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
        "water_recovery": np.linspace(0.35, 0.5, 4),
        "heat_price": np.linspace(0.00447, 0.011175, 4),  # $/kwh
        "grid_frac_heat": np.linspace(0.5, 0.9, 5),
        "dwi_lcow": np.linspace(4.2, 10.5, 4),  # $/m3 treated water
        "cst_cost_per_total_aperture_area": np.linspace(186.5, 466.25, 4),
        "cst_cost_per_storage_capital": np.linspace(31, 77.5, 4),
    }

    input_dict = {
        "Qin": 5,
        "tds": 130,
        "water_recovery": 0.5,
        "grid_frac_heat": 0.5,
        "heat_price": 0.00894,
        "electricity_price": 0.04346,
        "dwi_lcow": 8.4,
        "cst_cost_per_total_aperture_area": 373,
        "cst_cost_per_storage_capital": 62,
    }

    #############################################################################################
    # Select sweep type
    #############################################################################################
    
    sweep_type = "grid_frac_heat"

    sweep_type = "cst_cost_per_storage_capital"
    only_plot = False
    only_plot = True

    xcol_dict = {
        "water_recovery": "fs.water_recovery",
        "heat_price": "fs.costing.heat_cost_buy",
        "hours_storage": "fs.energy.FPC.hours_storage",
        "grid_frac_heat": "fs.costing.frac_heat_from_grid",
        "dwi_lcow": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        "cst_cost_per_total_aperture_area": "fs.energy.costing.trough_surrogate.cost_per_total_aperture_area",
        "cst_cost_per_storage_capital": "fs.energy.costing.trough_surrogate.cost_per_storage_capital",
    }

    ax_dict = {
        "water_recovery": "MD Water Recovery (%)",
        "heat_price": "Heat Price ($/kWh)",
        "hours_storage": "Hours Storage (h)",
        "grid_frac_heat": "Grid Fraction (Heat)",
        "dwi_lcow": "DWI LCOW (\$/m$^3$)",
        "cst_cost_per_total_aperture_area": "Cost per Total Aperture Area ($/m2)",
        "cst_cost_per_storage_capital": "Cost per Thermal Storage Capacity ($/kWh)",
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

    if only_plot == False:
        xcol = xcol_dict[sweep_type]
        m = run_permian_st1_md(
            Qin=input_dict["Qin"],
            tds=input_dict["tds"],
            grid_frac_heat=input_dict["grid_frac_heat"],
            water_recovery=input_dict["water_recovery"],
            heat_price=input_dict["heat_price"],
            electricity_price=input_dict["electricity_price"],
            dwi_lcow=input_dict["dwi_lcow"],
            cst_cost_per_total_aperture_area=input_dict[
                "cst_cost_per_total_aperture_area"
            ],
            cst_cost_per_storage_capital=input_dict["cst_cost_per_storage_capital"],
        )

        results_dict_test = build_results_dict(m, skips=skips)

        for i in sweep_dict[sweep_type]:
            input_dict[sweep_type] = i
            print(input_dict)
            m = run_permian_st1_md(
                Qin=input_dict["Qin"],
                tds=input_dict["tds"],
                grid_frac_heat=input_dict["grid_frac_heat"],
                water_recovery=input_dict["water_recovery"],
                heat_price=input_dict["heat_price"],
                electricity_price=input_dict["electricity_price"],
                dwi_lcow=input_dict["dwi_lcow"],
                cst_cost_per_total_aperture_area=input_dict[
                    "cst_cost_per_total_aperture_area"
                ],
                cst_cost_per_storage_capital=input_dict["cst_cost_per_storage_capital"],
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

    #     filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/permian-case-study-md/ST1_MD_sweep_results//permian_RPT1_MD_"+ sweep_type + "_grid_frac_" + grid_frac_var + "_recovery_" + rec_var + ".csv"
    #     df.to_csv(filename)
    #     # df_T= pd.DataFrame.from_dict(results_dict_test, orient='index')
    #     # df_T.to_csv("/Users/mhardika/Documents/watertap-seto/Mukta-Work//permian-case-study-md/ST1_MD_sweep_results//"+"grid_frac_heat_0.5"+ "_T.csv")

    # if sweep_type != "grid_frac_heat":
    #     grid_frac_var = str(input_dict["grid_frac_heat"])
    # else:
    #     grid_frac_var = "var"

    # if sweep_type != "water_recovery":
    #     rec_var = str(input_dict["water_recovery"])
    # else:
    #     rec_var = "var"

    # filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work//permian-case-study-md/ST1_MD_sweep_results//permian_RPT1_MD_" + sweep_type + "_grid_frac_" + grid_frac_var + "_recovery_" + rec_var + ".csv"
    # df = pd.read_csv(filename).drop(columns="Unnamed: 0")
    # plot_case_study(df, xcol=xcol_dict[sweep_type],ax_dict=ax_dict[sweep_type])
