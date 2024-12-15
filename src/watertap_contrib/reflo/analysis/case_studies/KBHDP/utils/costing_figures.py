from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from psPlotKit.data_plotter.ps_line_plotter import linePlotter
from psPlotKit.data_manager.ps_data_manager import psDataManager, psData
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.figure_device_groups import (
    figure_device_groups,
)

pal = {
    "Baseline": "#1f78b4",
    "FFRRO": "#33a02c",
    # "#b2df8a",
    # "#33a02c",
    # "#fb9a99",
    # "#e31a1c",
    # "#fdbf6f",
    # "#ff7f00",
    # "#cab2d6",
    # "#6a3d9a",
    # "#ffff99",
    # "#b15928",
}

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(os.path.dirname(parent_dir), "sweep_yamls")
sweep_results_dir = os.path.join(os.path.dirname(parent_dir), "sweep_results", "output")

# def create_sweep_cost_breakdown(
#     sweep_file=None, x_data=None, save_name=None, xrange=None
# ):
#     if sweep_file is None:
#         RuntimeError("Please provide a sweep file")
#     # if device_groups is None:
#     #     RuntimeError("Please provide device groups")
#     if save_name is None:
#         RuntimeError("Please provide a save name")

#     filepath = os.path.abspath(__file__)
#     parent_dir = os.path.dirname(filepath)
#     top_level_dir = os.path.dirname(parent_dir)
#     save_path = os.path.join(top_level_dir, "figures")
#     data_path = os.path.join(top_level_dir, "sweep_results/output/")

#     """ import data"""
#     costing_data = psDataManager(
#         os.path.join(
#             data_path,
#             sweep_file,
#         )
#     )

#     """ here we define the costing groups, these should be either units on the block
#     another example is in bgw_ref/costing_plotting_trains_rkt.py that shows different ways to isolate units for cost plotting
#     or explicit keys that point to capex or opex for each units.
#     IF you provide exact keys only they will be used (If you duplicate a key it will be used twice!)"""
#     device_groups = {
#         "Heat": {
#             "OPEX": {
#                 "units": {
#                     "fs.costing.total_heat_operating_cost",
#                     #   "fs.costing.aggregate_flow_electricity_sold"
#                 },
#             },
#         },
#         "Electricity": {
#             "OPEX": {
#                 "units": {
#                     "fs.costing.aggregate_flow_electricity_purchased",
#                     #   "fs.costing.aggregate_flow_electricity_sold"
#                 },
#             },
#         },
#         "Injection": {
#             "OPEX": {
#                 "units": {"fs.treatment.DWI.unit.costing.variable_operating_cost"},
#             },
#         },
#         "Pumps": {
#             "CAPEX": {
#                 "units": {"fs.treatment.pump.costing.capital_cost"},
#             },
#         },
#         "PV": {
#             "CAPEX": {
#                 "units": {"fs.energy.pv.costing.capital_cost"},
#             },
#             "OPEX": {
#                 "units": {
#                     "fs.energy.pv.costing.fixed_operating_cost",
#                 }
#             },
#         },
#         "FPC": {
#             "CAPEX": {
#                 "units": {"fs.energy.FPC.costing.capital_cost"},
#             },
#             "OPEX": {
#                 "units": {
#                     "fs.energy.FPC.costing.fixed_operating_cost",
#                 }
#             },
#         },
#         "RO": {
#             "CAPEX": {
#                 "units": {
#                     "fs.treatment.RO.stage[1].module.costing.capital_cost",
#                 }
#             },
#             "OPEX": {
#                 "units": {
#                     "fs.treatment.RO.stage[1].module.costing.fixed_operating_cost",
#                 }
#             },
#         },
#         "UF": {
#             "CAPEX": {
#                 "units": {
#                     "fs.treatment.UF.unit.costing.capital_cost",
#                 },
#             },
#         },
#         # "EC": {
#         #     "CAPEX": {
#         #         "units": {
#         #             "fs.treatment.EC.ec.costing.capital_cost",
#         #         },
#         #     },
#         #     "OPEX": {
#         #         "units": {
#         #             "fs.treatment.EC.ec.costing.fixed_operating_cost",
#         #             "fs.treatment.costing.aggregate_flow_costs[aluminum]",
#         #         },
#         #     },
#         # },
#         "LTMED": {
#             "CAPEX": {
#                 "units": {
#                     "fs.treatment.LTMED.unit.costing.capital_cost",
#                 },
#             },
#             "OPEX": {
#                 "units": {
#                     "fs.treatment.LTMED.unit.costing.fixed_operating_cost",
#                 },
#             },
#         },
#     }

#     costing_data.load_data(
#         [
#             {
#                 "filekey": "fs.costing.frac_heat_from_grid",
#                 "return_key": "Grid Frac Heat",
#                 "units": "%",
#             },
#             {
#                 "filekey": "fs.costing.heat_cost_buy",
#                 "return_key": "Heat Cost",
#                 "units": "USD/kWh",
#             },
#             {
#                 "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
#                 "return_key": "FPC Cost",
#                 "units": "USD/a/kW",
#             },
#             {
#                 "filekey": "fs.water_recovery",
#                 "return_key": "Water Recovery",
#                 "units": "%",
#             },
#             {
#                 "filekey": "fs.treatment.costing.deep_well_injection.dwi_lcow",
#                 "return_key": "Disposal Cost",
#                 "units": "USD_2021/m**3",
#             },
#         ],
#     )

#     """ define the base costing block and flow (This is used to normalize LCOW)
#     the costing block is used to pull out default values for cost of
#     electricity, costing factors (TIC etc) and deice costs if applicable, so in theory if you have multiple costing
#     blocks they should all share same factors, so it matters not which one you pick
#     The CAPEX and OPEX will be aggregated from keys in supplied costing groups, the
#     capital and opex costs that are attached to the costing packages will not be used!
#     """
#     costing_data.get_costing(
#         device_groups,
#         costing_block="fs.costing",
#         default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
#     )
#     """ note how you have access to total and levelized cost. You can plot either one"""
#     costing_data.display()

#     """lets plot our data"""
#     cost_plotter = breakDownPlotter(
#         costing_data,
#         save_location=save_path,
#         save_name=save_name,
#     )
#     print(cost_plotter)
#     """ define the costing groups, this will be order they are plotted in"""
#     cost_plotter.define_area_groups(
#         [
#             "Heat",
#             "Electricity",
#             "Injection",
#             "FPC",
#             "UF",
#             "Pumps",
#             "RO",
#             "LTMED",
#         ]
#     )
#     """ define if you want to plot specific groups, for example CAPEX, OPEX or TOTAl isstead"""
#     cost_plotter.define_hatch_groups({"CAPEX": {"hatch": ""}, "OPEX": {"hatch": "//"}})
#     cost_plotter.plotbreakdown(
#         xdata=x_data,
#         ydata=[
#             "cost_breakdown",
#             "levelized",
#         ],  # remove "levelized to plot absolute capex/opex
#         axis_options={
#             "yticks": [0, 1, 2, 3, 4, 5, 6],  # adjust as needed
#             # "yticks": [0, 0.25 ,0.5],  # adjust as needed
#             "xticks": costing_data[costing_data.directory_keys[0], x_data].data,
#         },
#         legend_loc="upper right",
#         generate_figure=True,
#     )

#     return costing_data, device_groups, cost_plotter


def create_sweep_cost_breakdown(
    costing_data, device_groups=None, save_name=None, xrange=None
):
    """here we define the costing groups, these should be either units on the block
    another example is in bgw_ref/costing_plotting_trains_rkt.py that shows different ways to isolate units for cost plotting
    or explicit keys that point to capex or opex for each units.
    IF you provide exact keys only they will be used (If you duplicate a key it will be used twice!)
    """
    if device_groups == None:
        device_groups = {
            "Heat": {
                "OPEX": {
                    "units": {
                        "fs.costing.total_heat_operating_cost",
                    },
                },
            },
            "Electricity": {
                "OPEX": {
                    "units": {
                        "fs.costing.aggregate_flow_electricity_purchased",
                    },
                },
            },
        }

    costing_data.load_data(
        [
            {
                "filekey": "fs.costing.frac_heat_from_grid",
                "return_key": "Grid Frac Heat",
                # "units": "%",
            },
            {
                "filekey": "fs.costing.heat_cost_buy",
                "return_key": "Heat Cost",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                "return_key": "FPC Cost",
                "units": "USD/a/kW",
            },
            {
                "filekey": "fs.water_recovery",
                "return_key": "fs.water_recovery",
                "units": "%",
            },
            {
                "filekey": "fs.treatment.costing.LCOW",
                "return_key": "LCOW",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.energy.pv.annual_energy",
                "return_key": "fs.energy.pv.annual_energy",
                # "units": "kWh",
            },
            {
                "filekey": "fs.treatment.costing.deep_well_injection.dwi_lcow",
                "return_key": "fs.treatment.costing.deep_well_injection.dwi_lcow",
                # "units": "kWh",
            },
            {
                "filekey": "fs.treatment.costing.deep_well_injection.dwi_lcow",
                "return_key": "Disposal Cost",
                "units": "USD/m**3",
            },
            {
                "filekey": "fs.costing.frac_heat_from_grid",
                "return_key": "fs.costing.frac_heat_from_grid",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                "return_key": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                # "units": "USD/m**3",
            },
        ],
    )

    x_var = costing_data.psDataImportInstances[0].file_index[
        costing_data.directory_keys[0]
    ]["sweep_params"]["data_keys"][0]

    """ define the base costing block and flow (This is used to normalize LCOW)
    the costing block is used to pull out default values for cost of 
    electricity, costing factors (TIC etc) and deice costs if applicable, so in theory if you have multiple costing 
    blocks they should all share same factors, so it matters not which one you pick 
    The CAPEX and OPEX will be aggregated from keys in supplied costing groups, the 
    capital and opex costs that are attached to the costing packages will not be used!
    """
    try:
        costing_data.get_costing(
            device_groups,
            costing_block="fs.costing",
            default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
        )
    except:
        costing_data.get_costing(
            device_groups,
            costing_block="fs.costing",
            default_flow="fs.product.properties[0.0].flow_vol_phase[Liq]",
        )

    """ note how you have access to total and levelized cost. You can plot either one"""
    costing_data.display()

    """lets plot our data"""
    cost_plotter = breakDownPlotter(
        costing_data,
        save_name=save_name,
    )
    print(cost_plotter)
    """ define the costing groups, this will be order they are plotted in"""
    cost_plotter.define_area_groups(
        [
            "Heat",
            "Electricity",
            "Injection",
            "FPC",
            "UF",
            "Pumps",
            "RO",
            "LTMED",
        ]
    )
    """ define if you want to plot specific groups, for example CAPEX, OPEX or TOTAl isstead"""
    cost_plotter.define_hatch_groups({"CAPEX": {"hatch": ""}, "OPEX": {"hatch": "//"}})

    y_max = (
        np.array(costing_data[costing_data.directory_keys[0], "LCOW"].data)
        .max()
        .round()
    )
    if np.isnan(y_max):
        print("No figure generated\nLCOW returned NaN\nMoving on to next figure")
    else:
        y_axis_lims = np.linspace(0, y_max, 5)
        cost_plotter.plotbreakdown(
            xdata=x_var,
            ydata=[
                "cost_breakdown",
                "levelized",
            ],
            axis_options={
                "yticks": y_axis_lims,
                # "yticks": [0, 0.25 ,0.5],  # adjust as needed
                "xticks": costing_data[costing_data.directory_keys[0], x_var].data,
            },
            legend_loc="upper right",
            generate_figure=True,
        )

    return costing_data, device_groups, cost_plotter


def create_case_figures(case_name=None, sweep_file=None, device_groups=None):
    if sweep_file is None:
        print(sweep_results_dir)
        for root, dirs, files in os.walk(sweep_results_dir):
            for file in files:
                file_id = file.split("_")
                case_id = case_name.split("_")
                if file_id[:3] == case_id[:3]:
                    costing_data = psDataManager(os.path.join(sweep_results_dir, file))
                    create_sweep_cost_breakdown(
                        costing_data, device_groups=device_groups
                    )


def create_all_figures():
    create_case_figures(
        case_name="KBHDP_SOA_1", device_groups=figure_device_groups["KBHDP_SOA_1"]
    )
    # create_case_figures(
    #     case_name="KBHDP_RPT_1", device_groups=figure_device_groups["KBHDP_RPT_1"]
    # )
    # create_case_figures(
    #     case_name="KBHDP_RPT_2", device_groups=figure_device_groups["KBHDP_RPT_2"]
    # )


if __name__ == "__main__":
    create_all_figures()
