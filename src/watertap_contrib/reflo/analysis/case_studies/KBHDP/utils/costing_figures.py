from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from psPlotKit.data_plotter.ps_line_plotter import linePlotter
from psPlotKit.data_manager.ps_data_manager import psDataManager, psData
import numpy as np
import os
import pandas as pd

# import seaborn as sns
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


def create_sweep_cost_breakdown(
    costing_data, device_groups=None, save_name=None, xrange=None
):
    """here we define the costing groups, these should be either units on the block
    another example is in bgw_ref/costing_plotting_trains_rkt.py that shows different ways to isolate units for cost plotting
    or explicit keys that point to capex or opex for each units.
    IF you provide exact keys only they will be used (If you duplicate a key it will be used twice!)
    """
    filepath = os.path.abspath(__file__)
    util_dir = os.path.dirname(filepath)
    parent_dir = os.path.dirname(util_dir)
    save_path = os.path.join(parent_dir, "figures")

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
                "return_key": "fs.costing.heat_cost_buy",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.cost_per_area_collector",
                "return_key": "fs.energy.costing.flat_plate.cost_per_area_collector",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.costing.electricity_cost_buy",
                "return_key": "fs.costing.electricity_cost_buy",
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
                "filekey": "fs.costing.LCOT",
                "return_key": "LCOT",
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
    print(costing_data.directory_keys)
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
    if save_name == "KBHDP_SOA_1":
        costing_data.get_costing(
            device_groups,
            costing_block="fs.treatment.costing",
            default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
        )
    else:
        costing_data.get_costing(
            device_groups,
            costing_block="fs.costing",
            default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
        )

    """ note how you have access to total and levelized cost. You can plot either one"""
    costing_data.display()

    """lets plot our data"""
    cost_plotter = breakDownPlotter(
        costing_data,
        save_location=save_path,
        save_name=save_name + "_" + x_var,
    )
    print(save_name + "_" + x_var)

    """ define the costing groups, this will be order they are plotted in"""
    ## This is why the rest of the device groups aren't showing up
    cost_plotter.define_area_groups(list(device_groups.keys()))
    """ define if you want to plot specific groups, for example CAPEX, OPEX or TOTAl isstead"""
    cost_plotter.define_hatch_groups({"CAPEX": {"hatch": ""}, "OPEX": {"hatch": "//"}})

    # print(list(device_groups.keys()))
    # assert False
    y_max = np.ceil(
        np.array(costing_data[costing_data.directory_keys[0], "LCOT"].data).max()
    )
    print(np.array(costing_data[costing_data.directory_keys[0], "LCOT"].data))
    print(f"Y Max {y_max}")
    # assert False
    if np.isnan(y_max):
        print("No figure generated\nLCOW returned NaN\nMoving on to next figure")
        y_axis_lims = np.linspace(0, 10, 5)
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
                    print(f"\n\nCreating Figures for {file} sweep\n\n")
                    costing_data = psDataManager(os.path.join(sweep_results_dir, file))
                    create_sweep_cost_breakdown(
                        costing_data, device_groups=device_groups, save_name=case_name
                    )


def create_all_figures():
    # create_case_figures(
    #     case_name="KBHDP_SOA_1", device_groups=figure_device_groups["KBHDP_SOA_1"]
    # )
    # create_case_figures(
    #     case_name="KBHDP_RPT_1", device_groups=figure_device_groups["KBHDP_RPT_1"]
    # )
    create_case_figures(
        case_name="KBHDP_RPT_2", device_groups=figure_device_groups["KBHDP_RPT_2"]
    )


if __name__ == "__main__":
    create_all_figures()
