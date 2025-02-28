from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from psPlotKit.data_plotter.ps_line_plotter import linePlotter
from psPlotKit.data_plotter import fig_generator
from psPlotKit.data_manager.ps_data_manager import psDataManager, psData
import numpy as np
import os
import pandas as pd

# import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, CenteredNorm, TwoSlopeNorm
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

def import_data(data_manager):
    import_keys = [
            {
                "filekey": "fs.costing.frac_heat_from_grid",
                "return_key": "Grid Frac Heat",
                # "units": "%",
            },
            {
                "filekey": "fs.costing.frac_elec_from_grid",
                "return_key": "Grid Electricity Fraction",
                # "units": "%",
            },
            {
                "filekey": "fs.costing.total_heat_operating_cost",
                "return_key": "fs.costing.total_heat_operating_cost",
                # "units": "%",
            },
            {
                "filekey": "heat_cost_buy",
                "return_key": "Heat Cost",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.cost_per_area_collector",
                "return_key": "FPC Cost Per Area",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.energy.costing.pv_surrogate.cost_per_watt_installed",
                "return_key": "PV Cost Per Watt Installed",
                # "units": "USD/kWh",
            },
            {
                "filekey": "fs.costing.electricity_cost_buy",
                "return_key": "Electricity Cost",
                "units": "USD/kWh",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                "return_key": "FPC Cost",
                "units": "USD/a/kW",
            },
            {
                "filekey": "fs.water_recovery",
                "return_key": "Water Recovery",
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
                "return_key": "Brine Injection Cost",
                "units": "USD/m**3",
            },
            {
                "filekey": "fs.costing.frac_heat_from_grid",
                "return_key": "Grid Heat Fraction",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                "return_key": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.electrocoagulation.sludge_handling_cost[kbhdp]",
                "return_key": "EC Sludge Disposal Cost",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.aluminum_cost",
                "return_key": "Aluminum Cost",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.lime.cost",
                "return_key": "Lime Cost",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.soda_ash.cost",
                "return_key": "Soda Ash Cost",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.co2.cost",
                "return_key": "CO2 Cost",
                # "units": "USD/m**3",
            },
            {
                "filekey": "fs.treatment.costing.mgcl2.cost",
                "return_key": "MgCl2 Cost",
                # "units": "USD/m**3",
            }
        ]

    data_manager.load_data(
            import_keys,
            exact_keys=False,
        )

    data_manager.display()
    return data_manager


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

    import_data(costing_data)
    print(costing_data.directory_keys)
    print(costing_data.data_keys)

    # assert False
    x_var = costing_data.directory_keys[0].split('/')[-1]

    # assert False

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
        save_name=save_name + "_" + x_var + '_final',
    )
    print(save_name + "_" + x_var + '_final')
    print(x_var)
    """ define the costing groups, this will be order they are plotted in"""
    ## This is why the rest of the device groups aren't showing up
    cost_plotter.define_area_groups(list(device_groups.keys()))
    """ define if you want to plot specific groups, for example CAPEX, OPEX or TOTAl isstead"""
    cost_plotter.define_hatch_groups({"CAPEX": {"hatch": ""}, "OPEX": {"hatch": "//"}})

    print(cost_plotter.area_groups)
    print(cost_plotter.hatch_groups)

    cost_plotter._select_data(x_var, "levelized")
    cost_plotter.selected_data = cost_plotter.psData.get_selected_data()
    for group, items in cost_plotter.hatch_groups.items():
        print(group, items)
    # assert False
    y_max = np.ceil(
        np.array(costing_data[costing_data.directory_keys[0], "LCOT"].data).max()
    )
    xmin = np.array(costing_data[costing_data.directory_keys[0], x_var].data).min()
    xmax = np.array(costing_data[costing_data.directory_keys[0], x_var].data).max()
    print(np.array(costing_data[costing_data.directory_keys[0], "LCOT"].data))
    print(f"Y Max {y_max}")
    # assert False
    if np.isnan(y_max):
        data = np.array(costing_data[costing_data.directory_keys[0], "LCOT"].data)
        data = data[~np.isnan(data)]
        y_max = np.ceil(data.max())
        y_axis_lims = np.linspace(0, y_max, 5)
        cost_plotter.plotbreakdown(
            xdata=x_var,
            ydata=[
                # "cost_breakdown",
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
        x_axis_lims = np.linspace(xmin, xmax, 5)
        cost_plotter.plotbreakdown(
            xdata=x_var,
            ydata="levelized",
            axis_options={
                "yticks": y_axis_lims,
                # "yticks": [0, 0.25 ,0.5],  # adjust as needed
                "xticks": x_axis_lims,
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
                print('\n')
                print(file_id, case_id, file)
                
                if file_id[:3] == case_id[:3]:
                    print('    ',file_id[:3], case_id[:3])
                    print(f"\n\nCreating Figures for {file} sweep\n\n")
                    if file_id[-1] == 'map.h5':
                        costing_data = psDataManager(os.path.join(sweep_results_dir, file))
                        create_map_figure(costing_data,
                                          x_data="Brine Injection Cost",
                                          y_data="Electricity Cost",
                                          z_data="Grid Electricity Fraction",
                                        #   xticks=[50,60,70,80],
                                        #   yticks=[0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0],
                                        #   zticks=[0.9,1,1.1,1.2, 1.3, 1.4, 1.5],
                                          save_name = case_name ,
                                          show=True)
                    else:
                        costing_data = psDataManager(os.path.join(sweep_results_dir, file))
                        print(device_groups, case_name)
                        create_sweep_cost_breakdown(
                            costing_data, device_groups=device_groups, save_name=case_name
                        )
                #         # pass

# def create_map_figures(file, file_id = None, case_id = None, device_groups=None, save_name=None):
#     costing_data = psDataManager(os.path.join(sweep_results_dir, file))

def create_map_figure(
    costing_data,
    x_data=None,
    y_data=None,
    z_data=None,
    x_units="",
    y_units="",
    save=True,
    show=False,
    save_name=None,
    xticks=None,
    yticks=None,
    zticks=None,
    zlabel="LCOW ($\$$/m$^3$)",
    norm=None,
    digitize_levels=None,
    digitize_colors=None,
):
    filepath = os.path.abspath(__file__)
    util_dir = os.path.dirname(filepath)
    parent_dir = os.path.dirname(util_dir)
    save_path = os.path.join(parent_dir, "figures/")
    save_name = save_name + z_data +'_vs_' + x_data +'_vs_' + y_data
    fig = fig_generator.figureGenerator()
    fig.init_figure()
    zformat_prec = 2
    text = True

    # Define the colors: center (white), and edges (#065a82, #942911)
    colors = ['#005AA0', '#DBDBDB', '#d62728']
    diverging_color_map = LinearSegmentedColormap.from_list('custom_diverging', colors)

    if digitize_levels is not None:
        text = False
        digitize_colors = LinearSegmentedColormap.from_list('custom_diverging', colors, N=len(digitize_levels))
    
    import_data(costing_data)
    num_sweep_vars = (len(costing_data.psDataImportInstances[0].file_index[costing_data.directory_keys[0]]['sweep_params']))

    if (x_data is not None) & (y_data is not None) & (z_data is not None):
        xdata = costing_data[costing_data.directory_keys[0], x_data].data
        ydata = costing_data[costing_data.directory_keys[0], y_data].data
        zdata = costing_data[costing_data.directory_keys[0], z_data].data

        if norm is not None:
            if norm == 'max':
                zdata = (zdata - np.nanmax(zdata)) / np.nanmax(zdata) * 100
                zlabel = "% Change in " + zlabel
                save_name = save_name + "_norm"
                zformat_prec = 0
            elif norm == 'min':
                zdata = (zdata - np.nanmin(zdata)) / np.nanmin(zdata) * 100
                zlabel = "% Change in " + zlabel
                save_name = save_name + "_norm"
                zformat_prec = 0
            # elif type(norm) == np.ndarray:
            #     zdata = (zdata - norm) / norm * 100
            else:
                # print(f"Norm Max: {np.nanmax(zdata)}")
                # print(f"Norm Min: {np.nanmin(zdata)}")
                zdata = ((zdata - norm) / norm) * 100
                zlabel = "% Change in " + zlabel
                save_name = save_name + "_norm"
                zformat_prec = 0
        if zticks is None:
            zticks = np.linspace(np.nanmin(zdata), np.nanmax(zdata), 5)
    #     #TODO Add new normalization dynamic range

        # divnorm = TwoSlopeNorm(vmin=min(zticks), vcenter=0, vmax=max(zticks))

        # assert False
        fig.plot_map(
            xdata=xdata,
            ydata=ydata,
            zdata=zdata,
            vmin=min(zticks),
            vmax=max(zticks),
            digitize_levels=digitize_levels,
            digitize_colors=digitize_colors,
            text_color="auto",
            text=text,
            # cmap="coolwarm",
            cmap=diverging_color_map,
            # # plot_contour_lines=[0],
            # plot_contour= False
        )

    #     # zdata = np.nan_to_num(zdata)

        # print(f"Colorar Range {z_range}")
        # cbar = plt.colorbar(fig.colorFig, cmap=diverging_color_map, norm=divnorm)
        # cbar.ax.set_yscale("linear")
        # fig.add_colorbar(fig.colorFig)
        # print(zticks)
        # fig.add_colorbar(zlabel, zticks=zticks)
        fig.add_colorbar(zlabel, zticks=zticks, zformat=zformat_prec, labelpad=17)
        # fig.cbar.ax.set_yscale("linear")
    # else:
    #     print("No data provided")

    if xticks is None:
        xticks = np.linspace(min(np.unique(xdata)),max(np.unique(xdata)), int(len(np.unique(xdata))/2)).round(2)
    if yticks is None:
        yticks = np.linspace(min(np.unique(ydata)),max(np.unique(ydata)), len(np.unique(ydata))).round(2)
    # yticks = np.unique(ydata)

    fig.set_axis_ticklabels(
        xlabel=x_data + x_units,
        ylabel=y_data + y_units,
        xticklabels=xticks,
        yticklabels=yticks,
        angle=0,
        ha="center",
        va="center",
    )

    print(save_name)
    if save == True:
        if save_name is None:
            print(f'Saving Figure as {os.path.join(save_path, "temp")}')
            fig.save(save_location=os.path.join(save_path, "temp"))
        else:
            print(f'Saving Figure as {os.path.join(save_path, save_name)}')
            fig.save(save_location=save_path, file_name=save_name)

    if show == True:
        fig.show()



def create_all_figures():
    # create_case_figures(
    #     case_name="KBHDP_SOA_1", device_groups=figure_device_groups["KBHDP_SOA_1"]
    # )
    create_case_figures(
        case_name="KBHDP_RPT_1", device_groups=figure_device_groups["KBHDP_RPT_1"]
    )
    # create_case_figures(
    #     case_name="KBHDP_RPT_2", device_groups=figure_device_groups["KBHDP_RPT_2"]
    # )
    # create_case_figures(
    #     case_name="KBHDP_RPT_3", device_groups=figure_device_groups["KBHDP_RPT_3"]
    # )


if __name__ == "__main__":
    create_all_figures()