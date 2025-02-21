from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from psPlotKit.data_plotter.ps_line_plotter import linePlotter
from psPlotKit.data_manager.ps_data_manager import psDataManager, psData
import numpy as np
import os
import pandas as pd

# import seaborn as sns
import matplotlib.pyplot as plt
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.plotting_data_inputs import (
    figure_device_groups,
    costing_data_keys,
    default_device_groups,
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
fig_dir = os.path.join(os.path.dirname(parent_dir), "figures")
sweep_yaml_dir = os.path.join(os.path.dirname(parent_dir), "sweep_yamls")
sweep_results_dir = os.path.join(os.path.dirname(parent_dir), "sweep_results", "output")


def create_sweep_cost_breakdown(
    costing_data,
    device_groups=None,
    save_name=None,
    xrange=None,
    y_var="LCOT", 
    fig_options=dict(),
    generate_fig=False,
):
    global x_var, fig_save_file
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
        device_groups = default_device_groups

    costing_data.load_data(costing_data_keys)

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

    costing_data.get_costing(
        device_groups,
        costing_block="fs.treatment.costing",
        default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    )

    """ note how you have access to total and levelized cost. You can plot either one"""
    costing_data.display()
    # assert False

    """lets plot our data"""
    fig_save_file = save_name + "_" + x_var
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
    # assert False
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
            generate_figure=generate_fig,
            fig_options=fig_options,
        )
    else:
        y_axis_lims = np.linspace(0, y_max, 5)
        cost_plotter.plotbreakdown(
            xdata=x_var,
            ydata=[
                "cost_breakdown",
                "levelized",
            ],
            # axis_options={
            #     "yticks": y_axis_lims,
                # "yticks": [0, 0.25 ,0.5],  # adjust as needed
                # "xticks": costing_data[costing_data.directory_keys[0], x_var].data,
            # },
            # legend_loc="upper right",
            generate_figure=generate_fig,
            fig_options=fig_options,
        )

    return costing_data, device_groups, cost_plotter


def create_case_figures(case_name=None, sweep_file=None, device_groups=None, **kwargs):
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
                        costing_data,
                        device_groups=device_groups,
                        save_name=case_name,
                        **kwargs,
                    )
    else:
        costing_data = psDataManager(sweep_file)
        create_sweep_cost_breakdown(
            costing_data, device_groups=device_groups, save_name=case_name, **kwargs
        )


def create_all_figures(**kwargs):
    # create_case_figures(
    #     case_name="KBHDP_SOA_1", device_groups=figure_device_groups["KBHDP_SOA_1"]
    # )
    # create_case_figures(
    #     case_name="KBHDP_RPT_1", device_groups=figure_device_groups["KBHDP_RPT_1"]
    # )
    # create_case_figures(
    #     case_name="KBHDP_RPT_2", device_groups=figure_device_groups["KBHDP_RPT_2"], **kwargs
    # )
    create_case_figures(
        case_name="KBHDP_RPT_3", device_groups=figure_device_groups["KBHDP_RPT_3"]
    )

    pass


def create_fig_legend(cost_plotter, leg1_kwargs=dict(), leg2_kwargs=dict()):

    fig = cost_plotter.fig
    ax = fig.ax[0]
    fig = fig.fig
    leg = ax.legend()
    # add_patch(leg)
    leg.remove()
    handles, labels = ax.get_legend_handles_labels()
    # leg1 is OPEX/CAPEX
    leg1 = ax.legend(handles[2:], labels[2:], **leg1_kwargs)
    ax.add_artist(leg1)
    # leg2 is units
    leg2 = ax.legend(handles[:2], labels[:2], **leg2_kwargs)
    ax.add_artist(leg2)

    return fig, ax


def create_axis_labels_etc(fig, ax, label_kwargs=dict()):
    ax.set(**label_kwargs)
    return fig, ax
    # costing_data.load_data(costing_data_keys)
    # x_var = costing_data.psDataImportInstances[0].file_index[
    #     costing_data.directory_keys[0]
    # ]["sweep_params"]["data_keys"][0]
    # # costing_data[costing_data.directory_keys[0], x_var].data
    # costing_data.get_costing(
    #     figure_device_groups[case_name],
    #     costing_block="fs.treatment.costing",
    #     default_flow="fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    # )
    # # costing_data.display()
    # print(x_var)
    # print(costing_data[costing_data.directory_keys[0], x_var])
    # # costing_data[costing_data.directory_keys[0], x_var]
    # assert False

if __name__ == "__main__":
    fig_options = dict(width=5, height=5)

    f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_RPT_3_20241220-085750_analysisType_KBHDP_RPT_3_heat_price_sweep.h5"

    # create_all_figures(fig_options=fig_options)
    case_name = "KBHDP_RPT_3"
    xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    print(xvar)

    costing_data = psDataManager(f)


    costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
        costing_data,
        save_name=case_name,
        device_groups=figure_device_groups[case_name],
        fig_options=fig_options,
    )

    leg1_kwargs = dict(
        bbox_to_anchor=(0.0, 1.0),
        loc="upper left",
        ncols=2,
        fontsize=8,
        # mode="expand"
    )

    leg2_kwargs = dict(
        # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
        loc="best",
        fontsize=8,
    )

    fig, ax = create_fig_legend(
        cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    )

    label_kwargs = dict(
        # ylabel="LCOW (\$/m$^3$)",
        # xlabel="Fraction Heat From Grid",
        # xlabel="Brine Disposal Cost (\$/m$^3$)", 
        # xlabel="Water Recovery (%)",
        # ylim=(ax.get_ylim()[0], 7),
        # title=f"{case_name.replace('_', ' ')}\nLCOW vs. Brine Disposal Cost",
        yscale="log"
    )

    fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # plt.tight_layout()
    plt.show()

    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.svg", dpi=500)
    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.png", dpi=500)





    # ######## SOA ################ SOA ################ SOA ################ SOA ################ SOA ################ SOA ########
    # ######## SOA ################ SOA ################ SOA ################ SOA ################ SOA ################ SOA ########
    # ######## SOA ################ SOA ################ SOA ################ SOA ################ SOA ################ SOA ########

    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-203726_analysisType_KBHDP_SOA_1_water_recovery_sweep.h5"

    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)

    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name="KBHDP_RPT_2-test",
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     # xlabel="Deep Well Injection Cost (\$/m$^3$)",
    #     # xlabel="Fraction Heat From Grid",
    #     xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 4.5),
    #     title="KBHDP SOA\nLCOW vs. Water Recovery",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()

    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.svg", dpi=500)
    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.png", dpi=500)

    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-203726_analysisType_KBHDP_SOA_1_disposal_cost_sweep.h5"
    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)

    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name="KBHDP_RPT_2-test",
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     xlabel="Brine Disposal Cost (\$/m$^3$)",
    #     # xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 4.5),
    #     title="KBHDP SOA\nLCOW vs. Brine Disposal Costs",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()

    # # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.svg", dpi=500)
    # # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.png", dpi=500)

    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-204718_analysisType_KBHDP_SOA_1_soda_ash_cost_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)

    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     xlabel="Soda Ash Cost ($/kg)",
    #     # xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     # ylim=(ax.get_ylim()[0], 4.5),
    #     title="KBHDP SOA\nLCOW vs. Soda Ash Costs",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()


    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-204718_analysisType_KBHDP_SOA_1_co2_cost_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)


    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     xlabel="Carbon Dioxide Cost ($/kg)",
    #     # xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 3.5),
    #     title="KBHDP SOA\nLCOW vs. Carbon Dioxide Costs",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()



    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-204718_analysisType_KBHDP_SOA_1_mgcl2_cost_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)


    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     xlabel="MgCl$_2$ Cost (\$/kg)",
    #     # xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 3.5),
    #     title="KBHDP SOA\nLCOW vs. MgCl$_2$ Costs",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()

    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_SOA_1_20241219-204718_analysisType_KBHDP_SOA_1_lime_cost_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_SOA_1"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)


    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     xlabel="Lime Cost (\$/kg)",
    #     # xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 3.5),
    #     title="KBHDP SOA\nLCOW vs. Lime Costs",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()

######## RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ###########
######## RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ###########
######## RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ################### RPT1 ###########







# ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ##########
# ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ##########
# ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ########### ######## RPT2 ##########


    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_RPT_2_20241219-213759_analysisType_KBHDP_RPT_2_frac_heat_from_grid_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_RPT_2"
    # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # print(xvar)

    # costing_data = psDataManager(f)


    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     # xlabel="Lime Cost (\$/kg)",
    #     xlabel="Fraction Heat From Grid",
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 8),
    #     title=f"{case_name.replace('_', ' ')}\nLCOW vs. Fraction Heat From Grid",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()



    # f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/sweep_results/output/KBHDP_RPT_2_20241219-213759_analysisType_KBHDP_SOA_2_disposal_cost_sweep.h5"

    # # create_all_figures(fig_options=fig_options)
    # case_name = "KBHDP_RPT_2"
    # # xvar = f.split(case_name + "_")[-1].replace(".h5", "").replace("_sweep", "")
    # xvar = f.split(case_name + "_")
    # print(xvar)
    # assert False
    # costing_data = psDataManager(f)


    # costing_data, device_groups, cost_plotter = create_sweep_cost_breakdown(
    #     costing_data,
    #     save_name=case_name,
    #     device_groups=figure_device_groups[case_name],
    #     fig_options=fig_options,
    # )

    # leg1_kwargs = dict(
    #     bbox_to_anchor=(0.0, 1.0),
    #     loc="upper left",
    #     ncols=2,
    #     fontsize=8,
    #     # mode="expand"
    # )

    # leg2_kwargs = dict(
    #     # bbox_to_anchor=(0.28, 1.0, 1.0, 0.102),
    #     loc="best",
    #     fontsize=8,
    # )

    # fig, ax = create_fig_legend(
    #     cost_plotter, leg1_kwargs=leg1_kwargs, leg2_kwargs=leg2_kwargs
    # )

    # label_kwargs = dict(
    #     ylabel="LCOW (\$/m$^3$)",
    #     # xlabel="Fraction Heat From Grid",
    #     xlabel="Brine Disposal Cost (\$/m$^3$)", 
    #     # xlabel="Water Recovery (%)",
    #     ylim=(ax.get_ylim()[0], 7),
    #     title=f"{case_name.replace('_', ' ')}\nLCOW vs. Brine Disposal Cost",
    # )

    # fig, ax = create_axis_labels_etc(fig, ax, label_kwargs=label_kwargs)

    # # plt.tight_layout()
    # plt.show()

    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.svg", dpi=500)
    # fig.savefig(f"{fig_dir}/{case_name}_{xvar}.png", dpi=500)