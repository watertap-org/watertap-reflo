import os
import pprint
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import OrderedDict
import os
import math
from idaes.core.base.costing_base import register_idaes_currency_units

# from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.h5_to_csv import (
#     convert_h5_to_csv,
# )

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.dirname(__location__)
sweep_csv_dir = os.path.join(parent_dir, "sweep_csvs")
fig_save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/figures"
unit_list = [
    "Chem. Softening",
    "CST",
    "H$_2$O$_2$ Addition",
    "EC",
    "UF",
    "RO",
    "LT-MED",
    "MEC",
    "DWI",
    "FO",
    "MD",
    "Pump",
    "PV",
    "FPC",
    "CF",
    "MVC",
]

flow_list = [
    "H$_2$O$_2$",
    "Heat",
    "Electricity",
    "Aluminum",
    "Lime",
    "CO$_2$",
    "Soda Ash",
]

unit_colors = plt.cm.tab20(np.arange(len(unit_list)).astype(int))
unit_color_dict_default = dict(zip(unit_list, unit_colors))
unit_color_dict_default = {
    "CF": np.array([0.49803922, 0.49803922, 0.49803922, 1.0]),
    #  'CST': np.array([0.68235294, 0.78039216, 0.90980392, 1.        ]),
    "CST": "goldenrod",
    "Chem. Softening": np.array([0.12156863, 0.46666667, 0.70588235, 1.0]),
    "DWI": np.array([0.58039216, 0.40392157, 0.74117647, 1.0]),
    "EC": np.array([1.0, 0.73333333, 0.47058824, 1.0]),
    "FO": np.array([0.77254902, 0.69019608, 0.83529412, 1.0]),
    # "FPC": np.array([0.96862745, 0.71372549, 0.82352941, 1.0]),
    "FPC": "lightsteelblue",
    "H$_2$O$_2$ Addition": np.array([1.0, 0.49803922, 0.05490196, 1.0]),
    "LT-MED": np.array([0.83921569, 0.15294118, 0.15686275, 1.0]),
    "MD": np.array([0.54901961, 0.3372549, 0.29411765, 1.0]),
    "MEC": np.array([1.0, 0.59607843, 0.58823529, 1.0]),
    "MVC": np.array([0.78039216, 0.78039216, 0.78039216, 1.0]),
    "PV": np.array([0.89019608, 0.46666667, 0.76078431, 1.0]),
    "Pump": np.array([0.76862745, 0.61176471, 0.58039216, 1.0]),
    "RO": np.array([0.59607843, 0.8745098, 0.54117647, 1.0]),
    "UF": np.array([0.17254902, 0.62745098, 0.17254902, 1.0]),
}
flow_colors = plt.cm.Dark2(np.arange(len(flow_list)).astype(int))
flow_color_dict_default = dict(zip(flow_list, flow_colors))
flow_color_dict_default["electric"] = flow_color_dict_default["Electricity"]

register_idaes_currency_units()


def heat_func_row(row):
    heat_cost = 0.00894 * pyunits.USD_2023 / pyunits.kWh  # $ / kWh, KBHDP
    heat_cost = 0.0166 * pyunits.USD_2023 / pyunits.kWh  # $ / kWh, Permian
    agg_flow_heat = (
        row.loc[f"fs.treatment.costing.aggregate_flow_heat"] * pyunits.kW
    )  # kW
    agg_flow_heat_annual = pyunits.convert(
        agg_flow_heat, to_units=pyunits.kWh / pyunits.year
    )  # kWh / year
    return pyunits.convert(
        agg_flow_heat_annual * heat_cost, to_units=pyunits.USD_2023 / pyunits.year
    )()


def elec_func_row(row):
    elec_cost = 0.066 * pyunits.USD_2023 / pyunits.kWh  # $ / kWh KBHDP
    elec_cost = 0.0575 * pyunits.USD_2023 / pyunits.kWh  # $ / kWh Permian
    agg_flow_elec = (
        row.loc[f"fs.treatment.costing.aggregate_flow_electricity"] * pyunits.kW
    )  # kW
    agg_flow_elec_annual = pyunits.convert(
        agg_flow_elec, to_units=pyunits.kWh / pyunits.year
    )  # kWh / year
    return pyunits.convert(
        agg_flow_elec_annual * elec_cost, to_units=pyunits.USD_2023 / pyunits.year
    )()


def case_study_stacked_plot(
    df,
    fig=None,
    ax=None,
    fig_rel=None,
    ax_rel=None,
    global_costing_blk="fs.treatment.costing",
    costing_blk="fs.costing",
    xcol=None,
    flow_col=None,  # column to be used as denominator in LCOW calculations, assumed to be in m3/s
    unit_dict=dict(),  # (unit name: unit location)
    agg_flows=dict(),  # dict of aggregated flows
    figsize=(6, 4),
    unit_color_dict=unit_color_dict_default,
    flow_color_dict=flow_color_dict_default,
    capex_hatch="",
    opex_hatch="\\\\\\",
    flow_hatch="..",
    ax_dict=dict(),
    check_calc=True,
    label_fontsize=16,
    tick_fontsize=14,
    add_legend=True,
    leg_kwargs=dict(
        # loc="lower left",
        frameon=False,
        ncol=3,
        handlelength=1,
        handleheight=1,
        labelspacing=0.2,
        columnspacing=0.9,
        # bbox_to_anchor=(0., 1.02, 1., .102),
        # mode="expand",
    ),
    leg_kwargs_rel=dict(
        loc="lower left",
        frameon=False,
        ncol=3,
        handlelength=1,
        handleheight=1,
        labelspacing=0.2,
        columnspacing=0.9,
        bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
        mode="expand",
    ),
    ylims=None,
    xlims=None,
    save=False,
    save_name="temp",
    rel=False,
    actual_lcow_row=None,
    heat_func=None,
    heat_cost_col=None,
    elec_func=None,
    elec_cost_col=None,
    use_calc_for_rel=False,
    **kwargs,
):
    global figure_csv_rel, figure_csv

    if flow_hatch is None:
        flow_hatch = capex_hatch

    if global_costing_blk is None:
        global_costing_blk = costing_blk

    global_params = [
        "maintenance_labor_chemical_factor",
        "utilization_factor",
        "capital_recovery_factor",
        "total_investment_factor",
    ]

    opex_lcow = defaultdict(list)
    capex_lcow = defaultdict(list)
    agg_flow_lcow = defaultdict(list)

    opex_lcow_rel = defaultdict(list)
    capex_lcow_rel = defaultdict(list)
    agg_flow_lcow_rel = defaultdict(list)

    actual_lcow = list()
    calc_lcow = list()
    calc_to_actual_lcow = list()
    capex = list()
    opex = list()

    costing_params = {
        "maintenance_labor_chemical_factor": 0.03,
        "utilization_factor": 1.0,
        "capital_recovery_factor": 0.1119559492025644,
        "total_investment_factor": 1.0,
    }

    df.set_index(xcol, inplace=True)
    df.sort_values(by=[xcol], inplace=True)

    for x, row in df.iterrows():

        total_lcow = 0
        total_capex = 0
        total_opex = 0
        total_flow_cost = 0
        total_annualized_cost = 0

        if actual_lcow_row is None:
            row_lcow = row[f"{costing_blk}.LCOW"]
            actual_lcow_row = f"{costing_blk}.LCOW"
        else:
            row_lcow = row[actual_lcow_row]

        actual_lcow.append(row_lcow)
        denominator = (
            value(
                pyunits.convert(
                    row.loc[flow_col] * pyunits.m**3 / pyunits.s,
                    to_units=pyunits.m**3 / pyunits.year,
                )
            )
            * costing_params["utilization_factor"]
        )

        for u, b in unit_dict.items():

            print(f"Getting CAPEX and OPEX for {u} in {b}")
            unit_capex = 0

            try:
                unit_capex += (
                    row.loc[f"{b}.capital_cost"]
                    * costing_params["total_investment_factor"]
                )  # USD2023

            except KeyError:
                print(f"No CAPEX for {u} found.")
                pass

            total_capex += unit_capex  # $
            total_annualized_cost += (
                unit_capex * costing_params["capital_recovery_factor"]
            )  # $ / year

            unit_capex_lcow = (
                unit_capex * costing_params["capital_recovery_factor"]
            ) / denominator  # $ / m3

            unit_capex_lcow_rel = unit_capex_lcow / row_lcow
            total_lcow += unit_capex_lcow  # $ / m3

            if unit_capex_lcow > 1e-13:
                capex_lcow[u].append(unit_capex_lcow)

            capex_lcow_rel[u].append(unit_capex_lcow_rel)  # $ / m3

            ### OPEX
            unit_opex_total = 0
            unit_opex_total += (
                unit_capex * costing_params["maintenance_labor_chemical_factor"]
            )  # $ / year

            try:
                unit_opex_total += row.loc[f"{b}.fixed_operating_cost"]
            except KeyError:
                print(f"No Fixed OPEX for {u} found.")
                pass

            try:
                unit_opex_total += row.loc[f"{b}.variable_operating_cost"]
            except KeyError:
                print(f"No Variable OPEX for {u} found at {b}.variable_operating_cost.")
                pass

            total_opex += unit_opex_total  # $ / year
            unit_opex_lcow = unit_opex_total / denominator  # $ / m3

            unit_opex_lcow_rel = unit_opex_lcow / row_lcow

            if unit_opex_lcow > 1e-12:
                opex_lcow[u].append(unit_opex_lcow)
            if unit_opex_lcow_rel < 1e-8:
                unit_opex_lcow_rel = 0

            opex_lcow_rel[u].append(unit_opex_lcow_rel)
            total_lcow += unit_opex_lcow  # $ / m3
            total_annualized_cost += unit_opex_total  # $ / year

        for flow_label, flow_name in agg_flows.items():
            print(f"Calculating LCOW for {flow_label} at {xcol} = {x}")

            flow_lcow = 0
            # Electricity and heat flows can get tricky depending on how the model was run
            if flow_name == "electricity":
                if elec_func is not None:
                    # Use external function to calculate electricity cost
                    # This is useful for cases where the electricity cost is not a simple column in the dataframe
                    flow_lcow = elec_func(row) / denominator
                elif elec_cost_col is not None:
                    # Use a specific column for electricity cost
                    flow_lcow = row.loc[elec_cost_col] / denominator
                else:
                    try:
                        # Default case, use the aggregate flow costs from the costing block
                        flow_lcow += (
                            row.loc[f"{costing_blk}.aggregate_flow_costs[{flow_name}]"]
                            / denominator
                        )
                    except:
                        raise ValueError(
                            f"No Elecricity OPEX for {flow_name} found in {costing_blk}.aggregate_flow_costs[{flow_name}]"
                        )

                agg_flow_lcow[flow_name].append(flow_lcow)
                flow_lcow_rel = flow_lcow / row_lcow
                agg_flow_lcow_rel[flow_name].append(flow_lcow_rel)
                total_lcow += flow_lcow
                continue

            if flow_name == "heat":
                if heat_func is not None:
                    # Use external function to calculate heat cost
                    # This is useful for cases where the heat cost is not a simple column in the dataframe
                    flow_lcow = heat_func(row) / denominator
                elif heat_cost_col is not None:
                    # Use a specific column for heat cost
                    flow_lcow = row.loc[heat_cost_col] / denominator
                else:
                    try:
                        # Default case, use the aggregate flow costs from the costing block
                        # This is where we might expect to find them
                        flow_lcow += (
                            row.loc[f"{costing_blk}.aggregate_flow_costs[{flow_name}]"]
                            / denominator
                        )
                    except:
                        raise ValueError(
                            f"No Heat OPEX for {flow_name} found in {costing_blk}.aggregate_flow_costs[{flow_name}]"
                        )

                agg_flow_lcow[flow_name].append(flow_lcow)
                flow_lcow_rel = flow_lcow / row_lcow
                agg_flow_lcow_rel[flow_name].append(flow_lcow_rel)
                total_lcow += flow_lcow
                continue

            flow_cost = row.loc[
                f"{global_costing_blk}.aggregate_flow_costs[{flow_name}]"
            ]  # $ / year
            total_annualized_cost = row.loc[
                f"{global_costing_blk}.aggregate_flow_costs[{flow_name}]"
            ]
            flow_lcow = flow_cost / denominator  # $ / m3

            agg_flow_lcow[flow_name].append(flow_lcow)
            flow_lcow_rel = flow_lcow / row_lcow
            agg_flow_lcow_rel[flow_name].append(flow_lcow_rel)
            total_lcow += flow_lcow

        calc_lcow.append(total_lcow)
        calc_to_actual_lcow.append(total_lcow / row_lcow)

        if check_calc:
            print(f"\nFor {xcol} = {x}:")
            print(f"\tActual LCOW: {row_lcow:.6f}")
            print(f"\tCalculated LCOW: {total_lcow:.6f}")
            print(f"\tRel Diff: {row_lcow / total_lcow:.6f}")

        capex.append(total_capex)
        opex.append(total_opex)

    if use_calc_for_rel:
        # Recalculate the relative LCOW values using the calculated LCOW
        # This is helpful if the actual LCOW is not available in data due to how the model was set up
        opex_lcow_rel = defaultdict(list)
        capex_lcow_rel = defaultdict(list)
        agg_flow_lcow_rel = defaultdict(list)
        for i, lcow in enumerate(calc_lcow):
            for u in unit_dict.keys():
                if u in capex_lcow.keys():
                    x = capex_lcow[u][i] / lcow
                    capex_lcow_rel[u].append(capex_lcow[u][i] / lcow)
                if u in opex_lcow.keys():
                    opex_lcow_rel[u].append(opex_lcow[u][i] / lcow)
            for flow, y in agg_flow_lcow.items():
                # print(flow_label)
                agg_flow_lcow_rel[flow].append(y[i] / lcow)

    ###################################
    ###################################
    # Create the inputs needed for stacked plot

    stacked_cols = list()
    stacked_labels = list()
    stacked_hatch = list()
    stacked_colors = list()

    stacked_cols_rel = list()
    stacked_labels_rel = list()
    stacked_hatch_rel = list()
    stacked_colors_rel = list()

    legend_elements = [
        Patch(facecolor="white", hatch=capex_hatch, label="CAPEX", edgecolor="k"),
        Patch(facecolor="white", hatch=opex_hatch, label="OPEX", edgecolor="k"),
    ]

    for flow_label, flow_name in agg_flows.items():

        stacked_cols.append(agg_flow_lcow[flow_name])
        stacked_cols_rel.append(agg_flow_lcow_rel[flow_name])

        stacked_labels.append(flow_label)
        stacked_labels_rel.append(flow_label)

        stacked_colors.append(flow_color_dict[flow_label])
        stacked_colors_rel.append(flow_color_dict[flow_label])

        stacked_hatch.append(flow_hatch)
        stacked_hatch_rel.append(flow_hatch)

        legend_elements.append(
            Patch(
                facecolor=flow_color_dict[flow_label],
                label=flow_label,
                hatch=flow_hatch,
                edgecolor="k",
            )
        )

    for u in unit_dict.keys():
        if u in capex_lcow.keys():
            stacked_cols.append(capex_lcow[u])
            stacked_cols_rel.append(capex_lcow_rel[u])

            stacked_labels.append(f"{u} CAPEX")
            stacked_labels_rel.append(f"{u} CAPEX")

            stacked_hatch.append(capex_hatch)
            stacked_hatch_rel.append(capex_hatch)

            stacked_colors.append(unit_color_dict[u])
            stacked_colors_rel.append(unit_color_dict[u])
        if u in opex_lcow.keys():
            stacked_cols.append(opex_lcow[u])
            stacked_cols_rel.append(opex_lcow_rel[u])

            stacked_labels.append(f"{u} OPEX")
            stacked_labels_rel.append(f"{u} OPEX")

            stacked_hatch.append(opex_hatch)
            stacked_hatch_rel.append(opex_hatch)

            stacked_colors.append(unit_color_dict[u])
            stacked_colors_rel.append(unit_color_dict[u])

        legend_elements.append(
            Patch(facecolor=unit_color_dict[u], label=u, edgecolor="k")
        )

    ###################################
    ###################################
    # Create stacked plot with absolute values

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(figsize=figsize, layout="constrained")
        fig.set_size_inches(5, 5, forward=True)

    ax.stackplot(
        df.index,
        stacked_cols,
        # labels=stacked_labels,
        hatch=stacked_hatch,
        colors=stacked_colors,
        edgecolor="black",
    )

    if add_legend:
        # legend = ax.legend(handles=legend_elements, **leg_kwargs)
        legend = ax.legend(handles=legend_elements, **leg_kwargs_rel)

    ax.set(**ax_dict)
    ax.set_xlabel(ax_dict["xlabel"], fontsize=label_fontsize)
    ax.set_ylabel(ax_dict["ylabel"], fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.set_xlim(df.index.min(), df.index.max())
    if any(x in xcol for x in ["soda_ash", "cost_per_watt"]):
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"${x:.2f}"))
    elif xcol == "heat_cost":
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"¢{x*100:.1f}"))
    else:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%"))
    # ax2 = ax.twinx()
    # ax2.plot(df.index, df[actual_lcow_row], color="k", linewidth=5)
    if ylims is not None:
        ax.set_ylim(ylims)
        # ax2.set_ylim(ylims)
    else:
        ax.set_ylim(0, np.ceil(max(actual_lcow)))
        # ax2.set_ylim(0, np.ceil(max(actual_lcow)))

    if xlims is not None:
        ax.set_xlim(xlims)
    else:
        ax.set_xlim(0, 1)

    fig.tight_layout()

    ###################################
    ###################################
    # Create stacked plot with relative values

    if (fig_rel, ax_rel) == (None, None):
        fig_rel, ax_rel = plt.subplots(figsize=figsize, layout="constrained")
        fig_rel.set_size_inches(5, 5, forward=True)

    ax_rel.stackplot(
        df.index,
        stacked_cols_rel,
        hatch=stacked_hatch_rel,
        colors=stacked_colors_rel,
        edgecolor="black",
    )

    if add_legend:
        # legend_rel = ax_rel.legend(handles=legend_elements, **leg_kwargs)
        legend_rel = ax_rel.legend(handles=legend_elements, **leg_kwargs_rel)

    ax_rel.set(**ax_dict)
    ax_rel.set_xlabel(ax_dict["xlabel"], fontsize=label_fontsize)
    ax_rel.set_ylabel("% LCOW", fontsize=label_fontsize)
    ax_rel.tick_params(axis="both", labelsize=tick_fontsize)
    ax_rel.set_xlim(df.index.min(), df.index.max())

    if any(x in xcol for x in ["soda_ash", "cost_per_watt"]):
        ax_rel.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"${x:.2f}"))
    elif xcol == "heat_cost":
        ax_rel.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"¢{x*100:.1f}"))
    else:
        ax_rel.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%")
        )
    if ylims is not None:
        ax_rel.set_ylim(ylims)
    else:
        ax_rel.set_ylim(0, np.ceil(max(actual_lcow)))
    ax_rel.set_ylim((0, 1.05))
    ax_rel.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%"))

    ax_rel2 = ax_rel.twinx()
    ax_rel2.plot(df.index, df[actual_lcow_row], color="black", linewidth=3, alpha=0.5)
    ax_rel2.set_ylim((0, 1.05 * df[actual_lcow_row].max()))
    ax_rel2.set_ylabel("", fontsize=label_fontsize)

    ax_rel2_ticks = np.linspace(0, df[actual_lcow_row].max(), 5)
    ax_rel2.tick_params(axis="both", labelsize=tick_fontsize)
    ax_rel2.set_yticks(ax_rel2_ticks)
    ax_rel2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.1f}"))

    if xlims is not None:
        ax_rel.set_xlim(xlims)
        # ax_rel2.set_xlim(xlims)
    else:
        ax_rel.set_xlim(0, 1)

    fig_rel.tight_layout()

    if global_save:
        fig.savefig(f"{fig_save_path}/{save_name}", dpi=300, bbox_inches="tight")
        fig_rel.savefig(
            f"{fig_save_path}/{save_name.replace('.png', '_rel.png')}",
            dpi=300,
            bbox_inches="tight",
        )

    ###################################
    ###################################
    # Create DataFrame for both figures

    unit_lcow = dict()
    unit_lcow_rel = dict()
    for unit, key in unit_dict.items():
        unit_lcow_rel[f"{unit} OPEX LCOW rel"] = opex_lcow_rel[unit]
        if len(np.array(capex_lcow[unit])) > 1:
            unit_lcow[unit] = np.array(capex_lcow[unit]) + np.array(opex_lcow[unit])
            unit_lcow_rel[f"{unit} CAPEX LCOW rel"] = capex_lcow_rel[unit]
        else:
            unit_lcow[unit] = np.array(opex_lcow[unit])

    for flow, key in agg_flows.items():
        unit_lcow[key] = np.array(agg_flow_lcow[key.lower()])
        unit_lcow_rel[key] = np.array(agg_flow_lcow_rel[key.lower()])

    # for k, v in unit_lcow.items():
    #     print(k, len(v))
    # assert False
    # figure_csv = pd.DataFrame.from_dict(unit_lcow)
    # figure_csv["LCOW"] = figure_csv.sum(axis=1)
    # figure_csv["LCOW_data"] = actual_lcow
    # figure_csv.index = df.index

    # # if global_save:
    # #     figure_csv.to_csv(f"{fig_save_path}/{save_name.replace('.png', '.csv')}", index=True)

    # figure_csv_rel = pd.DataFrame.from_dict(unit_lcow_rel)
    # figure_csv_rel["LCOW"] = figure_csv_rel.sum(axis=1)
    # figure_csv_rel["LCOW_data"] = actual_lcow
    # figure_csv_rel["LCOW_calc"] = calc_lcow
    # figure_csv_rel["LCOW_diff_rel"] = calc_to_actual_lcow
    # figure_csv_rel.index = df.index

    # # if global_save:
    # #     figure_csv_rel.to_csv(f"{fig_save_path}/{save_name.replace('.png', '_rel.csv')}", index=True)
    # # print(figure_csv.head(30))
    # print(figure_csv_rel.head(30))

    # if "electricity" in figure_csv_rel.columns:
    #     print(figure_csv_rel.electricity)
    # if "heat" in figure_csv_rel.columns:
    #     print(figure_csv_rel.heat)
    # assert False
    # print(figure_csv_rel["DWI OPEX LCOW rel"])
    return fig, ax, fig_rel, ax_rel, legend


cases_kbhdp_soa = {
    "KBHDP": {
        "KBHDP_SOA_1": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_SOA_1_water_recovery.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "Chem. Softening": "fs.treatment.softener.unit.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "Pump": "fs.treatment.pump.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Soda Ash": "soda_ash",
                    "Lime": "lime",
                    "CO$_2$": "co2",
                },
                "xcol": "fs.water_recovery",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 6),
                "xlims": (0.3, 0.8),
                "save_name": "kbhdp_soa_water_recovery_stacked_plot.png",
                "save": True,
            },
            "soda_ash": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_SOA_1_soda_ash_cost.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "Chem. Softening": "fs.treatment.softener.unit.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "Pump": "fs.treatment.pump.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Soda Ash": "soda_ash",
                    "Lime": "lime",
                    "CO$_2$": "co2",
                },
                "xcol": "fs.treatment.costing.soda_ash.cost",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Soda Ash Cost ($/kg)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 3),
                "xlims": (0.12, 0.36),
                "save_name": "kbhdp_soa_soda_ash_sweep_stacked_plot.png",
                # "save": True,
            },
        }
    }
}

kbhdp_opt = {
    "KBHDP": {
        "KBHDP_RPT_1": {
            "cost_per_watt_installed": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/test.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "Pump": "fs.treatment.pump.costing",
                    "PV": "fs.energy.pv.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.energy.costing.pv_surrogate.cost_per_watt_installed",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol",
                "ax_dict": dict(
                    xlabel="PV Cost Per Watt Installed ($/W)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 1.25),
                "xlims": (0.25, 1.5),
                # "save_name": "kbhdp_rpt1_water_recovery_stacked_plot.png",  # Change this
                "save_name": "kbhdp_optimization_cost_per_watt_installed.png",  # Change this
                "save": False,
            },
        },
    }
}


permian_opt = {
    "Permian": {
        "Permian_RPT2_FO": {
            "grid_fraction_optimize": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/sweep_results/permian_RPT2_FO_DWI_RPT_heat_cost_sweep-OPTIMIZE-2.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOW",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "heat_cost",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Heat Cost (¢/kWh)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 10),
                "xlims": (0.017, 0.133986),  # Change this
                "save_name": "permian_heat_price_stacked_plot-OPTIMIZE.png",  # Change this
                "save": False,
                "use_calc_for_rel": True,
            },
        },
        "Permian_ZLD2_FO_Cryst": {
            "grid_fraction_optimize": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/sweep_results/permian_ZLD2_FO_cryst_RPT_heat_cost_sweep-OPTIMIZE-3.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOW",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                # "xcol": "fs.costing.RE Fraction",  # Change this
                "xcol": "heat_cost",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Heat Cost (¢/kWh)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 10),
                "xlims": (0.017, 0.133986),  # Change this
                "save_name": "permian_heat_price_stacked_plot-OPTIMIZE.png",  # Change this
                "save": False,
                "use_calc_for_rel": True,
            },
        }
    }
}

kbhdp_wr = {
    "KBHDP": {
        "KBHDP_RPT_1": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_RPT_1_water_recovery.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "Pump": "fs.treatment.pump.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.water_recovery",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 5),
                "xlims": (0.3, 0.8),
                # "save_name": "kbhdp_rpt1_water_recovery_stacked_plot.png",  # Change this
                "save_name": "kbhdp_rpt_water_recovery_stacked_plot.png",  # Change this
                "save": False,
            },
        },
        "KBHDP_RPT_2": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_RPT_2_water_recovery.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                # "actual_lcow_row": "fs.treatment.costing.LCOW",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_func": heat_func_row,
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "LT-MED": "fs.treatment.LTMED.unit.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                    "Heat": "heat",
                },
                "xcol": "fs.water_recovery",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 5),
                "xlims": (0.3, 0.8),
                # "save_name": "kbhdp_rpt2_water_recovery_stacked_plot.png",  # Change this
                "save_name": "kbhdp_rpt_water_recovery_stacked_plot.png",  # Change this
                # "save": True,
                "use_calc_for_rel": True,
            },
        },
        "KBHDP_RPT_3": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/kbhdp_RPT3_water_recovery.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                # "actual_lcow_row": "fs.costing.LCOT",
                "unit_dict": {
                    "MD": "fs.treatment.md.unit",
                    "DWI": "fs.treatment.dwi.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                },
                "xcol": "fs.water_recovery",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 5),
                "xlims": (0.3, 0.8),
                # "save_name": "kbhdp_rpt3_water_recovery_stacked_plot.png",  # Change this
                "save_name": "kbhdp_rpt_water_recovery_stacked_plot.png",  # Change this
                # "save": False,
            },
        },
    }
}

kbhdp_grid_frac = {
    "KBHDP": {
        "KBHDP_RPT_1": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_RPT_1_grid_fraction.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "Pump": "fs.treatment.pump.costing",
                    "PV": "fs.energy.pv.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 2),
                "xlims": (0.1, 0.8),
                # "save_name": "kbhdp_rpt1_grid_frac_stacked_plot.png",  # Change this
                "save_name": "kbhdp_grid_frac_stacked_plot.png",  # Change this
                "save": False,
            },
        },
        "KBHDP_RPT_2": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/sweep_data_KBHDP_RPT_2_frac_heat_from_grid.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_func": heat_func_row,
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "LT-MED": "fs.treatment.LTMED.unit.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "FPC": "fs.energy.FPC.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                    "Heat": "heat",
                },
                "xcol": "fs.costing.RE Fraction",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "xlims": (0.1, 0.8),
                "ylims": (0, 16.5),
                "save_name": "kbhdp_grid_frac_stacked_plot.png",  # Change this
                "save": True,
                "use_calc_for_rel": True,
            },
        },
        "KBHDP_RPT_3": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/kbhdp_RPT3_grid_frac_heat_grid_frac_var_recovery_0.8.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "MD": "fs.treatment.md.unit",
                    "DWI": "fs.treatment.dwi.unit.costing",
                    "FPC": "fs.energy.FPC.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                },
                "xcol": "fs.costing.RE Fraction",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "xlims": (0.1, 0.8),
                "ylims": (0, 16.5),
                "save_name": "kbhdp_grid_frac_stacked_plot.png",  # Change this
                "save": True,
            },
        },
    }
}


kbhdp_zld = {
    "KBHDP": {
        "KBHDP_ZLD": {
            "cost_per_aperture_area": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/kbhdp/kbhdp_ZLD_cst_cost_per_total_aperture_area.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "EC": "fs.treatment.EC.ec.costing",
                    "UF": "fs.treatment.UF.unit.costing",
                    "RO": "fs.treatment.RO.stage[1].module.costing",
                    "MD": "fs.treatment.md.unit",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "PV": "fs.energy.pv.costing",
                    "CST": "fs.energy.cst.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.energy.costing.trough_surrogate.cost_per_total_aperture_area",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Cost Per CST Area (\$/m$^2$)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 2),
                "xlims": (145, 445),
                "save_name": "kbhdp_zld_cost_per_aperture_area_stacked_plot.png",  # Change this
                "save": False,
            },
        },
    }
}


permian_wr = {
    "Permian": {
        "Permian_RPT1_MD": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT1_MD_water_recovery_grid_frac_1_recovery_var.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.EC.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "MD": "fs.treatment.md.unit",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.water_recovery",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 25),
                "xlims": (0.35, 0.5),
                "save_name": "permian_water_recovery_stacked_plot.png",  # Change this
                "save": False,
            }
        },
        "Permian_RPT2_FO": {
            "water_recovery": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT2_FO_DWI_no_CST_sweep_recovery_ratio.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.treatment.FO.fs.fo.recovery_ratio",
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 25),
                "xlims": (0.35, 0.5),
                "save_name": "permian_water_recovery_stacked_plot.png",  # Change this
                "save": False,
            }
        },
    }
}
#         f"No Heat OPEX for {flow_name} found in in fs.costing.total_heat_operating_cost."

permian_grid_frac = {
    "Permian": {
        "Permian_RPT1_MD": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.EC.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "MD": "fs.treatment.md.unit",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 35),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "permian_grid_frac_stacked_plot.png",  # Change this
                "save": False,
            },
        },
        "Permian_RPT2_FO": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT2_FO_DWI_RPT_grid_frac.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOW",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 20),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "permian_grid_frac_stacked_plot.png",  # Change this
                "save": False,
            },
        },
        "Permian_ZLD1_MD_Cryst": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_ZLD1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOT",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.EC.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "MD": "fs.treatment.md.unit",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 20),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "permian_grid_frac_stacked_plot.png",  # Change this
                "save": False,
            },
        },
        "Permian_ZLD2_FO_Cryst": {
            "grid_fraction": {
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_ZLD2_FO_cryst_RPT_grid_frac.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "actual_lcow_row": "fs.costing.LCOW",
                "heat_cost_col": "fs.costing.total_heat_operating_cost",
                "elec_cost_col": "fs.costing.total_electric_operating_cost",
                "unit_dict": {
                    "H$_2$O$_2$ Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H$_2$O$_2$": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(xlabel="Solar Energy (%)", ylabel="LCOW (\$/m$^3$)"),
                "ylims": (0, 20),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "permian_grid_frac_stacked_plot.png",  # Change this
                "save": False,
                "use_calc_for_rel": True,
            },
        },
    }
}


def plot_case(case, fig=None, ax=None, fig_rel=None, ax_rel=None):
    print(case)
    df = pd.read_csv(case["file"])
    df.dropna(inplace=True)

    try:
        try:
            df["fs.costing.RE Fraction"] = 1 - df["fs.costing.frac_heat_from_grid"]
        except:
            df["fs.costing.RE Fraction"] = 1 - df["fs.costing.frac_elec_from_grid"]
    except:
        print("NO GRID FRACTION FOUND")
        pass
    # else:
    #     print('NO GRID FRACTION FOUND')
    #     pass
    # print(df['fs.costing.RE Fraction'])
    if "actual_lcow_row" not in case.keys():
        case["actual_lcow_row"] = None
    if "heat_func" not in case.keys():
        case["heat_func"] = None
    case["save"] = global_save

    fig_, ax_, fig_rel_, ax_rel_, legend_ = case_study_stacked_plot(
        df,
        fig=fig,
        ax=ax,
        fig_rel=fig_rel,
        ax_rel=ax_rel,
        **case,
        # global_costing_blk=case["global_costing_blk"],
        # costing_blk=case["costing_blk"],
        # unit_dict=case["unit_dict"],
        # agg_flows=case["agg_flows"],
        # xcol=case["xcol"],
        # flow_col=case["flow_col"],
        # ax_dict=case["ax_dict"],
        # opex_hatch="\\\\\\",
        # flow_hatch="..",
        # ylims=case["ylims"],
        # xlims=case["xlims"],
        # # save=case["save"],
        # save=True,
        # # save=False,
        # save_name=case["save_name"],
        # actual_lcow_row=case["actual_lcow_row"],
        # heat_func=case["heat_func"],
    )

    return legend_


def plot_all_cases(cases, xdim=4, ydim=3.5, legend_rows=2, bboxy=1.1):

    figs = {}
    axs = {}
    figs_rel = {}
    axs_rel = {}
    legends = {}

    for idx0, study in enumerate(cases):

        figs[study] = {}
        axs[study] = {}
        figs_rel[study] = {}
        axs_rel[study] = {}
        legends[study] = {}

        for key in cases[study][list(cases[study].keys())[0]].keys():

            n_cases = len(cases[study])
            n_rows = int(n_cases // 4 + 1)
            n_cols = int(n_cases / n_rows)

            figs[study][key], axs[study][key] = plt.subplots(
                n_rows,
                n_cols,
                figsize=(xdim * n_cols, ydim * n_rows),
                # constrained_layout=True,
            )

            figs_rel[study][key], axs_rel[study][key] = plt.subplots(
                n_rows,
                n_cols,
                figsize=((xdim + 0.6) * n_cols, (ydim - 0.5) * n_rows),
                # constrained_layout=True,
            )
            legends[study][key] = []
            axs[study][key] = np.atleast_1d(axs[study][key]).ravel()
            axs_rel[study][key] = np.atleast_1d(axs_rel[study][key]).ravel()

    fig_letters = ["A", "B", "C", "D"]

    for idx0, study in enumerate(cases):
        for key in cases[study][list(cases[study].keys())[0]].keys():
            for idx1, case in enumerate(cases[study]):
                legend_ = plot_case(
                    cases[study][case][key],
                    fig=figs[study][key],
                    ax=axs[study][key][idx1],
                    fig_rel=figs_rel[study][key],
                    ax_rel=axs_rel[study][key][idx1],
                )
                legends[study][key].append(legend_)

                axs[study][key][idx1].text(
                    -0.25,
                    0.95,
                    f"{fig_letters[idx1]})",
                    fontsize=20,
                    ha="center",
                    va="center",
                    transform=axs[study][key][idx1].transAxes,
                )

                axs_rel[study][key][idx1].text(
                    -0.33,
                    0.95,
                    f"{fig_letters[idx1]})",
                    fontsize=20,
                    ha="center",
                    va="center",
                    transform=axs_rel[study][key][idx1].transAxes,
                )

            # figs[study][key].tight_layout()
            # figs_rel[study][key].tight_layout()

            if len(axs[study][key]) < 4:
                figs[study][key].subplots_adjust(top=0.85)
                figs_rel[study][key].subplots_adjust(top=0.85)
            else:
                figs[study][key].subplots_adjust(top=0.925)
                figs_rel[study][key].subplots_adjust(top=0.925)

            # Collect handles and labels
            handles, labels = [], []
            for legend in legends[study][key]:
                h = legend.legend_handles
                l = [text.get_text() for text in legend.get_texts()]
                handles.extend(h)
                labels.extend(l)

            # Remove duplicates while preserving order
            unique = list(OrderedDict(zip(labels, handles)).items())
            # Ensure order to be:
            # - flows
            # - units
            unique_flows = list()
            unique_units = list()
            for u in unique:
                if u[0] in ["CAPEX", "OPEX"]:
                    continue
                if u[0] in flow_color_dict_default.keys():
                    unique_flows.append(u)
                else:
                    unique_units.append(u)

            unique = [unique[0], unique[1]] + unique_flows + unique_units

            # Add the unified legend
            # legend_rows = 2
            # legend_rows = 1
            # if len(axs[study][key]) > 3:
            #     legend_rows = 2
            legend_cols = len(unique) // legend_rows + 1
            # legend_cols = 6

            master_legend = figs[study][key].legend(
                [h for _, h in unique],
                [l for l, _ in unique],
                loc="upper center",  # or 'lower center' if you prefer
                ncol=legend_cols,
                bbox_to_anchor=(0, bboxy, 1, 0.05),  # above the plots
                frameon=True,
                # handlelength=1.4,
                # handleheight=1.4,
                # labelspacing=0.25,
                # handletextpad=0.4,
                # columnspacing=0.9,
                fontsize=10,
                mode="expand",
            )
            master_legend_rel = figs_rel[study][key].legend(
                [h for _, h in unique],
                [l for l, _ in unique],
                loc="upper center",  # or 'lower center' if you prefer
                ncol=legend_cols,
                bbox_to_anchor=(0, bboxy, 1, 0.05),  # above the plots
                frameon=True,
                # handlelength=1.4,
                # handleheight=1.4,
                # labelspacing=0.25,
                # handletextpad=0.4,
                # columnspacing=0.9,
                fontsize=10,
                mode="expand",
            )
            master_legend.get_frame().set_edgecolor("black")  # Frame line color
            master_legend.get_frame().set_linewidth(1.5)  # Thickness of the frame

            master_legend_rel.get_frame().set_edgecolor("black")  # Frame line color
            master_legend_rel.get_frame().set_linewidth(1.5)  # Thickness of the frame

            for axis in axs[study][key]:
                axis.legend_.remove()  # This works if the legend exists

            for axis in axs_rel[study][key]:
                axis.legend_.remove()  # This works if the legend exists

            # Get the bounding box of the original legend
            bb = master_legend.get_bbox_to_anchor().transformed(
                axs[study][key][idx1].transAxes.inverted()
            )
            bb_rel = master_legend_rel.get_bbox_to_anchor().transformed(
                axs_rel[study][key][idx1].transAxes.inverted()
            )

            # # Change to location of the legend.
            # yOffset = 0.025
            # bb.y0 += yOffset
            # bb.y1 += yOffset
            # bb_rel.y0 += yOffset
            # bb_rel.y1 += yOffset

            figs[study][key].tight_layout()
            figs_rel[study][key].tight_layout()

            # master_legend.set_bbox_to_anchor(
            #     bb, transform=axs[study][key][idx1].transAxes
            # )
            # master_legend_rel.set_bbox_to_anchor(
            #     bb_rel, transform=axs_rel[study][key][idx1].transAxes
            # )

            print("Saving figure...")
            print(f'{fig_save_path}/{cases[study][case][key]["save_name"]}')

            if global_save:
                figs[study][key].savefig(
                    f'{fig_save_path}/{cases[study][case][key]["save_name"]}',
                    dpi=300,
                    bbox_inches="tight",
                )
                figs_rel[study][key].savefig(
                    f'{fig_save_path}/{cases[study][case][key]["save_name"].replace(".png", "_rel.png")}',
                    dpi=300,
                    bbox_inches="tight",
                )


# global_save =  False
global_save = True


if __name__ == "__main__":
    # pprint.pprint(unit_color_dict_default)

    # plot_case(cases_kbhdp_soa["KBHDP"]["KBHDP_SOA_1"]["soda_ash"])
    # plot_case(cases_kbhdp_soa["KBHDP"]["KBHDP_SOA_1"]["water_recovery"])

    # fig_csv = pd.DataFrame()

    # plot_all_cases(kbhdp_wr, bboxy=1.03, legend_rows=1)

    # fig_csv.to_csv(f"{fig_save_path}/kbhdp_water_recovery_fig.csv", index=True)
    # plot_all_cases(kbhdp_grid_frac, bboxy=1.03, legend_rows=1)

    # plot_all_cases(permian_wr, legend_rows=2)
    # plot_all_cases(permian_grid_frac, bboxy=1.03)
    plot_all_cases(permian_opt, bboxy=1.1, xdim=5, ydim=4)

    # plot_case(permian_opt["Permian"]["Permian_ZLD2_FO_Cryst"]["grid_fraction_optimize"])
    # plot_case(kbhdp_opt["KBHDP"]["KBHDP_RPT_1"]["cost_per_watt_installed"])
    # plot_case(kbhdp_wr["KBHDP"]["KBHDP_RPT_3"]["water_recovery"])
    # plot_case(kbhdp_grid_frac["KBHDP"]["KBHDP_RPT_2"]["grid_fraction"])
    # plot_case(kbhdp_zld["KBHDP"]["KBHDP_ZLD"]["cost_per_aperture_area"])
    # plot_case(kbhdp_grid_frac["KBHDP"]["KBHDP_RPT_1"]["grid_fraction"])
    # plt.tight_layout()
    # print(figure_csv_rel.columns)
    # print(figure_csv_rel.heat)
    # print(figure_csv_rel.aluminum)
    # print(figure_csv_rel["FPC CAPEX LCOW rel"])
    # print(figure_csv_rel["FPC OPEX LCOW rel"])
    plt.show()

    # e = 150959.708 * pyunits.kilowatt
    # q = 5 * pyunits.Mgallons / pyunits.day
    # r = 0.5

    # sec = e / (q * r)
    # x = pyunits.convert(sec, to_units=pyunits.kilowatt * pyunits.hour / pyunits.meter**3)
    # print(x())

    # print(3.3 + 0.45 * x())
