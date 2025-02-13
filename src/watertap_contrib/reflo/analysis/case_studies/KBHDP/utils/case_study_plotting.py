import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


color_dict = {
    "Oxidation": '#1f77b4',
    "Chem. Soft": '#ff7f0e',
    "EC": '#2ca02c',
    "UF": '#d62728',
    "RO": '#9467bd',
    "Pump": "#e377c2",
    "MD": '#8c564b',
    "LT-MED": 'lightsteelblue',
    "MEC": "cyan",
    "DWI": "tan",
    "FO": "teal", 
    "PV": "yellow",
    "CST": "xkcd:melon",
    "FPC": "aquamarine",
    "aluminum": 'lightgray',
    "electricity": "deepskyblue",
    "heat": "deeppink",
    "NaCl_recovered": "lawngreen"
}


def case_study_stacked_plot(
    df,
    fig=None,
    ax=None,
    costing_blk="fs.costing",
    xcol=None,
    flow_col=None,  # column to be used as denominator in LCOW calculations, assumed to be in m3/s
    unit_dict=dict(),  # (unit name: unit location)
    agg_flows=list(), # list of aggregated flows
    figsize=(6, 4),
    color_dict=None,
    capex_hatch="",
    opex_hatch="\\\\\\",
    flow_hatch="..",
    ax_dict=dict(),
    check_calc=True,
    label_fontsize=16,
    tick_fontsize=14,
    add_legend=True,
    leg_kwargs=dict(bbox_to_anchor=(1, 1), loc="upper left"),
):

    if flow_hatch is None:
        flow_hatch = capex_hatch

    if color_dict is None:
        raise ValueError("Must provide color_dict")

    global_params = [
        "maintenance_labor_chemical_factor",
        "utilization_factor",
        "capital_recovery_factor",
        "total_investment_factor",
    ]

    opex_lcow = defaultdict(list)
    capex_lcow = defaultdict(list)
    agg_flow_lcow = defaultdict(list)
    actual_lcow = list()
    capex = list()
    opex = list()
    flow_opex = list()

    costing_params = dict()

    for gp in global_params:
        # make sure all the global parameters have the same value
        if not len(df[f"{costing_blk}.{gp}"].unique()) == 1:
            print(df[f"{costing_blk}.{gp}"].unique())
            raise ValueError(f"Global parameter {gp} does not have uniform value.")
        costing_params[gp] = df[f"{costing_blk}.{gp}"].iloc[0]

    # assert len(df[flow_col].unique()) == 1

    df.set_index(xcol, inplace=True)

    for x, row in df.iterrows():

        total_lcow = 0
        total_capex = 0
        total_opex = 0
        total_flow_cost = 0
        total_annualized_cost = 0
        row_lcow = row[f"{costing_blk}.LCOW"]
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

            ## CAPEX
            unit_capex = 0

            try:
                unit_capex += (
                    row.loc[f"{b}.costing.capital_cost"]
                    * costing_params["total_investment_factor"]
                )  # USD2023

            except KeyError:
                # print(f"No CAPEX for {u} found.")
                pass

            total_capex += unit_capex  # $
            total_annualized_cost += (
                unit_capex * costing_params["capital_recovery_factor"]
            )  # $ / year

            unit_capex_lcow = (
                unit_capex * costing_params["capital_recovery_factor"]
            ) / denominator  # $ / m3
            total_lcow += unit_capex_lcow  # $ / m3

            if unit_capex_lcow != 0:
                capex_lcow[u].append(unit_capex_lcow)

            ### OPEX
            unit_opex_total = 0
            unit_opex_total += (
                unit_capex * costing_params["maintenance_labor_chemical_factor"]
            )  # $ / year

            try:
                unit_opex_total += row.loc[f"{b}.costing.fixed_operating_cost"]
            except KeyError:
                # print(f"No Fixed OPEX for {u} found.")
                pass

            try:
                unit_opex_total += row.loc[f"{b}.costing.variable_operating_cost"]
            except KeyError:
                # print(f"No Variable OPEX for {u} found.")
                pass

            total_opex += unit_opex_total  # $ / year
            unit_opex_lcow = unit_opex_total / denominator  # $ / m3

            if unit_opex_lcow != 0:
                opex_lcow[u].append(unit_opex_lcow)

            total_lcow += unit_opex_lcow  # $ / m3
            total_annualized_cost += unit_opex_total  # $ / year

        for flow in agg_flows:
            flow_lcow = 0
            try:
                total_flow_cost += row.loc[
                    f"{costing_blk}.aggregate_flow_costs[{flow}]"
                ]  # $ / year
                flow_lcow += (
                    row.loc[f"{costing_blk}.aggregate_flow_costs[{flow}]"] / denominator
                )  # $ / m3
                agg_flow_lcow[flow].append(flow_lcow)
            except KeyError:
                # print(f"No aggregate cost for {flow} found.")
                pass

            total_lcow += flow_lcow
            total_annualized_cost += row.loc[
                f"{costing_blk}.aggregate_flow_costs[{flow}]"
            ]
        if check_calc:
            print(f"\nFor {xcol} = {x}:")
            print(f"\tActual LCOW: {row_lcow:.6f}")
            print(f"\tCalculated LCOW: {total_lcow:.6f}")

        capex.append(total_capex)
        opex.append(total_opex)
        flow_opex.append(total_flow_cost)

    stacked_cols = list()
    stacked_labels = list()
    stacked_hatch = list()
    stacked_colors = list()

    legend_elements = [
        Patch(facecolor="white", hatch=capex_hatch, label="CAPEX", edgecolor="k"),
        Patch(facecolor="white", hatch=opex_hatch, label="OPEX", edgecolor="k"),
    ]

    for flow in agg_flows:
        stacked_cols.append(agg_flow_lcow[flow])
        stacked_labels.append(flow.replace("_", " ").title())
        stacked_colors.append(color_dict[flow])
        stacked_hatch.append(flow_hatch)
        legend_elements.append(
            Patch(
                facecolor=color_dict[flow],
                label=flow.replace("_", " ").title(),
                hatch=flow_hatch,
                edgecolor="k",
            )
        )

    for u in unit_dict.keys():
        if u in capex_lcow.keys():
            stacked_cols.append(capex_lcow[u])
            stacked_labels.append(f"{u} CAPEX")
            stacked_hatch.append(capex_hatch)
            stacked_colors.append(color_dict[u])
        if u in opex_lcow.keys():
            stacked_cols.append(opex_lcow[u])
            stacked_labels.append(f"{u} OPEX")
            stacked_hatch.append(opex_hatch)
            stacked_colors.append(color_dict[u])
        legend_elements.append(Patch(facecolor=color_dict[u], label=u, edgecolor="k"))

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(figsize=figsize)

    ax.stackplot(
        df.index,
        stacked_cols,
        labels=stacked_labels,
        hatch=stacked_hatch,
        colors=stacked_colors,
        edgecolor="black",
    )

    if add_legend:
        ax.legend(handles=legend_elements, **leg_kwargs)

    ax.set(**ax_dict)
    ax.set_xlabel(ax_dict["xlabel"], fontsize=label_fontsize)
    ax.set_ylabel(ax_dict["ylabel"], fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)

    plt.tight_layout()

    return fig, ax


if __name__ == "__main__":
        
    f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/KBHDP/ZLD/test_stacked_plot.csv"
    df = pd.read_csv(f)
    df = df[df["fs.costing.electricity_cost"] == 0.049886703].copy()
    
    xcol = "fs.costing.electricity_cost"
    xcol = "fs.water_recovery"

    flow_col = "fs.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "UF": "fs.UF.unit",
        "EC": "fs.EC.ec",
        "Pump": "fs.pump",
        "RO": "fs.RO.stage[1].module",
    }
    agg_flows = [
        "aluminum",
        "electricity",
        # "heat",
    ]
    
    ax_dict = dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        unit_dict=unit_dict,
        agg_flows=agg_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
        color_dict=color_dict,
    )

    plt.show()