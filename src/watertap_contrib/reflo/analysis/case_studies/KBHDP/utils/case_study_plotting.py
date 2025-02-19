import os
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

unit_list = [
    "Chem. Softening",
    "H2O2 Addition",
    "EC",
    "UF",
    "RO",
    "Pump",
    "MD",
    "LT-MED",
    "MEC",
    "DWI",
    "FO",
    "PV",
    "CST",
    "FPC",
    "CF",
    "MVC",
]

flow_list = [
    "electricity",
    "heat",
    "hydrogen_peroxide",
    "aluminum",
    "soda_ash",
    "lime", 
    "CO2",
    "mgcl2",
]

unit_colors = plt.cm.tab20(np.arange(len(unit_list)).astype(int))
unit_color_dict_default = dict(zip(unit_list, unit_colors))

flow_colors = plt.cm.Dark2(np.arange(len(flow_list)).astype(int))
flow_color_dict_default = dict(zip(flow_list, flow_colors))
flow_color_dict_default["electric"] = flow_color_dict_default["electricity"]

def case_study_stacked_plot(
    df,
    fig=None,
    ax=None,
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
        loc="best",
        frameon=False,
        ncol=4,
        handlelength=1,
        handleheight=1,
        labelspacing=0.2,
        columnspacing=0.9,
    ),
):

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
    actual_lcow = list()
    capex = list()
    opex = list()
    flow_opex = list()

    costing_params = dict()

    for gp in global_params:
        # make sure all the global parameters have the same value
        if not len(df[f"{global_costing_blk}.{gp}"].unique()) == 1:
            print(df[f"{global_costing_blk}.{gp}"].unique())
            raise ValueError(f"Global parameter {gp} does not have uniform value.")
        costing_params[gp] = df[f"{global_costing_blk}.{gp}"].iloc[0]

    # assert len(df[flow_col].unique()) == 1

    df.set_index(xcol, inplace=True)
    df.sort_values(by=[xcol], inplace=True)

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
                print(f"{b}.capital_cost")
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
            total_lcow += unit_capex_lcow  # $ / m3

            if unit_capex_lcow != 0:
                capex_lcow[u].append(unit_capex_lcow)

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
                print(f"No Variable OPEX for {u} found.")
                pass

            total_opex += unit_opex_total  # $ / year
            unit_opex_lcow = unit_opex_total / denominator  # $ / m3

            if unit_opex_lcow != 0:
                opex_lcow[u].append(unit_opex_lcow)

            total_lcow += unit_opex_lcow  # $ / m3
            total_annualized_cost += unit_opex_total  # $ / year
        
        for flow_label, flow_name in agg_flows.items():
            flow_lcow = 0
            try:
                # First try to find it via aggregate_flow_costs
                total_flow_cost += row.loc[
                    f"{costing_blk}.aggregate_flow_costs[{flow_name}]"
                ]  # $ / year
                total_annualized_cost += row.loc[
                    f"{costing_blk}.aggregate_flow_costs[{flow_name}]"
                ]
                flow_lcow += (
                    row.loc[f"{costing_blk}.aggregate_flow_costs[{flow_name}]"] / denominator
                )  # $ / m3
                agg_flow_lcow[flow_name].append(flow_lcow)
            except KeyError:
                # print(f"No aggregate cost for {flow_name} found.")
                # pass

                try:
                    total_flow_cost += row.loc[
                        f"{costing_blk}.total_{flow_name}_operating_cost"
                    ]  # $ / year
                    total_annualized_cost += row.loc[
                        f"{costing_blk}.total_{flow_name}_operating_cost"
                    ]  # $ / year
                    flow_lcow += (
                        row.loc[f"{costing_blk}.total_{flow_name}_operating_cost"] / denominator
                    )  # $ / m3
                    agg_flow_lcow[flow_name].append(flow_lcow)
                except KeyError:
                    # print(f"No aggregate cost for {flow_name} found.")
                    # pass
                    raise ValueError(f"No cost for {flow_name} found.")

            total_lcow += flow_lcow

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

    for flow_label, flow_name in agg_flows.items():
        stacked_cols.append(agg_flow_lcow[flow_name])
        stacked_labels.append(flow_label)
        stacked_colors.append(flow_color_dict[flow_name])
        stacked_hatch.append(flow_hatch)
        legend_elements.append(
            Patch(
                facecolor=flow_color_dict[flow_name],
                label=flow_label,
                hatch=flow_hatch,
                edgecolor="k",
            )
        )

    for u in unit_dict.keys():
        if u in capex_lcow.keys():
            stacked_cols.append(capex_lcow[u])
            stacked_labels.append(f"{u} CAPEX")
            stacked_hatch.append(capex_hatch)
            stacked_colors.append(unit_color_dict[u])
        if u in opex_lcow.keys():
            stacked_cols.append(opex_lcow[u])
            stacked_labels.append(f"{u} OPEX")
            stacked_hatch.append(opex_hatch)
            stacked_colors.append(unit_color_dict[u])
        legend_elements.append(Patch(facecolor=unit_color_dict[u], label=u, edgecolor="k"))

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(figsize=figsize)
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
        ax.legend(handles=legend_elements, **leg_kwargs)

    ax.set(**ax_dict)
    ax.set_xlabel(ax_dict["xlabel"], fontsize=label_fontsize)
    ax.set_ylabel(ax_dict["ylabel"], fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.set_xlim(df.index.min(), df.index.max())
    # ax.set_ylim(0, np.ceil(max(actual_lcow)))

    plt.tight_layout()

    return fig, ax


if __name__ == "__main__":

    test_file = f"{__location__}/test_stacked_plot.csv"
    df = pd.read_csv(test_file)
    df = df[df["fs.costing.electricity_cost"] == 0.049886703].copy()

    xcol = "fs.costing.electricity_cost"
    xcol = "fs.water_recovery"

    flow_col = "fs.product.properties[0.0].flow_vol_phase[Liq]"
    # df = df[df["fs.water_recovery"] == 0.8].copy()

    unit_dict = {
        "UF": "fs.UF.unit.costing",
        "EC": "fs.EC.ec.costing",
        "Pump": "fs.pump.costing",
        "RO": "fs.RO.stage[1].module.costing",
    }
    agg_flows = [
        "aluminum",
        "electricity",
        # "heat",
    ]

    agg_flows = {
        "Aluminum": "aluminum",
        "Electricity": "electricity"
    }

    ax_dict = dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)")

    fig, ax = case_study_stacked_plot(
        df,
        global_costing_blk=None,
        unit_dict=unit_dict,
        agg_flows=agg_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
    )

    # plt.show()
    f = "/Users/ksitterl/Downloads/water_recovery.csv"
    df = pd.read_csv(f)

    unit_dict = {
        "FPC": "fs.energy.FPC.costing",
        "MD": "fs.treatment.md.unit",
        "DWI": "fs.treatment.dwi.unit.costing",
    }
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    
    agg_flows = {
        "Heat": "heat",
        "Electricity": "electric"
    }

    fig, ax = case_study_stacked_plot(
        df,
        # global_costing_blk=None,
        unit_dict=unit_dict,
        agg_flows=agg_flows,
        xcol=xcol,
        flow_col=flow_col,
        ax_dict=ax_dict,
        opex_hatch="\\\\\\",
        flow_hatch="..",
    )

    plt.show()



