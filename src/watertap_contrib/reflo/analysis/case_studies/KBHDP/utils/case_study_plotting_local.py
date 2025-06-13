import os
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import OrderedDict
import os
import math

# from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.h5_to_csv import (
#     convert_h5_to_csv,
# )

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.dirname(__location__)
sweep_csv_dir = os.path.join(parent_dir, "sweep_csvs")

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
    "Electricity",
    "Heat",
    "hydrogen_peroxide",
    "Aluminum",
    "Soda Ash",
    "Lime",
    "CO2",
    "MgCl2",
    "H2O2",
]

unit_colors = plt.cm.tab20(np.arange(len(unit_list)).astype(int))
unit_color_dict_default = dict(zip(unit_list, unit_colors))

flow_colors = plt.cm.Dark2(np.arange(len(flow_list)).astype(int))
flow_color_dict_default = dict(zip(flow_list, flow_colors))
flow_color_dict_default["electric"] = flow_color_dict_default["Electricity"]


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
        ncol=3,
        handlelength=1,
        handleheight=1,
        labelspacing=0.2,
        columnspacing=0.9,
    ),
    ylims=None,
    xlims=None,
    save=False,
    save_name="temp",
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
            print(f"{global_costing_blk}.{gp}")
            print(df[f"{global_costing_blk}.{gp}"].unique())
            print(df[f"{global_costing_blk}.{gp}"])
            # raise ValueError(f"Global parameter {gp} does not have uniform value.")
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
            # print(f"Calculating LCOW for {u} at {xcol} = {x}")
            ## CAPEX
            print(f"Getting CAPEX and OPEX for {u} in {b}")
            unit_capex = 0

            try:
                # print(f"{b}.capital_cost")
                unit_capex += (
                    row.loc[f"{b}.capital_cost"]
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

            if unit_capex_lcow > 1e-13:
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
                print(f"No Variable OPEX for {u} found at {b}.variable_operating_cost.")
                pass

            # print(unit_opex_total)
            total_opex += unit_opex_total  # $ / year
            unit_opex_lcow = unit_opex_total / denominator  # $ / m3

            if unit_opex_lcow > 1e-12:
                opex_lcow[u].append(unit_opex_lcow)

            total_lcow += unit_opex_lcow  # $ / m3
            total_annualized_cost += unit_opex_total  # $ / year

        print(total_lcow)
        # assert False
        print("\n\n")
        for flow_label, flow_name in agg_flows.items():
            print(f"Calculating LCOW for {flow_label} at {xcol} = {x}")
            flow_lcow = 0
            if flow_name == "electricity":
                try:
                    print(
                        "Looking for Electricity OPEX in fs.costing.total_electric_operating_cost"
                    )
                    print("Found:")
                    flow_lcow += (
                        row.loc[f"fs.costing.total_electric_operating_cost"]
                        / denominator
                    )
                    agg_flow_lcow[flow_name].append(flow_lcow)
                    continue
                except:
                    print(
                        f"No Electricity OPEX for {flow_name} found in in fs.costing.total_electric_operating_cost."
                    )
                    pass
            if flow_name == "heat":
                try:
                    print(
                        "Looking for Heat OPEX in fs.costing.total_heat_operating_cost"
                    )
                    print("Found:")
                    flow_lcow += (
                        row.loc[f"fs.costing.total_heat_operating_cost"] / denominator
                    )
                    agg_flow_lcow[flow_name].append(flow_lcow)
                    continue
                except:
                    print(
                        f"No Electricity OPEX for {flow_name} found in in fs.costing.total_electric_operating_cost."
                    )
                    pass

            print("trying to find electricity in aggregate_flow_costs")
            print(
                f"Looking for Electricity OPEX in {global_costing_blk}.aggregate_flow_costs[{flow_name}]"
            )
            print(f"Flow Label {flow_label}, Flow Name {flow_name}")
            print(f"Searching in {global_costing_blk}")
            total_flow_cost = row.loc[
                f"{global_costing_blk}.aggregate_flow_costs[{flow_name}]"
            ]  # $ / year
            total_annualized_cost = row.loc[
                f"{global_costing_blk}.aggregate_flow_costs[{flow_name}]"
            ]
            flow_lcow = total_flow_cost / denominator  # $ / m3
            print(row["fs.treatment.costing.aggregate_flow_costs[electricity]"])
            "fs.treatment.costing.aggregate_flow_costs[electricity]"
            "fs.treatment.costing.aggregate_flow_costs[electricity]"
            print(total_flow_cost)
            print(total_annualized_cost)
            print(
                f"Flow Name {flow_name} Flow LCOW {flow_lcow} Flow Cost {total_flow_cost} Denominator {denominator}"
            )

            agg_flow_lcow[flow_name].append(flow_lcow)
            # print(f"Flow Name {flow_name} Flow LCOW: {agg_flow_lcow[flow_name]}")

            # raise ValueError(f"No cost for {flow_name} found.")
            # print(f"Flow OPEX: {total_flow_cost}\n")
            # print(f"Flow LCOW: {flow_lcow}\n")
            total_lcow = flow_lcow

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
        stacked_colors.append(flow_color_dict[flow_label])
        stacked_hatch.append(flow_hatch)
        legend_elements.append(
            Patch(
                facecolor=flow_color_dict[flow_label],
                label=flow_label,
                hatch=flow_hatch,
                edgecolor="k",
            )
        )

    for u in unit_dict.keys():
        print(u)
        if u in capex_lcow.keys():
            stacked_cols.append(capex_lcow[u])
            stacked_labels.append(f"{u} CAPEX")
            stacked_hatch.append(capex_hatch)
            stacked_colors.append(unit_color_dict[u])
            print(capex_lcow[u])
        if u in opex_lcow.keys():
            stacked_cols.append(opex_lcow[u])
            stacked_labels.append(f"{u} OPEX")
            stacked_hatch.append(opex_hatch)
            stacked_colors.append(unit_color_dict[u])
            print(opex_lcow[u])
        legend_elements.append(
            Patch(facecolor=unit_color_dict[u], label=u, edgecolor="k")
        )

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(figsize=figsize)
        fig.set_size_inches(5, 5, forward=True)

    # print(stacked_cols)
    for item in stacked_cols:
        print(len(item), item)
    print("\n\n")
    ax.stackplot(
        df.index,
        stacked_cols,
        # labels=stacked_labels,
        hatch=stacked_hatch,
        colors=stacked_colors,
        edgecolor="black",
    )

    if add_legend:
        legend = ax.legend(handles=legend_elements, **leg_kwargs)

    ax.set(**ax_dict)
    ax.set_xlabel(ax_dict["xlabel"], fontsize=label_fontsize)
    ax.set_ylabel(ax_dict["ylabel"], fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.set_xlim(df.index.min(), df.index.max())
    # ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"${x:.2f}"))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}"))
    if ylims is not None:
        ax.set_ylim(ylims)
    else:
        ax.set_ylim(0, np.ceil(max(actual_lcow)))
    if xlims is not None:
        ax.set_xlim(xlims)
    else:
        ax.set_xlim(0, 1)

    plt.tight_layout()

    if save is True:
        plt.savefig(
            os.path.join(
                "/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/figures",
                save_name,
            ),
            dpi=300,
        )

    unit_lcow = dict()
    print(capex_lcow)
    print(opex_lcow)
    for unit, key in unit_dict.items():
        print(f"\nUnit: {unit}")
        print((np.array(capex_lcow[unit])))
        print((np.array(opex_lcow[unit])))
        if len(np.array(capex_lcow[unit])) > 1:
            unit_lcow[unit] = np.array(capex_lcow[unit]) + np.array(opex_lcow[unit])
        else:
            unit_lcow[unit] = np.array(opex_lcow[unit])

    print(f"Agg Flow Dict {agg_flows}")
    print(f"Agg Flow LCOW Dict {agg_flow_lcow}")
    for flow, key in agg_flows.items():
        print(f"\nFlow: {key.lower()}")
        print(f"Len {len(agg_flow_lcow[key.lower()])}")
        unit_lcow[key] = np.array(agg_flow_lcow[key.lower()])
        print(unit_lcow[key])

    for unit in unit_lcow.keys():
        print(unit, len(unit_lcow[unit]))
    figure_csv = pd.DataFrame.from_dict(unit_lcow)
    figure_csv["LCOW"] = figure_csv.sum(axis=1)
    figure_csv.index = df.index
    figure_csv.to_csv(
        os.path.join(
            "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/figures",
            save_name + "_data.csv",
        )
    )

    print(figure_csv.head(20))

    return fig, ax, legend


cases = {
    # "KBHDP": {
    #     "KBHDP_SOA_1": {
    #         # "water_recovery": {
    #         #     "file": os.path.join(
    #         #         # sweep_csv_dir, "sweep_data_analysisType_KBHDP_SOA_1_water_recovery.csv"
    #         #         sweep_csv_dir, "sweep_data_KBHDP_SOA_1_water_recovery.csv"
    #         #     ),
    #         #     "global_costing_blk": "fs.treatment.costing",
    #         #     "costing_blk": "fs.treatment.costing",
    #         #     "unit_dict": {
    #         #         "Chem. Softening": "fs.treatment.softener.unit.costing",
    #         #         "UF": "fs.treatment.UF.unit.costing",
    #         #         "RO": "fs.treatment.RO.stage[1].module.costing",
    #         #         "DWI": "fs.treatment.DWI.unit.costing",
    #         #         "Pump": "fs.treatment.pump.costing",
    #         #     },
    #         #     "agg_flows": {
    #         #         "Electricity": "electricity",
    #         #         "Soda Ash": "soda_ash",
    #         #         "Lime": "lime",
    #         #         "CO2": "co2",
    #         #         "MgCl2": "mgcl2",
    #         #     },
    #         #     "xcol": "fs.water_recovery",
    #         #     "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #         #     "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
    #         #     "ylims": (0, 6),
    #         #     "xlims": (0.3, 0.8),
    #         #     "save_name": "KBHDP_SOA_1_stacked_plot_water_recovery_sweep.png",
    #         #     "save": True,
    #         # },
    #         "water_recovery": {
    #             "file": os.path.join(
    #                 # sweep_csv_dir, "sweep_data_analysisType_KBHDP_SOA_1_water_recovery.csv"
    #                 sweep_csv_dir,
    #                 "sweep_data_KBHDP_SOA_1_soda_ash_cost.csv",
    #             ),
    #             "global_costing_blk": "fs.treatment.costing",
    #             "costing_blk": "fs.treatment.costing",
    #             "unit_dict": {
    #                 "Chem. Softening": "fs.treatment.softener.unit.costing",
    #                 "UF": "fs.treatment.UF.unit.costing",
    #                 "RO": "fs.treatment.RO.stage[1].module.costing",
    #                 "DWI": "fs.treatment.DWI.unit.costing",
    #                 "Pump": "fs.treatment.pump.costing",
    #             },
    #             "agg_flows": {
    #                 "Electricity": "electricity",
    #                 "Soda Ash": "soda_ash",
    #                 "Lime": "lime",
    #                 "CO2": "co2",
    #                 "MgCl2": "mgcl2",
    #             },
    #             "xcol": "fs.treatment.costing.soda_ash.cost",
    #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #             "ax_dict": dict(
    #                 xlabel="Soda Ash Cost ($/kg)", ylabel="LCOW (\$/m$^3$)"
    #             ),
    #             "ylims": (0, 3),
    #             "xlims": (0.12, 0.36),
    #             "save_name": "KBHDP_SOA_1_stacked_plot_soda_ash_sweep.png",
    #             "save": True,
    #         }
    #     },
    #     #     "KBHDP_RPT_1": {
    #     #         "water_recovery": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir, "sweep_data_KBHDP_RPT_1_water_recovery.csv"
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "EC": "fs.treatment.EC.ec.costing",
    #     #                 "UF": "fs.treatment.UF.unit.costing",
    #     #                 "RO": "fs.treatment.RO.stage[1].module.costing",
    #     #                 "DWI": "fs.treatment.DWI.unit.costing",
    #     #                 "Pump": "fs.treatment.pump.costing",
    #     #                 # "PV": "fs.energy.pv.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Aluminum": "aluminum",
    #     #             },
    #     #             "xcol": "fs.water_recovery",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
    #     #             "ylims": (0, 5),
    #     #             "xlims": (0.3, 0.8),
    #     #             "save_name": "KBHDP_RPT_1_stacked_plot_water_recovery_sweep.png",
    #     #             "save": False,
    #     #         },
    #     #         "grid_fraction": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir, "sweep_data_KBHDP_RPT_1_grid_fraction.csv"
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "EC": "fs.treatment.EC.ec.costing",
    #     #                 "UF": "fs.treatment.UF.unit.costing",
    #     #                 "RO": "fs.treatment.RO.stage[1].module.costing",
    #     #                 "DWI": "fs.treatment.DWI.unit.costing",
    #     #                 "Pump": "fs.treatment.pump.costing",
    #     #                 "PV": "fs.energy.pv.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Aluminum": "aluminum",
    #     #             },
    #     #             "xcol": "fs.costing.RE Fraction",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(
    #     #                 xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
    #     #             ),
    #     #             "ylims": (0, 2),
    #     #             "xlims": (0.1, 0.8),
    #     #             "save_name": "KBHDP_RPT_1_stacked_plot_grid_fraction_sweep.png",
    #     #             "save": True,
    #     #         },
    #     #     },
    #     #     "KBHDP_RPT_2": {
    #     #         "water_recovery": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir, "sweep_data_KBHDP_RPT_2_water_recovery.csv"
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "EC": "fs.treatment.EC.ec.costing",
    #     #                 "UF": "fs.treatment.UF.unit.costing",
    #     #                 "LT-MED": "fs.treatment.LTMED.unit.costing",
    #     #                 "DWI": "fs.treatment.DWI.unit.costing",
    #     #                 # "Pump": "fs.treatment.pump.costing",
    #     #                 # "FPC": "fs.energy.FPC.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Aluminum": "aluminum",
    #     #                 "Heat": "heat",
    #     #             },
    #     #             "xcol": "fs.water_recovery",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
    #     #             "ylims": (0, 5),
    #     #             "xlims": (0.3, 0.8),
    #     #             "save_name": "KBHDP_RPT_2_stacked_plot_water_recovery_sweep.png",
    #     #             "save": True,
    #     #         },
    #     #         "grid_fraction": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir, "sweep_data_KBHDP_RPT_2_frac_heat_from_grid.csv"
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "EC": "fs.treatment.EC.ec.costing",
    #     #                 "UF": "fs.treatment.UF.unit.costing",
    #     #                 "LT-MED": "fs.treatment.LTMED.unit.costing",
    #     #                 "DWI": "fs.treatment.DWI.unit.costing",
    #     #                 # "Pump": "fs.treatment.pump.costing",
    #     #                 "FPC": "fs.energy.FPC.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Aluminum": "aluminum",
    #     #                 "Heat": "heat",
    #     #             },
    #     #             "xcol": "fs.costing.RE Fraction",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(
    #     #                 xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
    #     #             ),
    #     #             "xlims": (0.1, 0.8),
    #     #             "ylims": (0, 20),
    #     #             "save_name": "KBHDP_RPT_2_stacked_plot_grid_fraction_sweep.png",
    #     #             "save": True,
    #     #         },
    #     #     },
    #     #     "KBHDP_RPT_3": {
    #     #         "water_recovery": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir,
    #     #                 "water_recovery_grid_frac_0.5_recovery_var.csv"
    #     #                 # "kbhdp_RPT3_water_recovery_grid_frac_0.5_recovery_varcheck.csv",
    #     #                 # "kbhdp_RPT3_water_recovery_grid_frac_1_recovery_var.csv",
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "MD": "fs.treatment.md.unit",
    #     #                 "DWI": "fs.treatment.dwi.unit.costing",
    #     #                 # "Pump": "fs.treatment.pump.costing",
    #     #                 # "FPC": "fs.energy.FPC.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Heat": "heat",
    #     #             },
    #     #             "xcol": "fs.water_recovery",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
    #     #             "ylims": (0, 5),
    #     #             "xlims": (0.3, 0.8),
    #     #             "save_name": "KBHDP_RPT_3_stacked_plot_water_recovery_sweep.png",
    #     #             "save": False,
    #     #         },
    #     #         "grid_fraction": {
    #     #             "file": os.path.join(
    #     #                 sweep_csv_dir,
    #     #                 # "kbhdp_RPT3_grid_frac_heat_grid_frac_var_recovery_0.8 1.csv",
    #     #                 "grid_frac_heat_grid_frac_var_recovery_0.7.csv",
    #     #             ),
    #     #             "global_costing_blk": "fs.treatment.costing",
    #     #             "costing_blk": "fs.treatment.costing",
    #     #             "unit_dict": {
    #     #                 "MD": "fs.treatment.md.unit",
    #     #                 "DWI": "fs.treatment.dwi.unit.costing",
    #     #                 # "Pump": "fs.treatment.pump.costing",
    #     #                 "FPC": "fs.energy.FPC.costing",
    #     #             },
    #     #             "agg_flows": {
    #     #                 "Electricity": "electricity",
    #     #                 "Heat": "heat",
    #     #             },
    #     #             "xcol": "fs.costing.RE Fraction",
    #     #             "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
    #     #             "ax_dict": dict(
    #     #                 xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
    #     #             ),
    #     #             "xlims": (0.1, 0.8),
    #     #             "ylims": (0, 20),
    #     #             "save_name": "KBHDP_RPT_3_stacked_plot_grid_fraction_sweep.png",
    #     #             "save": True,
    #     #         },
    #     #     },
    # },
    "Permian": {
        "Permian_RPT1_MD": {
            # "water_recovery": {
            #     # "file": os.path.join(
            #     #     sweep_csv_dir,
            #     #     "permian_RPT1_MD_water_recovery_grid_frac_1_recovery_var.csv",
            #     # ),
            #     "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT1_MD_water_recovery_grid_frac_1_recovery_var.csv",
            #     "global_costing_blk": "fs.treatment.costing",
            #     "costing_blk": "fs.treatment.costing",
            #     "unit_dict": {
            #         "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
            #         "EC": "fs.treatment.EC.unit.costing",
            #         "CF": "fs.treatment.cart_filt.unit.costing",
            #         "MD": "fs.treatment.md.unit",
            #         "DWI": "fs.treatment.DWI.unit.costing",
            #     },
            #     "agg_flows": {
            #         "Electricity": "electricity",
            #         "Heat": "heat",
            #         "H2O2": "hydrogen_peroxide",
            #         "Aluminum": "aluminum",
            #     },
            #     "xcol": "fs.water_recovery",
            #     "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
            #     "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
            #     "ylims": (0, 35),
            #     "xlims": (0.35, 0.5),
            #     "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",
            #     "save": False,
            # },
            "grid_fraction": {
                # "file": os.path.join(
                #     sweep_csv_dir,
                #     "permian_RPT1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv", # Change this
                # ),
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.EC.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "MD": "fs.treatment.md.unit",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H2O2": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 40),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",  # Change this
                "save": False,
            },
        },
        "Permian_RPT1_FO": {
            # "water_recovery": {
            #     # "file": os.path.join(
            #     #     # sweep_csv_dir, "FO_DWI_fo_recovery_ratio_sweep.csv"
            #     #     sweep_csv_dir,
            #     #     "permian_RPT2_FO_DWI_no_CST_sweep_recovery_ratio.csv",
            #     # ),
            #     "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT2_FO_DWI_no_CST_sweep_recovery_ratio.csv",
            #     "global_costing_blk": "fs.treatment.costing",
            #     "costing_blk": "fs.treatment.costing",
            #     "unit_dict": {
            #         "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
            #         "EC": "fs.treatment.ec.unit.costing",
            #         "CF": "fs.treatment.cart_filt.unit.costing",
            #         "FO": "fs.treatment.FO.fs.fo.costing",
            #         "DWI": "fs.treatment.DWI.unit.costing",
            #         # "CST": "fs.energy.cst.unit.costing",
            #     },
            #     "agg_flows": {
            #         "Electricity": "electricity",
            #         "Heat": "heat",
            #         "H2O2": "hydrogen_peroxide",
            #         "Aluminum": "aluminum",
            #     },
            #     "xcol": "fs.treatment.FO.fs.fo.recovery_ratio",
            #     "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
            #     "ax_dict": dict(xlabel="Water Recovery (%)", ylabel="LCOW (\$/m$^3$)"),
            #     "ylims": (0, 35),
            #     "xlims": (0.35, 0.5),
            #     "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",
            #     "save": False,
            # },
            "grid_fraction": {
                # "file": os.path.join(
                #     sweep_csv_dir,
                #     "permian_RPT2_FO_DWI_RPT_grid_frac.csv", # Change this
                # ),
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_RPT2_FO_DWI_RPT_grid_frac.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "DWI": "fs.treatment.DWI.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H2O2": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 40),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",  # Change this
                "save": False,
            },
        },
        "Permian_RPT2_MD_Cryst": {
            "grid_fraction": {
                # "file": os.path.join(
                #     sweep_csv_dir,
                #     "permian_ZLD1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv", # Change this
                # ),
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_ZLD1_MD_grid_frac_heat_grid_frac_var_recovery_0.5.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.EC.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "MD": "fs.treatment.md.unit",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H2O2": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 40),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",  # Change this
                "save": False,
            },
        },
        "Permian_RPT2_FO_Cryst": {
            "grid_fraction": {
                # "file": os.path.join(
                #     sweep_csv_dir,
                #     "permian_ZLD2_FO_cryst_RPT_grid_frac.csv",  # Change this
                # ),
                "file": "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/finalized results/permian/permian_ZLD2_FO_cryst_RPT_grid_frac.csv",
                "global_costing_blk": "fs.treatment.costing",
                "costing_blk": "fs.treatment.costing",
                "unit_dict": {
                    "H2O2 Addition": "fs.treatment.chem_addition.unit.costing",
                    "EC": "fs.treatment.ec.unit.costing",
                    "CF": "fs.treatment.cart_filt.unit.costing",
                    "FO": "fs.treatment.FO.fs.fo.costing",
                    "MEC": "fs.treatment.mec.unit.costing",
                    "CST": "fs.energy.cst.unit.costing",
                },
                "agg_flows": {
                    "Electricity": "electricity",
                    "Heat": "heat",
                    "H2O2": "hydrogen_peroxide",
                    "Aluminum": "aluminum",
                },
                "xcol": "fs.costing.RE Fraction",  # Change this
                "flow_col": "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]",
                "ax_dict": dict(
                    xlabel="Renewable Energy (%)", ylabel="LCOW (\$/m$^3$)"
                ),
                "ylims": (0, 40),
                "xlims": (0.1, 0.5),  # Change this
                "save_name": "Permian_ST1_stacked_plot_water_recovery_sweep.png",  # Change this
                "save": False,
            },
        },
    },
}


def plot_case(case, fig=None, ax=None):
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

    fig_, ax_, legend_ = case_study_stacked_plot(
        df,
        fig=fig,
        ax=ax,
        global_costing_blk=case["global_costing_blk"],
        costing_blk=case["costing_blk"],
        unit_dict=case["unit_dict"],
        agg_flows=case["agg_flows"],
        xcol=case["xcol"],
        flow_col=case["flow_col"],
        ax_dict=case["ax_dict"],
        opex_hatch="\\\\\\",
        flow_hatch="..",
        ylims=case["ylims"],
        xlims=case["xlims"],
        save=case["save"],
        save_name=case["save_name"],
    )

    return legend_


def plot_all_cases():
    figs = {}
    axs = {}
    legends = {}
    for idx0, study in enumerate(cases):
        figs[study] = {}
        axs[study] = {}
        legends[study] = {}
        for key in cases[study][list(cases[study].keys())[0]].keys():
            n_cases = len(cases[study])
            n_rows = int(n_cases // 4 + 1)
            n_cols = int(n_cases / n_rows)

            print(n_cols, n_rows)
            print(f"{study} {key}")
            figs[study][key], axs[study][key] = plt.subplots(
                n_rows,
                n_cols,
                figsize=(4.0 * n_cols, 4 * n_rows),
                # 1, len(cases[study]), figsize=(4.0 * len(cases[study]), 4)
            )
            legends[study][key] = []
            axs[study][key] = np.atleast_1d(axs[study][key]).ravel()

        for axis_id, axis in enumerate(axs[study][key]):
            print(f"Axis {axis_id}:{axis}")

    # assert False
    fig_letters = ["A", "B", "C", "D"]
    for idx0, study in enumerate(cases):
        for key in cases[study][list(cases[study].keys())[0]].keys():
            print(f"{study} {key}")
            for idx1, case in enumerate(cases[study]):
                print("CASE:", case)
                legend_ = plot_case(
                    cases[study][case][key],
                    fig=figs[study][key],
                    ax=axs[study][key][idx1],
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

            figs[study][key].tight_layout()
            if len(axs[study][key]) < 4:
                figs[study][key].subplots_adjust(top=0.85)
            else:
                figs[study][key].subplots_adjust(top=0.925)

            # Collect handles and labels
            handles, labels = [], []
            for legend in legends[study][key]:
                # for legend in [legends[1], legends[0], legends[2]]:
                h = legend.legend_handles
                l = [text.get_text() for text in legend.get_texts()]
                handles.extend(h)
                labels.extend(l)
            # Remove duplicates while preserving order
            unique = list(OrderedDict(zip(labels, handles)).items())
            # Add the unified legend
            legend_rows = 2
            if len(axs[study][key]) > 3:
                legend_rows = 2
            legend_cols = len(unique) // legend_rows + 1
            print(legend_cols)
            # assert False
            master_legend = figs[study][key].legend(
                [h for _, h in unique],
                [l for l, _ in unique],
                loc="upper center",  # or 'lower center' if you prefer
                ncol=legend_cols,
                # bbox_to_anchor=(0.1, 1.0, 0.8, 0.03),  # above the plots
                frameon=True,
                handlelength=1.4,
                handleheight=1.4,
                labelspacing=0.25,
                handletextpad=0.4,
                # columnspacing=0.9,
                fontsize=10,
                mode="expand",
            )
            master_legend.get_frame().set_edgecolor("black")  # Frame line color
            master_legend.get_frame().set_linewidth(1.5)  # Thickness of the frame

            for axis in axs[study][key]:
                axis.legend_.remove()  # This works if the legend exists

            # Get the bounding box of the original legend
            bb = master_legend.get_bbox_to_anchor().transformed(
                axs[study][key][idx1].transAxes.inverted()
            )

            # Change to location of the legend.
            yOffset = 0.025
            bb.y0 += yOffset
            bb.y1 += yOffset
            master_legend.set_bbox_to_anchor(
                bb, transform=axs[study][key][idx1].transAxes
            )

            # if len(axs[study][key]) > 3:
            #     print("Adjusting layout for 2x2 plot...")
            #     figs[study][key].tight_layout()

            # figs[study][key].subplots_adjust(wspace=0.3, bottom=0.15, top=0.95)

            print("Saving figure...")
            figs[study][key].savefig(
                os.path.join(
                    "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/figures",
                    study + "_" + case + "_" + key + "_stacked_plot.png",
                ),
                dpi=300,
            )


if __name__ == "__main__":

    plot_all_cases()
    # plot_case(cases["KBHDP"]["KBHDP_SOA_1"]["water_recovery"])
    # plot_case(cases["Permian"]["Permian_RPT1_MD"]["water_recovery"])
    # plot_case(cases["Permian_RPT1_MD"]["water_recovery"])
    # plot_case(cases["Permian"]["Permian_RPT1_FO"]["water_recovery"])
    # plot_case(cases["Permian"]["Permian_RPT1_MD_Cryst"]["grid_fraction"])
    # plot_case(cases["Permian"]["Permian_RPT1_FO_Cryst"]["grid_fraction"])
    # plot_case(cases["Permian"]["Permian_RPT1_MD"]["grid_fraction"])
    # plot_case(cases["Permian_RPT1_FO"]["water_recovery"])
    # plot_case(cases["Permian_ST2"]["water_recovery"])
    plt.show()
