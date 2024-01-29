#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import numpy as np
from pyomo.environ import value, units as pyunits
import pandas as pd
import matplotlib.pyplot as plt


def get_data(mp):
    n = 24
    hour = [i for i in range(1, n + 1)]
    battery_state = np.array(
        [value(mp.blocks[i].process.fs.battery.state_of_charge[0]) for i in range(n)]
    )
    pv_gen = np.array(
        [value(mp.blocks[i].process.fs.elec_generation) for i in range(n)]
    )
    pv_curtail = np.array(
        [value(mp.blocks[i].process.fs.curtailment) for i in range(n)]
    )
    electric_price = np.array(
        [value(mp.blocks[i].process.fs.elec_price) for i in range(n)]
    )
    ro_demand = np.array([value(mp.blocks[i].process.fs.elec_price) for i in range(n)])
    grid_cost = [
        value(mp.blocks[i].process.grid_cost)
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    lcow = [
        0.4
        + (
            (value(mp.annualized_capital_cost) / 365 / 24)
            + value(mp.blocks[i].process.grid_cost)
        )
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    pv_to_ro = np.array([value(mp.blocks[i].process.fs.pv_to_ro) for i in range(n)])
    pv_to_battery = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_in[0]) for i in range(n)]
    )
    battery_to_ro = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_out[0]) for i in range(n)]
    )
    grid_to_ro = np.array([value(mp.blocks[i].process.fs.grid_to_ro) for i in range(n)])
    labels = ["PV to RO", "Battery to RO", "Grid to RO", "PV to Battery"]

    frames = []
    for label in labels:
        frames.append(
            pd.DataFrame([hour, pv_to_ro, [label] * len(pv_to_ro)]).transpose()
        )
    df2 = pd.concat(frames)
    df2.columns = ["Hour", "State", "Type"]

    df = pd.DataFrame(
        [hour, pv_to_ro, pv_to_battery, battery_to_ro, grid_to_ro]
    ).transpose()
    df.columns = ["Hour", "PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    features = ["PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    df["Total"] = (
        df["PV to RO"] + df["Battery to RO"] + df["Grid to RO"] + df["PV to Battery"]
    )

    return df


def create_plot(mp):
    fig, ax = plt.subplots(figsize=(20, 4))
    colors = ["#235789", "#4A7C59", "#F1A208"]
    color1 = "#3971ad"
    color2 = "#c07432"
    color3 = "#8c8b8b"
    n = 24
    hour = [i for i in range(1, n + 1)]
    battery_state = np.array(
        [value(mp.blocks[i].process.fs.battery.state_of_charge[0]) for i in range(n)]
    )
    pv_gen = np.array(
        [value(mp.blocks[i].process.fs.elec_generation) for i in range(n)]
    )
    pv_curtail = np.array(
        [value(mp.blocks[i].process.fs.curtailment) for i in range(n)]
    )
    electric_price = np.array(
        [value(mp.blocks[i].process.fs.electricity_price) for i in range(n)]
    )
    ro_demand = np.array(
        [value(mp.blocks[i].process.fs.RO.power_demand) for i in range(n)]
    )
    grid_cost = [
        value(mp.blocks[i].process.grid_cost)
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    lcow = [
        0.4
        + (
            (value(mp.annualized_capital_cost) / 365 / 24)
            + value(mp.blocks[i].process.grid_cost)
        )
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    pv_to_ro = np.array([value(mp.blocks[i].process.fs.pv_to_ro) for i in range(n)])
    pv_to_battery = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_in[0]) for i in range(n)]
    )
    battery_to_ro = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_out[0]) for i in range(n)]
    )
    grid_to_ro = np.array([value(mp.blocks[i].process.fs.grid_to_ro) for i in range(n)])

    labels = ["PV to RO", "Battery to RO", "Grid to RO", "PV to Battery"]
    norm_labels = ["PV to RO", "Battery to RO", "PV to Battery", "Grid to RO"]

    frames = []
    for label in labels:
        frames.append(
            pd.DataFrame([hour, pv_to_ro, [label] * len(pv_to_ro)]).transpose()
        )
    df2 = pd.concat(frames)
    df2.columns = ["Hour", "State", "Type"]

    df = pd.DataFrame(
        [hour, pv_to_ro, pv_to_battery, battery_to_ro, grid_to_ro]
    ).transpose()
    df.columns = ["Hour", "PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    features = ["PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    df["Total"] = (
        df["PV to RO"] + df["Battery to RO"] + df["Grid to RO"] + df["PV to Battery"]
    )

    ax.stackplot(
        hour,
        df["PV to RO"],
        df["Battery to RO"],
        df["Grid to RO"],
        baseline="zero",
        colors=["#1f77b4", "#ff7f0e", "#d62728"],
        labels=labels,
        alpha=1,
        ec="white",
    )
    ax.plot(
        hour, df["PV to Battery"], label="PV to Battery", color="#2ca02c", linewidth=2
    )
    ax.fill_between(
        hour,
        df["PV to Battery"],
        color="#2ca02c",
        hatch="////",
        edgecolor="#515251",
        linewidth=2,
        alpha=0.5,
    )
    ax.set_ylabel("  Power (kW)", loc="center", fontsize=16)

    ax.set_xlabel("Operation Hours", fontsize=16)
    # leg3 = ax.legend(loc="lower left", frameon = True, bbox_to_anchor=(0, 1.0, 0.65, 1),
    #                 ncols=4, mode="expand", fontsize=14, borderaxespad=0.,
    #                 facecolor='#8f8f8f', edgecolor=None, framealpha=1)
    ax.set_xlim([1, n])
    ax.set_ylim([0, 1000])
    # ax.set_title(titles[idx], loc='center', x=-0.08, y=0.5, rotation=90, fontweight='bold', ha='center', va='center', fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    ax4 = ax.twinx()
    ax5 = ax.twinx()

    ax4.spines.right.set_position(("axes", 1.1))

    line1 = ax4.plot(
        hour, electric_price, dashes=[4, 4], lw=2.5, color="white", label="Grid Price"
    )
    # line1 = ax4.plot(hour, electric_price, dashes=[6, 4], color="white", label='Grid Price')
    ax4.set_ylabel(
        "Grid Price ($/kWh)", ha="center", va="center", fontsize=16, labelpad=20
    )
    ax4.set_ylim([0, 0.6])
    ax4.yaxis.set_major_formatter("${x:1.2f}")
    # leg4 = ax4.legend(loc="lower left", frameon = True, bbox_to_anchor=(0.8, 1.0, 0.15, 1),
    #     ncols=1, mode="expand", fontsize=14, borderaxespad=0.,
    #     facecolor='#8f8f8f', edgecolor=None, framealpha=1)

    line2 = ax5.plot(hour, lcow, linestyle="dashed", color="k", lw=2.5, label="LCOW")
    ax5.set_ylabel(
        f"LCOW ({str('$')}" + f"/m3)",
        ha="center",
        va="center",
        fontsize=16,
        labelpad=20,
    )
    ax5.set_ylim([0, 1.5])
    ax5.yaxis.set_major_formatter("${x:1.2f}")
    # leg5 = ax5.legend(loc="lower left", frameon = True, bbox_to_anchor=(0.67, 1.0, 0.15, 1),
    #     ncols=1, mode="expand", fontsize=14, borderaxespad=0.,
    #     facecolor='#8f8f8f', edgecolor=None, framealpha=1)

    ax4.tick_params(axis="y", labelsize=16)
    ax5.tick_params(axis="y", labelsize=16)

    ax4.spines["right"].set_color("#6e6e6e")
    ax5.spines["right"].set_color("k")
    ax4.tick_params(axis="y", colors="#6e6e6e")
    ax5.tick_params(axis="y", colors="k")
    ax4.yaxis.label.set_color("#6e6e6e")
    ax5.yaxis.label.set_color("k")

    ab = ax5.annotate(
        f"LCOW=${mp.LCOW():1.2f}",
        (0.01, 0.9),
        xycoords="axes fraction",
        fontsize=16,
        color="k",
        bbox=dict(boxstyle="square", fc="white", ec="k", lw=1),
    )

    ac = ax5.annotate(
        f"Battery Size={value(mp.blocks[0].process.fs.battery.nameplate_energy):1.0f} kWh",
        (0.82, 0.9),
        xycoords="axes fraction",
        fontsize=16,
        color="k",
        bbox=dict(boxstyle="square", fc="white", ec="k", lw=1),
    )

    ab.set_zorder(100)

    plt.figlegend(
        loc="upper left",
        fontsize=16,
        frameon=True,
        bbox_to_anchor=(0.12, 0.51, 1.0, 0.51),
        ncols=6,
        facecolor="#d4d2d2",
        edgecolor=None,
        framealpha=1,
    )


def create_long_plot(mp):
    fig, ax = plt.subplots(figsize=(20, 4))
    colors = ["#235789", "#4A7C59", "#F1A208"]
    color1 = "#3971ad"
    color2 = "#c07432"
    color3 = "#8c8b8b"
    n = len(mp.blocks)
    hour = [i for i in range(1, n + 1)]
    battery_state = np.array(
        [value(mp.blocks[i].process.fs.battery.state_of_charge[0]) for i in range(n)]
    )
    pv_gen = np.array(
        [value(mp.blocks[i].process.fs.elec_generation) for i in range(n)]
    )
    pv_curtail = np.array(
        [value(mp.blocks[i].process.fs.curtailment) for i in range(n)]
    )
    electric_price = np.array(
        [value(mp.blocks[i].process.fs.electricity_price) for i in range(n)]
    )
    ro_demand = np.array(
        [value(mp.blocks[i].process.fs.RO.power_demand) for i in range(n)]
    )
    grid_cost = [
        value(mp.blocks[i].process.grid_cost)
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    lcow = [
        0.4
        + (
            (value(mp.annualized_capital_cost) / 365 / 24)
            + value(mp.blocks[i].process.grid_cost)
        )
        / value(
            pyunits.convert(
                6000 * pyunits.m**3 / pyunits.day,
                to_units=pyunits.m**3 / pyunits.hour,
            )
        )
        for i in range(n)
    ]
    pv_to_ro = np.array([value(mp.blocks[i].process.fs.pv_to_ro) for i in range(n)])
    pv_to_battery = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_in[0]) for i in range(n)]
    )
    battery_to_ro = np.array(
        [value(mp.blocks[i].process.fs.battery.elec_out[0]) for i in range(n)]
    )
    grid_to_ro = np.array([value(mp.blocks[i].process.fs.grid_to_ro) for i in range(n)])

    labels = ["PV to RO", "Battery to RO", "Grid to RO", "PV to Battery"]
    norm_labels = ["PV to RO", "Battery to RO", "PV to Battery", "Grid to RO"]

    frames = []
    for label in labels:
        frames.append(
            pd.DataFrame([hour, pv_to_ro, [label] * len(pv_to_ro)]).transpose()
        )
    df2 = pd.concat(frames)
    df2.columns = ["Hour", "State", "Type"]

    df = pd.DataFrame(
        [hour, pv_to_ro, pv_to_battery, battery_to_ro, grid_to_ro]
    ).transpose()
    df.columns = ["Hour", "PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    features = ["PV to RO", "PV to Battery", "Battery to RO", "Grid to RO"]
    df["Total"] = (
        df["PV to RO"] + df["Battery to RO"] + df["Grid to RO"] + df["PV to Battery"]
    )

    ax.stackplot(
        hour,
        df["PV to RO"],
        df["Battery to RO"],
        df["Grid to RO"],
        baseline="zero",
        colors=["#1f77b4", "#ff7f0e", "#d62728"],
        labels=labels,
        alpha=1,
        ec="white",
    )
    ax.plot(
        hour, df["PV to Battery"], label="PV to Battery", color="#2ca02c", linewidth=2
    )
    ax.fill_between(
        hour,
        df["PV to Battery"],
        color="#2ca02c",
        hatch="////",
        edgecolor="#515251",
        linewidth=2,
        alpha=0.5,
    )
    ax.set_ylabel("  Power (kW)", loc="center", fontsize=16)

    ax.set_xlabel("Operation Hours", fontsize=16)
    # leg3 = ax.legend(loc="lower left", frameon = True, bbox_to_anchor=(0, 1.0, 0.65, 1),
    #                 ncols=4, mode="expand", fontsize=14, borderaxespad=0.,
    #                 facecolor='#8f8f8f', edgecolor=None, framealpha=1)
    ax.set_xlim([1, n])
    ax.set_ylim([0, 1000])
    # ax.set_title(titles[idx], loc='center', x=-0.08, y=0.5, rotation=90, fontweight='bold', ha='center', va='center', fontsize=16)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    ax4 = ax.twinx()
    ax5 = ax.twinx()

    ax4.spines.right.set_position(("axes", 1.1))

    line1 = ax4.plot(
        hour, electric_price, dashes=[4, 4], lw=2.5, color="white", label="Grid Price"
    )
    # line1 = ax4.plot(hour, electric_price, dashes=[6, 4], color="white", label='Grid Price')
    ax4.set_ylabel(
        "Grid Price ($/kWh)", ha="center", va="center", fontsize=16, labelpad=20
    )
    ax4.set_ylim([0, 0.6])
    ax4.yaxis.set_major_formatter("${x:1.2f}")
    # leg4 = ax4.legend(loc="lower left", frameon = True, bbox_to_anchor=(0.8, 1.0, 0.15, 1),
    #     ncols=1, mode="expand", fontsize=14, borderaxespad=0.,
    #     facecolor='#8f8f8f', edgecolor=None, framealpha=1)

    line2 = ax5.plot(hour, lcow, linestyle="dashed", color="k", lw=2.5, label="LCOW")
    ax5.set_ylabel(
        f"LCOW ({str('$')}" + f"/m3)",
        ha="center",
        va="center",
        fontsize=16,
        labelpad=20,
    )
    ax5.set_ylim([0, 1.5])
    ax5.yaxis.set_major_formatter("${x:1.2f}")
    # leg5 = ax5.legend(loc="lower left", frameon = True, bbox_to_anchor=(0.67, 1.0, 0.15, 1),
    #     ncols=1, mode="expand", fontsize=14, borderaxespad=0.,
    #     facecolor='#8f8f8f', edgecolor=None, framealpha=1)

    ax4.tick_params(axis="y", labelsize=16)
    ax5.tick_params(axis="y", labelsize=16)

    ax4.spines["right"].set_color("#6e6e6e")
    ax5.spines["right"].set_color("k")
    ax4.tick_params(axis="y", colors="#6e6e6e")
    ax5.tick_params(axis="y", colors="k")
    ax4.yaxis.label.set_color("#6e6e6e")
    ax5.yaxis.label.set_color("k")

    ab = ax5.annotate(
        f"LCOW=${mp.LCOW():1.2f}",
        (0.01, 0.9),
        xycoords="axes fraction",
        fontsize=16,
        color="k",
        bbox=dict(boxstyle="square", fc="white", ec="k", lw=1),
    )

    ac = ax5.annotate(
        f"Battery Size={value(mp.blocks[0].process.fs.battery.nameplate_energy):1.0f} kWh",
        (0.82, 0.9),
        xycoords="axes fraction",
        fontsize=16,
        color="k",
        bbox=dict(boxstyle="square", fc="white", ec="k", lw=1),
    )

    ab.set_zorder(100)

    plt.figlegend(
        loc="upper left",
        fontsize=16,
        frameon=True,
        bbox_to_anchor=(0.12, 0.51, 1.0, 0.51),
        ncols=6,
        facecolor="#d4d2d2",
        edgecolor=None,
        framealpha=1,
    )
