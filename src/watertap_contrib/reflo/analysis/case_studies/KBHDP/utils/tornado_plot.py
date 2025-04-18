import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import math


def tornado_plot(
    baseline_lcow,
    df,
    title="<Low> VS <High> values",
    offset=1,
    xlim=None,
    ylim=None,
    y_offset=0.3,
    ylab_offset=0.1,
    xticks=None,
    textsize=13,
    ticksize=13,
    fig=None,
    ax=None,
):
    """
    Parameters
    ----------
    baseline_lcow: The LCOW assuming baseline costing parameters
    df: Dataframe containing-
        1. Sensitivity parameter labels
        2. Sensitivity parameter baseline price
        2. Low values
        3. High values

    tite: Title of the plot
    """

    color_low = "royalblue"  #'#e1ceff'
    color_high = "darkorange"  # ff6262'

    ys = range(len(df["labels"]))[::1]  # iterate through # of labels

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 5, forward=True)

    # make vertical dashed line at baseline LCOW
    ax.axvline(baseline_lcow, ymax=10, color="black", linewidth=2, ls=":", zorder=0)

    for y, low_value, high_value in zip(
        ys, df["lcow_low_values"], df["lcow_high_values"]
    ):

        low_width = baseline_lcow - low_value
        high_width = high_value - baseline_lcow

        ax.broken_barh(
            [(low_value, low_width), (baseline_lcow, high_width)],
            (y - y_offset, y_offset * 2),  # thickness of bars and their offset
            facecolors=[color_low, color_high],
            edgecolors=["black", "black"],
            linewidth=1,
        )

        # offset = 1  # offset value labels from end of bar

        if high_value > low_value:
            x_high = baseline_lcow + high_width + offset
            x_low = baseline_lcow - low_width - offset
        else:
            x_high = baseline_lcow + high_width - offset
            x_low = baseline_lcow - low_width + offset

        x_low_change = (baseline_lcow - low_value) / baseline_lcow * 100
        x_high_change = (high_value - baseline_lcow) / baseline_lcow * 100

        ax.text(
            x_high,
            y + ylab_offset,
            f"{high_value:0.2f}",
            va="center",
            ha="center",
            fontsize=textsize,
        )
        ax.text(
            x_low,
            y + ylab_offset,
            f"{low_value:0.2f}",
            va="center",
            ha="center",
            fontsize=textsize,
        )

        print(x_low_change, x_high_change)
        ax.text(
            x_low,
            y - ylab_offset,
            "(-" + str(f"{x_low_change:0.1f}") + "%)",
            va="center",
            ha="center",
            fontsize=textsize,
        )
        if x_high_change == 0:
            ax.text(
                x_high,
                y - ylab_offset,
                "(0.0%)",
                va="center",
                ha="center",
                fontsize=textsize,
            )
        else:
            ax.text(
                x_high,
                y - ylab_offset,
                "(+" + str(f"{x_high_change:0.1f}") + "%)",
                va="center",
                ha="center",
                fontsize=textsize,
            )

    # ax.text(
    #     (
    #         math.floor(min(df["lcow_low_values"]))
    #         + math.ceil(max(df["lcow_high_values"]))
    #     )
    #     / 2,
    #     len(df["labels"]) + 0.2,
    #     title,
    #     va="center",
    #     ha="center",
    #     # fontdict={"fontsize": 20},
    #     fontsize=textsize
    # )

    # ax.set_xlabel("LCOW (\$/m$^3$)")
    # fig.supxlabel("LCOW (\$/m$^3$)", fontsize=14)
    # ax.set_title(title, fontdict=dict(fontsize=ticksize + 2))
    ax.set_yticks(ys, df["labels"], verticalalignment="center")
    if xticks is None:
        ax.set_xticks(np.arange(0, max(df["lcow_high_values"]) + 3, 1))

    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(
            math.floor(min(df["lcow_low_values"])) - 2,
            math.ceil(max(df["lcow_high_values"])) + 2,
        )
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(-0.5, len(df["labels"]) - 0.5)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    # ax.text(xmax, ymax, f"{baseline_lcow:0.2f} \$/m$^3$", fontsize=10)
    # props = dict(boxstyle="round", facecolor=None, alpha=0.25, bbox=None)
    props = dict(edgecolor="none", facecolor="white", alpha=0)
    ax.text(
        0.05,
        0.85,
        # f"{title}\n{baseline_lcow:0.2f} \$/m$^3$",
        f"{title}",
        fontsize=ticksize + 10,
        transform=ax.transAxes,
        bbox=props,
        ha="center",
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=ticksize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.grid(visible=True)

    # fig = plt.gcf()
    # fig.tight_layout()

    # plt.show()
    # plt.tight_layout()

    return fig, ax


def tornado_plot_rel(
    baseline_lcow,
    df,
    title="<Low> VS <High> values",
    offset=1,
    xlim=None,
    ylim=None,
    y_offset=0.3,
    ylab_offset=0.1,
    xticks=None,
    textsize=13,
    ticksize=13,
    fig=None,
    ax=None,
    title_x=0.05,
    title_y=0.8,
):
    """
    Parameters
    ----------
    baseline_lcow: The LCOW assuming baseline costing parameters
    df: Dataframe containing-
        1. Sensitivity parameter labels
        2. Sensitivity parameter baseline price
        2. Low values
        3. High values

    tite: Title of the plot
    """

    color_low = "royalblue"  #'#e1ceff'
    color_high = "darkorange"  # ff6262'

    ys = range(len(df["labels"]))[::1]  # iterate through # of labels

    if (fig, ax) == (None, None):
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 5, forward=True)

    # make vertical dashed line at baseline LCOW
    ax.axvline(0, ymax=10, color="black", linewidth=2, ls=":", zorder=0)

    for y, low_value, high_value in zip(
        ys, df["lcow_low_values"], df["lcow_high_values"]
    ):

        low_width = (baseline_lcow - low_value) / baseline_lcow
        high_width = (high_value - baseline_lcow) / baseline_lcow

        ax.broken_barh(
            [(0 - low_width, low_width), (0, high_width)],
            (y - y_offset, y_offset * 2),  # thickness of bars and their offset
            facecolors=[color_low, color_high],
            edgecolors=["black", "black"],
            linewidth=1,
        )

        print(f"xmin low = {1 - low_value / baseline_lcow}")
        print(f"low_width = {low_width}")
        print(f"xmin high = {1}")
        print(f"low_width = {high_width}")

        # offset = 1  # offset value labels from end of bar

        if high_value > low_value:
            x_high = 0 + high_width + offset
            x_low = 0 - low_width - offset
        else:
            x_high = 0 + high_width - offset
            x_low = 0 - low_width + offset

        x_low_change = (baseline_lcow - low_value) / baseline_lcow * 100
        x_high_change = (high_value - baseline_lcow) / baseline_lcow * 100

        ax.text(
            x_high,
            y + ylab_offset,
            f"{high_value:0.2f}",
            va="center",
            ha="center",
            fontsize=textsize,
        )
        ax.text(
            x_low,
            y + ylab_offset,
            f"{low_value:0.2f}",
            va="center",
            ha="center",
            fontsize=textsize,
        )

        print(x_low_change, x_high_change)
        ax.text(
            x_low,
            y - ylab_offset,
            "(-" + str(f"{x_low_change:0.1f}") + "%)",
            va="center",
            ha="center",
            fontsize=textsize,
        )
        if x_high_change == 0:
            ax.text(
                x_high,
                y - ylab_offset,
                "(0.0%)",
                va="center",
                ha="center",
                fontsize=textsize,
            )
        else:
            ax.text(
                x_high,
                y - ylab_offset,
                "(+" + str(f"{x_high_change:0.1f}") + "%)",
                va="center",
                ha="center",
                fontsize=textsize,
            )

    # ax.text(
    #     (
    #         math.floor(min(df["lcow_low_values"]))
    #         + math.ceil(max(df["lcow_high_values"]))
    #     )
    #     / 2,
    #     len(df["labels"]) + 0.2,
    #     title,
    #     va="center",
    #     ha="center",
    #     # fontdict={"fontsize": 20},
    #     fontsize=textsize
    # )

    # ax.set_xlabel("LCOW (\$/m$^3$)")
    # fig.supxlabel("LCOW (\$/m$^3$)", fontsize=14)
    # ax.set_title(title, fontdict=dict(fontsize=ticksize + 2))
    ax.set_yticks(ys, df["labels"], verticalalignment="center")
    if xticks is None:
        # ax.set_xticks(np.arange(0, max(df["lcow_high_values"]) + 3, 1))
        pass

    if xlim is not None:
        ax.set_xlim(xlim)
    # else:
    #     ax.set_xlim(
    #         # math.floor(min(df["lcow_low_values"])) - 2,
    #         # math.ceil(max(df["lcow_high_values"])) + 2,
    #     )
    if ylim is not None:
        ax.set_ylim(ylim)
    # else:
    #     ax.set_ylim(-0.5, len(df["labels"]) - 0.5)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    # ax.text(xmax, ymax, f"{baseline_lcow:0.2f} \$/m$^3$", fontsize=10)
    # props = dict(boxstyle="round", facecolor=None, alpha=0.25, bbox=None)
    props = dict(edgecolor="none", facecolor="white", alpha=0)
    ax.text(
        title_x,
        title_y,
        # f"{title}\n{baseline_lcow:0.2f} \$/m$^3$",
        f"{title}",
        fontsize=ticksize + 10,
        transform=ax.transAxes,
        bbox=props,
        ha="center",
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=ticksize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.0f}%"))
    # ax.yaxis.grid(True)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.grid(visible=True)

    # fig = plt.gcf()
    # fig.tight_layout()

    # plt.show()
    # plt.tight_layout()

    return fig, ax


if __name__ == "__main__":

    labels = pd.Series(
        {
            "fpc_collector": "FPC Collector Cost (\$/m$^2$)\n 600 [300,1200]",
            "fpc_storage": "Thermal Storage Cost (\$/m$^3$)\n 2000 [300,1200]",
            "dwi_lcow": "Injection Cost (\$/m$^3$)\n 0.0587 [300,1200]",
            "heat_price": "Heat Price ($/kWh)\n 0.0166 [300,1200]",
        }
    )

    # value order corresponds to label order
    lcow_low_values = pd.Series(
        {
            "fpc_collector": 7.46,
            "fpc_storage": 8.21,
            "dwi_lcow": 9.63,
            "heat_price": 9.35,
        }
    )
    lcow_high_values = pd.Series(
        {
            "fpc_collector": 10.73,
            "fpc_storage": 10.35,
            "dwi_lcow": 9.65,
            "heat_price": 9.78,
        }
    )

    df = pd.DataFrame(columns=["labels", "lcow_low_values", "lcow_high_values"])
    df["labels"] = labels
    df["lcow_low_values"] = lcow_low_values
    df["lcow_high_values"] = lcow_high_values

    baseline_lcow = 9.64

    var_effect = np.abs(lcow_high_values - lcow_low_values) / baseline_lcow

    df["range"] = var_effect

    df = df.sort_values(
        "range", ascending=True, inplace=False, ignore_index=False, key=None
    )

    tornado_plot(baseline_lcow, df, title="RPT3")
