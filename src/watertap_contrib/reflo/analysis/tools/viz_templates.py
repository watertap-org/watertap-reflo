import matplotlib.pyplot as plt

# from matplotlib.ticker import PercentFormatter
import matplotlib.colors as colors
from matplotlib import ticker
import numpy as np


def contour(
    x,
    y,
    z,
    levels,
    x_label="",
    y_label="",
    z_label="",
    x_scale="linear",
    y_scale="linear",
    low=-1,
    mid=0,
    high=1,
    xlimits=[0.001, 0.1],
    ylimits=[30, 150],
    cmap_pallete="RdBu_r",
    **kwargs
):

    divnorm = colors.TwoSlopeNorm(vmin=low, vcenter=mid, vmax=high)
    divnorm = colors.CenteredNorm(vcenter=mid)

    fig, ax = plt.subplots(figsize=(5, 4))
    # fig.figsize(6,6)
    # fig.subplots_adjust()
    # levels = np.linspace(low, high, resolution)
    if levels != None:
        cs2 = ax.contourf(x, y, z, 100, cmap=cmap_pallete, norm=divnorm)
        cs3 = ax.contour(
            x, y, z, levels, colors="k", linewidths=2, linestyles="dashed", norm=divnorm
        )
        ax.clabel(cs3, fmt="$%1.2f  ", colors="k", inline_spacing=60, fontsize=12)
        cbar = fig.colorbar(cs2)
    else:
        cs1 = ax.contourf(x, y, z, 100, cmap=cmap_pallete, norm=divnorm)
        cbar = fig.colorbar(cs1)

    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    # cbar.ax.set_yticklabels

    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)

    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel(z_label, rotation=90, fontsize=16)

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_xlim(xlimits[0], xlimits[1])
    ax.set_ylim(ylimits[0], ylimits[1])
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    fig.tight_layout()
    plt.locator_params(axis="y", nbins=5)
    return fig, ax, cbar


# ex.
# contour([1,2], [1,2], [[1,3],[2,3]], levels=np.array([1,2,3]))
# plt.show()
