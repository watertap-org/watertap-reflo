import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mtick
import matplotlib.patches as patches
from matplotlib import ticker

from pyomo.environ import (Var,
                           value,
                           units as pyunits)

def mgd_to_m3s(Q):
    return 0.0438*Q

def get_data(data_dict, metric, gridshape: tuple):
    x, y, z = np.array([]), np.array([]), []
    grid = np.empty(gridshape, dtype=float, order='C')
    for key, val in data_dict.items():
        x = np.append(x, val['var2'])
        y = np.append(y, val['var1'])
        if val['results'] != np.nan:
            try:
                z.append(getattr(val['model'].fs.sys_costing, metric)())
                grid[val['idx1'],val['idx2']] = getattr(val['model'].fs.sys_costing, metric)()
            except:
                print('fail')
                z.append(np.nan)
        else:
            z.append(np.nan)
            grid[val['idx1'],val['idx2']] = getattr(val['model'].fs.sys_costing, metric)()
    return np.unique(x), np.unique(y), z, grid

def print_coag_floc_results(m):
    print(f'\n{"=======> COAG FLOC RESULTS <=======":^60}\n')
    print(f'{"PARAM":<25s}{"VALUE":<25s}{"UNITS":<10s}')
    print(
        f'{"Capital Cost":<24s}{f"${value(pyunits.convert(m.fs.coag_and_floc.costing.capital_cost, to_units=pyunits.USD_2021)):<25,.3f}"}{"$":<10s}'
    )
    print(
        f'{"LCOW":<24s}{f"${m.fs.coag_and_floc.costing.LCOW():<25,.3f}"}{"$/m3":<10s}'
    )

def contour_LCOW(x, y, z, levels=None, x_label='', y_label='', z_label='', x_scale='linear', y_scale='linear', low=-1, mid=0, high=1, xlimits = [0.001, 0.1], ylimits = [30, 150], cmap_pallete="RdBu_r", contour_x_pos=1, contour_label_space=0, auto_ticks=None, **kwargs):

    divnorm = colors.TwoSlopeNorm(vmin=low, vcenter=mid, vmax=high)
    
    if levels == None:
        levels = np.arange(round(low,1), round(high,1), (round(high,1)-round(low,1))/4)
    
    fig, ax = plt.subplots(figsize=(5,4))
    
    cs2 = ax.contourf(x, y, z, 100, cmap=cmap_pallete, norm=divnorm)
    cs3 = ax.contour(x, y, z, levels, colors='k', linewidths=2, linestyles='dashed', norm=divnorm)
    ax.clabel(cs3, fmt='$%1.2f ', colors='k', fontsize=12, inline_spacing=contour_label_space)
    cbar = fig.colorbar(cs2, ticks = auto_ticks)

    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)

    cbar.ax.yaxis.set_major_formatter('${x:1.2f}')
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel(z_label, rotation=90, fontsize=16)

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_xlim(xlimits[0],xlimits[1])
    ax.set_ylim(ylimits[0],ylimits[1])
    ax.tick_params(axis='x', labelsize = 18)
    ax.tick_params(axis='y', labelsize = 18)
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    fig.tight_layout()
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='x', nbins=2)
    plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xy = (xmin,ymin)
    width = xmax - xmin
    height = ymax - ymin
    p = patches.Rectangle(xy, width, height, fc='#393E41', fill=True, zorder=-10)
    ax.add_patch(p)

    return fig, ax, cbar