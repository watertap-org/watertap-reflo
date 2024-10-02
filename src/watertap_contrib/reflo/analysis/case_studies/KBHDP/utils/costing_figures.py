from psPlotKit.data_plotter.ps_break_down_plotter import breakDownPlotter
from psPlotKit.data_plotter.ps_line_plotter import linePlotter
from psPlotKit.data_manager.ps_data_manager import psDataManager, psData
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pal = {
    "Baseline":"#1f78b4",
    "FFRRO":"#33a02c",
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

def create_sweep_cost_breakdown():
    filepath = os.path.abspath(__file__)
    parent_dir = os.path.dirname(filepath)
    top_level_dir = os.path.dirname(parent_dir)
    save_path = os.path.join(top_level_dir,'figures')
    data_path = os.path.join(top_level_dir,'sweep_results/output/')

    """ import data"""
    costing_data = psDataManager(
        os.path.join(
            data_path,
            "KBHDP_SOA_1_analysisType_KBHDP_SOA_1_sweep.h5",
        )
    )

    """ here we define the costing groups, these should be either units on the block 
    another example is in bgw_ref/costing_plotting_trains_rkt.py that shows different ways to isolate units for cost plotting
    or explicit keys that point to capex or opex for each units. 
    IF you provide exact keys only they will be used (If you duplicate a key it will be used twice!)"""
    device_groups = {
        "Pumps": {
            "CAPEX": {
                "units": {
                    "fs.pump.costing.capital_cost"
                },
            },
            "OPEX": {
                "units": {
                    "fs.pump.control_volume.work[0.0]",
                }
            },
        },
        "RO": {
            "CAPEX": {
                "units": {
                    "fs.RO.stage[1].module.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.RO.stage[1].module.costing.fixed_operating_cost",
                }
            },
        },
        "UF": {
            "CAPEX": {
                "units": {
                    "fs.UF.unit.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.UF.unit.electricity[0.0]",
                },
            }
        },
        # "Deep Well Injection": {
        #     "OPEX": {
        #         "units": {
        #             "fs.disposal.costing.fixed_operating_cost",
        #         },
        #     }
        # },
        "Brine disposal": {
            "OPEX": {
                "units": {
                    "fs.disposal.costing.fixed_operating_cost",
                },
            }
        },
        "Softening": {
            "CAPEX": {
                "units": {
                    "fs.softener.unit.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    'fs.softener.unit.costing.fixed_operating_cost',
                    # "fs.costing.aggregate_flow_costs[co2]",
                    # "fs.costing.aggregate_flow_costs[soda_ash]",
                    # "fs.costing.aggregate_flow_costs[mgcl2]"
                    "fs.costing.aggregate_flow_costs[electricity]",
                },
            }
        },
        "Lime Dosing": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[lime]",
                },
            }
        },
        "Soda Ash Dosing": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[soda_ash]",
                },
            }
        },
        "CO2 Dosing": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[co2]",
                },
            }
        },
        "MgCl2 Dosing": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_costs[mgcl2]",
                },
            }
        },
    }
    costing_data.load_data(
        [
            {
                "filekey": "fs.water_recovery",
                "return_key": "Water recovery",
                "units": "%",
            },
        ],
    )

    """ define the base costing block and flow (This is used to normalize LCOW)
    the costing block is used to pull out default values for cost of 
    electricity, costing factors (TIC etc) and deice costs if applicable, so in theory if you have multiple costing 
    blocks they should all share same factors, so it matters not which one you pick 
    The CAPEX and OPEX will be aggregated from keys in supplied costing groups, the 
    capital and opex costs that are attached to the costing packages will not be used!
    """
    costing_data.get_costing(
        device_groups,
        costing_block="fs.costing",
        default_flow="fs.product.properties[0.0].flow_vol_phase[Liq]",
    )
    """ note how you have access to total and levelized cost. You can plot either one"""
    costing_data.display()

    # # """ get only data for when membrane cost is 1000"""
    # # d = costing_data["ffrro_multiperiod_sweep/map_sweep", "Brine disposal cost"].data
    # # # mask = np.where(d == 1000)[0]
    # # # print(mask)

    # # # costing_data.add_mask("ffrro_multiperiod_sweep/map_sweep", mask)

    """lets plot our data"""
    cost_plotter = breakDownPlotter(
        costing_data,
        save_location=save_path,
        save_name="cost_breakdown_blank",
    )
    print(cost_plotter)
    """ define the costing groups, this will be order they are plotted in"""
    cost_plotter.define_area_groups(
        [
            "Brine disposal",
            "UF",
            "Pumps",
            "RO",
            "Lime Dosing",
            "Soda Ash Dosing",
            "CO2 Dosing",
            "MgCl2 Dosing",
            "Softening",
        ]
    )
    """ define if you want to plot specific groups, for example CAPEX, OPEX or TOTAl isstead"""
    cost_plotter.define_hatch_groups({"CAPEX": {"hatch": ""}, "OPEX": {"hatch": "//"}})
    cost_plotter.plotbreakdown(
        xdata="Water recovery",
        ydata=[
            "cost_breakdown",
            "levelized",
        ],  # remove "levelized to plot absolute capex/opex
        axis_options={
            "yticks": [0, 1, 2, 3, 4, 5, 6, 7],  # adjust as needed
            "xticks": [20, 30, 40, 50, 60, 70],
        },
        legend_loc="upper right",
        generate_figure=True,
    )

    # # cost_plotter.fig.plot_line(
    # #                 xdata=[80, 95],
    # #                 ydata=[0.223 , 0.223],
    # #                 ls='--',
    # #                 color='#e31a1c',
    # #                 label="2-Stage Baseline LCOW",
    # #                 lw=2,
    # #             )
    
    # # cost_plotter.fig.plot_line(
    # #                 xdata=[84, 84],
    # #                 ydata=[0 , 0.223],
    # #                 ls='--',
    # #                 color='k',
    # #                 label="Breakeven Recovery",
    # #                 lw=2,
    # #             )
    
    # cost_plotter.generate_figure(loc="upper left")

def create_system_comparison():
    filepath = os.path.abspath(__file__)
    parent_dir = os.path.dirname(filepath)
    save_path = os.path.join(parent_dir,'figures')
    data_path = os.path.join(parent_dir,'presentations/data/system_comp.csv')

    # default_font = {
    #         "family": "serif",
    #         "serif": "Arial",
    #         "weight": "normal",
    #         "size": 10,
    #     }

    # df = pd.read_csv(data_path)

    # sns.axes_style({'font.serif': 'Arial'})
    # g = sns.FacetGrid(df, col="Variable", hue="System", col_wrap=3, height=4.5, aspect=.55, palette=pal, sharey=False)

    # g.map(sns.barplot, "System", "Value", alpha=.99, edgecolor='k')

    # g.axes[0].set_ylabel("MGD",fontsize=12)
    # g.axes[1].set_ylabel("MGD",fontsize=12)
    # g.axes[2].set_ylabel("LCOW $\$$/m$^3$",fontsize=12)

    # g.axes[0].set_title("Water Production",fontsize=12)
    # g.axes[1].set_title("Brine Disposal",fontsize=12)
    # g.axes[2].set_title("LCOW",fontsize=12)

    # g.axes[0].set_ylim(0, 20)
    # g.axes[1].set_ylim(0, 4)
    # g.axes[2].set_ylim(0, 0.25)

    # g.axes[0].set_xlabel("")
    # g.axes[1].set_xlabel("")
    # g.axes[2].set_xlabel("")

    # for ax in g.axes:
    #     ax.set_xticklabels(ax.get_xticklabels(),fontsize=12)
    #     ax.set_yticklabels(ax.get_yticklabels(),fontsize=12)

    # g.figure.subplots_adjust(wspace=.7)
    # plt.tight_layout()
    # # plt.savefig(os.path.join(save_path, "system_comparison.png"), dpi=900)
    # plt.show()

if __name__ == "__main__":
    # create_system_comparison()
    create_sweep_cost_breakdown()