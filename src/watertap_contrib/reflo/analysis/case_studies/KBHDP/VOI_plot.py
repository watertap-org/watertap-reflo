import os

from analysis_plot_kit.core import (
    # fig_generator,
    data_import,
    # data_collator,
)

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import yaml
from scipy import stats

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(os.path.dirname(parent_dir), "sweep_yamls")
sweep_results_dir = os.path.join(os.path.dirname(parent_dir), "sweep_results", "output")
save_path = os.path.join(os.path.dirname(parent_dir), "figures/")

default_font = {
            "family": "serif",
            "serif": "Arial",
            "weight": "normal",
            "size": 10,
        }

default_label_size = {
            "labelsize": 12,
        }

matplotlib.rc("font", **default_font)
matplotlib.rc("axes", **default_label_size)

"""Set global math text to regular"""
default_math_text = {"mathtext.default": "regular"}
plt.rcParams.update(default_math_text)
plt.rcParams.update({"svg.fonttype": "none"})

diff_files = ["diff_sweep_analysisType_Permian_ST_1_diff_analysis.h5"]

def fetch_file(file=None, yaml_file = None, case = None):
    yaml_file_path = os.path.join(sweep_yaml_dir, yaml_file)
    yaml_file = yaml.safe_load(open(yaml_file_path, "r"))
    tree_top = list(yaml_file.keys())[0]
    sweeps = list(yaml_file[tree_top]['diff_param_loop'].keys())
    groups = {}
    x_keys = {}
    conversions = {}
    print(sweeps)
    print(tree_top)
    for sweep in sweeps[:-1]:
        groups[sweep] = yaml_file[tree_top]['diff_param_loop'][sweep]['group']
        x_keys[sweep] = yaml_file[tree_top]['diff_param_loop'][sweep]['param']
        conversions[sweep] = yaml_file[tree_top]['diff_param_loop'][sweep]['y_convert']

    data_types = []
    for sweep in sweeps[:-1]:
        data = data_import.waterTAP_dataImport(
            file
            # os.path.join(sweep_results_dir, file)
        )

        data.set_data_keys(
            [
                list(yaml_file.keys())[0],
                sweep,
            ]
        )

        data_types.append(
            {
                "data_manager": data,
                "label": "VOI",
                "sweep": sweep,
                "case": case,
                "color": "#a6cee3",
                "ls": "-",
                "marker": "",
            },
        )

    grouped_figures = [
        {
            "y_key": "fs.costing.LCOT",
            "x_convert": 1,
            "y_convert": 1,
            "xlabel": "VOI",
            "position": 0,
        },
    ]

    x_data = []
    y_data = []
    x_diff = []
    y_diff = []
    VOI = []
    sweep = []
    group = []
    case = []

    for ft in grouped_figures:
        for dt in data_types:
            # print(dt['case'])
            # print(' ',x_keys[dt['case']])
            print(dt['sweep'])
            print(conversions[dt["sweep"]])

            x = (
                dt["data_manager"].retrive_loop_data(
                    x_keys[dt['sweep']], [1], x_keys[dt['sweep']], stack_data=False
                )
                * ft["x_convert"]
            )[0]
            y = (
                dt["data_manager"].retrive_loop_data(
                    x_keys[dt['sweep']], [1], ft["y_key"], stack_data=False
                )
                * ft["y_convert"] * conversions[dt["sweep"]]
            )[0]
            dx = (
                dt["data_manager"].retrive_loop_data(
                    x_keys[dt['sweep']], [1], x_keys[dt['sweep']], stack_data=False, get_diff=True
                )
                * ft["x_convert"]
            )[0]
            dy = (
                dt["data_manager"].retrive_loop_data(
                    x_keys[dt['sweep']], [1], ft["y_key"], stack_data=False, get_diff=True
                )
                * ft["y_convert"] * conversions[dt["sweep"]]
            )[0]

            # print(x)
            x_data.extend(x)
            y_data.extend(y)
            x_diff.extend(dx)
            y_diff.extend(dy)
            if dt['sweep'] == 'Water Permeability':
                print(dx)
                print(dy)
                # assert False
            sweep.extend([dt["sweep"] for i in range(len(dx))])
            group.extend([groups[dt["sweep"]] for i in range(len(dx))])
            case.extend([dt["case"] for i in range(len(dx))])

    voi = np.array(y_diff)/np.array(x_diff)
    VOI.extend(voi)

    # print(x_diff)
    # print(y_diff)
    # print(VOI)
    # print(sweep)
    # print(group)
    # print(case)

    df = pd.DataFrame({'del_perform': x_diff,'del_LCOW': y_diff,'VOI': VOI, 'Sweep': sweep, 'Group': group, 'Case': case})
    df = df.sort_values(by='Group')
    # df.loc[df['Sweep'] == 'EC Current Efficiency', 'VOI'] = df.loc[df['Sweep'] == 'EC Current Efficiency', 'VOI'] * -1
    # df.loc[df['Sweep'] == 'Water Permeability', 'VOI'] = df.loc[df['Sweep'] == 'Water Permeability', 'VOI'] * -1
    # df.loc[df['Sweep'] == 'Salt Permeability', 'VOI'] = df.loc[df['Sweep'] == 'Salt Permeability', 'VOI'] * -1
    # df.loc[df['Sweep'] == 'Pump Efficiency', 'VOI'] = df.loc[df['Sweep'] == 'Pump Efficiency', 'VOI'] * -1
    print(df.head(21))
    print('\n\n')
    return df



def plot_voi(df, filename=None):
    fig, ax = plt.subplots(figsize=(7,9))
    sns.boxplot(
        df,
        x='VOI',
        y='Sweep',
        hue='Case',
        gap=0.1,
        fliersize=0,
        # whis=[20, 80],
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(df['Sweep'].unique())-0.5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylabel('')
    ax.set_xlabel("VOI (%$_{\Delta LCOW}$/%$_{\Delta performance}$)", fontsize=16)
    color_bands = ["#d1d1d1", "white"]
    for i in range(0, len(ax.get_yticks())):
        ax.axhspan(i - 0.5, i + 0.5, facecolor=color_bands[i%2], alpha=0.75)

    # Create a annotated label for each group
    subgroup = 0
    for group in df['Group'].unique():
        subgroups = len(df.loc[df['Group'] == group, "Sweep"].unique())
        subgroup += subgroups
        print(group, subgroup)
        ax.plot([-1.5,1],[subgroup-0.5,subgroup-0.5], color="k", clip_on=False)

    plt.legend(
        bbox_to_anchor=(-0.3, 1.1), 
        loc='upper left',
        ncol=len(df['Case'].unique()), 
        frameon=False, 
        fontsize=14,
        alignment='left',
        handlelength=1, 
        handleheight=1)
    
    plt.tight_layout()
    # plt.savefig(f'{filename}.svg', format='svg')
    # plt.savefig(f'{filename}.png', format='png', dpi=300)
    plt.show()

def create_KBHDP():
    df1 = fetch_file(file = "diff_sweep_analysisType_RPT1_diff_analysis.h5", yaml_file="KBHDP_RPT_1_diff.yaml", case = "RPT1")
    df2 = fetch_file(file = "diff_sweep_analysisType_RPT2_diff_analysis.h5", yaml_file="KBHDP_RPT_2_diff.yaml", case = "RPT2")
    df3 = fetch_file(file = "diff_sweep_analysisType_RPT3_diff_analysis.h5", yaml_file="KBHDP_RPT_3_diff.yaml", case = "RPT3")
    # df4 = fetch_file(file = "diff_sweep_analysisType_RPT1_diff_analysis.h5", yaml_file="KBHDP_RPT_1_diff.yaml", case = "ZLD")
    df = pd.concat([df1, df2, df3])
    # df = pd.concat([df1, df2, df3, df4])
    # df = df3
    df = df.sort_values(by=['Group', 'Case', 'Sweep'])


    df['z-score'] = 0
    for sweep in df['Sweep'].unique():
        for case in df['Case'].unique():
            df.loc[(df['Sweep'] == sweep) & (df['Case'] == case), 'z-score'] = np.abs(stats.zscore(df.loc[(df['Sweep'] == sweep) & (df['Case'] == case), 'VOI']))

    df.drop(df[df['z-score'] > 2].index, inplace=True)
    print(df.head(10))

    df.to_csv('KBHDP_VOI.csv')
    plot_voi(df, filename='KBHDP_VOI')

def create_Permian():
    # df1 = fetch_file(file = "diff_sweep_analysisType_Permian_ST1_diff_analysis.h5", yaml_file="Permian_ST_1_diff.yaml", case = "ST1")
    df2 = fetch_file(file = "/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/sweep_results/archive/diff_sweep_analysisType_Permian_ST1_diff_analysis_20250427-110116.h5", yaml_file="Permian_ST_2_diff.yaml", case = "RPT2")
    # df2 = fetch_file(file = "diff_sweep_analysisType_RPT2_diff_analysis.h5", yaml_file="KBHDP_RPT_2_diff.yaml", case = "RPT2")
    # df3 = fetch_file(file = "diff_sweep_analysisType_RPT3_diff_analysis.h5", yaml_file="KBHDP_RPT_3_diff.yaml", case = "RPT3")
    # df4 = fetch_file(file = "diff_sweep_analysisType_RPT1_diff_analysis.h5", yaml_file="KBHDP_RPT_1_diff.yaml", case = "ZLD")
    # df = pd.concat([df1])
    df = pd.concat([df2])
    # df = pd.concat([df1, df2, df3, df4])
    # df = df3
    df = df.sort_values(by=['Group', 'Case', 'Sweep'])


    df['z-score'] = 0
    for sweep in df['Sweep'].unique():
        for case in df['Case'].unique():
            df.loc[(df['Sweep'] == sweep) & (df['Case'] == case), 'z-score'] = np.abs(stats.zscore(df.loc[(df['Sweep'] == sweep) & (df['Case'] == case), 'VOI']))

    df.drop(df[df['z-score'] > 2].index, inplace=True)
    print(df.head(10))

    df.to_csv('Permian_VOI.csv')
    plot_voi(df, filename='Permian_VOI')

def main():
    # create_Permian()
    create_KBHDP()

if __name__ == "__main__":
    main()