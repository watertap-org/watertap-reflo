import os
from psPlotKit.data_manager.ps_data_manager import psDataManager
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(os.path.dirname(parent_dir), "sweep_yamls")
sweep_results_dir = os.path.join(os.path.dirname(parent_dir), "sweep_results", "output")
save_path = os.path.join(os.path.dirname(parent_dir), "figures/")

print("sweep_yaml_dir: ", sweep_yaml_dir)
print("sweep_results_dir: ", sweep_results_dir)
print("save_path: ", save_path)

def plot_voi(df):
    fig, ax = plt.subplots(figsize=(7, 9))
    sns.boxplot(
        df,
        x="VOI",
        y="Sweep",
        hue='Case',
        gap=0.1,
        fliersize=0,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(df["Sweep"].unique()) - 0.5)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_ylabel("")
    ax.set_xlabel("VOI (%$_{\Delta LCOW}$/%$_{\Delta performance}$)", fontsize=16)
    color_bands = ["#d1d1d1", "white"]
    for i in range(0, len(ax.get_yticks())):
        ax.axhspan(i - 0.5, i + 0.5, facecolor=color_bands[i % 2], alpha=0.75)

    plt.legend(
        bbox_to_anchor=(-0.3, 1.1),
        loc="upper left",
        ncol=len(df["Sweep"].unique()),
        frameon=False,
        fontsize=14,
        alignment="left",
        handlelength=1,
        handleheight=1,
    )

    plt.tight_layout()
    plt.show()

def get_case_VOIs(yaml_file, h5_file):
    """
    Extracts the VOI from the data manager based on the provided YAML file.
    """
    yaml_file = os.path.join(sweep_yaml_dir, "Permian_RPT_2_diff.yaml")
    h5_file = os.path.join(sweep_results_dir, "diff_sweep_analysisType_Permian_RPT2_diff_analysis.h5")
    data_manager = psDataManager(h5_file)

    VOI = data_manager.get_voi(yaml_file)

    frames = [pd.DataFrame({"VOI": VOI[key]["VOI"], "Sweep": key}) for key in VOI.keys()]
    df = pd.concat(frames, ignore_index=True)

    return df

def plot_permian_voi():
    """
    Plots the VOI for the Permian case.
    """

    cases = {
        'RPT1': {
            'yaml_file': os.path.join(sweep_yaml_dir, "Permian_RPT_1_diff.yaml"),
            'h5_file': os.path.join(sweep_results_dir, "diff_sweep_analysisType_Permian_ST1_diff_analysis.h5")
        },
        'RPT2': {
            'yaml_file': os.path.join(sweep_yaml_dir, "Permian_RPT_2_diff.yaml"),
            'h5_file': os.path.join(sweep_results_dir, "diff_sweep_analysisType_Permian_RPT2_diff_analysis.h5")
        },
        # 'ZLD1': {},
        # 'ZLD2': {},
    }

    frames = []
    for case in cases:
        frame = get_case_VOIs(yaml_file, h5_file)
        frame["Case"] = case
        frames.append(frame)
    df = pd.concat(frames)

    return df



yaml_file = os.path.join(sweep_yaml_dir, "Permian_RPT_2_diff.yaml")
h5_file = os.path.join(sweep_results_dir, "diff_sweep_analysisType_Permian_RPT2_diff_analysis.h5")

df = plot_permian_voi()
plot_voi(df)