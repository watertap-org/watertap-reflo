import json
import os
from math import floor, ceil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
import time
# import PySAM.Pvsamv1 as pv
# import PySAM.Grid as grid
# import PySAM.Utilityrate5 as utilityrate
# import PySAM.Singleowner as singleowner
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP

absolute_path = os.path.dirname(__file__)
print(absolute_path)

tech_config_file = "/pvsamv1.json"
tech_config_file = absolute_path + tech_config_file
grid_config_file = "/grid.json"
grid_config_file = absolute_path + grid_config_file
rate_config_file = "/utilityrate5.json"
rate_config_file = absolute_path + rate_config_file
cash_config_file = "/singleowner.json"
cash_config_file = absolute_path + cash_config_file
weather_file = "/phoenix_az_33.450495_-111.983688_psmv3_60_tmy.csv"
weather_file = absolute_path + weather_file

dataset_filename = os.path.join(
        os.path.dirname(__file__), "dataset.pkl"
    )  # output dataset for surrogate training

pysam = PySAMWaterTAP(
        pysam_model="pv",
        tech_config_file=tech_config_file,
        grid_config_file=grid_config_file,
        rate_config_file=rate_config_file,
        cash_config_file=cash_config_file,
        weather_file=weather_file,
    )

def run(desired_size):  
    pysam.run_pv_single_owner(desired_size=desired_size)
    results = {
        "design_size": desired_size, # [kWh]
        "annual_energy": pysam.annual_energy,  # [kWh] annual net thermal energy in year 1
        # "hourly_energy": np.mean(pysam.hourly_energy),  # [kWht]
        # 'capital_cost': pysam.cash_model.Outputs.capital_cost
        }
    return pd.DataFrame.from_dict([results])

def plot_2D(frame, x_var='design_size', y_var ='annual_energy', x_label='Desired Size', y_label='Annual Energy'):
    fig, ax = plt.subplots(figsize=(6,6))
    plt.scatter(x=frame[x_var], y=frame[y_var])

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.tick_params(axis='x', labelsize = 18)
    ax.tick_params(axis='y', labelsize = 18)
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.show()

if __name__ == "__main__":
    # Run parametrics via multiprocessing
    data = []
    pv_size = np.linspace(10,10000,1000)
    df = pd.DataFrame()

    time_start = time.process_time()
    with multiprocessing.Pool(processes=8) as pool:
        results = pool.map(run, pv_size)
    time_stop = time.process_time()
    print("Multiprocessing time:", time_stop - time_start, "\n")

    df = pd.concat(results)

    df.to_pickle(dataset_filename)
    print(f'Data written to {dataset_filename}')
    plot_2D(df)