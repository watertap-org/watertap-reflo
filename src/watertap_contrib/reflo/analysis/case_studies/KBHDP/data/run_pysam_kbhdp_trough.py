#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import os
import json
import time
import math
import multiprocessing
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import PySAM.TroughPhysicalProcessHeat as iph

__author__ = "Kurban Sitterley"

__all__ = [
    "read_trough_module_datafile",
    "setup_pysam_trough_model",
    "run_pysam_trough_model",
    "setup_and_run_trough",
    "run_pysam_kbhdp_trough_sweep",
]

model_name = "PhysicalTroughIPHLCOHCalculator"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
param_file = os.path.join(__location__, "cst/trough_physical_process_heat-kbhdp.json")
weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")


def read_trough_module_datafile(file_name):
    with open(file_name, "r") as file:
        data = json.load(file)
    return data


def setup_pysam_trough_model(weather_file=None, config_data=None):

    tech_model = iph.new()

    for k, v in config_data.items():
        try:
            tech_model.value(k, v)
        except:
            pass

    if weather_file is not None:
        tech_model.Weather.file_name = weather_file

    return tech_model


def run_pysam_trough_model(
    tech_model,
    heat_load=None,
    hours_storage=None,
    return_tech_model=False,
    temperature_hot=200,
    verbosity=0,  # 0 - no PySAM log; 1 - include PySAM log
):
    # Trough Physical SAM Documentation:
    # https://nrel-pysam.readthedocs.io/en/main/modules/TroughPhysicalIph.html

    if heat_load is None:
        raise RuntimeError("heat_load must be specified.")
    # From SAM inputs browser: q_pb_design = Heat sink thermal power
    tech_model.value("q_pb_design", heat_load)
    if hours_storage is None:
        raise RuntimeError("hours_storage must be specified.")
    # From SAM inputs browswer: tshours = Hours of storage at design point
    tech_model.value("tshours", hours_storage)
    if temperature_hot is None:
        raise RuntimeError("temperature_hot must be specified.")
    # Target loop outlet temperature [C]
    # tech_model.value("T_loop_out", temperature_hot)

    print(f"Running:")
    print(f"\tHeat Load = {heat_load} MW")
    print(f"\tHours Storage = {hours_storage} hr")
    tech_model.execute(verbosity)

    # NOTE: freeze_protection_field can sometimes be nan (when it should be 0) and this causes other nan's
    #  Thus, freeze_protection, annual_energy and capacity_factor must be calculated manually
    # annual_energy = tech_model.Outputs.annual_energy                            # [kWht] net, does not include that used for freeze protection
    # freeze_protection = tech_model.Outputs.annual_thermal_consumption           # [kWht]
    # capacity_factor = tech_model.Outputs.capacity_factor                        # [%]
    freeze_protection_field = tech_model.Outputs.annual_field_freeze_protection
    freeze_protection_field = (
        0 if math.isnan(freeze_protection_field) else freeze_protection_field
    )  # occasionally seen to be nan
    freeze_protection_tes = tech_model.Outputs.annual_tes_freeze_protection
    freeze_protection_tes = (
        0 if math.isnan(freeze_protection_tes) else freeze_protection_tes
    )
    freeze_protection = freeze_protection_field + freeze_protection_tes
    annual_energy = (
        sum(tech_model.Outputs.gen) - freeze_protection
    )  # [kWht] net, does not include that used for freeze protection
    annual_energy_2 = (
        tech_model.Outputs.annual_energy - freeze_protection
    )  # [kWht] net, does not include that used for freeze protection
    capacity_factor = (
        annual_energy / (tech_model.value("q_pb_design") * 1e3 * 8760) * 100
    )  # [%]
    electrical_load = tech_model.Outputs.annual_electricity_consumption  # [kWhe]
    print(f"\tAnnual Heat Delivered {annual_energy_2:.2f} kWh")
    print(f"\tAnnual Heat Delivered {annual_energy:.2f} kWh")
    results = {
        "annual_energy": annual_energy,  # [kWh] annual net thermal energy in year 1
        "freeze_protection": freeze_protection,  # [kWht]
        "capacity_factor": capacity_factor,  # [%] capacity factor
        "electrical_load": electrical_load,  # [kWhe]
    }
    if return_tech_model:
        return results, tech_model
    else:
        return results


def setup_and_run_trough(weather_file, config_data, heat_load, hours_storage):
    model_name = "PhysicalTroughIPHLCOHCalculator"
    tech_model = setup_pysam_trough_model(
        weather_file=weather_file, config_data=config_data
    )
    result = run_pysam_trough_model(tech_model, heat_load, hours_storage)
    return result


def run_pysam_kbhdp_trough_sweep(
    heat_loads=np.linspace(1, 25, 25),
    hours_storages=np.linspace(1, 24, 24),
    processes=6,
    save_data=True,
    dataset_filename="",
):

    config_data = read_trough_module_datafile(param_file)

    # output dataset for surrogate training
    dataset_filename = os.path.join(os.path.dirname(__file__), dataset_filename)

    arguments = list(itertools.product(heat_loads, hours_storages))
    df = pd.DataFrame(arguments, columns=["heat_load", "hours_storage"])

    time_start = time.process_time()
    with multiprocessing.Pool(processes=processes) as pool:
        args = [(weather_file, config_data, *args) for args in arguments]
        results = pool.starmap(setup_and_run_trough, args)
    time_stop = time.process_time()
    print("Multiprocessing time:", time_stop - time_start, "\n")

    df_results = pd.DataFrame(results)

    df = pd.concat(
        [
            df,
            df_results[
                [
                    "annual_energy",
                    "freeze_protection",
                    "capacity_factor",
                    "electrical_load",
                ]
            ],
        ],
        axis=1,
    )
    if save_data:
        df.to_pickle(dataset_filename)


if __name__ == "__main__":

    config_data = read_trough_module_datafile(param_file)

    tech_model = setup_pysam_trough_model(
        weather_file=weather_file, config_data=config_data
    )

    result, tech_model = run_pysam_trough_model(
        tech_model, 8.12, 24, return_tech_model=True
    )
    result, tech_model = run_pysam_trough_model(
        tech_model, 1, 6, return_tech_model=True
    )
    # result, tech_model = run_pysam_trough_model(tech_model, 50, 12, return_tech_model=True)
    # result, tech_model = run_pysam_trough_model(tech_model, 100, 24, return_tech_model=True)

    # heat_loads = np.linspace(1, 10, 3)  # [MWt]
    # hours_storages = np.linspace(0, 24, 3)  # [hr]
    # # dataset_filename = f"trough_data_heat_load_{int(min(heat_loads))}_{int(max(heat_loads))}_hours_storage_{int(min(hours_storages))}_{int(max(hours_storages))}.pkl"
    # dataset_filename = "test.pkl"

    # run_pysam_kbhdp_trough_sweep(
    #     heat_loads=heat_loads,
    #     hours_storages=hours_storages,
    #     dataset_filename=dataset_filename,
    # )
