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

__all__ = [
    "read_module_datafile",
    "load_pysam_trough_config",
    "setup_pysam_trough_model",
    "run_pysam_trough_model",
    "setup_and_run_trough",
    "run_pysam_kbhdp_trough_sweep",
]

model_name = "PhysicalTroughIPHLCOHCalculator"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

config_files = [
    os.path.join(__location__, "cst/trough_physical_process_heat-reflo.json"),
]
weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")


def read_module_datafile(file_name):
    with open(file_name, "r") as file:
        data = json.load(file)
    return data


def load_pysam_trough_config(modules, file_names=None, module_data=None):
    """
    Loads parameter values into PySAM modules, either from files or supplied dicts

    :param modules: List of PySAM modules
    :param file_names: List of JSON file paths containing parameter values for respective modules
    :param module_data: List of dictionaries containing parameter values for respective modules

    :returns: no return value
    """
    for i in range(len(modules)):
        if file_names is not None:
            assert len(file_names) == len(modules)
            data = read_module_datafile(file_names[i])
        elif module_data is not None:
            assert len(module_data) == len(modules)
            data = module_data[i]
        else:
            raise Exception("Either file_names or module_data must be assigned.")

        missing_values = []  # for debugging
        for k, v in data.items():
            if k != "number_inputs":
                try:
                    modules[i].value(k, v)
                except:
                    missing_values.append(k)
        pass


def tes_cost(tech_model):
    storage_cost_specific = 62  # [$/kWht] borrowed from physical power trough
    tes_thermal_capacity = (
        tech_model.value("q_pb_design") * 1e3 * tech_model.value("tshours")
    )  # [kWht]
    return tes_thermal_capacity * storage_cost_specific


def system_capacity(tech_model):
    return (
        tech_model.value("q_pb_design")
        * tech_model.value("specified_solar_multiple")
        * 1e3
    )  # [kW]


def setup_pysam_trough_model(
    model_name,
    weather_file=None,
    weather_data=None,
    config_files=None,
    config_data=None,
):
    tech_model = iph.new()
    modules = [tech_model]

    load_pysam_trough_config(modules, config_files, config_data)
    if weather_file is not None:
        tech_model.Weather.file_name = weather_file
    elif weather_data is not None:
        tech_model.Weather.solar_resource_data = weather_data
    else:
        raise Exception("Either weather_file or weather_data must be specified.")

    return {
        "tech_model": tech_model,
    }


def run_pysam_trough_model(
    modules,
    heat_load=None,
    hours_storage=None,
    return_tech_model=False,
):
    tech_model = modules["tech_model"]

    if heat_load is not None:
        tech_model.value("q_pb_design", heat_load)
    if hours_storage is not None:
        tech_model.value("tshours", hours_storage)

    tech_model.execute()

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
        tech_model.Outputs.annual_energy - freeze_protection
    )  # [kWht] net, does not include that used for freeze protection
    capacity_factor = (
        annual_energy / (tech_model.value("q_pb_design") * 1e3 * 8760) * 100
    )  # [%]
    electrical_load = tech_model.Outputs.annual_electricity_consumption  # [kWhe]
    print(f"Running:")
    print(f"\tHeat Load = {heat_load} MW")
    print(f"\tHours Storage = {hours_storage} hr")
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


def setup_and_run_trough(model_name, weather_file, config_data, heat_load, hours_storage):
    modules = setup_pysam_trough_model(
        model_name, weather_file=weather_file, config_data=config_data
    )
    result = run_pysam_trough_model(modules, heat_load, hours_storage)
    return result


def run_pysam_kbhdp_trough_sweep(
    heat_loads=np.linspace(1, 25, 25),
    hours_storages=np.linspace(1, 24, 24),
    processes=6,
    save_data=True,
    dataset_filename="",
):

    config_data = [read_module_datafile(config_file) for config_file in config_files]
    del config_data[0]["file_name"]  # remove weather filename

    # output dataset for surrogate training
    dataset_filename = os.path.join(os.path.dirname(__file__), dataset_filename)

    arguments = list(itertools.product(heat_loads, hours_storages))
    df = pd.DataFrame(arguments, columns=["heat_load", "hours_storage"])

    time_start = time.process_time()
    with multiprocessing.Pool(processes=processes) as pool:
        args = [(model_name, weather_file, config_data, *args) for args in arguments]
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

    heat_loads = np.linspace(1, 10, 3)  # [MWt]
    hours_storages = np.linspace(0, 24, 3)  # [hr]
    # dataset_filename = f"trough_data_heat_load_{int(min(heat_loads))}_{int(max(heat_loads))}_hours_storage_{int(min(hours_storages))}_{int(max(hours_storages))}.pkl"
    dataset_filename = "test.pkl"

    run_pysam_kbhdp_trough_sweep(
        heat_loads=heat_loads,
        hours_storages=hours_storages,
        dataset_filename=dataset_filename,
    )
