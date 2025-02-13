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

import json
from os.path import join, dirname
from math import floor, ceil, isnan
import numpy as np
import pandas as pd
import time
import multiprocessing
from itertools import product
import matplotlib.pyplot as plt
import PySAM.TroughPhysicalIph as iph
import os

# TODO:
# Annual energy for year 1 is a little different than that calculated in SAM


def read_module_datafile(file_name):
    with open(file_name, "r") as file:
        data = json.load(file)
    return data


def load_config(modules, file_names=None, module_data=None):
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


def setup_model(
    model_name,
    weather_file=None,
    weather_data=None,
    config_files=None,
    config_data=None,
):
    tech_model = iph.new()
    modules = [tech_model]

    load_config(modules, config_files, config_data)
    if weather_file is not None:
        tech_model.Weather.file_name = weather_file
    elif weather_data is not None:
        tech_model.Weather.solar_resource_data = weather_data
    else:
        raise Exception("Either weather_file or weather_data must be specified.")

    return {
        "tech_model": tech_model,
    }


def run_model(modules, heat_load=None, hours_storage=None):
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
        0 if isnan(freeze_protection_field) else freeze_protection_field
    )  # occasionally seen to be nan
    freeze_protection_tes = tech_model.Outputs.annual_tes_freeze_protection
    freeze_protection_tes = 0 if isnan(freeze_protection_tes) else freeze_protection_tes
    freeze_protection = freeze_protection_field + freeze_protection_tes
    annual_energy = (
        tech_model.Outputs.annual_energy - freeze_protection
    )  # [kWht] net, does not include that used for freeze protection
    capacity_factor = (
        annual_energy / (tech_model.value("q_pb_design") * 1e3 * 8760) * 100
    )  # [%]

    electrical_load = tech_model.Outputs.annual_electricity_consumption  # [kWhe]
    solar_multiplier = tech_model.Outputs.solar_mult
    total_aperture_reflective_area = tech_model.Outputs.total_aperture  # [m2]
    nloops = tech_model.Outputs.nLoops

    return {
        "annual_energy": annual_energy,  # [kWh] annual net thermal energy in year 1
        "freeze_protection": freeze_protection,  # [kWht]
        "capacity_factor": capacity_factor,  # [%] capacity factor
        "electrical_load": electrical_load,  # [kWhe]
        "solar_multiplier": solar_multiplier,
        "total_aperture_area": total_aperture_reflective_area,
        "number_loops": nloops,
    }


def setup_and_run(
    model_name, weather_file, config_data, heat_load
):  # , hours_storage):
    modules = setup_model(
        model_name, weather_file=weather_file, config_data=config_data
    )
    result = run_model(modules, heat_load)  # S, hours_storage)
    return result


#########################################################################################################
if __name__ == "__main__":
    # Model name is not relevant in WaterTAP-REFLO package because cost is calculated using REFLO costing packages
    model_name = "TroughPhysicalIph_PhysicalTroughIPHLCOHCalculator"
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )

    config_files = [
        os.path.join(__location__, "cst/trough_physical_iph-reflo.json"),
    ]
    weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")
    dataset_filename = os.path.join(
        __location__,
        "cst/trough_kbhdp_heat_load_1_100_hours_storage_24_T_loop_out_300.pkl",
    )  # output dataset for surrogate training

    config_data = [read_module_datafile(config_file) for config_file in config_files]
    del config_data[0]["file_name"]  # remove weather filename

    # Run parametrics via multiprocessing
    data = []
    heat_loads = np.linspace(1, 100, 200)  # [MWt]
    # hours_storages = np.linspace(20, 24, 5)  # [hr]
    arguments = list(product(heat_loads))  # , hours_storages))
    df = pd.DataFrame(arguments, columns=["heat_load"])  # , "hours_storage"])

    time_start = time.process_time()
    with multiprocessing.Pool(processes=6) as pool:
        args = [(model_name, weather_file, config_data, *args) for args in arguments]
        results = pool.starmap(setup_and_run, args)
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
                    "solar_multiplier",
                    "total_aperture_area",
                    "number_loops",
                ]
            ],
        ],
        axis=1,
    )
    df.to_pickle(dataset_filename)

    pass
