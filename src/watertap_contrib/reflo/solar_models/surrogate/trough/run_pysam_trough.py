import os
import json
import time
import multiprocessing

import numpy as np
import pandas as pd
from math import isnan
from itertools import product
import PySAM.TroughPhysicalIph as iph

# Documentation for PySAM module used in this script:
# https://nrel-pysam.readthedocs.io/en/main/modules/TroughPhysicalIph.html

__all__ = [
    "setup_model_trough",
    "run_model_trough",
    "setup_and_run_trough",
    "generate_trough_data",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
# Tuscon, AZ
weather_file_default = os.path.join(__location__, "data/test_trough_weather_data.csv")
config_file_default = os.path.join(__location__, "data/trough_physical_iph-reflo.json")
# NOTE:
# The initial version of this PySAM run and the config file was created using nrel-pysam version 6.0.0;
# version 7.0.0 requires the parameter hs_htf_mdot_max_frac be specified but it is not included in the 6.0.0 config file.
# This parameter was manually added to trough_physical_iph-reflo.json with a default value of 1.1.
# If you are using a different version of nrel-pysam, you may need to adjust the config file accordingly.
# You could also use default configurations available here:
# https://nrel-pysam.readthedocs.io/en/latest/sam-configurations.html


def setup_model_trough(
    weather_file=None,
    config_file=None,
):
    """
    Create the PySAM technology model for Trough industrial process heat (IPH) system.
    """

    if weather_file is None:
        raise RuntimeError(
            "Weather file must be specified for trough PySAM model setup. "
            "Use the 'weather_file' argument to specify the path to the weather data file."
        )

    if config_file is None:
        raise RuntimeError(
            "Configuration file must be specified for trough PySAM model setup. "
            "Use the 'config_file' argument to specify the path to the configuration data file."
        )

    tech_model = iph.new()

    with open(config_file, "r") as file:
        config_data = json.load(file)

    for k, v in config_data.items():
        if k != "number_inputs":
            try:
                tech_model.value(k, v)
            except AttributeError:
                print(f"Warning: {k} not found in the technology model. Skipping.")

    tech_model.Weather.file_name = weather_file

    return tech_model


def run_model_trough(
    tech_model,
    system_capacity=None,
    hours_storage=None,
    temperature_loop=300,
    return_tech_model=False,
):
    """
    Runs the trough model with specified heat load and storage hours.

    :param tech_model: PySAM TroughPhysicalIph model object
    :param system_capacity: System capacity in MWt to be supplied by the trough system
    :param hours_storage: Number of hours of thermal storage (optional)
    :param temperature_loop: Loop outlet temperature in Celsius (default is 300 C)
    :param return_tech_model: If True, returns the tech_model object along with results
    """
    if system_capacity is None:
        raise ValueError("system_capacity must be specified for trough model run.")

    tech_model.value("q_pb_design", system_capacity)
    tech_model.value("T_loop_out", temperature_loop)  # [C] default is 300 C

    if hours_storage is not None:
        tech_model.value("tshours", hours_storage)

    print(
        f"\nRunning:"
        f"\n\tSystem Capacity = {system_capacity}"
        f"\n\tHours Storage = {hours_storage}"
        f"\n\tTemperature Hot = {temperature_loop}"
    )
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

    # [kWht] net, does not include that used for freeze protection
    heat_annual = tech_model.Outputs.annual_energy - freeze_protection
    capacity_factor = (
        heat_annual / (tech_model.value("q_pb_design") * 1e3 * 8760) * 100
    )  # [%]

    electricity_annual = tech_model.Outputs.annual_electricity_consumption  # [kWhe]
    solar_multiplier = tech_model.Outputs.solar_mult
    total_aperture_reflective_area = tech_model.Outputs.total_aperture  # [m2]
    nloops = tech_model.Outputs.nLoops

    results = {
        "heat_annual": heat_annual,  # [kWh] annual net thermal energy production in year 1
        "electricity_annual": electricity_annual,  # [kWhe] annual electricity consumption in year 1
        "freeze_protection": freeze_protection,  # [kWht]
        "capacity_factor": capacity_factor,  # [%] capacity factor
        "solar_multiplier": solar_multiplier,
        "total_aperture_area": total_aperture_reflective_area,
        "number_loops": nloops,
    }

    if return_tech_model:
        return results, tech_model
    else:
        return results


def setup_and_run_trough(
    weather_file, config_file, system_capacity, hours_storage, temperature_loop
):
    """
    Set up and run the PySAM TroughPhysicalIph model with specified parameters once.

    :param weather_file: Path to the weather data file
    :param config_file: Path to the configuration data file
    :param system_capacity: Heat load in MWt to be supplied by the trough system
    :param hours_storage: Number of hours of thermal storage (optional)
    :param temperature_loop: Loop outlet temperature in Celsius (default is 300 C)
    :return: Dictionary containing results from the model run
    """
    tech_model = setup_model_trough(weather_file=weather_file, config_file=config_file)

    results = run_model_trough(
        tech_model,
        system_capacity=system_capacity,
        hours_storage=hours_storage,
        temperature_loop=temperature_loop,
    )
    return results


def generate_trough_data(
    system_capacities=np.geomspace(1, 50, 3),
    hours_storages=[24],
    temperatures_loop=[300],
    weather_file=weather_file_default,
    config_file=config_file_default,
    save_data=True,
    use_multiprocessing=True,
    dataset_filename=None,
):
    """
    Run the PySAM TroughPhysicalIph model for a range of heat loads, storage hours, and loop temperatures,
    and save the results to a DataFrame.
    This data could then be used to train a surrogate model.

    :param system_capacities: List of system capacities in MWt to be supplied by the trough system
    :param hours_storages: List of number of hours of thermal storage (optional)
    :param temperatures_loop: List of loop outlet temperatures in Celsius (default is 300 C)
    :param weather_file: Path to the weather data file
    :param config_file: Path to the configuration data file
    :param save_data: If True, saves the results to a pickle file
    :param use_multiprocessing: If True, uses multiprocessing to run the model in parallel
    :param dataset_filename: Path to the output dataset file (if save_data is True)
    :return: DataFrame containing the results from the model runs
    """
    if dataset_filename is None:
        # assume it is run for testing purposes
        dataset_filename = os.path.join(__location__, "data/test_data.pkl")

    tech_model = setup_model_trough(weather_file=weather_file, config_file=config_file)

    if use_multiprocessing:
        combos = list(product(system_capacities, hours_storages, temperatures_loop))
        df = pd.DataFrame(
            combos, columns=["system_capacity", "hours_storage", "temperature_loop"]
        )
        time_start = time.process_time()
        with multiprocessing.Pool(processes=6) as pool:
            args_in = [(weather_file, config_file, *combo) for combo in combos]
            results = pool.starmap(setup_and_run_trough, args_in)
        time_stop = time.process_time()
        # print("Multiprocessing time:", time_stop - time_start, "\n")
        df_results = pd.DataFrame(results)
        df = pd.concat(
            [
                df,
                df_results[
                    [
                        "heat_annual",
                        "electricity_annual",
                        "freeze_protection",
                        "capacity_factor",
                        "solar_multiplier",
                        "total_aperture_area",
                        "number_loops",
                    ]
                ],
            ],
            axis=1,
        )
    else:
        data = []
        combos = list(product(system_capacities, hours_storages, temperatures_loop))
        for system_capacity, hours_storage, temperature_loop in combos:
            results = run_model_trough(
                tech_model,
                system_capacity=system_capacity,
                hours_storage=hours_storage,
                temperature_loop=temperature_loop,
            )
            data.append(
                [
                    system_capacity,
                    hours_storage,
                    temperature_loop,
                    results["heat_annual"],
                    results["electricity_annual"],
                    results["freeze_protection"],
                    results["capacity_factor"],
                    results["solar_multiplier"],
                    results["total_aperture_area"],
                    results["number_loops"],
                ]
            )
        df = pd.DataFrame(
            data,
            columns=[
                "system_capacity",
                "hours_storage",
                "temperature_loop",
                "heat_annual",
                "electricity_annual",
                "freeze_protection",
                "capacity_factor",
                "solar_multiplier",
                "total_aperture_area",
                "number_loops",
            ],
        )
    if save_data:
        df.to_pickle(dataset_filename)
    return df


if __name__ == "__main__":
    df = generate_trough_data(use_multiprocessing=True)
    print(df.head())
    os.remove(os.path.join(__location__, "data/test_data.pkl"))
