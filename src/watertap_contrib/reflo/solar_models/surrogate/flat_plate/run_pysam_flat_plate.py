import os
import json
import time
import multiprocessing

import numpy as np
import pandas as pd
from itertools import product
import PySAM.Swh as swh

__all__ = [
    "setup_model_fpc",
    "run_model_fpc",
    "setup_and_run_fpc",
    "generate_fpc_data",
]

# Parameters used in the model
# These are the same as in the SAM SolarWaterHeating model
cp_water = 4.181  # [kJ/kg-K]
density_water = 1000  # [kg/m3]
pump_power_per_collector = 45 / 2  # [W]
pipe_length_fixed = 9  # [m]
pipe_length_per_collector = 0.5  # [m]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
weather_file_default = os.path.join(__location__, "data/test_fpc_weather_data.csv")

# see: https://nrel-pysam.readthedocs.io/en/latest/modules/Swh.html


def system_capacity_computed(tech_model):
    """
    Computes the system capacity in kW

    Equation taken from SAM UI
    """
    system_capacity = (
        tech_model.value("area_coll")
        * tech_model.value("ncoll")
        * (tech_model.value("FRta") - tech_model.value("FRUL") * 30 / 1000)
    )
    return system_capacity


def setup_model_fpc(
    temperatures,
    weather_file=None,
    config_file=None,
):
    """
    Create the PySAM technology model for flat plate collectors (FPC).
    """

    if weather_file is None:
        raise RuntimeError(
            "Weather file must be specified for FPC model setup. "
            "Use the 'weather_file' argument to specify the path to the weather data file."
        )

    if config_file is not None:
        # Load the SWH system with custom parameters
        tech_model = swh.new()
        with open(config_file, "r") as file:
            config_data = json.load(file)

        for k, v in config_data.items():
            try:
                tech_model.value(k, v)
            except AttributeError:
                print(f"Warning: {k} not found in the technology model. Skipping.")
    else:
        # Load the default SWH system
        tech_model = swh.default("SolarWaterHeatingCommercial")
        tech_model.value("solar_resource_file", weather_file)
        tech_model.value("azimuth", 180)  # [deg] south facing
        tech_model.value("tilt", 32)  # [deg] south facing
        tech_model.value("use_custom_mains", 1)
        tech_model.value("sky_model", 2)

    tech_model.value("solar_resource_file", weather_file)

    # Set constant temperatures
    tech_model.value("custom_mains", 8760 * (temperatures["T_cold"],))
    tech_model.value("T_set", temperatures["T_hot"])
    tech_model.value("T_room", temperatures["T_amb"])

    # Set collector loop mass flow (this should be done automatically in SAM)
    tech_model.value(
        "mdot", tech_model.value("test_flow") * tech_model.value("ncoll")
    )  # [kg/s]

    # Ensure system capacity parameter agreement
    system_capacity_actual = system_capacity_computed(tech_model)
    tech_model.value("system_capacity", system_capacity_actual)

    return tech_model


def run_model_fpc(
    tech_model,
    system_capacity_mwt=None,
    hours_storage=24,
    temperature_hot=80,
    temp_lb_thresh=0.5,
    temp_ub_thresh=None,
    scaled_draw_min=1,
    num_scaled_draw_pts=200,
    temp_frac=0.05,
    return_tech_model=False,
):
    """
    Run the PySAM technology model for flat plate collectors (FPC).

    NOTE: In this approach, it is assumed that the user wants the delivered temperature of the
    of the FPC system to be within a certain range around the desired hot temperature to minimize
    the need for electric heating. The delivered temperature is adjusted by changing the scaled draw
    until the delivered temperature is within the bounds defined by the hot temperature and the
    temperature bounds (temp_lb_thresh and temp_ub_thresh).

    :param tech_model: PySAM technology model
    :param system_capacity_mwt: desired system capacity [MWt]
    :param hours_storage: hours of storage to use in the simulation [hr]
    :param temperature_hot: desired temperature setpoint for the simulation [C]
    :param temp_lb_thresh: [C] lower bound threshold for delivered temperature
    :param temp_ub_thresh: [C] upper bound threshold for delivered temperature
    :param scaled_draw_min: [kg/hr] minimum scaled draw to start the search for
                            the delivered temperature within bounds
    :param num_scaled_draw_pts: number of points to search for the delivered temperature
    :param return_tech_model: if True, returns the technology model along with the results
    :param temp_frac: fraction of the temperature setpoint to use for the maximum tank temperature
                        (default is 0.05, which means the maximum tank temperature will be
                        temperature_hot * (1 + temp_frac))
    :return: dictionary with results including heat_annual, electricity_annual,
             aux_power_electric_heating, system_capacity_actual, scaled_draw, and temperature_delivered
    :rtype: dict
    :raises RuntimeError: if the delivered temperature is not within the bounds after the search

    """

    if system_capacity_mwt is None:
        raise ValueError("system_capacity_mwt must be specified for FPC model run.")

    system_capacity = system_capacity_mwt * 1e3  # convert MW to kW
    T_cold = tech_model.value("custom_mains")[0]  # [C]

    if temp_ub_thresh is None:
        temp_ub_thresh = temp_lb_thresh

    print(
        f"\nRunning:"
        f"\n\tSystem Capacity = {system_capacity_mwt}"
        f"\n\tHours Storage = {hours_storage}"
        f"\n\tTemperature Hot = {temperature_hot}"
    )

    if system_capacity is not None:
        # Set system capacity
        n_collectors = round(
            system_capacity
            / (
                tech_model.value("area_coll")
                * (tech_model.value("FRta") - tech_model.value("FRUL") * 30 / 1000)
            )
        )  # from SAM UI
        tech_model.value("ncoll", n_collectors)
        system_capacity_actual = system_capacity_computed(tech_model)
        tech_model.value("system_capacity", system_capacity_actual)  # [kW]
    if temperature_hot is not None:
        # Set hot outlet temperature
        tech_model.value("T_set", temperature_hot)
    if hours_storage is not None:
        # Set hours of storage (tank volume)
        hours_storage = max(hours_storage, 1e-3)  # don't accept 0 hours
        system_capacity = tech_model.value("system_capacity")  # [kW]
        T_hot = tech_model.value("T_set")  # [C]
        mass_tank_water = (
            hours_storage * 3600 * system_capacity / (cp_water * (T_hot - T_cold))
        )  # [kg]
        volume_tank = mass_tank_water / density_water  # [m3]
        tech_model.value("V_tank", volume_tank)

    # Set collector loop and hot water mass flow rates
    mdot_collectors = tech_model.value("test_flow") * tech_model.value("ncoll")
    tech_model.value("mdot", mdot_collectors)  # [kg/s]
    T_hot = tech_model.value("T_set")  # [C]
    mdot_hot = (
        tech_model.value("system_capacity") / (cp_water * (T_hot - T_cold)) * 3600
    )  # [kg/hr]
    T_tank_max = T_hot * (1 + temp_frac)
    if T_tank_max > 99:
        T_tank_max = 99
    tech_model.value(
        "T_tank_max", T_tank_max
    )  # [C], max tank temperature is the temperature we want

    # Set pipe diameter and pump power
    pipe_length = pipe_length_fixed + pipe_length_per_collector * tech_model.value(
        "ncoll"
    )
    tech_model.value("pipe_length", pipe_length)  # [m] default is 0.019 m
    pumping_power = pump_power_per_collector * tech_model.value("ncoll")
    tech_model.value("pump_power", pumping_power)  # [W]

    assert tech_model.value("T_set") == temperature_hot

    sd = scaled_draw_min
    lb = temperature_hot - temp_lb_thresh
    ub = temperature_hot + temp_ub_thresh
    increment = mdot_hot / num_scaled_draw_pts
    num_runs = 0

    for sd in np.linspace(scaled_draw_min, mdot_hot, num_scaled_draw_pts):
        tech_model.value("scaled_draw", 8760 * (sd,))
        tech_model.execute()
        temp_delivered = np.mean(tech_model.Outputs.T_deliv)
        num_runs += 1

        if temp_delivered < lb:
            break
        if lb < temp_delivered < ub:
            break

    if not lb < temp_delivered < ub:
        # scaled draw interval was probably high
        # rerun again with 10x more points
        for sd in np.linspace(scaled_draw_min, mdot_hot, num_scaled_draw_pts * 100):

            tech_model.value("scaled_draw", 8760 * (sd,))
            tech_model.execute()
            temp_delivered = np.mean(tech_model.Outputs.T_deliv)
            num_runs += 1

            if temp_delivered < lb:
                break
            if lb < temp_delivered < ub:
                break

    if not lb < temp_delivered < ub:
        msg = f"Final design results in delivered temperature that is outside the bounds.\n"
        msg += f"For {system_capacity_mwt} MW, {hours_storage} hrs storage, {temperature_hot} C:"
        msg += f"Delivered temperature {temp_delivered:.2f} C is not between "
        msg += f"{lb:.2f} C and {ub:.2f} C with scaled_draw {sd:.2f} kg/hr.\n"
        msg += "Try setting more num_scaled_draw_pts and rerunning."
        raise RuntimeError(msg)

    # [kWh] does not include electric heat, includes losses
    heat_annual = tech_model.value("annual_Q_deliv")
    # [kWh] auxiliary power used for electric heating
    aux_power_annual = sum(tech_model.Outputs.Q_aux)
    # [kWh] includes pumping and auxiliary power
    electricity_annual = sum(tech_model.value("P_pump")) + sum(tech_model.Outputs.Q_aux)
    # [-] for analysis only, plant beneficial if < 1
    frac_electricity_annual = electricity_annual / heat_annual

    results = {
        "heat_annual": heat_annual,  # [kWh] annual net thermal energy in year 1
        "electricity_annual": electricity_annual,  # [kWhe]
        "aux_power_electric_heating": aux_power_annual,  # [kWhe]
        "system_capacity_actual": system_capacity_actual,
        "scaled_draw": np.mean(tech_model.Outputs.draw),
        "temperature_delivered": np.mean(tech_model.Outputs.T_deliv),
    }

    if return_tech_model:
        return results, tech_model
    else:
        return results


def setup_and_run_fpc(
    temperatures,
    weather_file,
    config_file,
    system_capacity,
    hours_storage,
    temperature_hot,
):
    """
    Setup and run the FPC model with the given parameters a single time.

    :param temperatures: dictionary with temperatures (T_cold, T_hot, T_amb)
    :param weather_file: path to the weather file
    :param config_file: configuration file for the model
    :param system_capacity: desired heat load in MWt
    :param hours_storage: hours of storage to use in the simulation
    :param temperature_hot: desired hot temperature in C
    :return: dictionary with results from the model run
    :rtype: dict
    """

    tech_model = setup_model_fpc(
        temperatures, weather_file=weather_file, config_file=config_file
    )
    result = run_model_fpc(
        tech_model,
        system_capacity_mwt=system_capacity,
        temperature_hot=temperature_hot,
        hours_storage=hours_storage,
    )

    return result


def generate_fpc_data(
    system_capacities=np.geomspace(1, 50, 3),
    hours_storages=[24],
    temperatures_hot=[80],
    temperature_cold=20,
    weather_file=weather_file_default,
    config_file=None,
    save_data=True,
    use_multiprocessing=False,
    number_processes=6,
    dataset_filename=None,
):
    """
    Run PySAM to collect data for FPC surrogate model.
    This function sets up the FPC model, runs it for various heat loads, hours of storage,
    and hot temperatures, and returns a DataFrame with the results.

    :param system_capacitys: list of heat loads in MWt
    :param hours_storages: list of hours of storage
    :param temperatures_hot: list of hot temperatures in C
    :param temperature_cold: cold temperature in C
    :param weather_file: path to the weather file
    :param config_file: configuration file for the model
    :param save_data: if True, saves the data to a pickle file
    :param use_multiprocessing: if True, uses multiprocessing to run the model
    :param dataset_filename: filename to save the data to, if None, uses default test data file
    :return: DataFrame with results
    :rtype: pd.DataFrame

    """

    temperatures = {
        "T_cold": temperature_cold,
        "T_hot": 70,  # this will be overwritten by temperature_hot value
        "T_amb": 18,
    }

    if dataset_filename is None:
        # assume it is run for testing purposes
        dataset_filename = os.path.join(__location__, "data/test_data.pkl")

    print(f"Saving data to {dataset_filename}")

    tech_model = setup_model_fpc(
        temperatures=temperatures,
        weather_file=weather_file,
        config_file=config_file,
    )

    # Run pysam
    if use_multiprocessing:
        combos = list(product(system_capacities, hours_storages, temperatures_hot))
        df = pd.DataFrame(
            combos, columns=["system_capacity", "hours_storage", "temperature_hot"]
        )

        time_start = time.process_time()
        with multiprocessing.Pool(processes=number_processes) as pool:
            args_in = [
                (temperatures, weather_file, config_file, *combo) for combo in combos
            ]
            results = pool.starmap(setup_and_run_fpc, args_in)
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
                        "aux_power_electric_heating",
                        "system_capacity_actual",
                        "scaled_draw",
                        "temperature_delivered",
                    ]
                ],
            ],
            axis=1,
        )

    else:
        # run model sequentially
        data = []
        combos = list(product(system_capacities, hours_storages, temperatures_hot))
        for system_capacity, hours_storage, temperature_hot in combos:
            result = run_model_fpc(
                tech_model, system_capacity, hours_storage, temperature_hot
            )
            data.append(
                [
                    system_capacity,
                    hours_storage,
                    temperature_hot,
                    result["heat_annual"],
                    result["electricity_annual"],
                    result["aux_power_electric_heating"],
                    result["system_capacity_actual"],
                    result["scaled_draw"],
                    result["temperature_delivered"],
                ]
            )
        df = pd.DataFrame(
            data,
            columns=[
                "system_capacity",
                "hours_storage",
                "temperature_hot",
                "heat_annual",
                "electricity_annual",
                "aux_power_electric_heating",
                "system_capacity_actual",
                "scaled_draw",
                "temperature_delivered",
            ],
        )

    if save_data:
        df.to_pickle(dataset_filename)

    return df


if __name__ == "__main__":

    # system_capacities = np.linspace(1, 50, 10)
    # hours_storages = [12, 24]
    # temperatures_hot = [60, 80]
    # df = generate_fpc_data(
    #     system_capacities=system_capacities,
    #     hours_storages=hours_storages,
    #     temperatures_hot=temperatures_hot,
    #     temperature_cold=20,
    #     weather_file=weather_file_default,
    #     save_data=True,
    #     use_multiprocessing=True,
    #     dataset_filename=os.path.join(__location__, "data/test_data-with-defaults.pkl"),
    # )

    df = generate_fpc_data()
    print(df.head(20))
    os.remove(os.path.join(__location__, "data/test_data.pkl"))
