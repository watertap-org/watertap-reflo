#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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

import numpy as np
import pandas as pd
import time
import multiprocessing
from itertools import product
import matplotlib.pyplot as plt
import PySAM.Swh as swh

__all__ = [
    "read_module_datafile",
    "load_config",
    "setup_model",
    "run_model",
    "setup_and_run",
    "run_pysam_kbhdp_fpc",
]


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(__location__, "swh-kbhdp.json")


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


def setup_model(
    temperatures,
    weather_file=None,
    weather_data=None,
    config_file=None,
    config_data=None,
):

    tech_model = swh.new()

    for k, v in config_data.items():
        tech_model.value(k, v)

    if weather_file is not None:
        tech_model.value("solar_resource_file", weather_file)
    elif weather_data is not None:
        tech_model.value("solar_resource_data", weather_data)
    else:
        raise Exception("Either weather_file or weather_data must be specified.")

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


def run_model(tech_model, heat_load_mwt=None, hours_storage=None, temperature_hot=None):
    """
    :param tech_model: PySAM technology model
    :param heat_load_mwt: [MWt]
    :param hours_storage: [hr]
    :param temperature_hot: [C]
    """
    cp_water = 4.181  # [kJ/kg-K]
    density_water = 1000  # [kg/m3]
    pump_power_per_collector = 45 / 2  # [W]
    pipe_length_fixed = 9  # [m]
    pipe_length_per_collector = 0.5  # [m]

    T_cold = tech_model.value("custom_mains")[0]  # [C]
    heat_load = heat_load_mwt * 1e3 if heat_load_mwt is not None else None  # [kWt]

    print(
        f"Running:\n\tHeat Load = {heat_load_mwt}\n\tHours Storage = {hours_storage}\n\tTemperature Hot = {temperature_hot}"
    )

    if heat_load is not None:
        # Set heat load (system capacity)
        n_collectors = round(
            heat_load
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
    tech_model.value("scaled_draw", 8760 * (mdot_hot,))  # [kg/hr]

    # Set pipe diameter and pump power
    pipe_length = pipe_length_fixed + pipe_length_per_collector * tech_model.value(
        "ncoll"
    )
    tech_model.value("pipe_length", pipe_length)  # [m] default is 0.019 m
    pumping_power = pump_power_per_collector * tech_model.value("ncoll")
    tech_model.value("pump_power", pumping_power)  # [W]

    tech_model.execute()

    heat_annual = tech_model.value(
        "annual_Q_deliv"
    )  # [kWh] does not include electric heat, includes losses
    electricity_annual = sum(tech_model.value("P_pump"))  # [kWh]
    frac_electricity_annual = (
        electricity_annual / heat_annual
    )  # [-] for analysis only, plant beneficial if < 1

    return {
        "heat_annual": heat_annual,  # [kWh] annual net thermal energy in year 1
        "electricity_annual": electricity_annual,  # [kWhe]
    }


def setup_and_run(
    temperatures, weather_file, config_data, heat_load, hours_storage, temperature_hot
):

    tech_model = setup_model(
        temperatures, weather_file=weather_file, config_data=config_data
    )
    result = run_model(tech_model, heat_load, hours_storage, temperature_hot)
    return result


def plot_contour(
    df, x_label, y_label, z_label, units=None, grid=False, countour_lines=False
):
    def _set_aspect(ax, aspect):
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * aspect)

    levels = 25
    df2 = df[[x_label, y_label, z_label]].pivot([y_label, x_label, z_label])
    y = df2.index.values
    x = df2.columns.values
    z = df2.values
    fig, ax = plt.subplots(1, 1)
    cs = ax.contourf(x, y, z, levels=levels)
    if countour_lines:
        cl = ax.contour(x, y, z, colors="black", levels=levels)
        ax.clabel(cl, colors="black", fmt="%#.4g")
    if grid:
        ax.grid(color="black")
    _set_aspect(ax, 0.5)
    fig.colorbar(cs)

    if units is None:
        x_units, y_units, z_units = "", "", ""
    else:
        x_units, y_units, z_units = ["  [" + unit + "]" for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    ax.set_title(z_label + z_units)
    plt.show()


def plot_3d(df, x_label, y_label, z_label, units=None):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    surf = ax.plot_trisurf(
        df[x_label], df[y_label], df[z_label], cmap=plt.cm.viridis, linewidth=0.2
    )

    if units is None:
        x_units, y_units, z_units = "", "", ""
    else:
        x_units, y_units, z_units = ["  [" + unit + "]" for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    ax.set_title(z_label + z_units)
    plt.show()


def plot_2d(df, x_label, y_label, units=None):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    surf = ax.plot(df[x_label], df[y_label])

    if units is None:
        x_units, y_units = "", ""
    else:
        x_units, y_units = ["  [" + unit + "]" for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    plt.show()


def plot_contours(df):
    plot_contour(
        df.query("temperature_hot == 70"),
        "heat_load",
        "hours_storage",
        "heat_annual",
        units=["MWt", "hr", "kWht"],
    )
    plot_contour(
        df.query("hours_storage == 12"),
        "heat_load",
        "temperature_hot",
        "heat_annual",
        units=["MWt", "C", "kWht"],
    )
    plot_contour(
        df.query("heat_load == 500"),
        "hours_storage",
        "temperature_hot",
        "heat_annual",
        units=["hr", "C", "kWht"],
    )
    plot_contour(
        df.query("temperature_hot == 70"),
        "heat_load",
        "hours_storage",
        "electricity_annual",
        units=["MWt", "hr", "kWhe"],
    )
    plot_contour(
        df.query("hours_storage == 12"),
        "heat_load",
        "temperature_hot",
        "electricity_annual",
        units=["MWt", "C", "kWhe"],
    )
    plot_contour(
        df.query("heat_load == 500"),
        "hours_storage",
        "temperature_hot",
        "electricity_annual",
        units=["hr", "C", "kWhe"],
    )


def plot_3ds(df):
    plot_3d(
        df.query("temperature_hot == 70"),
        "heat_load",
        "hours_storage",
        "heat_annual",
        units=["MWt", "hr", "kWht"],
    )
    plot_3d(
        df.query("hours_storage == 12"),
        "heat_load",
        "temperature_hot",
        "heat_annual",
        units=["MWt", "C", "kWht"],
    )
    plot_3d(
        df.query("heat_load == 6"),
        "hours_storage",
        "temperature_hot",
        "heat_annual",
        units=["hr", "C", "kWht"],
    )
    plot_3d(
        df.query("temperature_hot == 70"),
        "heat_load",
        "hours_storage",
        "electricity_annual",
        units=["MWt", "hr", "kWhe"],
    )
    plot_3d(
        df.query("hours_storage == 6"),
        "heat_load",
        "temperature_hot",
        "electricity_annual",
        units=["MWt", "C", "kWhe"],
    )
    plot_3d(
        df.query("heat_load == 6"),
        "hours_storage",
        "temperature_hot",
        "electricity_annual",
        units=["hr", "C", "kWhe"],
    )


def run_pysam_kbhdp_fpc(
    heat_loads=np.linspace(1, 25, 25),
    hours_storages=np.linspace(0, 12, 13),
    temperature_hots=np.arange(50, 102, 2),
    temperature_cold=20,
    plot_saved_dataset=False,
    run_pysam=True,
    save_data=True,
    use_multiprocessing=True,
    dataset_filename="FPC_KBHDP_el_paso.pkl",
):
    """
    Run PySAM to collect data for FPC surrogate model
    for KBHPD case study
    """

    pysam_model_name = "SolarWaterHeatingCommercial"

    temperatures = {
        "T_cold": temperature_cold,
        "T_hot": 70,  # this will be overwritten by temperature_hot value
        "T_amb": 18,
    }

    dataset_filename = os.path.join(os.path.dirname(__file__), dataset_filename)
    config_data = read_module_datafile(param_file)

    if "solar_resource_file" in config_data:
        del config_data["solar_resource_file"]
    tech_model = setup_model(
        temperatures=temperatures,
        weather_file=weather_file,
        config_data=config_data,
    )

    # Run pysam
    data = []
    if run_pysam:
        if use_multiprocessing:
            arguments = list(product(heat_loads, hours_storages, temperature_hots))
            df = pd.DataFrame(
                arguments, columns=["heat_load", "hours_storage", "temperature_hot"]
            )

            time_start = time.process_time()
            with multiprocessing.Pool(processes=6) as pool:
                args = [
                    (temperatures, weather_file, config_data, *args)
                    for args in arguments
                ]
                results = pool.starmap(setup_and_run, args)
            time_stop = time.process_time()
            print("Multiprocessing time:", time_stop - time_start, "\n")

            df_results = pd.DataFrame(results)
            df = pd.concat(
                [
                    df,
                    df_results[
                        [
                            "heat_annual",
                            "electricity_annual",
                        ]
                    ],
                ],
                axis=1,
            )
        else:
            comb = [
                (hl, hs, th)
                for hl in heat_loads
                for hs in hours_storages
                for th in temperature_hots
            ]
            for heat_load, hours_storage, temperature_hot in comb:
                result = run_model(
                    tech_model, heat_load, hours_storage, temperature_hot
                )
                data.append(
                    [
                        heat_load,
                        hours_storage,
                        temperature_hot,
                        result["heat_annual"],
                        result["electricity_annual"],
                    ]
                )
            df = pd.DataFrame(data, columns=["heat_annual", "electricity_annual"])

        if save_data:
            df.to_pickle(dataset_filename)


if __name__ == "__main__":

    run_pysam_kbhdp_fpc(dataset_filename="test.pkl")
