import json
from os.path import join, dirname
from math import floor, ceil, isnan
import numpy as np
import pandas as pd
import time
import multiprocessing
from itertools import product
import matplotlib.pyplot as plt
import PySAM.Swh as swh
import PySAM.Utilityrate5 as utilityrate5
import PySAM.Cashloan as cashloan


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
    model_name,
    temperatures,
    weather_file=None,
    weather_data=None,
    config_files=None,
    config_data=None,
):

    tech_model = swh.new()
    bill_model = utilityrate5.from_existing(tech_model, model_name)
    cash_model = cashloan.from_existing(tech_model, model_name)
    modules = [tech_model, bill_model, cash_model]
    load_config(modules, config_files, config_data)
    if weather_file is not None:
        tech_model.value("solar_resource_file", weather_file)
    elif weather_data is not None:
        tech_model.value("solar_resource_data", weather_data)
    else:
        raise Exception("Either weather_file or weather_data must be specified.")

    # Set constant temperatures
    tech_model.value("use_custom_mains", 1)
    tech_model.value("custom_mains", 8760 * (temperatures["T_cold"],))
    tech_model.value("use_custom_set", 0)  # ensure constant hot set temperature is used
    tech_model.value("T_set", temperatures["T_hot"])
    tech_model.value("T_room", temperatures["T_amb"])

    # Set collector loop mass flow (this should be done automatically in SAM)
    tech_model.value(
        "mdot", tech_model.value("test_flow") * tech_model.value("ncoll")
    )  # [kg/s]

    # Set pump and piping parameters
    tech_model.value("pipe_diam", 0.019)  # [m] default is 0.019 m
    tech_model.value("pump_eff", 0.85)  # [-] pump efficiency default is 0.85

    # Ensure system capacity parameter agreement
    system_capacity_actual = system_capacity_computed(tech_model)
    tech_model.value("system_capacity", system_capacity_actual)

    return {
        "tech_model": tech_model,
        "bill_model": bill_model,
        "cash_model": cash_model,
    }


def run_model(modules, heat_load_mwt=None, hours_storage=None):
    """
    :param modules: dict of PySAM modules with keys 'tech_model', 'bill_model' and 'cash_model'
    :param heat_load_mwt: [MWt]
    :param hours_storage: [hr]
    """
    CP_WATER = 4.181  # [kJ/kg-K]
    DENSITY_WATER = 1000  # [kg/m3]
    PUMP_POWER_PER_COLLECTOR = 45 / 2  # [W]
    PIPE_LENGTH_FIXED = 9  # [m]
    PIPE_LENGTH_PER_COLLECTOR = 0.5  # [m]

    tech_model = modules["tech_model"]
    bill_model = modules["bill_model"]
    cash_model = modules["cash_model"]

    T_cold = tech_model.value("custom_mains")[0]  # [C]
    T_hot = tech_model.value("T_set")  # [C]
    heat_load = heat_load_mwt * 1e3  # [kWt]

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
    if hours_storage is not None:
        # Set hours of storage (tank volume)
        hours_storage = max(hours_storage, 1e-3)  # don't accept 0 hours
        system_capacity = tech_model.value("system_capacity")  # [kW]
        mass_tank_water = (
            hours_storage * 3600 * system_capacity / (CP_WATER * (T_hot - T_cold))
        )  # [kg]
        volume_tank = mass_tank_water / DENSITY_WATER  # [m3]
        tech_model.value("V_tank", volume_tank)

    # Set collector loop and hot water mass flow rates
    tech_model.value(
        "mdot", tech_model.value("test_flow") * tech_model.value("ncoll")
    )  # [kg/s]
    mdot = tech_model.value("system_capacity") / (CP_WATER * (T_hot - T_cold))  # [kg/s]
    tech_model.value("scaled_draw", 8760 * (mdot * 3600,))  # [kg/hr]

    # Set pipe diameter and pump power
    tech_model.value(
        "pipe_length",
        PIPE_LENGTH_FIXED + PIPE_LENGTH_PER_COLLECTOR * tech_model.value("ncoll"),
    )  # [m] default is 0.019 m
    tech_model.value(
        "pump_power", PUMP_POWER_PER_COLLECTOR * tech_model.value("ncoll")
    )  # [W]

    tech_model.execute()

    annual_energy = (
        tech_model.Outputs.annual_Q_deliv
    )  # [kWh] does not include electric heat, includes losses
    electrical_load = sum(tech_model.Outputs.P_pump)  # [kWh]
    frac_electrical_load = (
        electrical_load / annual_energy
    )  # [-] for analysis only, plant beneficial if < 1

    # NOTE: running these are not required for generating the annual energy and electrical load
    bill_model.execute()
    cash_model.execute()

    return {
        "annual_energy": annual_energy,  # [kWh] annual net thermal energy in year 1
        "electrical_load": electrical_load,  # [kWhe]
    }


def setup_and_run(
    model_name, temperatures, weather_file, config_data, heat_load, hours_storage
):

    modules = setup_model(
        model_name, temperatures, weather_file=weather_file, config_data=config_data
    )
    result = run_model(modules, heat_load, hours_storage)
    return result


def plot_3d(df, x_index=0, y_index=1, z_index=2, grid=True, countour_lines=True):
    """
    index 0 = x axis
    index 1 = y axis
    index 2 = z axis
    """

    def _set_aspect(ax, aspect):
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * aspect)

    # CONTOUR PLOT
    levels = 25
    df2 = df.pivot(df.columns[y_index], df.columns[x_index], df.columns[z_index])
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
    ax.set_xlabel(df.columns[x_index])
    ax.set_ylabel(df.columns[y_index])
    ax.set_title(df.columns[z_index])
    plt.show()


#########################################################################################################
if __name__ == "__main__":
    PLOT_SAVED_DATASET = False  # plot previously run, saved data?
    USE_MULTIPROCESSING = True
    HEAT_LOADS = np.arange(100, 1100, 25)  # [MWt]
    HEAT_LOADS = np.arange(5, 115, 10)  # [MWt]
    HOURS_STORAGES = np.arange(0, 27, 1)  # [hr]
    TEMPERATURES = {
        "T_cold": 20,
        "T_hot": 70,
        "T_amb": 18,
    }
    MODEL_NAME = "SolarWaterHeatingCommercial"
    CONFIG_FILES = [
        join(dirname(__file__), "untitled_swh.json"),
        join(dirname(__file__), "untitled_utilityrate5.json"),
        join(dirname(__file__), "untitled_cashloan.json"),
    ]
    WEATHER_FILE = join(
        dirname(__file__), "tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
    )
    DATASET_FILENAME = join(
        dirname(__file__), "flat_plate_data.pkl"
    )  # output dataset for surrogate training

    config_data = [read_module_datafile(config_file) for config_file in CONFIG_FILES]
    if "solar_resource_file" in config_data[0]:
        del config_data[0]["solar_resource_file"]
    modules = setup_model(
        model_name=MODEL_NAME,
        temperatures=TEMPERATURES,
        weather_file=WEATHER_FILE,
        config_data=config_data,
    )

    # Run model for single parameter set
    result = run_model(modules, heat_load_mwt=1000, hours_storage=1)

    if PLOT_SAVED_DATASET:
        # Load and plot saved df
        df = pd.read_pickle(DATASET_FILENAME)
        plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)  # annual energy
        plot_3d(df, 0, 1, 3, grid=False, countour_lines=False)  # electrical load

    # Run parametrics
    data = []
    if USE_MULTIPROCESSING:
        arguments = list(product(HEAT_LOADS, HOURS_STORAGES))
        df = pd.DataFrame(arguments, columns=["heat_load", "hours_storage"])

        time_start = time.process_time()
        with multiprocessing.Pool(processes=6) as pool:
            args = [
                (MODEL_NAME, TEMPERATURES, WEATHER_FILE, config_data, *args)
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
                        "annual_energy",
                        "electrical_load",
                    ]
                ],
            ],
            axis=1,
        )
    else:
        comb = [(hl, hs) for hl in HEAT_LOADS for hs in HOURS_STORAGES]
        for heat_load, hours_storage in comb:
            result = run_model(modules, heat_load, hours_storage)
            data.append(
                [
                    heat_load,
                    hours_storage,
                    result["annual_energy"],
                    result["electrical_load"],
                ]
            )
        df = pd.DataFrame(data, columns=["annual_energy", "electrical_load"])

    df.to_pickle(DATASET_FILENAME)
    plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)  # annual energy
    plot_3d(df, 0, 1, 3, grid=False, countour_lines=False)  # electrical_load
