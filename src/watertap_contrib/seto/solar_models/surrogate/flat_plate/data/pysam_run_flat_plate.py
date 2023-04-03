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
    load_config([tech_model], config_file, [config_data])
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

    return tech_model


def run_model(tech_model, heat_load_mwt=None, hours_storage=None, temperature_hot=None):
    """
    :param tech_model: PySAM technology model
    :param heat_load_mwt: [MWt]
    :param hours_storage: [hr]
    :param temperature_hot: [C]
    """
    CP_WATER = 4.181  # [kJ/kg-K]
    DENSITY_WATER = 1000  # [kg/m3]
    PUMP_POWER_PER_COLLECTOR = 45 / 2  # [W]
    PIPE_LENGTH_FIXED = 9  # [m]
    PIPE_LENGTH_PER_COLLECTOR = 0.5  # [m]

    T_cold = tech_model.value("custom_mains")[0]  # [C]
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
    if temperature_hot is not None:
        # Set hot outlet temperature
        tech_model.value("T_set", temperature_hot)
    if hours_storage is not None:
        # Set hours of storage (tank volume)
        hours_storage = max(hours_storage, 1e-3)  # don't accept 0 hours
        system_capacity = tech_model.value("system_capacity")  # [kW]
        T_hot = tech_model.value("T_set")  # [C]
        mass_tank_water = (
            hours_storage * 3600 * system_capacity / (CP_WATER * (T_hot - T_cold))
        )  # [kg]
        volume_tank = mass_tank_water / DENSITY_WATER  # [m3]
        tech_model.value("V_tank", volume_tank)

    # Set collector loop and hot water mass flow rates
    mdot_collectors = tech_model.value("test_flow") * tech_model.value("ncoll")
    tech_model.value("mdot", mdot_collectors)  # [kg/s]
    T_hot = tech_model.value("T_set")  # [C]
    mdot_hot = tech_model.value("system_capacity") / (CP_WATER * (T_hot - T_cold)) * 3600  # [kg/hr]
    tech_model.value("scaled_draw", 8760 * (mdot_hot,))  # [kg/hr]

    # Set pipe diameter and pump power
    pipe_length = PIPE_LENGTH_FIXED + PIPE_LENGTH_PER_COLLECTOR * tech_model.value("ncoll")
    tech_model.value("pipe_length", pipe_length)  # [m] default is 0.019 m
    pumping_power = PUMP_POWER_PER_COLLECTOR * tech_model.value("ncoll")
    tech_model.value("pump_power", pumping_power)  # [W]

    tech_model.execute()

    annual_energy = (
        tech_model.value("annual_Q_deliv")
    )  # [kWh] does not include electric heat, includes losses
    electrical_load = sum(tech_model.value("P_pump"))  # [kWh]
    frac_electrical_load = (
        electrical_load / annual_energy
    )  # [-] for analysis only, plant beneficial if < 1

    return {
        "annual_energy": annual_energy,  # [kWh] annual net thermal energy in year 1
        "electrical_load": electrical_load,  # [kWhe]
    }


def setup_and_run(
    temperatures, weather_file, config_data, heat_load, hours_storage, temperature_hot
):

    tech_model = setup_model(
        temperatures, weather_file=weather_file, config_data=config_data
    )
    result = run_model(tech_model, heat_load, hours_storage, temperature_hot)
    return result


def plot_contour(df, x_label, y_label, z_label, units=None, grid=False, countour_lines=False):
    def _set_aspect(ax, aspect):
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * aspect)

    levels = 25
    df2 = df[[x_label, y_label, z_label]].pivot(y_label, x_label, z_label)
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
        x_units, y_units, z_units = '', '', ''
    else:
        x_units, y_units, z_units = ['  [' + unit + ']' for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    ax.set_title(z_label + z_units)
    plt.show()


def plot_3d(df, x_label, y_label, z_label, units=None):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    surf = ax.plot_trisurf(df[x_label], df[y_label], df[z_label], cmap=plt.cm.viridis, linewidth=0.2)

    if units is None:
        x_units, y_units, z_units = '', '', ''
    else:
        x_units, y_units, z_units = ['  [' + unit + ']' for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    ax.set_title(z_label + z_units)
    plt.show()


def plot_2d(df, x_label, y_label, units=None):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1)
    surf = ax.plot(df[x_label], df[y_label])

    if units is None:
        x_units, y_units = '', ''
    else:
        x_units, y_units = ['  [' + unit + ']' for unit in units]
    ax.set_xlabel(x_label + x_units)
    ax.set_ylabel(y_label + y_units)
    plt.show()


def plot_contours(df):
    plot_contour(df.query('temperature_hot == 70'), 'heat_load', 'hours_storage', 'annual_energy', units=['MWt', 'hr', 'kWht'])
    plot_contour(df.query('hours_storage == 12'), 'heat_load', 'temperature_hot', 'annual_energy', units=['MWt', 'C', 'kWht'])
    plot_contour(df.query('heat_load == 500'), 'hours_storage', 'temperature_hot', 'annual_energy', units=['hr', 'C', 'kWht'])
    plot_contour(df.query('temperature_hot == 70'), 'heat_load', 'hours_storage', 'electrical_load', units=['MWt', 'hr', 'kWhe'])
    plot_contour(df.query('hours_storage == 12'), 'heat_load', 'temperature_hot', 'electrical_load', units=['MWt', 'C', 'kWhe'])
    plot_contour(df.query('heat_load == 500'), 'hours_storage', 'temperature_hot', 'electrical_load', units=['hr', 'C', 'kWhe'])


def plot_3ds(df):
    plot_3d(df.query('temperature_hot == 70'), 'heat_load', 'hours_storage', 'annual_energy', units=['MWt', 'hr', 'kWht'])
    plot_3d(df.query('hours_storage == 12'), 'heat_load', 'temperature_hot', 'annual_energy', units=['MWt', 'C', 'kWht'])
    plot_3d(df.query('heat_load == 500'), 'hours_storage', 'temperature_hot', 'annual_energy', units=['hr', 'C', 'kWht'])
    plot_3d(df.query('temperature_hot == 70'), 'heat_load', 'hours_storage', 'electrical_load', units=['MWt', 'hr', 'kWhe'])
    plot_3d(df.query('hours_storage == 12'), 'heat_load', 'temperature_hot', 'electrical_load', units=['MWt', 'C', 'kWhe'])
    plot_3d(df.query('heat_load == 500'), 'hours_storage', 'temperature_hot', 'electrical_load', units=['hr', 'C', 'kWhe'])


def debug_t_hot(tech_model):
    # Running at 73.5 and 74 C
    result = run_model(tech_model, 500, 12, 73)
    result = run_model(tech_model, 500, 12, 73.5)
    result = run_model(tech_model, 500, 12, 74)


    dataset_filename = join(
        dirname(__file__), "debugging_t_hot.pkl"
    )  # output dataset for surrogate training

    if False:
        df = pd.read_pickle(dataset_filename)
        plot_2d(df.query('hours_storage == 12 & heat_load == 500'), 'temperature_hot', 'annual_energy', units=['C', 'kWht'])

    heat_loads = [500]                              # [MWt]
    hours_storages = [12]                           # [hr]
    temperature_hots = np.arange(50, 100.5, 0.5)    # [C]
    comb = [(hl, hs, th) for hl in heat_loads for hs in hours_storages for th in temperature_hots]
    data = []
    for heat_load, hours_storage, temperature_hot in comb:
        result = run_model(tech_model, heat_load, hours_storage, temperature_hot)
        data.append(
            [
                heat_load,
                hours_storage,
                temperature_hot,
                result["annual_energy"],
                result["electrical_load"],
            ]
        )
    df = pd.DataFrame(data, columns=["heat_load", "hours_storage", "temperature_hot", "annual_energy", "electrical_load"])
    df.to_pickle(dataset_filename)
    plot_2d(df.query('hours_storage == 12 & heat_load == 500'), 'temperature_hot', 'annual_energy', units=['C', 'kWht'])
    x = 1


#########################################################################################################
if __name__ == "__main__":
    DEBUG = False
    PLOT_SAVED_DATASET = True                   # plot previously run, saved data?
    RUN_PARAMETRICS = False
    USE_MULTIPROCESSING = True

    # HEAT_LOADS = np.arange(5, 115, 10)          # [MWt]
    HEAT_LOADS = np.arange(100, 1100, 25)       # [MWt]
    HOURS_STORAGES = np.arange(0, 27, 1)        # [hr]
    TEMPERATURE_HOTS = np.arange(50, 102, 2)    # [C]
    TEMPERATURES = {
        "T_cold": 20,
        "T_hot": 70,        # this will be overwritten by temperature_hot value
        "T_amb": 18,
    }
    MODEL_NAME = "SolarWaterHeatingCommercial"
    PARAM_FILE = join(dirname(__file__), "untitled_swh.json")
    WEATHER_FILE = join(
        dirname(__file__), "tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
    )
    DATASET_FILENAME = join(
        dirname(__file__), "flat_plate_data.pkl"
    )  # output dataset for surrogate training

    config_data = read_module_datafile(PARAM_FILE)
    if "solar_resource_file" in config_data:
        del config_data["solar_resource_file"]
    tech_model = setup_model(
        temperatures=TEMPERATURES,
        weather_file=WEATHER_FILE,
        config_data=config_data,
    )

    if DEBUG:
        debug_t_hot(tech_model)

    if PLOT_SAVED_DATASET:
        # Load and plot saved df (x, y z)
        df = pd.read_pickle(DATASET_FILENAME)
        plot_2d(df.query('hours_storage == 12 & heat_load == 500'), 'temperature_hot', 'annual_energy', units=['C', 'kWht'])
        # plot_3ds(df)
        # plot_contours(df)

    # Run model for single parameter set
    result = run_model(tech_model, heat_load_mwt=1000, hours_storage=1, temperature_hot=70)

    # Run parametrics
    data = []
    if RUN_PARAMETRICS:
        if USE_MULTIPROCESSING:
            arguments = list(product(HEAT_LOADS, HOURS_STORAGES, TEMPERATURE_HOTS))
            df = pd.DataFrame(arguments, columns=["heat_load", "hours_storage", "temperature_hot"])

            time_start = time.process_time()
            with multiprocessing.Pool(processes=6) as pool:
                args = [
                    (TEMPERATURES, WEATHER_FILE, config_data, *args)
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
            comb = [(hl, hs, th) for hl in HEAT_LOADS for hs in HOURS_STORAGES for th in TEMPERATURE_HOTS]
            for heat_load, hours_storage, temperature_hot in comb:
                result = run_model(tech_model, heat_load, hours_storage, temperature_hot)
                data.append(
                    [
                        heat_load,
                        hours_storage,
                        temperature_hot,
                        result["annual_energy"],
                        result["electrical_load"],
                    ]
                )
            df = pd.DataFrame(data, columns=["annual_energy", "electrical_load"])

        df.to_pickle(DATASET_FILENAME)
        plot_contours(df)
