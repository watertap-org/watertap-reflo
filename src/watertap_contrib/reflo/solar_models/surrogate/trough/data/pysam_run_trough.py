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
import PySAM.TroughPhysicalProcessHeat as iph
import PySAM.IphToLcoefcr as iph_to_lcoefcr
import PySAM.Lcoefcr as lcoefcr


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


def tes_cost(tech_model):
    STORAGE_COST_SPECIFIC = 62  # [$/kWht] borrowed from physical power trough
    tes_thermal_capacity = (
        tech_model.value("q_pb_design") * 1e3 * tech_model.value("tshours")
    )  # [kWht]
    return tes_thermal_capacity * STORAGE_COST_SPECIFIC


def system_capacity(tech_model):
    return (
        tech_model.value("q_pb_design")
        * tech_model.value("specified_solar_multiple")
        * 1e3
    )  # [kW]


def setup_model(
    model_name,
    weather_file=None,
    weather_data=None,
    config_files=None,
    config_data=None,
):
    tech_model = iph.new()
    post_model = iph_to_lcoefcr.from_existing(tech_model, model_name)
    cash_model = lcoefcr.from_existing(tech_model, model_name)
    modules = [tech_model, post_model, cash_model]
    load_config(modules, config_files, config_data)
    if weather_file is not None:
        tech_model.Weather.file_name = weather_file
    elif weather_data is not None:
        tech_model.Weather.solar_resource_data = weather_data
    else:
        raise Exception("Either weather_file or weather_data must be specified.")

    # Determine storage cost component
    capital_cost_orig = cash_model.value("capital_cost")
    storage_cost_orig = tes_cost(tech_model)
    capital_cost_minus_storage = capital_cost_orig - storage_cost_orig
    capital_cost_minus_storage_per_kW = capital_cost_minus_storage / system_capacity(
        tech_model
    )

    fixed_operating_cost_per_kW = cash_model.value(
        "fixed_operating_cost"
    ) / system_capacity(tech_model)

    return {
        "tech_model": tech_model,
        "post_model": post_model,
        "cash_model": cash_model,
        "capital_cost_minus_storage_per_kW": capital_cost_minus_storage_per_kW,
        "fixed_operating_cost_per_kW": fixed_operating_cost_per_kW,
    }


def run_model(modules, heat_load=None, hours_storage=None):
    tech_model = modules["tech_model"]
    post_model = modules["post_model"]
    cash_model = modules["cash_model"]

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
        tech_model.Outputs.annual_gross_energy - freeze_protection
    )  # [kWht] net, does not include that used for freeze protection
    capacity_factor = (
        annual_energy / (tech_model.value("q_pb_design") * 1e3 * 8760) * 100
    )  # [%]
    electrical_load = tech_model.Outputs.annual_electricity_consumption  # [kWhe]

    post_model.value(
        "fixed_operating_cost",
        modules["fixed_operating_cost_per_kW"] * system_capacity(tech_model),
    )
    post_model.execute()

    cash_model.value(
        "annual_energy", annual_energy
    )  # override the linked value from the tech model, which could be nan
    cash_model.value(
        "capital_cost",
        modules["capital_cost_minus_storage_per_kW"] * system_capacity(tech_model)
        + tes_cost(tech_model),
    )
    cash_model.execute()

    lcoh = cash_model.Outputs.lcoe_fcr  # [$/kWht]
    capital_cost = cash_model.value("capital_cost")  # [$]
    fixed_operating_cost = cash_model.value(
        "fixed_operating_cost"
    )  # [$] more than shown in UI, includes electricity purchases
    variable_operating_cost = cash_model.value("variable_operating_cost")  # [$]

    return {
        "annual_energy": annual_energy,  # [kWh] annual net thermal energy in year 1
        "freeze_protection": freeze_protection,  # [kWht]
        "capacity_factor": capacity_factor,  # [%] capacity factor
        "electrical_load": electrical_load,  # [kWhe]
        "lcoh": lcoh,  # [$/kWht] LCOH
        "capital_cost": capital_cost,  # [$]
        "fixed_operating_cost": fixed_operating_cost,  # [$]
        "variable_operating_cost": variable_operating_cost,  # [$]
    }


def setup_and_run(model_name, weather_file, config_data, heat_load, hours_storage):
    modules = setup_model(
        model_name, weather_file=weather_file, config_data=config_data
    )
    result = run_model(modules, heat_load, hours_storage)
    return result


def plot_3d(df, x_index=0, y_index=1, z_index=2, grid=True, countour_lines=True):
    """
    index 0 = x axis
    index 1 = y axis
    index 2 = z axis
    """
    # 3D PLOT
    # fig = plt.figure(figsize=(8,6))
    # ax = fig.add_subplot(1, 1, 1, projection='3d')
    # surf = ax.plot_trisurf(df.iloc[:,0], df.iloc[:,1], df.iloc[:,2], cmap=plt.cm.viridis, linewidth=0.2)
    # modld_pts = ax.scatter(df.iloc[:,0], df.iloc[:,1], df.iloc[:,2], c='black', s=15)
    # ax.set_xlabel(df.columns[0])
    # ax.set_ylabel(df.columns[1])
    # ax.set_zlabel(df.columns[2])
    # plt.show()

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
    model_name = "PhysicalTroughIPHLCOHCalculator"
    config_files = [
        join(dirname(__file__), "trough_physical_process_heat-reflo.json"),
        join(dirname(__file__), "iph_to_lcoefcr-reflo.json"),
        join(dirname(__file__), "lcoefcr-reflo.json"),
    ]
    weather_file = join(
        dirname(__file__), "tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
    )
    dataset_filename = join(
        dirname(__file__), "test_trough_data.pkl"
    )  # output dataset for surrogate training

    # modules = setup_model(model_name, weather_file, config_files=config_files)
    config_data = [read_module_datafile(config_file) for config_file in config_files]
    del config_data[0]["file_name"]  # remove weather filename
    # modules = setup_model(model_name, weather_file, config_data=config_data)

    # Run default model
    # result = run_model(modules, heat_load=None, hours_storage=None)

    # Run model at specific parameters
    # result = run_model(modules, heat_load=600, hours_storage=3)

    # Run model conducive to multiprocessing
    # weather_data = utils.read_weather_data(weather_file)    # passing of weather data not yet enabled
    # result_check = setup_and_run(model_name, weather_file, config_data, heat_load=600, hours_storage=3)

    # # Load and plot saved df
    # df = pd.read_pickle('pickle_multiproc2.pkl')
    # plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)      # annual energy
    # plot_3d(df, 0, 1, 3, grid=False, countour_lines=False)      # capacity factor
    # plot_3d(df, 0, 1, 4, grid=False, countour_lines=False)      # LCOH

    # Run parametrics via multiprocessing
    data = []
    # heat_loads =        np.arange(100, 1100, 100)
    # hours_storages =    np.arange(0, 27, 3)
    heat_loads = np.arange(100, 1100, 25)  # [MWt]
    hours_storages = np.arange(0, 27, 1)  # [hr]
    # heat_loads =        np.arange(100, 300, 100)
    # hours_storages =    np.arange(0, 6, 3)
    arguments = list(product(heat_loads, hours_storages))
    df = pd.DataFrame(arguments, columns=["heat_load", "hours_storage"])

    time_start = time.process_time()
    with multiprocessing.Pool(processes=6) as pool:
        args = [(model_name, weather_file, config_data, *args) for args in arguments]
        results = pool.starmap(setup_and_run, args)
    time_stop = time.process_time()
    print("Multiprocessing time:", time_stop - time_start, "\n")

    df_results = pd.DataFrame(results)
    df_results["total_cost"] = (
        df_results["capital_cost"]
        + df_results["fixed_operating_cost"]
        + df_results["variable_operating_cost"]
    )
    df = pd.concat(
        [
            df,
            df_results[
                [
                    "annual_energy",
                    "freeze_protection",
                    "capacity_factor",
                    "electrical_load",
                    "lcoh",
                    "capital_cost",
                    "fixed_operating_cost",
                    "variable_operating_cost",
                    "total_cost",
                ]
            ],
        ],
        axis=1,
    )
    df.to_pickle(dataset_filename)

    # # Run parametrics
    # data = []
    # # heat_loads =        np.arange(5, 115, 10)       # [MWt]
    # # heat_loads =        np.arange(100, 1100, 25)   # [MWt]
    # # hours_storages =    np.arange(0, 27, 1)         # [hr]
    # heat_loads =        np.arange(100, 1100, 100)   # [MWt]
    # hours_storages =    np.arange(0, 27, 3)         # [hr]
    # comb = [(hl, hs) for hl in heat_loads for hs in hours_storages]
    # for heat_load, hours_storage in comb:
    #     result = run_model(modules, heat_load, hours_storage)
    #     total_cost = result['capital_cost'] + result['fixed_operating_cost'] + result['variable_operating_cost']
    #     data.append([
    #         heat_load,
    #         hours_storage,
    #         result['annual_energy'],
    #         result['capacity_factor'],
    #         result['lcoh'],
    #         total_cost,
    #     ])
    # df = pd.DataFrame(data, columns=[
    #     'heat_load',
    #     'hours_storage',
    #     'annual_energy',
    #     'capacity_factor',
    #     'lcoh',
    #     'total_cost'])

    # df.to_pickle('pickle_filename5.pkl')
    plot_3d(df, 0, 1, 2, grid=False, countour_lines=False)  # annual energy
    plot_3d(df, 0, 1, 4, grid=False, countour_lines=False)  # capacity factor
    plot_3d(df, 0, 1, 6, grid=False, countour_lines=False)  # lcoh

    x = 1  # for breakpoint
    pass
