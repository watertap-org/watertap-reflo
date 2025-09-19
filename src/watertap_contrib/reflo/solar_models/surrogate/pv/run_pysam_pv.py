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
import multiprocessing

import numpy as np
import pandas as pd
from math import floor, ceil

import PySAM.Pvsamv1 as pv
import PySAM.Grid as grid
import PySAM.Utilityrate5 as utilityrate
import PySAM.Singleowner as singleowner

__author__ = "Kurban Sitterley"
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

default_tech_config_file = f"{__location__}/data/pv_pvsamv1_config_default_reflo.json"
default_grid_config_file = f"{__location__}/data/pv_grid_config_default_reflo.json"
default_rate_config_file = (
    f"{__location__}/data/pv_utilityrate5_config_default_reflo.json"
)
default_cash_config_file = (
    f"{__location__}/data/pv_singleowner_config_default_reflo.json"
)
default_weather_file = (
    f"{__location__}/data/test_pv_weather_data.csv"  # from El Paso, TX
)

default_config_files = [
    default_tech_config_file,
    default_grid_config_file,
    default_rate_config_file,
    default_cash_config_file,
]


def _flatten_dict(d):
    def get_key_values(d):
        for key, value in d.items():
            if isinstance(value, dict):
                yield from get_key_values(value)
            else:
                yield key, value

    return {key: value for (key, value) in get_key_values(d)}


def _load_config_files(modules, config_files):
    for file_name, module in zip(config_files, modules):
        with open(file_name, "r") as file:
            data = json.load(file)
            missing_values = []  # for debugging
            for k, v in data.items():
                if k != "number_inputs":
                    try:
                        module.value(k, v)
                    except:
                        missing_values.append(k)
            pass


def _spe_power(spe_eff_level, spe_rad_level, spe_area):
    return spe_eff_level / 100 * spe_rad_level * spe_area


def size_pv_array(tech_model, system_capacity=50, desired_dcac_ratio=1.2):

    # Sizing rules
    # 1. Voc < Vdcmax
    # 2. Vmp > Vmin
    # 3. Vmp < Vmax
    # 4. num series * num_parallel is about desired array size (num_parallel = desired / (num series * mod_power)
    # 5. num inverters is about desired array size (num_inv = num_series * num_parallel * mod_power) / inv_power

    m = _flatten_dict(tech_model.export())

    # module parameters
    module_model = int(m["module_model"])
    m["spe_imp"] = (
        _spe_power(m["spe_eff4"], m["spe_rad4"], m["spe_area"]) / m["spe_vmp"]
    )  # 4 = reference conditions
    mod_imp = [
        m["spe_imp"],
        m["cec_i_mp_ref"],
        m["sixpar_imp"],
        m["snl_impo"],
        m["sd11par_Imp0"],
    ][module_model]
    mod_vmp = [
        m["spe_vmp"],
        m["cec_v_mp_ref"],
        m["sixpar_vmp"],
        m["snl_vmpo"],
        m["sd11par_Vmp0"],
    ][module_model]
    mod_voc = [
        m["spe_voc"],
        m["cec_v_oc_ref"],
        m["sixpar_voc"],
        m["snl_voco"],
        m["sd11par_Voc0"],
    ][module_model]
    mod_power = mod_vmp * mod_imp

    # inverter parameters
    inv_vmin = m["mppt_low_inverter"]
    inv_vmax = m["mppt_hi_inverter"]
    inverter_model = int(m["inverter_model"])
    m["inv_ds_pdco"] = m["inv_ds_paco"] / (m["inv_ds_eff"] / 100)
    inv_vdcmax = [
        m["inv_snl_vdcmax"],
        m["inv_ds_vdcmax"],
        m["inv_pd_vdcmax"],
        m["inv_cec_cg_vdcmax"],
    ][
        inverter_model
    ]  # Vdcmax
    inv_power = [
        m["inv_snl_paco"],
        m["inv_ds_paco"],
        m["inv_pd_paco"],
        m["inv_cec_cg_paco"],
    ][
        inverter_model
    ]  # Paco
    inv_dc_power = [
        m["inv_snl_pdco"],
        m["inv_ds_pdco"],
        m["inv_pd_pdco"],
        m["inv_cec_cg_pdco"],
    ][
        inverter_model
    ]  # Pdco

    if mod_vmp > 0:
        num_series = 0.5 * (inv_vmin + inv_vmax) / mod_vmp
        if inv_vdcmax > 0:
            while num_series > 0 and (num_series * mod_voc) > inv_vdcmax:
                num_series -= 1

    num_series = max(1, round(num_series))
    num_parallel = system_capacity * 1000 / (num_series * mod_power)
    num_parallel = max(1, round(num_parallel))
    if desired_dcac_ratio > 0:
        inverters = ((num_series * num_parallel * mod_power)) / (
            desired_dcac_ratio * inv_power
        )
        # round inverters for best DC-AC ratio
        if inverters - floor(inverters) < 0.5:
            num_inverters = floor(inverters)
        else:
            num_inverters = ceil(inverters)
    else:
        num_inverters = ceil(((num_series * num_parallel * mod_power)) / inv_power)
    num_inverters = max(1, num_inverters)
    total_modules = num_series * num_parallel
    total_module_area = total_modules * tech_model.value("spe_area")  # [m2]
    total_module_area = total_module_area * 0.000247105  # m2 to [acre]
    land_area = total_module_area / tech_model.value(
        "subarray1_gcr"
    )  # assumes land_area_multiplier = 1
    total_ac_capacity = inv_power * num_inverters / 1000
    total_dc_inverter_capacity = inv_dc_power * num_inverters / 1000

    # check that the sizing was close to the desired sizes, otherwise error out
    nameplate_dc = total_modules * mod_power / 1000
    proposed_ratio = nameplate_dc / (num_inverters * inv_power / 1000)
    if abs(nameplate_dc - system_capacity) / system_capacity > 0.2:
        num_inverters = None
        num_series = None
        num_parallel = None
        total_modules = None
        nameplate_dc = None
        total_ac_capacity = None
        total_dc_inverter_capacity = None
    else:
        num_inverters = num_inverters
        num_series = num_series
        num_parallel = num_parallel
        total_modules = total_modules
        nameplate_dc = nameplate_dc
        total_ac_capacity = total_ac_capacity
        total_dc_inverter_capacity = total_dc_inverter_capacity

    pv_array_design = {
        "inverter_count": num_inverters,
        "number_modules_per_string": num_series,
        "number_strings": num_parallel,
        "total_modules": total_modules,
        "total_module_area": total_module_area,
        "land_req": land_area,
        "nameplate_dc": nameplate_dc,
        "total_inverter_capacity": total_ac_capacity,
        "total_dc_inverter_capacity": total_dc_inverter_capacity,
    }

    if not all([num_inverters, num_parallel]):
        raise ValueError("One of inverters, num_inverters, num_parallel is None.")

    return pv_array_design


def setup_model_pv(pysam_model_config, config_files=[], weather_file=None):

    tech_model = pv.new()
    grid_model = grid.from_existing(tech_model, pysam_model_config)
    rate_model = utilityrate.from_existing(tech_model, pysam_model_config)
    cash_model = singleowner.from_existing(tech_model, pysam_model_config)
    modules = [
        tech_model,
        grid_model,
        rate_model,
        cash_model,
    ]
    _load_config_files(modules, config_files=config_files)
    tech_model.SolarResource.solar_resource_file = weather_file

    return modules


def run_pv_single_owner(
    modules,
    pysam_model_config,
    system_capacity=None,  # kW
    desired_dcac_ratio=1.2,
    tech_model_kwargs={},
    cash_model_kwargs={},
):
    if system_capacity is None:
        raise ValueError("system_capacity input must be provided")

    tech_model = modules[0]
    cash_model = modules[3]

    pv_array_design = size_pv_array(
        tech_model,
        system_capacity=system_capacity,
        desired_dcac_ratio=desired_dcac_ratio,
    )

    print(
        f"\nRunning PySAM model {pysam_model_config} for system capacity = {system_capacity:.2f} kW...\n"
    )
    tech_model.value("inverter_count", pv_array_design["inverter_count"])
    tech_model.value("subarray1_nstrings", pv_array_design["number_strings"])

    for param, val in tech_model_kwargs.items():
        tech_model.value(param, val)
    for param, val in cash_model_kwargs.items():
        if not isinstance(val, list):
            val = [val]
        cash_model.value(param, val)

    for mod in modules:
        mod.execute()

    return tech_model, cash_model, pv_array_design


def setup_and_run_pv(
    system_capacity,
    pysam_model_config,
    config_files,
    weather_file,
    return_tech_model=False,
):
    if config_files is None:
        config_files = default_config_files

    modules = setup_model_pv(
        pysam_model_config, config_files=config_files, weather_file=weather_file
    )

    tech_model, cash_model, pv_array_design = run_pv_single_owner(
        modules, pysam_model_config, system_capacity=system_capacity
    )
    electricity_annual = tech_model.Outputs.annual_energy  # kWh

    cash_model.execute()

    result = pv_array_design
    result["electricity_annual"] = electricity_annual
    result["lcoe_real"] = cash_model.Outputs.lcoe_real
    result["ac_capacity_factor"] = tech_model.Outputs.capacity_factor_ac
    result["dc_capacity_factor"] = tech_model.Outputs.capacity_factor
    print(f"System Capacity {system_capacity:.1f} kW completed.")
    print(f"\tAnnual Energy = {electricity_annual:.2f} kWh")
    print(f"\tLand Required = {result['land_req']:.2f} acre\n")

    if return_tech_model:
        return result, tech_model
    else:
        return result


def generate_pv_data(
    system_capacities=np.linspace(100, 10000, 3),
    pysam_model_config="FlatPlatePVSingleOwner",
    save_data=True,
    use_multiprocessing=True,
    processes=8,
    dataset_filename=None,
    weather_file=None,
    tech_config_file=None,
    grid_config_file=None,
    rate_config_file=None,
    cash_config_file=None,
):

    if dataset_filename is None:
        # assume it is run for testing purposes
        dataset_filename = os.path.join(__location__, "data/test_data.pkl")
    if weather_file is None:
        weather_file = default_weather_file
    if tech_config_file is None:
        tech_config_file = default_tech_config_file
    if grid_config_file is None:
        grid_config_file = default_grid_config_file
    if rate_config_file is None:
        rate_config_file = default_rate_config_file
    if cash_config_file is None:
        cash_config_file = default_cash_config_file

    config_files = [
        tech_config_file,
        grid_config_file,
        rate_config_file,
        cash_config_file,
    ]

    print(f"Saving data to {dataset_filename}")

    df = pd.DataFrame(system_capacities, columns=["system_capacity"])

    if use_multiprocessing:

        with multiprocessing.Pool(processes=processes) as pool:
            args = [
                (ds, pysam_model_config, config_files, weather_file)
                for ds in system_capacities
            ]
            results = pool.starmap(setup_and_run_pv, args)

        df_results = pd.DataFrame(results)

    else:
        results = list()
        for ds in system_capacities:
            result = setup_and_run_pv(ds, pysam_model_config)
            results.append(result)
        df_results = pd.DataFrame(results)

    df = pd.concat(
        [df, df_results],
        axis=1,
    )
    if save_data:
        if dataset_filename[-3:] == "pkl":
            df.to_pickle(dataset_filename)
        elif dataset_filename[-3:] == "csv":
            df.to_csv(dataset_filename)
        else:
            raise ValueError("dataset_filename must end in .pkl or .csv")

    return df


if __name__ == "__main__":
    # To create data used in PV surrogate test file:
    # system_capacities = np.linspace(1000, 1000000, 50)

    df = generate_pv_data()
    print(df.head(20))
    os.remove(os.path.join(__location__, "data/test_data.pkl"))
