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
import time
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

pysam_model = "pv"
tech_config_file = f"{__location__}/pv/untitled_pvsamv1.json"
grid_config_file = f"{__location__}/pv/untitled_grid.json"
rate_config_file = f"{__location__}/pv/untitled_utilityrate5.json"
cash_config_file = f"{__location__}/pv/untitled_singleowner.json"
weather_file = f"{__location__}/el_paso_texas-KBHDP-weather.csv"
dcac_ratio = 1.2

if pysam_model is None:
    raise Exception("PySAM module is required.")
if not all([tech_config_file, grid_config_file, rate_config_file, cash_config_file]):
    raise Exception("One of the PySAM model configuration files is missing.")
if weather_file is None:
    raise Exception("Weather file is required.")

_pysam_model = pysam_model
_tech_config_file = tech_config_file
_grid_config_file = grid_config_file
_rate_config_file = rate_config_file
_cash_config_file = cash_config_file
_weather_file = weather_file
_dcac_ratio = dcac_ratio
_config_files = [
    _tech_config_file,
    _grid_config_file,
    _rate_config_file,
    _cash_config_file,
]


def setup_pv_single_owner(pysam_model_config):
    # global modules
    # print(f"\nBuilding PySAM model {pysam_model_config}...\n")
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
    _load_config_files(modules)
    tech_model.SolarResource.solar_resource_file = _weather_file

    return modules


def run_pv_single_owner(
    modules,
    pysam_model_config,
    design_size=None,  # kW
    desired_dcac_ratio=1.2,
    tech_model_kwargs={},
    cash_model_kwargs={},
):
    if design_size is None:
        raise ValueError("design_size input must be provided")
    tech_model = modules[0]
    cash_model = modules[3]
    size_pv_array = _size_pv_array(
        tech_model, design_size=design_size, desired_dcac_ratio=desired_dcac_ratio
    )
    # inverters, num_inverters, num_parallel = _size_pv_array(
    #     tech_model, design_size=design_size, desired_dcac_ratio=desired_dcac_ratio
    # )
    print(
        f"\nRunning PySAM model {pysam_model_config} for design size = {design_size:.2f} kW...\n"
    )
    tech_model.value("inverter_count", size_pv_array["inverter_count"])
    tech_model.value("subarray1_nstrings", size_pv_array["number_strings"])
    
    # tech_model.value("inverter_count", num_inverters)
    # tech_model.value("subarray1_nstrings", num_parallel)

    # BUG Dec 2024: cec_gamma_r was a required parameter but wasn't provided by default config
    # Courtesy of Google AI:
    # The "cec_gamma_r" maximum power point temperature coefficient, typically
    # found on a solar panel datasheet, represents the percentage change
    # in maximum power output per degree Celsius change in temperature, usually
    # expressed as a negative value ranging from around -0.3% to -0.5% per degree Celsius.
    tech_model.value("cec_gamma_r", -0.3)

    for param, val in tech_model_kwargs.items():
        # if not isinstance(val, list):
        #     val = [val]
        tech_model.value(param, val)
    for param, val in cash_model_kwargs.items():
        if not isinstance(val, list):
            val = [val]
        cash_model.value(param, val)
    for mod in modules:
        mod.execute()
    # print("PySAM run finished.\n")
    # tech_model = modules[0]
    # cash_model = modules[3]
    # annual_energy = tech_model.Outputs.annual_energy
    # hourly_energy = tech_model.Outputs.gen
    # hourly_energy = [x if x > 0 else 0 for x in hourly_energy]
    # dc_capacity_factor = tech_model.Outputs.capacity_factor
    # ac_capacity_factor = tech_model.Outputs.capacity_factor_ac
    # lcoe_real = cash_model.Outputs.lcoe_real
    return tech_model, cash_model, size_pv_array


def _size_pv_array(tech_model, design_size=50, desired_dcac_ratio=1.2):

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

    # DC-connected battery parameters (assumed to use common inverter)
    batt_max_power_dc = 0
    if m["en_batt"] and m["batt_ac_or_dc"] == 0:
        batt_max_power_dc = m["batt_max_power"]

    if mod_vmp > 0:
        num_series = 0.5 * (inv_vmin + inv_vmax) / mod_vmp
        if inv_vdcmax > 0:
            while num_series > 0 and (num_series * mod_voc) > inv_vdcmax:
                num_series -= 1

    num_series = max(1, round(num_series))
    num_parallel = design_size * 1000 / (num_series * mod_power)
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
    total_module_area = total_module_area * 0.000247105 # m2 to [acre]
    land_area = total_module_area / tech_model.value("subarray1_gcr")  # assumes land_area_multiplier = 1
    total_ac_capacity = inv_power * num_inverters / 1000
    total_dc_inverter_capacity = inv_dc_power * num_inverters / 1000

    # check that the sizing was close to the desired sizes, otherwise error out
    nameplate_dc = total_modules * mod_power / 1000
    proposed_ratio = nameplate_dc / (num_inverters * inv_power / 1000)
    if abs(nameplate_dc - design_size) / design_size > 0.2:
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

    size_pv_array = {
        "inverter_count": num_inverters,
        "number_modules_per_string": num_series,
        "number_strings": num_parallel,
        "total_modules": total_modules,
        "total_module_area": total_module_area,
        "land_req": land_area,
        "system_capacity": nameplate_dc,
        "total_inverter_capacity": total_ac_capacity,
        "total_dc_inverter_capacity": total_dc_inverter_capacity,
    }

    if not all([num_inverters, num_parallel]):
        raise ValueError("One of inverters, num_inverters, num_parallel is None.")

    # return inverters, num_inverters, num_parallel
    return size_pv_array


def _flatten_dict(d):
    def get_key_values(d):
        for key, value in d.items():
            if isinstance(value, dict):
                yield from get_key_values(value)
            else:
                yield key, value

    return {key: value for (key, value) in get_key_values(d)}


def _load_config_files(modules):
    for file_name, module in zip(_config_files, modules):
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


def setup_and_run_pv(design_size, pysam_model_config, return_tech_model=False):
    modules = setup_pv_single_owner(pysam_model_config)
    tech_model, cash_model, size_pv_array = run_pv_single_owner(
        modules, pysam_model_config, design_size=design_size
    )
    annual_energy = tech_model.Outputs.annual_energy  # kWh
    # electricity_annual = tech_model.Outputs.gen # kW
    # hourly_energy = [x if x > 0 else 0 for x in hourly_energy]
    # dc_capacity_factor = tech_model.Outputs.capacity_factor
    # ac_capacity_factor = tech_model.Outputs.capacity_factor_ac

    cash_model.execute()
    # lcoe_real = cash_model.Outputs.lcoe_real

    result = size_pv_array
    result["annual_energy"] = annual_energy
    result["lcoe_real"] = cash_model.Outputs.lcoe_real
    result["ac_capacity_factor"] = tech_model.Outputs.capacity_factor_ac
    result["dc_capacity_factor"] = tech_model.Outputs.capacity_factor
    # result = {"electricity_annual": electricity_annual, "land_req": land_req}
    print(f"Design Size {design_size:.1f} kW completed.")
    print(f"\tAnnual Energy = {annual_energy:.2f} kWh")
    print(f"\tLand Required = {result['land_req']:.2f} acre\n")

    if return_tech_model:
        return result, tech_model
    else:
        return result


def run_pysam_kbhdp_pv(
    design_sizes=np.linspace(10, 10000, 10),
    dataset_filename=None,
    save_data=True,
    use_multiprocessing=True,
    processes=8,
    pysam_model_config="FlatPlatePVSingleOwner",
):
    if dataset_filename is None:
        raise ValueError("Must provide dataset_filename")
    dataset_filename = os.path.join(os.path.dirname(__file__), dataset_filename)
    data = list()
    if use_multiprocessing:
        df = pd.DataFrame(design_sizes, columns=["design_size"])
        # args = [(ds, *a) for ds in design_sizes]
        # print(args)
        # assert False
        with multiprocessing.Pool(processes=processes) as pool:
            args = [(ds, pysam_model_config) for ds in design_sizes]
            results = pool.starmap(setup_and_run_pv, args)
        df_results = pd.DataFrame(results)
        df = pd.concat(
            [
                df,
                df_results
            ],
            axis=1,
        )
    else:
        data = list()
        for ds in design_sizes:
            result = setup_and_run_pv(ds, pysam_model_config)
            data.append([ds] + [result.values()])
        df = pd.DataFrame(
            data, columns=["design_size"] + result.keys()
        )
    if save_data:
        df.to_pickle(dataset_filename)


if __name__ == "__main__":
    design_sizes = np.linspace(100000, 150000, 3)
    print(design_sizes)
    run_pysam_kbhdp_pv(
        design_sizes=design_sizes, dataset_filename="test.pkl", use_multiprocessing=False
    )
    