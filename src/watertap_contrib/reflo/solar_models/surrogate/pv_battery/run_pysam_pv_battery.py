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
import multiprocessing
import pandas as pd
from itertools import product
from math import floor, ceil

import PySAM.Pvsamv1 as pvsam
import PySAM.Grid as grid
import PySAM.ResourceTools as rtools
from PySAM.BatteryTools import battery_model_sizing


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
default_weather_file = os.path.join(
    __location__, "data/test_pv_battery_weather_data.csv"
)
# TODO: need to set the actual model dispatch schedule defaults here;
# the ones below are not the defaults
# default dispatch schedule:
# - charging from 8 am 8 pm
# - discharging from 8 pm to 8 am
# rows are for months of year
# columns are for hours of day
# each value represents the dispatch period assigned to that (month, hour)
# e.g., if 2 is in (month, hour) = (6, 12), the 12 hour of all days in June
# is assigned to dispatch period 2
default_dispatch_manual_sched = (
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
)
# of the 6 dispatch periods,
# this will set to charge during period 1 (8am-8pm for default)
dispatch_manual_charge_default = (1, 0, 0, 0, 0, 0)
# and discharge during period 2 (8pm-8am for default)
dispatch_manual_discharge_default = (0, 1, 0, 0, 0, 0)
# this will enable grid charging for battery only during period 1
dispatch_manual_gridcharge_default = (1, 0, 0, 0, 0, 0)
# allow 100 percent charge from grid during period 1
dispatch_manual_percent_gridcharge_default = (100,)
# allow 100 percent discharge to grid during period 2
dispatch_manual_percent_discharge_default = (100,)


def size_pv_array(tech_model, system_capacity=50, desired_dcac_ratio=1.2):

    # Sizing rules
    # 1. Voc < Vdcmax
    # 2. Vmp > Vmin
    # 3. Vmp < Vmax
    # 4. num series * num_parallel is about desired array size (num_parallel = desired / (num series * mod_power))
    # 5. num inverters is about desired array size (num_inv = num_series * num_parallel * mod_power) / inv_power)

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


def _flatten_dict(d):
    def get_key_values(d):
        for key, value in d.items():
            if isinstance(value, dict):
                yield from get_key_values(value)
            else:
                yield key, value

    return {key: value for (key, value) in get_key_values(d)}


def _spe_power(spe_eff_level, spe_rad_level, spe_area):
    return spe_eff_level / 100 * spe_rad_level * spe_area


def create_pv_batt_modules(weather_file=None):
    if weather_file is None:
        weather_file = default_weather_file

    tech_model = pvsam.default("PVBatterySingleOwner")
    # set the meter position
    # 0 = Behind-the-meter
    # 1 = Front-of-meter
    # tech_model.BatterySystem.batt_meter_position = 0
    grid_model = grid.from_existing(tech_model, "PVBatterySingleOwner")

    weather_data = rtools.SAM_CSV_to_solar_data(weather_file)
    tech_model.SolarResource.solar_resource_data = weather_data
    # tech_model.SolarResource.solar_resource_file = weather_file

    modules = [
        tech_model,
        grid_model,
    ]

    return modules


def set_battery_options(
    tech_model,
    battery_power,
    hours_storage,
    battery_volts=500,
    dispatch_manual_sched=None,
    dispatch_manual_charge=None,
    dispatch_manual_discharge=None,
    dispatch_manual_gridcharge=None,
    dispatch_manual_percent_gridcharge=None,
    dispatch_manual_percent_discharge=None,
    **kwargs,
):
    # TODO: This does not currently get called;
    # need to fix the defaults at top of file
    # see: https://nrel-pysam.readthedocs.io/en/main/modules/Battery.html#PySAM.Battery.Battery.BatteryDispatch

    # if dispatch_manual_sched is None:
    #     dispatch_manual_sched = default_dispatch_manual_sched
    # if dispatch_manual_charge is None:
    #     dispatch_manual_charge = dispatch_manual_charge_default
    # if dispatch_manual_discharge is None:
    #     dispatch_manual_discharge = dispatch_manual_discharge_default
    # if dispatch_manual_gridcharge is None:
    #     dispatch_manual_gridcharge = dispatch_manual_gridcharge_default
    # if dispatch_manual_percent_gridcharge is None:
    #     dispatch_manual_percent_gridcharge = dispatch_manual_percent_gridcharge_default
    # if dispatch_manual_percent_discharge is None:
    #     dispatch_manual_percent_discharge = dispatch_manual_percent_discharge_default

    batt = tech_model.BatteryDispatch
    tech_model.BatteryDispatch.batt_dispatch_auto_can_clipcharge = 1

    tech_model.BatteryDispatch.dispatch_manual_charge = dispatch_manual_charge
    tech_model.BatteryDispatch.dispatch_manual_discharge = dispatch_manual_discharge
    tech_model.BatteryDispatch.dispatch_manual_gridcharge = dispatch_manual_gridcharge
    tech_model.BatteryDispatch.dispatch_manual_percent_gridcharge = (
        dispatch_manual_percent_gridcharge
    )
    tech_model.BatteryDispatch.dispatch_manual_percent_discharge = (
        dispatch_manual_percent_discharge
    )
    # set the same dispatch schedule for weekday and weekend
    tech_model.BatteryDispatch.dispatch_manual_sched = dispatch_manual_sched
    tech_model.BatteryDispatch.dispatch_manual_sched_weekend = dispatch_manual_sched

    battery_kwh = battery_power * hours_storage
    battery_model_sizing(
        tech_model, battery_power, battery_kwh, battery_volts, **kwargs
    )


def run_pysam_pv_battery(
    system_capacity, battery_power, hours_storage, weather_file=None, **kwargs
):

    modules = create_pv_batt_modules(weather_file=weather_file, **kwargs)

    tech_model = modules[0]
    pv_array_design = size_pv_array(tech_model, system_capacity=system_capacity)

    tech_model.value("inverter_count", pv_array_design["inverter_count"])
    tech_model.value("subarray1_nstrings", pv_array_design["number_strings"])
    # tech_model.value("cec_gamma_r", -0.3)

    # TODO: need to set the actual model defaults here;
    # those at the top of this file are not the defaults.
    # set_battery_options(tech_model, battery_power, hours_storage, **kwargs)

    print(
        f"Running:\n\tPV Design Size {system_capacity:.1f} kW\n\tBattery Size {battery_power:.1f} kW\n\tBattery Storage {hours_storage:.1f} hours"
    )
    for mod in modules:
        mod.execute()
    print(f"Completed:\tAnnual Energy = {tech_model.Outputs.annual_energy:.2f}")

    return modules, pv_array_design


def setup_and_run_pv_battery(
    system_capacity,
    battery_power,
    hours_storage,
    weather_file=None,
    return_tech_model=False,
    **kwargs,
):
    """
    Setup and run the PV + Battery model with the given parameters.
    """
    if weather_file is None:
        weather_file = default_weather_file

    modules, pv_array_design = run_pysam_pv_battery(
        system_capacity,
        battery_power,
        hours_storage,
        weather_file,
        **kwargs,
    )

    tech_model = modules[0]

    result = pv_array_design
    gen_1yr = tech_model.Outputs.gen[:8760]
    gen_without_batt_1yr = tech_model.Outputs.gen_without_battery[:8760]
    net_electricity_produced = tech_model.Outputs.annual_energy - sum(
        tech_model.Outputs.batt_annual_charge_from_grid
    )
    result["electricity_annual"] = net_electricity_produced
    result["batt_annual_charge_from_grid"] = sum(
        tech_model.Outputs.batt_annual_charge_from_grid
    )
    result["batt_annual_charge_from_system"] = sum(
        tech_model.Outputs.batt_annual_charge_from_system
    )
    result["batt_annual_discharge_energy"] = sum(
        tech_model.Outputs.batt_annual_discharge_energy
    )
    result["ac_capacity_factor"] = tech_model.Outputs.capacity_factor_ac
    result["dc_capacity_factor"] = tech_model.Outputs.capacity_factor
    result["annual_gen_without_battery"] = sum(gen_without_batt_1yr)

    if return_tech_model:
        return result, tech_model
    else:
        return result


def generate_pv_battery_data(
    system_capacities=[100000, 200000],
    battery_powers=[10000, 60000],
    hours_storages=[6, 12],
    weather_file=None,
    save_data=True,
    use_multiprocessing=True,
    processes=8,
    dataset_filename=None,
    **kwargs,
):
    """
    Generate PV + Battery data for a range of conditions.
    """

    if dataset_filename is None:
        # assume it is run for testing purposes
        dataset_filename = os.path.join(__location__, "data/test_data.pkl")

    if weather_file is None:
        weather_file = default_weather_file

    combos = list(product(system_capacities, battery_powers, hours_storages))
    df = pd.DataFrame(
        combos, columns=["system_capacity", "battery_power", "hours_storage"]
    )
    if use_multiprocessing:

        with multiprocessing.Pool(processes=processes) as pool:
            args_in = [(*combo,) for combo in combos]
            results = pool.starmap(setup_and_run_pv_battery, args_in)
        df_results = pd.DataFrame(results)
    else:
        results = []
        for system_capacity, battery_power, hours_storage in combos:
            result = setup_and_run_pv_battery(
                system_capacity,
                battery_power,
                hours_storage,
                weather_file=weather_file,
                **kwargs,
            )
            results.append(result)
        df_results = pd.DataFrame(results)
    df = pd.concat([df, df_results], axis=1)

    if save_data:
        df.to_pickle(dataset_filename)

    return df


if __name__ == "__main__":

    df = generate_pv_battery_data()
    print(df.head(20))
    os.remove(os.path.join(__location__, "data/test_data.pkl"))

    # For generating test data:
    # design_sizes = np.linspace(100000, 500000, 5)
    # battery_powers = np.linspace(10000, 60000, 5)
    # hours_storages = [3, 6, 8, 12, 24]
    # df = generate_pv_battery_data(
    #     design_sizes=design_sizes,
    #     battery_powers=battery_powers,
    #     hours_storages=hours_storages,
    # )
