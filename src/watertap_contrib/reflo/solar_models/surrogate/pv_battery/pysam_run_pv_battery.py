import os
import PySAM
import PySAM.Pvsamv1 as pvsam
import PySAM.Grid as grid
import PySAM.Utilityrate5 as ur
import PySAM.Singleowner as so

from math import floor, ceil
import PySAM.ResourceTools as tools
from PySAM.BatteryTools import battery_model_sizing


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
default_weather_file = os.path.join(
    __location__, "data/tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
)
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


def _size_pv_array(tech_model, design_size=50, desired_dcac_ratio=1.2):

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
    total_module_area = total_module_area * 0.000247105  # m2 to [acre]
    land_area = total_module_area / tech_model.value(
        "subarray1_gcr"
    )  # assumes land_area_multiplier = 1
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

    return size_pv_array


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
    grid_model = grid.from_existing(tech_model, "PVBatterySingleOwner")
    rate_model = ur.from_existing(tech_model, "PVBatterySingleOwner")
    cash_model = so.from_existing(tech_model, "PVBatterySingleOwner")

    weather_data = tools.SAM_CSV_to_solar_data(weather_file)
    tech_model.SolarResource.solar_resource_data = weather_data

    modules = [
        tech_model,
        grid_model,
        rate_model,
        cash_model,
    ]

    return modules


def set_battery_options(
    tech_model,
    battery_kw,
    number_hrs,
    battery_volts=500,
    dispatch_manual_sched=None,
    **kwargs,
):
    # see: https://nrel-pysam.readthedocs.io/en/main/modules/Battery.html#PySAM.Battery.Battery.BatteryDispatch

    if dispatch_manual_sched is None:
        dispatch_manual_sched = default_dispatch_manual_sched

    batt = tech_model.BatteryDispatch
    tech_model.BatteryDispatch.batt_dispatch_auto_can_clipcharge = 1
    # of the 6 dispatch periods,
    # this will set to charge during period 1 (8am-8pm for default)
    # and discharge during period 2 (8pm-8am for default)
    tech_model.BatteryDispatch.dispatch_manual_charge = (1, 0, 0, 0, 0, 0)
    tech_model.BatteryDispatch.dispatch_manual_discharge = (0, 1, 0, 0, 0, 0)
    # this will enable grid charging for battery only during period 1
    tech_model.BatteryDispatch.dispatch_manual_gridcharge = (1, 0, 0, 0, 0, 0)
    # allow 100 percent charge from grid during period 1
    tech_model.BatteryDispatch.dispatch_manual_percent_gridcharge = (100,)
    # allow 100 percent discgarge to grid during period 2
    tech_model.BatteryDispatch.dispatch_manual_percent_discharge = (100,)
    # set the same dispatch schedule for weekday and weekend
    tech_model.BatteryDispatch.dispatch_manual_sched = dispatch_manual_sched
    tech_model.BatteryDispatch.dispatch_manual_sched_weekend = dispatch_manual_sched

    battery_kwh = battery_kw * number_hrs
    battery_model_sizing(tech_model, battery_kw, battery_kwh, battery_volts, **kwargs)


def run_pysam_pv_battery(design_size, battery_kw, number_hrs, weather_file=None, **kwargs):
    
    modules = create_pv_batt_modules(weather_file=weather_file, **kwargs)

    tech_model = modules[0]
    size_pv_array = _size_pv_array(tech_model, design_size=design_size)

    tech_model.value("inverter_count", size_pv_array["inverter_count"])
    tech_model.value("subarray1_nstrings", size_pv_array["number_strings"])
    tech_model.value("cec_gamma_r", -0.3)

    set_battery_options(tech_model, battery_kw, number_hrs, **kwargs)

    for mod in modules:
        mod.execute()

    tech_model = modules[0]

    print(f"Annual Energy = {tech_model.Outputs.annual_energy:.2f}")
    return modules, size_pv_array


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    
    design_size = 100000
    battery_kw = 50000
    number_hrs = 12

    modules, size_pv_array = run_pysam_pv_battery(design_size, battery_kw, number_hrs, weather_file=None)
    tech_model = modules[0]
    gen = tech_model.Outputs.gen

    start_day = 180
    start_i = start_day * 24
    end_day = 185
    end_i = end_day * 24

    fig, ax = plt.subplots()
    ax.plot(gen[start_i:end_i], label="gen")
    ax.plot(tech_model.Outputs.gen_without_battery[start_i:end_i], label="gen_without_battery")
    ax.plot(tech_model.Outputs.grid_to_batt[start_i:end_i], label="grid_to_batt")
    ax.plot(tech_model.Outputs.batt_to_grid[start_i:end_i], label="batt_to_grid")
    ax.legend()
    plt.show()