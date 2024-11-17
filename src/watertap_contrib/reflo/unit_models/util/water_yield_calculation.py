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

import math
from pandas import read_csv
import numpy as np

from idaes.core.util.exceptions import ConfigurationError
from watertap_contrib.reflo.unit_models.util.sw_props import (
    calculate_density,
    calculate_viscosity,
    calculate_specific_heat,
    calculate_thermal_conductivity,
)

days_in_year = 365
hours_per_day = 24
seconds_per_day = 86400
stefan_boltzmann = 5.6697e-8  # W / m2 / K4
thickness_insulation = 0.005  # Thickness of SS Insulation (m)
conductivity_insulation = 0.033  # Conductivity of SS Insulation (W/m.K)
thickness_glass = 0.004  # Thickness of glass m
conductivity_glass = 1.03  # Conductivity of glass W/m.K
gravity = 9.81  # Accelaration due to gravity (m/s^2)
density_nacl = 2165  # density of sodium chloride (kg/m^3)
maximum_solubility = 365  # maximum solubility of salt in water (g/l)

# Radiative properties
absorp_glass = 0.047  # Absorptivity of glass (-)
absorp_water = 0.20  # Absorptivity of water (-)
absorp_basin = 0.65  # Absorptivity of basin (-)

reflectivity_glass = 0.047  # Reflectivity of glass (-)
reflectivity_water = 0.08  # Reflectivity of Water (-)

emissivity_glass = 0.94  # Glass Emissivity (-)
emissivity_water = 0.95 * 1.0  # Water Emissivity (-)

# Effective emissivity of water to glass (-)
emissivity_effective = 1 / ((1 / emissivity_water) + (1 / emissivity_glass) - 1)

# No Attenuation factor considered
# Fraction of solar radiation absorbed by water (-)
absorp_effective_water = (
    (absorp_water)
    * (1 - absorp_glass)
    * (1 - reflectivity_glass)
    * (1 - reflectivity_water)
)
# Fraction of solar radiation absorbed by basin liner (-)
absorp_effective_basin = (
    (absorp_basin)
    * (1 - absorp_glass)
    * (1 - reflectivity_glass)
    * (1 - absorp_water)
    * (1 - reflectivity_water)
)
# Fraction of solar radiation absorbed by a glass cover (-)
absorp_effective_glass = (1 - reflectivity_glass) * absorp_glass

# Adaptive coeff heat transfer coeff with buoyancy
AA = 0.54
# Power of nondimensional numbers for heat transfer coeff with buoyancy
BB = 0.25


def create_input_arrays(
    weather_data,
    irradiance_threshold=0,  # irradiance values < threshold assumed to have negligible impact on calculation; W/m2
    irradiance_col=None,  # column from weather_data to use as irradiance input data
    temperature_col=None,  # column from weather_data to use as temperature input data
    wind_velocity_col=None,  # column from weather_data to use as wind velocity input data
):

    def generate_continuous_day_series():
        # Days in each month for non-leap year
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        continuous_day_series = []

        for day in range(1, 32):  # up to 31 days to cover all months
            for month_index, days in enumerate(days_in_month):
                if day <= days:  # check if the day exists in this month
                    # Calculate the day of the year for the nth day of each month
                    day_of_year = sum(days_in_month[:month_index]) + day
                    continuous_day_series.append(day_of_year)

        return continuous_day_series

    if irradiance_col is None:
        irradiance_col = "GHI"
    if temperature_col is None:
        temperature_col = "Tdry"
    if wind_velocity_col is None:
        wind_velocity_col = "Wspd"

    continuous_day_series = generate_continuous_day_series()

    ambient_temp_by_hr = []
    irradiance_by_hr = []
    wind_vel_by_hr = []

    # Collecting hourly weather data
    for day in continuous_day_series:
        if day == 365:
            day = 364
        start_row = day * 24
        for kk in range(24):
            # irr = float(weather_data.iloc[start_row + kk, 4]) + 10
            # temp = float(weather_data.iloc[start_row + kk, 7])
            # wind = float(weather_data.iloc[start_row + kk, 11])
            irr = float(weather_data[irradiance_col].loc[start_row + kk]) + 10
            temp = float(weather_data[temperature_col].loc[start_row + kk])
            wind = float(weather_data[wind_velocity_col].loc[start_row + kk])

            if irr > irradiance_threshold:
                ambient_temp_by_hr.append(temp)
                irradiance_by_hr.append(irr)
                wind_vel_by_hr.append(wind)

            else:
                ambient_temp_by_hr.append(temp)
                irradiance_by_hr.append(irr)
                wind_vel_by_hr.append(wind)

    ambient_temp_by_hr = np.array(ambient_temp_by_hr)
    irradiance_by_hr = np.array(irradiance_by_hr)
    wind_vel_by_hr = np.array(wind_vel_by_hr)

    return ambient_temp_by_hr, irradiance_by_hr, wind_vel_by_hr


def get_solar_still_daily_water_yield(
    blk,
    input_weather_file_path=None,  # path to input weather file
    initial_salinity=200,  # initial salinity of influent water; g/L
    initial_water_depth=0.01,  # initial depth of water in solar still basin; m
    length_basin=0.6,  # length of each side of basin (length=width); m
    **kwargs,  # kwargs can inlude column name arguments
):
    area_bottom_basin = length_basin**2  # Area of square basin (m^2)

    weather_data = read_csv(input_weather_file_path, skiprows=2)

    if not len(weather_data) >= 8760:
        err_msg = f"Water yield calculation for {blk.name} requires at least "
        err_msg += f"one year of hourly weather data, but the input dataset is "
        err_msg += f"{len(weather_data)} hours long."
        raise ConfigurationError(err_msg)

    if len(weather_data) > 8760:
        weather_data = weather_data.loc[:8760]

    ambient_temp_by_hr, irradiance_by_hr, wind_vel_by_hr = create_input_arrays(
        weather_data, **kwargs
    )

    len_data_hr = len(irradiance_by_hr)

    # Create arrays
    # Second variables
    depth = np.zeros(len_data_hr * 3600)
    salinity = np.zeros(len_data_hr * 3600)
    excess_salinity = np.zeros(len_data_hr * 3600)
    volume_scale_formation = np.zeros(len_data_hr * 3600)
    thickness_scale_formation = np.zeros(len_data_hr * 3600)
    evap_sw_mass = np.zeros(len_data_hr * 3600)
    irradiance = np.zeros(len_data_hr * 3600)
    salt_precipitated = np.zeros(len_data_hr * 3600)
    sw_mass = np.zeros(len_data_hr * 3600)
    fw_mass = np.zeros(len_data_hr * 3600)
    time = np.zeros(len_data_hr * 3600)
    wind_velocity = np.zeros(len_data_hr * 3600)
    ambient_temp = np.zeros(len_data_hr * 3600)
    basin_temp = np.zeros(len_data_hr * 3600)
    glass_temp = np.zeros(len_data_hr * 3600)
    sky_temp = np.zeros(len_data_hr * 3600)
    saltwater_temp = np.zeros(len_data_hr * 3600)

    # Converting hourly data into per second
    for hour in range(len_data_hr):
        start_idx = 3600 * hour
        end_idx = start_idx + 3600
        irradiance[start_idx:end_idx] = irradiance_by_hr[hour]
        wind_velocity[start_idx:end_idx] = wind_vel_by_hr[hour]
        ambient_temp[start_idx:end_idx] = ambient_temp_by_hr[hour]

    # Initializing Temperatures
    # Initial system is assumed to be in thermal equilibrium with ambient

    # Initial water temperature (°C)
    saltwater_temp[1] = ambient_temp_by_hr[1]
    # Initial basin temperature (°C)
    basin_temp[1] = ambient_temp_by_hr[1]
    # Initial glass temperature (°C)
    glass_temp[1] = ambient_temp_by_hr[1]

    initial_density = calculate_density(initial_salinity, saltwater_temp[1])

    salt_precipitated[1] = 0
    salt_precipitated[0] = salt_precipitated[1]

    salinity[1] = initial_salinity
    salinity[0] = initial_salinity
    depth[1] = initial_water_depth
    depth[0] = initial_water_depth
    sw_mass[1] = depth[1] * initial_density * area_bottom_basin  # kg
    fw_mass[1] = sw_mass[1] / (1 + salinity[1] / 1000)  # kg
    initial_mass_fw = fw_mass[1]
    salt_mass = (salinity[1] * fw_mass[1]) / 1000  # Mass of Sodium Chloride (kg)
    excess_salinity[1] = salinity[1]  # Salinity without maximum solublity (g/l)

    # Initial effective radiation temperature of the sky (deg C)
    sky_temp[1] = 0.0552 * ((ambient_temp[1]) ** 1.5)

    time[1] = 1

    for i in range(2, len(irradiance), 1):

        if depth[i - 1] <= 0 or fw_mass[i - 1] <= 0:
            depth[i - 1] = initial_water_depth
            salinity[i - 1] = initial_salinity
            sw_mass[i - 1] = depth[i - 1] * initial_density * area_bottom_basin
            fw_mass[i - 1] = sw_mass[i - 1] / (1 + salinity[i] / 1000)

        time[i] = time[i - 1] + 1

        # Avoiding singularities
        temp_diff_inside_basin = saltwater_temp[i - 1] - glass_temp[i - 1]
        if temp_diff_inside_basin <= 0:
            temp_diff_inside_basin = 0.01
        temp_diff_outside_basin = glass_temp[i - 1] - ambient_temp[i - 1]
        if temp_diff_outside_basin <= 0:
            temp_diff_outside_basin = 0.01

        # Effective radiation temperature of the sky
        sky_temp[i] = 0.0552 * ((ambient_temp[i]) ** 1.5)

        # Perimeter x depth of water (m^2)
        area_side_water = (2 * (2 * length_basin)) * depth[i - 1]

        density = calculate_density(salinity[i - 1], saltwater_temp[i - 1])

        dynamic_visc = calculate_viscosity(salinity[i - 1], saltwater_temp[i - 1])

        specific_heat = calculate_specific_heat(salinity[i - 1], saltwater_temp[i - 1])

        thermal_conductivity = calculate_thermal_conductivity(
            salinity[i - 1], saltwater_temp[i - 1]
        )

        kinem_visc_sw = dynamic_visc / density
        # Prandtl number for water (-)
        Pr = (specific_heat * dynamic_visc) / thermal_conductivity
        # Latent heat of vaporization of pure water J/kg
        freshwater_vap_latent_heat = (2501.67 - 2.389 * saltwater_temp[i - 1]) * 1000

        # Calculation of partial saturated vapor pressure of saltwater
        # According to parametric analysis and available literature,
        # the partial vapor pressure plays a major role in the evaporation of water.

        # A coeff obtained from the water molar fraction in salt solutions from 0-350 g/l Paper: Kokya and Kokya
        water_activity = (-0.000566 * salinity[i - 1]) + 0.99853070
        # Partial saturated vapor pressure at a saltwater temperature (N/m^2)
        sw_partial_vap_press = water_activity * math.exp(
            25.317 - (5144 / (saltwater_temp[i - 1] + 273))
        )
        # Partial saturated vapor pressure at glass cover temperature (N/m^2)
        partial_vap_press_at_glass = math.exp(
            25.317 - (5144 / (glass_temp[i - 1] + 273))
        )
        # Coefficient of volume expansion (1/°C) correlation gotten from Zhutovsky and Kovler (2015)
        beta = 1e-6 * (
            -0.000006 * saltwater_temp[i - 1] ** 4
            + 0.001667 * saltwater_temp[i - 1] ** 3
            - 0.197796 * saltwater_temp[i - 1] ** 2
            + 16.862446 * saltwater_temp[i - 1]
            - 64.319951
        )
        # Grashof number
        Gr = abs(
            (
                gravity
                * beta
                * (basin_temp[i - 1] - saltwater_temp[i - 1])
                * (depth[i - 1] ** 3)
            )
            / (kinem_visc_sw**2)
        )
        # Heat transfer coeff of water layer
        water_heat_trans_coeff = abs(
            (thermal_conductivity / depth[i - 1]) * AA * (Gr * Pr) ** BB
        )
        # Convective heat transfer coeff (W/m^2.°C) (Dunkel)
        conv_heat_trans_coeff_water_glass = 0.884 * (
            (
                abs(
                    (
                        (saltwater_temp[i - 1] - glass_temp[i - 1])
                        + (
                            (
                                (sw_partial_vap_press - partial_vap_press_at_glass)
                                * (saltwater_temp[i - 1] + 273.15)
                            )
                            / (268900 - sw_partial_vap_press)
                        )
                    )
                )
            )
            ** (1 / 3)
        )
        # Radiative heat transfer coeff from basin water to glass cover (W/m^2.°C)
        rad_heat_transf_coeff_water_glass = abs(
            emissivity_water
            * stefan_boltzmann
            * (
                (
                    ((saltwater_temp[i - 1] + 273) ** 2)
                    + ((glass_temp[i - 1] + 273) ** 2)
                )
                * (saltwater_temp[i - 1] + glass_temp[i - 1] + 546)
            )
        )
        # Evaporative heat transfer coeff from basin water to glass cover (W/m2.°C)
        evap_heat_trans_coeff_water_glass = abs(
            0.01628
            * conv_heat_trans_coeff_water_glass
            * (
                (sw_partial_vap_press - partial_vap_press_at_glass)
                / (temp_diff_inside_basin)
            )
        )
        # Total heat transfer coeff from basin water to glass cover (W/m2°C)
        tot_heat_trans_coeff_water_glass = (
            conv_heat_trans_coeff_water_glass
            + rad_heat_transf_coeff_water_glass
            + evap_heat_trans_coeff_water_glass
        )
        # Radiative heat transfer coeff from glass cover to ambient (W/m^2.°C)
        rad_heat_trans_coeff_glass_ambient = (
            stefan_boltzmann
            * emissivity_glass
            * (
                (((glass_temp[i - 1] + 273) ** 4) - ((sky_temp[i - 1] + 273) ** 4))
                / (temp_diff_outside_basin)
            )
        )
        if wind_velocity[i] > 5:
            # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            conv_heat_trans_coeff_basin_ambient = 2.8 + (3.0 * wind_velocity[i])
            # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
            conv_heat_trans_coeff_glass_ambient = 2.8 + (3.0 * wind_velocity[i])
        else:
            # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            conv_heat_trans_coeff_basin_ambient = 2.8 + (3.8 * wind_velocity[i])
            # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
            conv_heat_trans_coeff_glass_ambient = 2.8 + (3.8 * wind_velocity[i])
        # Total heat loss coeff from the glass cover to the outer atmosphere
        tot_heat_trans_coeff_glass_ambient = (
            conv_heat_trans_coeff_glass_ambient + rad_heat_trans_coeff_glass_ambient
        )
        # Heat loss coefficient from basin liner to the atmosphere
        tot_heat_trans_coeff_basin_ambient = 1 / (
            (thickness_insulation / conductivity_insulation)
            + (1 / (conv_heat_trans_coeff_basin_ambient))
        )

        # Effective overall absorptivity for energy balance equation
        effective_absorp = (
            (
                absorp_effective_basin
                * (
                    water_heat_trans_coeff
                    / (
                        water_heat_trans_coeff
                        + tot_heat_trans_coeff_basin_ambient
                        + conv_heat_trans_coeff_basin_ambient
                    )
                )
            )
            + (absorp_effective_water)
            + (
                (absorp_effective_glass)
                * (
                    tot_heat_trans_coeff_water_glass
                    / (
                        tot_heat_trans_coeff_water_glass
                        + tot_heat_trans_coeff_glass_ambient
                    )
                )
            )
        )

        # Calculation of overall heat transfer coefficients
        # Overall heat loss coefficient (W/m^2.°C)
        overall_heat_loss_coeff_glass_surr = (
            (conductivity_glass / thickness_glass)
            * (tot_heat_trans_coeff_glass_ambient)
        ) / (
            (conductivity_glass / thickness_glass) + tot_heat_trans_coeff_glass_ambient
        )
        # Overall bottom heat transfer coefficient between the water mass and the surroundings (W/m^2.°C)
        overall_bottom_heat_trans_coeff_water_mass_surr = (
            tot_heat_trans_coeff_water_glass * overall_heat_loss_coeff_glass_surr
        ) / (tot_heat_trans_coeff_water_glass + overall_heat_loss_coeff_glass_surr)
        # Overall bottom heat transfer coefficient from bottom to ambient (W/m^2.°C)
        overall_bottom_heat_loss_coeff_water_mass_surr = (
            water_heat_trans_coeff * tot_heat_trans_coeff_basin_ambient
        ) / (water_heat_trans_coeff + tot_heat_trans_coeff_basin_ambient)
        overall_side_heat_loss_coefficient = (
            area_side_water / area_bottom_basin
        ) * overall_bottom_heat_loss_coeff_water_mass_surr

        overall_heat_trans_coeff_basin_surr = (
            overall_bottom_heat_loss_coeff_water_mass_surr
            + overall_side_heat_loss_coefficient
        )
        overall_external_heat_trans_loss_coeff = (
            overall_bottom_heat_trans_coeff_water_mass_surr
            + overall_heat_trans_coeff_basin_surr
        )

        # Present [i] temperature calculation
        grouping_term = overall_external_heat_trans_loss_coeff / (
            sw_mass[i - 1] * specific_heat
        )

        time_dependent_term = (
            (effective_absorp * irradiance[i])
            + (overall_external_heat_trans_loss_coeff * ambient_temp[i])
        ) / (sw_mass[i - 1] * specific_heat)

        glass_temp[i] = (
            (absorp_effective_glass * irradiance[i])
            + (tot_heat_trans_coeff_water_glass * saltwater_temp[i - 1])
            + (overall_heat_loss_coeff_glass_surr * ambient_temp[i])
        ) / (tot_heat_trans_coeff_water_glass + overall_heat_loss_coeff_glass_surr)

        saltwater_temp[i] = (time_dependent_term / grouping_term) * (
            1 - np.exp(-grouping_term * time[i])
        ) + (saltwater_temp[i] * np.exp(-grouping_term * time[i]))

        basin_temp[i] = (
            (absorp_effective_basin * irradiance[i])
            + (water_heat_trans_coeff * saltwater_temp[i - 1])
            + (
                (
                    tot_heat_trans_coeff_basin_ambient
                    + conv_heat_trans_coeff_basin_ambient
                )
                * basin_temp[i - 1]
            )
        ) / (
            water_heat_trans_coeff
            + tot_heat_trans_coeff_basin_ambient
            + conv_heat_trans_coeff_basin_ambient
        )

        # Evaporation estimation of freshwater and saltwater
        # evap_fw_mass = (
        #     area_bottom_basin
        #     * evap_heat_trans_coeff_water_glass
        #     * (saltwater_temp[i] - glass_temp[i])
        #     * 3600
        # ) / freshwater_vap_latent_heat  # Distillated (kg)

        evap_fw_mass = (
            area_bottom_basin
            * evap_heat_trans_coeff_water_glass
            * (temp_diff_inside_basin)
        ) / freshwater_vap_latent_heat  # Distillated (kg)

        # Distillated saltwater conversion (Morton) (kg)
        # evap_sw_mass[i] = evap_fw_mass / (1 + salinity[i - 1] / 1e3)
        evap_sw_mass[i] = evap_fw_mass / (1 + salinity[i - 1] / 1e3)

        # Freshwater (kg) this iteration
        fw_mass[i] = fw_mass[i - 1] - evap_sw_mass[i]
        # Saltwater (kg) this iteration
        sw_mass[i] = fw_mass[i] + salt_mass
        # Water depth this iteration
        depth[i] = sw_mass[i] / (density * area_bottom_basin)  # Current depth (m)

        if salinity[i - 1] >= maximum_solubility:
            salinity[i] = maximum_solubility

            # Excess salinity (assuming no saturation possible) (g/l)
            excess_salinity[i] = (salt_mass * 1000) / fw_mass[i]

            if excess_salinity[i] < maximum_solubility:
                salt_precipitated[i] = 0
                volume_scale_formation[i] = 0
                thickness_scale_formation[i] = 0
            else:
                # Mass precipitated (kg)
                salt_precipitated[i] = (
                    fw_mass[i] * (excess_salinity[i] - maximum_solubility) / 1000
                )
                # Volume of scale formation (m3)
                volume_scale_formation[i] = salt_precipitated[i] / density_nacl
                # Thickness of scale formation (m)
                thickness_scale_formation[i] = (
                    volume_scale_formation[i] / area_bottom_basin
                )

        else:
            salinity[i] = (salt_mass * 1000) / fw_mass[i]
            excess_salinity[i] = salinity[i]

        if depth[i] <= 0 or fw_mass[i] <= 0:
            # At this point either the depth is negative
            # or the amount of freshwater available is negative
            # signaling we have reached the time required for one ZLD cycle for one solar still

            # Number of seconds for a single ZLD cycle
            num_seconds_for_zld_cycle = i
            break

    num_zld_cycles_per_year = (
        len(weather_data) * 3600
    ) / num_seconds_for_zld_cycle  # (s / year) / s = year**-1

    # Total water productivity in year [kg water per m2 area per year]
    annual_water_yield = (initial_mass_fw * num_zld_cycles_per_year) / area_bottom_basin
    # Daily water yield [kg water per m2 area per day]
    daily_water_yield = annual_water_yield / days_in_year

    # Volume_salt = zld_counter * max(
    #     volume_scale_formation
    # )  # total volume of salt precipitated in the year (m3)

    print(annual_water_yield, daily_water_yield, num_zld_cycles_per_year)
    return daily_water_yield


if __name__ == "__main__":
    get_solar_still_daily_water_yield_zld(
        input_weather_file_path="/Users/ksitterl/Documents/SETO/models/solar_still_zld/SS_ZLD_Model 2/Data/TMY2 SAM CSV/data.csv"
    )
