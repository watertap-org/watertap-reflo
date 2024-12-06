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

__author__ = "Kurban Sitterley, Nikhil Dani, Erick Moreno Resendiz"

days_in_year = 365
hours_per_day = 24
seconds_per_day = 86400
stefan_boltzmann = 5.6697e-8  # W / m2 / K4
thickness_insulation = 0.005  # Thickness of SS Insulation (m)
conductivity_insulation = 0.033  # Conductivity of SS Insulation (W/m.K)
thickness_glass = 0.004  # Thickness of glass m
conductivity_glass = 1.03  # Conductivity of glass W/m.K
gravity = 9.81  # Accelaration due to gravity (m/s^2)
density_nacl = 2165  # blk.density[i] of sodium chloride (kg/m^3)
maximum_solubility = 365  # maximum solubility of salt in water (g/L)

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
    blk,
    irradiance_threshold=0,  # blk.irradiance values < threshold assumed to have negligible impact on calculation; W/m2
    irradiance_col=None,  # column from weather_data to use as blk.irradiance input data
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

    # Default column names are those in files as downloaded from https://sam.nrel.gov/weather-data.html
    # Requires an irradiance column (W/m2), temperature column (°C), and wind velocity column (m/s)
    if irradiance_col is None:
        irradiance_col = "GHI"
    if temperature_col is None:
        temperature_col = "Tdry"
    if wind_velocity_col is None:
        wind_velocity_col = "Wspd"

    continuous_day_series = generate_continuous_day_series()

    blk.ambient_temp_by_hr = []
    blk.irradiance_by_hr = []
    blk.wind_vel_by_hr = []

    # Collecting hourly weather data
    for day in continuous_day_series:
        if day == 365:
            day = 364
        start_row = day * 24
        for kk in range(24):
            irr = float(blk.weather_data[irradiance_col].loc[start_row + kk]) + 10
            temp = float(blk.weather_data[temperature_col].loc[start_row + kk])
            wind = float(blk.weather_data[wind_velocity_col].loc[start_row + kk])

            if irr > irradiance_threshold:
                blk.ambient_temp_by_hr.append(temp)
                blk.irradiance_by_hr.append(irr)
                blk.wind_vel_by_hr.append(wind)

            else:
                blk.ambient_temp_by_hr.append(temp)
                blk.irradiance_by_hr.append(irr)
                blk.wind_vel_by_hr.append(wind)

    blk.ambient_temp_by_hr = np.array(blk.ambient_temp_by_hr)
    blk.irradiance_by_hr = np.array(blk.irradiance_by_hr)
    blk.wind_vel_by_hr = np.array(blk.wind_vel_by_hr)

    return blk.ambient_temp_by_hr, blk.irradiance_by_hr, blk.wind_vel_by_hr


def get_solar_still_daily_water_yield(
    blk,
    input_weather_file_path=None,  # path to input weather file
    initial_salinity=200,  # initial salinity of influent water; g/L
    initial_water_depth=0.01,  # initial depth of water in solar still basin; m
    length_basin=0.6,  # length of each side of basin (length=width); m
    **kwargs,  # kwargs can inlude column name arguments
):
    area_bottom_basin = length_basin**2  # Area of square basin (m^2)
    blk.initial_salinity = initial_salinity
    blk.initial_water_depth = initial_water_depth

    blk.weather_data = read_csv(input_weather_file_path, skiprows=2)

    if not len(blk.weather_data) >= 8760:
        err_msg = f"Water yield calculation for {blk.name} requires at least "
        err_msg += f"one year of hourly weather data, but the input dataset is "
        err_msg += f"only {len(blk.weather_data)} hours long."
        raise ConfigurationError(err_msg)

    if len(blk.weather_data) > 8760:
        blk.weather_data = blk.weather_data.loc[:8760]

    blk.ambient_temp_by_hr, blk.irradiance_by_hr, blk.wind_vel_by_hr = (
        create_input_arrays(blk, **kwargs)
    )

    len_data_hr = len(blk.irradiance_by_hr)

    # Create arrays
    # Values calculated per second
    blk.time = np.zeros(len_data_hr * 3600)
    blk.area_side_water = np.zeros(len_data_hr * 3600)
    blk.depth = np.zeros(len_data_hr * 3600)
    blk.salinity = np.zeros(len_data_hr * 3600)
    blk.excess_salinity = np.zeros(len_data_hr * 3600)
    blk.salt_precipitated = np.zeros(len_data_hr * 3600)
    blk.sw_mass = np.zeros(len_data_hr * 3600)
    blk.fw_mass = np.zeros(len_data_hr * 3600)
    blk.evap_sw_mass = np.zeros(len_data_hr * 3600)
    blk.evap_fw_mass = np.zeros(len_data_hr * 3600)
    blk.volume_scale_formation = np.zeros(len_data_hr * 3600)
    blk.thickness_scale_formation = np.zeros(len_data_hr * 3600)

    blk.irradiance = np.zeros(len_data_hr * 3600)
    blk.wind_velocity = np.zeros(len_data_hr * 3600)

    blk.ambient_temp = np.zeros(len_data_hr * 3600)
    blk.basin_temp = np.zeros(len_data_hr * 3600)
    blk.glass_temp = np.zeros(len_data_hr * 3600)
    blk.sky_temp = np.zeros(len_data_hr * 3600)
    blk.saltwater_temp = np.zeros(len_data_hr * 3600)

    blk.density = np.zeros(len_data_hr * 3600)
    blk.dynamic_visc = np.zeros(len_data_hr * 3600)
    blk.kinem_visc_sw = np.zeros(len_data_hr * 3600)
    blk.specific_heat = np.zeros(len_data_hr * 3600)
    blk.thermal_conductivity = np.zeros(len_data_hr * 3600)
    blk.Pr = np.zeros(len_data_hr * 3600)
    blk.Gr = np.zeros(len_data_hr * 3600)
    blk.freshwater_vap_latent_heat = np.zeros(len_data_hr * 3600)
    blk.water_activity = np.zeros(len_data_hr * 3600)
    blk.sw_partial_vap_press = np.zeros(len_data_hr * 3600)
    blk.partial_vap_press_at_glass = np.zeros(len_data_hr * 3600)

    blk.effective_absorp = np.zeros(len_data_hr * 3600)
    blk.water_heat_trans_coeff = np.zeros(len_data_hr * 3600)
    blk.conv_heat_trans_coeff_water_glass = np.zeros(len_data_hr * 3600)
    blk.rad_heat_transf_coeff_water_glass = np.zeros(len_data_hr * 3600)
    blk.evap_heat_trans_coeff_water_glass = np.zeros(len_data_hr * 3600)
    blk.tot_heat_trans_coeff_water_glass = np.zeros(len_data_hr * 3600)
    blk.rad_heat_trans_coeff_glass_ambient = np.zeros(len_data_hr * 3600)
    blk.conv_heat_trans_coeff_basin_ambient = np.zeros(len_data_hr * 3600)
    blk.conv_heat_trans_coeff_glass_ambient = np.zeros(len_data_hr * 3600)
    blk.tot_heat_trans_coeff_glass_ambient = np.zeros(len_data_hr * 3600)
    blk.tot_heat_trans_coeff_basin_ambient = np.zeros(len_data_hr * 3600)
    blk.overall_external_heat_trans_loss_coeff = np.zeros(len_data_hr * 3600)
    blk.overall_heat_loss_coeff_glass_surr = np.zeros(len_data_hr * 3600)
    blk.overall_bottom_heat_trans_coeff_water_mass_surr = np.zeros(len_data_hr * 3600)
    blk.overall_bottom_heat_loss_coeff_water_mass_surr = np.zeros(len_data_hr * 3600)
    blk.overall_side_heat_loss_coefficient = np.zeros(len_data_hr * 3600)
    blk.overall_heat_trans_coeff_basin_surr = np.zeros(len_data_hr * 3600)

    # Initializing Temperatures
    # Initial system is assumed to be in thermal equilibrium with ambient

    # Initial water temperature (°C)
    blk.saltwater_temp[1] = blk.ambient_temp_by_hr[1]
    # Initial basin temperature (°C)
    blk.basin_temp[1] = blk.ambient_temp_by_hr[1]
    # Initial glass temperature (°C)
    blk.glass_temp[1] = blk.ambient_temp_by_hr[1]

    blk.area_side_water[0] = (2 * (2 * length_basin)) * blk.initial_water_depth

    blk.density[0] = calculate_density(blk.initial_salinity, blk.saltwater_temp[1])

    blk.dynamic_visc[0] = calculate_viscosity(
        blk.initial_salinity, blk.saltwater_temp[1]
    )

    blk.specific_heat[0] = calculate_specific_heat(
        blk.initial_salinity, blk.saltwater_temp[1]
    )

    blk.thermal_conductivity[0] = calculate_thermal_conductivity(
        blk.initial_salinity, blk.saltwater_temp[1]
    )

    blk.salt_precipitated[1] = 0
    blk.salt_precipitated[0] = blk.salt_precipitated[1]

    blk.salinity[1] = blk.initial_salinity
    blk.salinity[0] = blk.initial_salinity
    blk.depth[1] = blk.initial_water_depth
    blk.depth[0] = blk.initial_water_depth
    blk.sw_mass[1] = blk.depth[1] * blk.density[0] * area_bottom_basin  # kg
    blk.fw_mass[1] = blk.sw_mass[1] / (1 + blk.salinity[1] / 1000)  # kg
    blk.initial_mass_fw = blk.fw_mass[1]
    blk.salt_mass = (
        blk.salinity[1] * blk.fw_mass[1]
    ) / 1000  # Mass of Sodium Chloride (kg)
    blk.excess_salinity[1] = blk.salinity[1]  # Salinity without maximum solublity (g/l)

    # Initial effective radiation temperature of the sky (°C)
    if blk.ambient_temp[1] <= 0:
        blk.sky_temp[1] = blk.ambient_temp[1]
    else:
        blk.sky_temp[1] = 0.0552 * ((blk.ambient_temp[1]) ** 1.5)

    blk.time[1] = 1

    # Converting hourly data into per second
    for hour in range(len_data_hr):
        start_idx = 3600 * hour
        end_idx = start_idx + 3600
        blk.irradiance[start_idx:end_idx] = blk.irradiance_by_hr[hour]
        blk.wind_velocity[start_idx:end_idx] = blk.wind_vel_by_hr[hour]
        blk.ambient_temp[start_idx:end_idx] = blk.ambient_temp_by_hr[hour]

    for i in range(2, len(blk.irradiance), 1):

        if blk.depth[i - 1] <= 0 or blk.fw_mass[i - 1] <= 0:
            blk.depth[i - 1] = blk.initial_water_depth
            blk.salinity[i - 1] = blk.initial_salinity
            blk.sw_mass[i - 1] = blk.depth[i - 1] * blk.density[0] * area_bottom_basin
            blk.fw_mass[i - 1] = blk.sw_mass[i - 1] / (1 + blk.salinity[i] / 1000)

        blk.time[i] = blk.time[i - 1] + 1

        # Avoiding singularities
        blk.temp_diff_inside_basin = blk.saltwater_temp[i - 1] - blk.glass_temp[i - 1]
        if blk.temp_diff_inside_basin <= 0:
            blk.temp_diff_inside_basin = 0.01
        blk.temp_diff_outside_basin = blk.glass_temp[i - 1] - blk.ambient_temp[i - 1]
        if blk.temp_diff_outside_basin <= 0:
            blk.temp_diff_outside_basin = 0.01

        # Effective radiation temperature of the sky
        if blk.ambient_temp[i] <= 0:
            blk.sky_temp[i] == blk.ambient_temp[i]
        else:
            blk.sky_temp[i] = 0.0552 * ((blk.ambient_temp[i]) ** 1.5)

        # Perimeter x depth of water (m^2)
        blk.area_side_water[i] = (2 * (2 * length_basin)) * blk.depth[i - 1]

        blk.density[i] = calculate_density(
            blk.salinity[i - 1], blk.saltwater_temp[i - 1]
        )

        blk.dynamic_visc[i] = calculate_viscosity(
            blk.salinity[i - 1], blk.saltwater_temp[i - 1]
        )

        blk.specific_heat[i] = calculate_specific_heat(
            blk.salinity[i - 1], blk.saltwater_temp[i - 1]
        )

        blk.thermal_conductivity[i] = calculate_thermal_conductivity(
            blk.salinity[i - 1], blk.saltwater_temp[i - 1]
        )

        # Kinemtic viscosity of salt water
        blk.kinem_visc_sw[i] = blk.dynamic_visc[i] / blk.density[i]

        # Prandtl number for water (-)
        blk.Pr[i] = (
            blk.specific_heat[i] * blk.dynamic_visc[i]
        ) / blk.thermal_conductivity[i]

        # Latent heat of vaporization of pure water J/kg
        blk.freshwater_vap_latent_heat[i] = (
            2501.67 - 2.389 * blk.saltwater_temp[i - 1]
        ) * 1000

        # Calculation of partial saturated vapor pressure of saltwater
        # According to parametric analysis and available literature,
        # the partial vapor pressure plays a major role in the evaporation of water.

        # A coeff obtained from the water molar fraction in salt solutions from 0-350 g/l Paper: Kokya and Kokya
        blk.water_activity[i] = (-0.000566 * blk.salinity[i - 1]) + 0.99853070

        # Partial saturated vapor pressure at a saltwater temperature (N/m^2)
        blk.sw_partial_vap_press[i] = blk.water_activity[i] * math.exp(
            25.317 - (5144 / (blk.saltwater_temp[i - 1] + 273))
        )

        # Partial saturated vapor pressure at glass cover temperature (N/m^2)
        blk.partial_vap_press_at_glass[i] = math.exp(
            25.317 - (5144 / (blk.glass_temp[i - 1] + 273))
        )

        # Coefficient of volume expansion (1/°C) correlation from Zhutovsky and Kovler (2015)
        beta = 1e-6 * (
            -0.000006 * blk.saltwater_temp[i - 1] ** 4
            + 0.001667 * blk.saltwater_temp[i - 1] ** 3
            - 0.197796 * blk.saltwater_temp[i - 1] ** 2
            + 16.862446 * blk.saltwater_temp[i - 1]
            - 64.319951
        )

        # Grashof number
        blk.Gr[i] = abs(
            (
                gravity
                * beta
                * (blk.basin_temp[i - 1] - blk.saltwater_temp[i - 1])
                * (blk.depth[i - 1] ** 3)
            )
            / (blk.kinem_visc_sw[i] ** 2)
        )

        # Heat transfer coeff of water layer
        blk.water_heat_trans_coeff[i] = abs(
            (blk.thermal_conductivity[i] / blk.depth[i - 1])
            * AA
            * (blk.Gr[i] * blk.Pr[i]) ** BB
        )

        # Convective heat transfer coeff (W/m^2.°C) (Dunkel)
        blk.conv_heat_trans_coeff_water_glass[i] = 0.884 * (
            (
                abs(
                    (
                        (blk.saltwater_temp[i - 1] - blk.glass_temp[i - 1])
                        + (
                            (
                                (
                                    blk.sw_partial_vap_press[i]
                                    - blk.partial_vap_press_at_glass[i]
                                )
                                * (blk.saltwater_temp[i - 1] + 273.15)
                            )
                            / (268900 - blk.sw_partial_vap_press[i])
                        )
                    )
                )
            )
            ** (1 / 3)
        )

        # Radiative heat transfer coeff from basin water to glass cover (W/m^2.°C)
        blk.rad_heat_transf_coeff_water_glass[i] = abs(
            emissivity_water
            * stefan_boltzmann
            * (
                (
                    ((blk.saltwater_temp[i - 1] + 273) ** 2)
                    + ((blk.glass_temp[i - 1] + 273) ** 2)
                )
                * (blk.saltwater_temp[i - 1] + blk.glass_temp[i - 1] + 546)
            )
        )

        # Evaporative heat transfer coeff from basin water to glass cover (W/m2.°C)
        blk.evap_heat_trans_coeff_water_glass[i] = abs(
            0.01628
            * blk.conv_heat_trans_coeff_water_glass[i]
            * (
                (blk.sw_partial_vap_press[i] - blk.partial_vap_press_at_glass[i])
                / (blk.temp_diff_inside_basin)
            )
        )

        # Total heat transfer coeff from basin water to glass cover (W/m2°C)
        blk.tot_heat_trans_coeff_water_glass[i] = (
            blk.conv_heat_trans_coeff_water_glass[i]
            + blk.rad_heat_transf_coeff_water_glass[i]
            + blk.evap_heat_trans_coeff_water_glass[i]
        )

        # Radiative heat transfer coeff from glass cover to ambient (W/m^2.°C)
        blk.rad_heat_trans_coeff_glass_ambient[i] = (
            stefan_boltzmann
            * emissivity_glass
            * (
                (
                    ((blk.glass_temp[i - 1] + 273) ** 4)
                    - ((blk.sky_temp[i - 1] + 273) ** 4)
                )
                / (blk.temp_diff_outside_basin)
            )
        )
        if blk.wind_velocity[i] > 5:
            # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            blk.conv_heat_trans_coeff_basin_ambient[i] = 2.8 + (
                3.0 * blk.wind_velocity[i]
            )
            # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
            blk.conv_heat_trans_coeff_glass_ambient[i] = 2.8 + (
                3.0 * blk.wind_velocity[i]
            )
        else:
            # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            blk.conv_heat_trans_coeff_basin_ambient[i] = 2.8 + (
                3.8 * blk.wind_velocity[i]
            )
            # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
            blk.conv_heat_trans_coeff_glass_ambient[i] = 2.8 + (
                3.8 * blk.wind_velocity[i]
            )

        # Total heat loss coeff from the glass cover to the outer atmosphere
        blk.tot_heat_trans_coeff_glass_ambient[i] = (
            blk.conv_heat_trans_coeff_glass_ambient[i]
            + blk.rad_heat_trans_coeff_glass_ambient[i]
        )

        # Heat loss coefficient from basin liner to the atmosphere
        blk.tot_heat_trans_coeff_basin_ambient[i] = 1 / (
            (thickness_insulation / conductivity_insulation)
            + (1 / (blk.conv_heat_trans_coeff_basin_ambient[i]))
        )

        # Effective overall absorptivity for energy balance equation
        blk.effective_absorp[i] = (
            (
                absorp_effective_basin
                * (
                    blk.water_heat_trans_coeff[i]
                    / (
                        blk.water_heat_trans_coeff[i]
                        + blk.tot_heat_trans_coeff_basin_ambient[i]
                        + blk.conv_heat_trans_coeff_basin_ambient[i]
                    )
                )
            )
            + (absorp_effective_water)
            + (
                (absorp_effective_glass)
                * (
                    blk.tot_heat_trans_coeff_water_glass[i]
                    / (
                        blk.tot_heat_trans_coeff_water_glass[i]
                        + blk.tot_heat_trans_coeff_glass_ambient[i]
                    )
                )
            )
        )

        # Calculation of overall heat transfer coefficients
        # Overall heat loss coefficient (W/m^2.°C)
        blk.overall_heat_loss_coeff_glass_surr[i] = (
            (conductivity_glass / thickness_glass)
            * (blk.tot_heat_trans_coeff_glass_ambient[i])
        ) / (
            (conductivity_glass / thickness_glass)
            + blk.tot_heat_trans_coeff_glass_ambient[i]
        )

        # Overall bottom heat transfer coefficient between the water mass and the surroundings (W/m^2.°C)
        blk.overall_bottom_heat_trans_coeff_water_mass_surr[i] = (
            blk.tot_heat_trans_coeff_water_glass[i]
            * blk.overall_heat_loss_coeff_glass_surr[i]
        ) / (
            blk.tot_heat_trans_coeff_water_glass[i]
            + blk.overall_heat_loss_coeff_glass_surr[i]
        )

        # Overall bottom heat transfer coefficient from bottom to ambient (W/m^2.°C)
        blk.overall_bottom_heat_loss_coeff_water_mass_surr[i] = (
            blk.water_heat_trans_coeff[i] * blk.tot_heat_trans_coeff_basin_ambient[i]
        ) / (blk.water_heat_trans_coeff[i] + blk.tot_heat_trans_coeff_basin_ambient[i])

        blk.overall_side_heat_loss_coefficient[i] = (
            blk.area_side_water[i] / area_bottom_basin
        ) * blk.overall_bottom_heat_loss_coeff_water_mass_surr[i]

        blk.overall_heat_trans_coeff_basin_surr[i] = (
            blk.overall_bottom_heat_loss_coeff_water_mass_surr[i]
            + blk.overall_side_heat_loss_coefficient[i]
        )
        blk.overall_external_heat_trans_loss_coeff[i] = (
            blk.overall_bottom_heat_trans_coeff_water_mass_surr[i]
            + blk.overall_heat_trans_coeff_basin_surr[i]
        )

        # Present [i] temperature calculation
        grouping_term = blk.overall_external_heat_trans_loss_coeff[i] / (
            blk.sw_mass[i - 1] * blk.specific_heat[i]
        )

        time_dependent_term = (
            (blk.effective_absorp[i] * blk.irradiance[i])
            + (blk.overall_external_heat_trans_loss_coeff[i] * blk.ambient_temp[i])
        ) / (blk.sw_mass[i - 1] * blk.specific_heat[i])

        blk.glass_temp[i] = (
            (absorp_effective_glass * blk.irradiance[i])
            + (blk.tot_heat_trans_coeff_water_glass[i] * blk.saltwater_temp[i - 1])
            + (blk.overall_heat_loss_coeff_glass_surr[i] * blk.ambient_temp[i])
        ) / (
            blk.tot_heat_trans_coeff_water_glass[i]
            + blk.overall_heat_loss_coeff_glass_surr[i]
        )

        blk.saltwater_temp[i] = (time_dependent_term / grouping_term) * (
            1 - np.exp(-grouping_term * blk.time[i])
        ) + (blk.saltwater_temp[i] * np.exp(-grouping_term * blk.time[i]))

        blk.basin_temp[i] = (
            (absorp_effective_basin * blk.irradiance[i])
            + (blk.water_heat_trans_coeff[i] * blk.saltwater_temp[i - 1])
            + (
                (
                    blk.tot_heat_trans_coeff_basin_ambient[i]
                    + blk.conv_heat_trans_coeff_basin_ambient[i]
                )
                * blk.basin_temp[i - 1]
            )
        ) / (
            blk.water_heat_trans_coeff[i]
            + blk.tot_heat_trans_coeff_basin_ambient[i]
            + blk.conv_heat_trans_coeff_basin_ambient[i]
        )

        # Evaporation estimation of freshwater and saltwater
        # Distilled (kg)
        blk.evap_fw_mass[i] = (
            area_bottom_basin
            * blk.evap_heat_trans_coeff_water_glass[i]
            * (blk.temp_diff_inside_basin)
        ) / blk.freshwater_vap_latent_heat[i]

        # Distilled saltwater conversion (Morton) (kg)
        blk.evap_sw_mass[i] = blk.evap_fw_mass[i] / (1 + blk.salinity[i - 1] / 1e3)

        # Freshwater (kg) this iteration
        blk.fw_mass[i] = blk.fw_mass[i - 1] - blk.evap_sw_mass[i]
        # Saltwater (kg) this iteration
        blk.sw_mass[i] = blk.fw_mass[i] + blk.salt_mass
        # Water blk.depth this iteration
        blk.depth[i] = blk.sw_mass[i] / (
            blk.density[i] * area_bottom_basin
        )  # Current blk.depth (m)

        if blk.salinity[i - 1] >= maximum_solubility:
            blk.salinity[i] = maximum_solubility

            # Excess blk.salinity (assuming no saturation possible) (g/l)
            blk.excess_salinity[i] = (blk.salt_mass * 1000) / blk.fw_mass[i]

            if blk.excess_salinity[i] < maximum_solubility:
                blk.salt_precipitated[i] = 0
                blk.volume_scale_formation[i] = 0
                blk.thickness_scale_formation[i] = 0
            else:
                # Mass precipitated (kg)
                blk.salt_precipitated[i] = (
                    blk.fw_mass[i]
                    * (blk.excess_salinity[i] - maximum_solubility)
                    / 1000
                )
                # Volume of scale formation (m3)
                blk.volume_scale_formation[i] = blk.salt_precipitated[i] / density_nacl
                # Thickness of scale formation (m)
                blk.thickness_scale_formation[i] = (
                    blk.volume_scale_formation[i] / area_bottom_basin
                )

        else:
            blk.salinity[i] = (blk.salt_mass * 1000) / blk.fw_mass[i]
            blk.excess_salinity[i] = blk.salinity[i]
        if blk.depth[i] <= 0 or blk.fw_mass[i] <= 0:
            # At this point either the blk.depth is negative
            # or the amount of freshwater available is negative
            # signaling we have reached the blk.time required for one ZLD cycle for one solar still

            # Number of seconds for a single ZLD cycle
            blk.num_seconds_for_zld_cycle = i
            blk.num_days_for_zld_cycle = i / 3600 / 24
            break

    num_zld_cycles_per_year = (
        len(blk.weather_data) * 3600
    ) / blk.num_seconds_for_zld_cycle  # (s / year) / s = year**-1

    # Total water productivity in year [kg water per m2 area per year]
    annual_water_yield = (
        blk.initial_mass_fw * num_zld_cycles_per_year
    ) / area_bottom_basin
    # Daily water yield [kg water per m2 area per day]
    daily_water_yield = annual_water_yield / days_in_year

    return daily_water_yield, num_zld_cycles_per_year
