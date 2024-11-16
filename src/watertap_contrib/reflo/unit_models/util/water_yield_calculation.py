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
import pandas as pd
import numpy as np

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

# Radiative properties
absorp_glass = 0.047  # Absorptivity of glass (-)
absorp_water = 0.20  # Absorptivity of water (-)
absorp_basin = 0.65  # Absorptivity of basin (-)

reflectivity_glass = 0.047  # Reflectivity of glass (-)
reflectivity_water = 0.08  # Reflectivity of Water (-)

emissivity_glass = 0.94  # Glass Emissivity (-)
emissivity_water = 0.95 * 1.0  # Water Emissivity (-)

emissivity_effective = 1 / (
    (1 / emissivity_water) + (1 / emissivity_glass) - 1
)  # Effective emissivity of water to glass (-)

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

# Adaptive coefficient heat transfer coefficient with buoyancy
AA = 0.54
# Power of nondimensional numbers for heat transfer coefficient with buoyancy
BB = 0.25


def get_solar_still_daily_water_yield(
    input_weather_file_path=None,  # path to input weather file
    interval_day=15,  # interval at which to calculate water yield (e.g., every 15 days)
    salinity=20,  # salinity of influent water; g/L
    water_depth_basin=0.02,  # depth of water in solar still basin; m
    length_basin=0.6,  # length of each side of basin (length=width); m
):
    """
    Uses input weather data to calculate daily water yield for solar still (m3/m2/day)
    """

    weather_data = pd.read_csv(input_weather_file_path, skiprows=2)
    daterow = np.zeros(days_in_year)

    for tt in range(0, days_in_year, interval_day + 1):
        daterow[tt] = tt * 24
        operational_matrix = [
            daterow[i] for i in range(len(daterow)) if daterow[i] != 0
        ]
        cumulative_yearly_water_yield = np.zeros(len(operational_matrix))

    for tt in range(len(operational_matrix)):
        u = int(operational_matrix[tt])

        # Allocating the variables
        # Second variables
        evap_sw_mass = np.zeros(seconds_per_day)
        irradiance = np.zeros(seconds_per_day)
        time = np.zeros(seconds_per_day)
        wind_velocity = np.zeros(seconds_per_day)
        ambient_temp = np.zeros(seconds_per_day)
        basin_temp = np.zeros(seconds_per_day)
        glass_temp = np.zeros(seconds_per_day)
        sky_temp = np.zeros(seconds_per_day)
        saltwater_temp = np.zeros(seconds_per_day)

        # Hour variables
        total_evaporated_saltwater_in_year = np.zeros(hours_per_day)
        progressive_evaporation_water_day = np.zeros(hours_per_day)
        hourly_irradiance = np.zeros(hours_per_day)
        hourly_wind_velocity = np.zeros(hours_per_day)
        hourly_ambient_temp = np.zeros(hours_per_day)

        # Collecting hourly weather data
        for kk in range(0, hours_per_day, 1):
            # Solar irradiation at given hour (W/m2) + 10 to avoid divisions by zero in eqns.
            hourly_irradiance[kk] = float(weather_data.iloc[u + kk, 4]) + 10
            # Ambient temperature at given hour (°C)
            hourly_ambient_temp[kk] = float(weather_data.iloc[u + kk, 7])
            # Windspeed at given hour (m/s)
            hourly_wind_velocity[kk] = float(weather_data.iloc[u + kk, 11])

        # Converting hourly data into per second
        for hour in range(hours_per_day):
            start_index = 3600 * hour
            end_index = start_index + 3600
            irradiance[start_index:end_index] = hourly_irradiance[hour]
            wind_velocity[start_index:end_index] = hourly_wind_velocity[hour]
            ambient_temp[start_index:end_index] = hourly_ambient_temp[hour]

        # Initializing Temperatures
        # Initial system is assumed to be in thermal equilibrium with ambient

        # Initial water temperature (°C)
        saltwater_temp[1] = float(weather_data.iloc[u, 7])
        # Initial basin temperature (°C)
        basin_temp[1] = float(weather_data.iloc[u, 7])
        # Initial glass temperature (°C)
        glass_temp[1] = float(weather_data.iloc[u, 7])

        # Geometrical properties of squared basin
        # Area of square basin (m^2)
        area_bottom_basin = length_basin**2
        # Perimeter x water_depth_basin of water (m^2)
        area_side_water = (2 * (2 * length_basin)) * water_depth_basin

        # Initial effective radiation temperature of the sky (deg C)
        sky_temp[1] = 0.0552 * ((ambient_temp[1]) ** 1.5)

        time[1] = 1

        for i in range(2, seconds_per_day, 1):
            time[i] = time[i - 1] + 1

            # Avoiding singularities
            temp_diff_inside_basin = saltwater_temp[i - 1] - glass_temp[i - 1]
            if temp_diff_inside_basin <= 0:
                temp_diff_inside_basin = 0.01

            # Avoiding singularities
            temp_diff_outside_basin = glass_temp[i - 1] - ambient_temp[i - 1]
            if temp_diff_outside_basin <= 0:
                temp_diff_outside_basin = 0.01
            # Effective radiation temperature of the sky
            sky_temp[i] = 0.0552 * ((ambient_temp[i]) ** 1.5)

            # Calculating thermophysical properties of seawater
            density = calculate_density(salinity, saltwater_temp[i - 1])

            dynamic_visc = calculate_viscosity(salinity, saltwater_temp[i - 1])

            specific_heat = calculate_specific_heat(salinity, saltwater_temp[i - 1])

            thermal_conductivity = calculate_thermal_conductivity(
                salinity, saltwater_temp[i - 1]
            )

            # Water Mass (kg)
            sw_mass = density * (area_bottom_basin * water_depth_basin)
            kinematic_visc = dynamic_visc / density
            # Prandtl number for water (-)
            Pr = (specific_heat * dynamic_visc) / thermal_conductivity
            # Latent heat of vaporization of pure water J/kg
            freshwater_vap_latent_heat = (
                2501.67 - 2.389 * saltwater_temp[i - 1]
            ) * 1000

            # Calculation of partial saturated vapor pressure of saltwater
            # According to parametric analysis and available literature,
            # the partial vapor pressure plays a major role in the evaporation of water.

            # A coefficient obtained from the water molar fraction in salt solutions from 0-350 g/l Paper: Kokya and Kokya
            water_activity = (-0.000566 * salinity) + 0.99853070
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
                    * (water_depth_basin**3)
                )
                / (kinematic_visc**2)
            )
            # Heat transfer coefficient of water layer
            water_heat_trans_coeff = abs(
                (thermal_conductivity / water_depth_basin) * AA * (Gr * Pr) ** BB
            )
            # Convective heat transfer coefficient (W/m^2.°C) (Dunkel)
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
            # Radiative heat transfer coefficient from basin water to glass cover (W/m^2.°C)
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
            # Evaporative heat transfer coefficient from basin water to glass cover (W/m2.°C)
            evaporative_heat_transfer_coefficient_water_glass = abs(
                0.01628
                * conv_heat_trans_coeff_water_glass
                * (
                    (sw_partial_vap_press - partial_vap_press_at_glass)
                    / (temp_diff_inside_basin)
                )
            )
            # Total heat transfer coefficient from basin water to glass cover (W/m2°C)
            tot_heat_trans_coeff_water_glass = (
                conv_heat_trans_coeff_water_glass
                + rad_heat_transf_coeff_water_glass
                + evaporative_heat_transfer_coefficient_water_glass
            )
            # Radiative heat transfer coefficient from glass cover to ambient (W/m^2.°C)
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
            # Total heat loss coefficient from the glass cover to the outer atmosphere
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
                (conductivity_glass / thickness_glass)
                + tot_heat_trans_coeff_glass_ambient
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
                sw_mass * specific_heat
            )
            time_dependent_term = (
                (effective_absorp * irradiance[i])
                + (overall_external_heat_trans_loss_coeff * ambient_temp[i])
            ) / (sw_mass * specific_heat)
            glass_temp[i] = (
                (absorp_effective_glass * irradiance[i])
                + (tot_heat_trans_coeff_water_glass * saltwater_temp[i - 1])
                + (overall_heat_loss_coeff_glass_surr * ambient_temp[i])
            ) / (tot_heat_trans_coeff_water_glass + overall_heat_loss_coeff_glass_surr)
            saltwater_temp[i] = (time_dependent_term / grouping_term) * (
                1 - np.exp(-grouping_term * time[i])
            ) + (saltwater_temp[1] * np.exp(-grouping_term * time[i]))
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
            # Distillated (kg)
            evap_fw_mass = (
                area_bottom_basin
                * evaporative_heat_transfer_coefficient_water_glass
                * (saltwater_temp[i] - glass_temp[i])
                * 3600
            ) / freshwater_vap_latent_heat
            # Distillated saltwater conversion (Morton) (kg)
            evap_sw_mass[i] = evap_fw_mass / (1 + salinity / 1e3)

        # Reshaping Mean hourly mass of evaporated water array for productivity and cumulative calculation
        reshaping_mass_evaporated_water = np.reshape(
            evap_sw_mass, (hours_per_day, 3600)
        )
        reshaping_mass_evaporated_water[reshaping_mass_evaporated_water < 0] = 0
        reshaping_mass_evaporated_water[reshaping_mass_evaporated_water > 1] = 0
        df_reshaping_mass_evaporated_water = pd.DataFrame(
            reshaping_mass_evaporated_water
        )
        total_evaporated_saltwater_in_year = df_reshaping_mass_evaporated_water.mean(
            axis=1
        )

        # Water distillated per hour kg/m^2
        hourly_productivity = total_evaporated_saltwater_in_year / area_bottom_basin
        # Total Water Productivity kg/m^2
        cumulative_daily_water_yield = np.nansum(hourly_productivity)
        # Total Water Productivity kg/m^2 in year
        cumulative_yearly_water_yield[tt] = cumulative_daily_water_yield

    # Total Water Productivity in year [m3 water per m2 area per year]
    annual_water_yield = (
        (days_in_year / (len(operational_matrix)))
        * np.nansum(cumulative_yearly_water_yield)
    ) / 1000

    print(annual_water_yield)

    # Daily water yield on volumetric basis [m3 water per m2 area per year]
    daily_water_yield_vol = annual_water_yield / days_in_year
    # Daily water yield on mass basis
    daily_water_yield_mass = daily_water_yield_vol * 1000

    return daily_water_yield_mass


if __name__ == "__main__":

    f = "/Users/ksitterl/Documents/SETO/models/solar_still/SS_model-Sept2024/SS_Model/Data/TMY2 SAM CSV/data.csv"
    water_yield = get_solar_still_daily_water_yield(
        input_weather_file_path=f, interval_day=100
    )
    print(f"water_yield = {water_yield}")
