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
from watertap_contrib.reflo.unit_models.util.water_yield_calculation import (
    calculate_density,
    calculate_viscosity,
    calculate_specific_heat,
    calculate_thermal_conductivity,
)

# reading the csv file
cvsDataframe = pd.read_csv(
    "/Users/ksitterl/Documents/SETO/models/solar_still_zld/SS_ZLD_Model 2/Data/TMY2 SAM CSV/data.csv",
    skiprows=2,
)

Input_data = cvsDataframe  # dataset

# Radiative properties
absorptivity_glass = 0.047  # Absorptivity of glass (-)
absorptivity_water = 0.20  # Absorptivity of water (-)
absorptivity_basin = 0.65  # Absorptivity of basin (-)

reflectivity_glass = 0.047  # Reflectivity of glass (-)
reflectivity_water = 0.08  # Reflectivity of Water (-)

emissivity_glass = 0.94  # Glass Emissivity (-)
emissivity_water = 0.95 * 1.0  # Water Emissivity (-)

emissivity_effective = 1 / (
    (1 / emissivity_water) + (1 / emissivity_glass) - 1
)  # Effective emissivity of water to glass (-)

# No Attenuation factor considered
absorptivity_effective_water = (
    (absorptivity_water)
    * (1 - absorptivity_glass)
    * (1 - reflectivity_glass)
    * (1 - reflectivity_water)
)  # Fraction of solar radiation absorbed by water (-)
absorptivity_effective_basin = (
    (absorptivity_basin)
    * (1 - absorptivity_glass)
    * (1 - reflectivity_glass)
    * (1 - absorptivity_water)
    * (1 - reflectivity_water)
)  # Fraction of solar radiation absorbed by basin liner (-)
absorptivity_effective_glass = (
    1 - reflectivity_glass
) * absorptivity_glass  # Fraction of solar radiation absorbed by a glass cover (-)
stefan_boltzmann = 5.6697 * (10**-8)  # Stephen-Boltzman Constant (W / (m^2.K^4))
thickness_insulation = 0.005  # Thickness of SS Insulation (m)
conductivity_insulation = 0.033  # Conductivity of SS Insulation (W/m.K)
thickness_glass = 0.004  # Thickness of glass m
conductivity_glass = 1.03  # Conductivity of glass W/m.K
gravity = 9.81  # Accelaration due to gravity (m/s^2)

No_months = 1
Hours = 24  # Determine how many hours are being considered
seconds = 3600 * Hours  # Converting hours to seconds
density_NaCl = 2165  # density of sodium chloride (kg/m^3)
K_NaCl = 3.5  # thermal conductivity of sodium chloride (W/m.K)
maximum_solubility = 365  # maximum solubility of salt in water (g/l)

# Adaptive coefficient heat transfer coefficient with buoyancy
AA = 0.54
# Power of nondimensional numbers for heat transfer coefficient with buoyancy
BB = 0.25
# User input from pyomo.ipynb
# initial_salinity = float(Input_data.iloc[-1, 0])  # Water salinity (g/l)
# initial_depth = float(Input_data.iloc[-1, 2])  # water depth in SS basin (m)
# side_basin = float(Input_data.iloc[-1, 1])  # Side of SS basin (m)


def create_input_arrays():

    def generate_continuous_day_series():
        # Days in each month for a non-leap year
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        continuous_day_series = []

        for day in range(1, 32):  # up to 31 days to cover all months
            for month_index, days in enumerate(days_in_month):
                if day <= days:  # check if the day exists in this month
                    # Calculate the day of the year for the nth day of each month
                    day_of_year = sum(days_in_month[:month_index]) + day
                    continuous_day_series.append(day_of_year)

        return continuous_day_series

    continuous_day_series = generate_continuous_day_series()

    ambient_temp_by_hr = []
    irradiance_by_hr = []
    wind_vel_by_hr = []

    # Collecting hourly weather data
    for day in continuous_day_series:
        if day == 365:
            day = 364
        start_row = day * 24
        # end_row = start_row + 24
        for kk in range(24):
            temp = float(Input_data.iloc[start_row + kk, 7])
            irr = float(Input_data.iloc[start_row + kk, 4]) + 10
            wind = float(Input_data.iloc[start_row + kk, 11])

            if irr > 20:
                ambient_temp_by_hr.append(temp)
                irradiance_by_hr.append(irr)
                wind_vel_by_hr.append(wind)

    ambient_temp_by_hr = pd.Series(ambient_temp_by_hr)
    irradiance_by_hr = pd.Series(irradiance_by_hr)
    wind_vel_by_hr = pd.Series(wind_vel_by_hr)

    return ambient_temp_by_hr, irradiance_by_hr, wind_vel_by_hr


def get_solar_still_daily_water_yield_zld(
    initial_salinity=200, initial_depth=0.01, side_basin=0.6
):
    area_bottom_basin = side_basin**2  # Area of square basin (m^2)

    ambient_temp_by_hr, irradiance_by_hr, wind_vel_by_hr = create_input_arrays()

    # Allocating the variables
    # Second variables
    depth = np.zeros(len(irradiance_by_hr) * 3600)
    salinity = np.zeros(len(irradiance_by_hr) * 3600)
    excess_salinity = np.zeros(len(irradiance_by_hr) * 3600)
    volume_scale_formation = np.zeros(len(irradiance_by_hr) * 3600)
    thickness_scale_formation = np.zeros(len(irradiance_by_hr) * 3600)
    Evap_mass_sw = np.zeros(len(irradiance_by_hr) * 3600)
    Evap_mass_sw_sec = np.zeros(len(irradiance_by_hr) * 3600)
    irradiance = np.zeros(len(irradiance_by_hr) * 3600)
    salt_precipitated = np.zeros(len(irradiance_by_hr) * 3600)
    sw_mass = np.zeros(len(irradiance_by_hr) * 3600)
    fw_mass = np.zeros(len(irradiance_by_hr) * 3600)
    time = np.zeros(len(irradiance_by_hr) * 3600)
    wind_velocity = np.zeros(len(irradiance_by_hr) * 3600)
    ambient_temp = np.zeros(len(irradiance_by_hr) * 3600)
    basin_temp = np.zeros(len(irradiance_by_hr) * 3600)
    glass_temp = np.zeros(len(irradiance_by_hr) * 3600)
    sky_temp = np.zeros(len(irradiance_by_hr) * 3600)
    saltwater_temp = np.zeros(len(irradiance_by_hr) * 3600)
    # Total_water_mass = np.zeros(len(irradiance_by_hr) * 3600)

    # Converting hourly data into per second
    for hour in range(len(irradiance_by_hr)):
        start_idx = 3600 * hour
        end_idx = start_idx + 3600
        irradiance[start_idx:end_idx] = irradiance_by_hr[hour]
        wind_velocity[start_idx:end_idx] = wind_vel_by_hr.loc[hour]
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
    depth[1] = initial_depth
    depth[0] = initial_depth
    sw_mass[1] = depth[1] * initial_density * area_bottom_basin
    fw_mass[1] = sw_mass[1] / (1 + salinity[1] / 1000)
    salt_mass = (salinity[1] * fw_mass[1]) / 1000  # Mass of Sodium Chloride (kg)
    excess_salinity[1] = salinity[1]  # salinity without maximum solublity (g/l)

    # Initial effective radiation temperature of the sky (deg AA)
    sky_temp[1] = 0.0552 * ((ambient_temp[1]) ** 1.5)

    time[1] = 1

    for i in range(2, len(irradiance), 1):

        if depth[i - 1] <= 0 or fw_mass[i - 1] <= 0:
            depth[i - 1] = initial_depth
            salinity[i - 1] = initial_salinity
            sw_mass[i - 1] = depth[i - 1] * initial_density * area_bottom_basin
            fw_mass[i - 1] = sw_mass[i - 1] / (1 + salinity[i] / 1000)

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

        # Perimeter x depth of water (m^2)
        area_side_water = (2 * (2 * side_basin)) * depth[i - 1]

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

        # A coefficient obtained from the water molar fraction in salt solutions from 0-350 g/l Paper: Kokya and Kokya
        water_activity = (-0.000566 * salinity[i - 1]) + 0.99853070
        # Partial saturated vapor pressure at a saltwater temperature (BB/m^2)
        sw_partial_vap_press = water_activity * math.exp(
            25.317 - (5144 / (saltwater_temp[i - 1] + 273))
        )
        # Partial saturated vapor pressure at glass cover temperature (BB/m^2)
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
        )  # Gr number
        # Heat transfer coefficient of water layer
        water_heat_transfer_coefficient = abs(
            (thermal_conductivity / depth[i - 1]) * AA * (Gr * Pr) ** BB
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
        evap_heat_trans_coeff_water_glass = abs(
            0.01628
            * conv_heat_trans_coeff_water_glass
            * (
                (sw_partial_vap_press - partial_vap_press_at_glass)
                / (temp_diff_inside_basin)
            )
        )  # Evaporative heat transfer coefficient from basin water to glass cover (W/m2.°C)
        hT_wgi = (
            conv_heat_trans_coeff_water_glass
            + rad_heat_transf_coeff_water_glass
            + evap_heat_trans_coeff_water_glass
        )  # Total heat transfer coefficient from basin water to glass cover (W/m2°C)
        hRgoa = (
            stefan_boltzmann
            * emissivity_glass
            * (
                (((glass_temp[i - 1] + 273) ** 4) - ((sky_temp[i - 1] + 273) ** 4))
                / (temp_diff_outside_basin)
            )
        )  # Radiative heat transfer coefficient from glass cover to ambient (W/m^2.°C)
        if wind_velocity[i] > 5:
            hcba = 2.8 + (
                3.0 * wind_velocity[i]
            )  # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            hCgoa = 2.8 + (
                3.0 * wind_velocity[i]
            )  # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
        else:
            hcba = 2.8 + (
                3.8 * wind_velocity[i]
            )  # Convective heat transfer coefficient from basin to ambient (W/m2°C)
            hCgoa = 2.8 + (
                3.8 * wind_velocity[i]
            )  # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)

        hTgoa = (
            hCgoa + hRgoa
        )  # total heat loss coefficient from the glass cover to the outer atmosphere
        htba = 1 / (
            (thickness_insulation / conductivity_insulation) + (1 / (hcba))
        )  # heat loss coefficient from basin liner to the atmosphere

        alpha_eff = (
            (
                absorptivity_effective_basin
                * (
                    water_heat_transfer_coefficient
                    / (water_heat_transfer_coefficient + htba + hcba)
                )
            )
            + (absorptivity_effective_water)
            + ((absorptivity_effective_glass) * (hT_wgi / (hT_wgi + hTgoa)))
        )  # effective overall absorptivity for energy balance equation

        # Calculation of overall heat transfer coefficients
        UTgi_a = ((conductivity_glass / thickness_glass) * (hTgoa)) / (
            (conductivity_glass / thickness_glass) + hTgoa
        )  # overall heat loss coefficient (W/m^2.°C)
        UT = (hT_wgi * UTgi_a) / (
            hT_wgi + UTgi_a
        )  # overall bottom heat transfer coefficient between the water mass and the surroundings (W/m^2.°C)
        Ub = (water_heat_transfer_coefficient * htba) / (
            water_heat_transfer_coefficient + htba
        )  # Overall bottom heat transfer coefficient from bottom to ambient (W/m^2.°C)
        Uss = (area_side_water / area_bottom_basin) * Ub
        Ubs = Ub + Uss
        ULs = UT + Ubs

        # Present [i] temperature calculation
        a = ULs / (sw_mass[i - 1] * specific_heat)
        f = ((alpha_eff * irradiance[i]) + (ULs * ambient_temp[i])) / (
            sw_mass[i - 1] * specific_heat
        )
        glass_temp[i] = (
            (absorptivity_effective_glass * irradiance[i])
            + (hT_wgi * saltwater_temp[i - 1])
            + (UTgi_a * ambient_temp[i])
        ) / (hT_wgi + UTgi_a)
        saltwater_temp[i] = (f / a) * (1 - np.exp(-a * time[i])) + (
            saltwater_temp[i] * np.exp(-a * time[i])
        )
        basin_temp[i] = (
            (absorptivity_effective_basin * irradiance[i])
            + (water_heat_transfer_coefficient * saltwater_temp[i - 1])
            + ((htba + hcba) * basin_temp[i - 1])
        ) / (water_heat_transfer_coefficient + htba + hcba)

        # Evaporation estimation of freshwater and saltwater
        Evap_mass_fw = (
            area_bottom_basin
            * evap_heat_trans_coeff_water_glass
            * (saltwater_temp[i] - glass_temp[i])
            * 3600
        ) / freshwater_vap_latent_heat  # Distillated (kg)

        Evap_mass_fw_sec = (
            area_bottom_basin
            * evap_heat_trans_coeff_water_glass
            * (temp_diff_inside_basin)
        ) / freshwater_vap_latent_heat  # Distillated (kg)

        Evap_mass_sw[i] = Evap_mass_fw / (
            1 + salinity[i - 1] / 1e3
        )  # Distillated saltwater conversion (Morton) (kg)
        Evap_mass_sw_sec[i] = Evap_mass_fw_sec / (
            1 + salinity[i - 1] / 1e3
        )  # Distillated saltwater conversion (Morton) (kg)

        fw_mass[i] = fw_mass[i - 1] - Evap_mass_sw_sec[i]  # Current Freshwater (kg)
        sw_mass[i] = fw_mass[i] + salt_mass  # Current saltwater (kg)
        depth[i] = sw_mass[i] / (density * area_bottom_basin)  # Current depth (m)

        if salinity[i - 1] >= maximum_solubility:
            salinity[i] = maximum_solubility

            excess_salinity[i] = (salt_mass * 1000) / fw_mass[
                i
            ]  # Excess salinity (assuming no saturation possible) (g/l)

            if excess_salinity[i] < maximum_solubility:
                salt_precipitated[i] = 0
                volume_scale_formation[i] = 0
                thickness_scale_formation[i] = 0
            else:
                salt_precipitated[i] = (
                    fw_mass[i] * (excess_salinity[i] - maximum_solubility) / 1000
                )  # Mass precipitated(kg)
                volume_scale_formation[i] = (
                    salt_precipitated[i] / density_NaCl
                )  # volume of scale formation(m3)
                thickness_scale_formation[i] = (
                    volume_scale_formation[i] / area_bottom_basin
                )  # thickness of scale formation(m)

        else:
            salinity[i] = (salt_mass * 1000) / fw_mass[i]
            excess_salinity[i] = salinity[i]

        if depth[i] <= 0 or fw_mass[i] <= 0:
            break

    ZLD_counter = len(irradiance) / i

    TotalYield = (
        fw_mass[1] * ZLD_counter / 1000
    )  # Yield in m3, assuming water density of 1000 kg/m3
    Volume_salt = ZLD_counter * max(
        volume_scale_formation
    )  # total volume of salt precipitated in the year (m3)

    print(TotalYield, ZLD_counter)


get_solar_still_daily_water_yield_zld()
