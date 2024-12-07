# Libraries
import math
import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# https://nsrdb.nrel.gov/data-viewer
f = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/development/evaporation_pond/evaporation_pond_test_data.csv"
df = pd.read_csv(f)


emissivity_water = 0.97
stefan_boltzmann = 4.903e-9  # Stefan-Boltzmann constant (MJ M-2 K-4)
long_wavelength_albedo = 0.03
short_wavelength_albedo = 0.05
# evap_rate
bowen_constant = 0.665
# Specific heat capacity MJ/kg degC
specific_heat_capacity = 1.013e-3
# salinity = int(df.iloc[-1, 0]) #salinity g/l
# evap_rate_method = int(df.iloc[-1, 1]) #salinity g/l
salinity = 100
evap_rate_method = 1
print(salinity, evap_rate_method)
# Create a new column to group by every 24 rows
df["group"] = np.arange(len(df)) // 24

# Group by the 'group' column and compute the mean for each group
daily_avg = df.groupby("group").mean().reset_index(drop=True)
daily_max = df.groupby("group").max().reset_index(drop=True)
daily_min = df.groupby("group").min().reset_index(drop=True)
daily_sum = df.groupby("group").sum().reset_index(drop=True)

# Drop the 'group' column as it's no longer needed
daily_avg = daily_avg.reset_index(drop=True)
daily_max = daily_max.reset_index(drop=True)
daily_min = daily_min.reset_index(drop=True)
daily_sum = daily_sum.reset_index(drop=True)

# Inputting characteristic data
pi = numpy.pi
days_in_year = 365  # Determine how many days_in_year are being considered

# Allocating the variables
air_temp_min = np.zeros(days_in_year)
air_temp_max = np.zeros(days_in_year)
air_temp_mean = np.zeros(days_in_year)
rad_shortwave = np.zeros(days_in_year)
rel_humidity_max = np.zeros(days_in_year)
rel_humidity_min = np.zeros(days_in_year)
pressure_atm = np.zeros(days_in_year)
rel_humidity_mean = np.zeros(days_in_year)
evap_rate_breb = np.zeros(days_in_year)
wind_velocity = np.zeros(days_in_year)
water_temp_mean = np.zeros(days_in_year)
rad_net = np.zeros(days_in_year)
# short_wavelength_albedo = np.zeros(days_in_year)
water_depth = np.zeros(days_in_year)
evap_rate_penman = np.zeros(days_in_year)
evap_rate = np.zeros(days_in_year)
wet_bulb_temp = np.zeros(days_in_year)
evap_rate_harbeck_waiv = np.zeros(days_in_year)

for kk in range(0, days_in_year, 1):
    # Daily
    rad_shortwave[kk] = float(daily_sum.iloc[kk, 8]) * 0.0036  # GHI
    air_temp_max[kk] = float(daily_max.iloc[kk, 5])  # Temperature
    air_temp_min[kk] = float(daily_min.iloc[kk, 5])  # Temperature
    air_temp_mean[kk] = float(daily_avg.iloc[kk, 5])  # Temperature
    rel_humidity_max[kk] = float(daily_max.iloc[kk, 10]) / 100  # Relative Humidity
    rel_humidity_min[kk] = float(daily_min.iloc[kk, 10]) / 100  # Relative Humidity
    rel_humidity_mean[kk] = float(daily_avg.iloc[kk, 10]) / 100  # Relative Humidity
    wind_velocity[kk] = float(daily_avg.iloc[kk, 17])  # Wind Speed
    pressure_atm[kk] = float(daily_avg.iloc[kk, 16]) / 10  # Pressure

if water_depth[0] <= 2:
    water_temp_mean[0] = 1.167 * air_temp_mean[0] - 0.175
else:
    water_temp_mean[0] = 0.955 * air_temp_mean[0] + 2.367

for i in range(1, days_in_year, 1):
    if water_depth[i - 1] <= 2:
        water_temp_mean[i] = 1.04 * air_temp_mean[i] + 0.22
    else:
        water_temp_mean[i] = 1.04 * air_temp_mean[i] + 0.22

    if salinity > 380:
        salinity = 380

    # Water latent heat of vaporization
    latent_heat_vap_water = (2501.67 - 2.389 * water_temp_mean[i - 1]) * 1000

    # Saturation vapor pressure at the min air temperature [kPa]
    pressure_sat_air_temp_min = 0.6108 * numpy.exp(
        (17.269 * air_temp_min[i]) / (air_temp_min[i] + 237.3)
    )
    # Saturation vapor pressure at the max air temperature [kPa]
    pressure_sat_air_temp_max = 0.6108 * numpy.exp(
        (17.269 * air_temp_max[i]) / (air_temp_max[i] + 237.3)
    )
    # Saturation vapor pressure at the mean air temperature [kPa]
    pressure_sat_water_temp = 0.6108 * numpy.exp(
        (17.269 * water_temp_mean[i]) / (water_temp_mean[i] + 237.3)
    )
    pressure_sat_actual = (
        (pressure_sat_air_temp_max * rel_humidity_min[i])
        + (pressure_sat_air_temp_min * rel_humidity_max[i])
    ) * 0.5

    if air_temp_mean[i] <= 0:
        emissivity_air = 0.99
    else:
        emissivity_air = 1.24 * ((pressure_sat_actual / air_temp_mean[i]) ** (1 / 7))
        if emissivity_air > 1:
            emissivity_air = 0.99
    # print(air_temp_mean[i], emissivity_air)

    rad_net_shortwave = (1 - short_wavelength_albedo) * rad_shortwave[i]
    rad_longwave_in = (
        emissivity_air * stefan_boltzmann * (air_temp_mean[i] + 273.15) ** 4
    )

    rad_net_longwave_in = (1 - long_wavelength_albedo) * rad_longwave_in
    rad_net_longwave_out = (
        emissivity_water * stefan_boltzmann * (water_temp_mean[i] + 273.15) ** 4
    )
    rad_net[i] = rad_net_shortwave + rad_net_longwave_in - rad_net_longwave_out

    # print(rad_shortwave[i], rad_net_shortwave, rad_longwave_in, rad_net_longwave_out, rad_net[i])
    psychometric_constant = (specific_heat_capacity * pressure_atm[i]) / (
        0.622 * (latent_heat_vap_water / 1e6)
    )

    # delta is the slope of the saturated vapor pressure-temperature curve at mean water temperature (kPa/degC)
    delta = (
        4098
        * (
            0.6108
            * numpy.exp((17.27 * water_temp_mean[i]) / (water_temp_mean[i] + 237.3))
        )
    ) / ((water_temp_mean[i] + 237.3) ** 2)

    wind_speed_func = 0.54 * wind_velocity[i]

    wet_bulb_temp[i] = (
        air_temp_mean[i]
        * np.arctan(0.151977 * (rel_humidity_mean[i] * 100 + 8.313659) ** (1 / 2))
        + 0.00391838
        * ((rel_humidity_mean[i] * 100) ** 3) ** (1 / 2)
        * numpy.arctan(0.023101 * rel_humidity_mean[i] * 100)
        - numpy.arctan(rel_humidity_mean[i] * 100 - 1.676331)
        + numpy.arctan(air_temp_mean[i] + rel_humidity_mean[i] * 100)
        - 4.686035
    )

    pressure_sat_wet_bulb_temp = 0.6108 * numpy.exp(
        (17.269 * wet_bulb_temp[i]) / (wet_bulb_temp[i] + 237.3)
    )
    wetted_surface_area = 5700

    area_coeff = 3.719e-9 * wetted_surface_area ** (-0.0459)
    water_activity = -0.00056678 * salinity + 0.9985307
    salinity_conversion = 1 / (1 + salinity / 1000)
    salinity_conversion = 1
    # Bowen ratio
    bowen_ratio = psychometric_constant * (
        (water_temp_mean[i] - air_temp_mean[i])
        / (water_activity * pressure_sat_water_temp - pressure_sat_actual)
    )
    print(i, bowen_ratio)
    # Daily evaporation Harbeck/WAIV model (mm/day)
    evap_rate_harbeck_waiv[i] = (
        area_coeff
        * (wind_velocity[i] * 8.64e7)
        * (
            (water_activity * pressure_sat_wet_bulb_temp)
            - (pressure_sat_actual * rel_humidity_mean[i])
        )
        * 10
    )

    if evap_rate_method == 0:

        evap_rate_breb[i] = 0
        evap_rate_penman[i] = 0
        evap_rate[i] = 0
        area_pond = 0

    elif evap_rate_method == 1:
        # Daily evaporation BREB model (mm/day)

        evap_rate_breb[i] = salinity_conversion * (
            (rad_net[i]) / ((latent_heat_vap_water / 1e6) * (1 + bowen_ratio))
        )
        # print(evap_rate_breb[i])
        evap_rate_penman[i] = 0
        evap_rate[i] = evap_rate_breb[i]
        area_pond = 10 * 4046.86

    elif evap_rate_method == 2:
        # Daily evaporation Penman model (mm/day)

        evap_rate_breb[i] = 0
        evap_rate_penman[i] = salinity_conversion * (
            (
                (delta / (delta + psychometric_constant))
                * (rad_net[i] / (latent_heat_vap_water / 1e6))
            )
            + (
                (psychometric_constant / (delta + psychometric_constant))
                * (
                    6.43
                    * wind_speed_func
                    * (water_activity * pressure_sat_water_temp - pressure_sat_actual)
                )
                / (latent_heat_vap_water / 1e6)
            )
        )
        evap_rate[i] = evap_rate_penman[i]
        area_pond = 10 * 4046.86

    elif evap_rate_method == 3:
        # Average daily evaporation rate of BREB and Penman (mm/day)

        # Daily evaporation BREB model (mm/day)
        evap_rate_breb[i] = salinity_conversion * (
            (rad_net[i]) / ((latent_heat_vap_water / 1e6) * (1 + bowen_ratio))
        )
        # Daily evaporation Penman model (mm/day)
        evap_rate_penman[i] = salinity_conversion * (
            (
                (delta / (delta + psychometric_constant))
                * (rad_net[i] / (latent_heat_vap_water / 1e6))
            )
            + (
                (psychometric_constant / (delta + psychometric_constant))
                * (
                    6.43
                    * wind_speed_func
                    * (water_activity * pressure_sat_water_temp - pressure_sat_actual)
                )
                / (latent_heat_vap_water / 1e6)
            )
        )
        # Average daily evaporation rate of two models (mm/day)
        evap_rate[i] = (evap_rate_breb[i] + evap_rate_penman[i]) / 2
        area_pond = 10 * 4046.86  # acre to m2


Tot_Evaporation = sum(evap_rate)  # Annual evaporation BREB model (mm)
Tot_E_Harbeck_WAIV = sum(
    evap_rate_harbeck_waiv
)  # Annual evaporation Harbeck/WAIV model (mm)

EP_evaporation_daily = area_pond * (
    Tot_Evaporation / 365000
)  # Daily evaporation BREB model(m^3/day)
Evap_day_WAIV = wetted_surface_area * (
    Tot_E_Harbeck_WAIV / 365000
)  # Daily evaporation Harbeck/WAIV model (m^3/day)


rad_net_mean = np.mean(rad_shortwave)


#   Literature
""" 
[1] Evaluating best evaporation estimate modelta for water surface evaporation in semi-arid region, India
[2] Estimation of lake evaporation from standard meteorological measurements: application to four Australian lakes in different climatic regions
[3] Crop evapotranspiration guideltaines for computing crop requirements. FAO Irrig. Drain. Report modeltaing and application
[4] Simulation of Lake evap_rate With Application to Modeltaing Lake Level Variations of Harney-Malheur Lake, Oregon
[5] A holistic water depth simulation modelta for small ponds
[6] Modeltaling historial lake levels and recent climate change at three closed lakes, Western Victoria, Australia (c. 1840-1990)
"""


def E_EPM():
    return EP_evaporation_daily, Evap_day_WAIV, rad_net_mean


print(E_EPM(), salinity_conversion, Tot_Evaporation)
