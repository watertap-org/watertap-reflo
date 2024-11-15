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
Stephen_Boltzman = 5.6697 * (10**-8)  # Stephen-Boltzman Constant (W / (m^2.K^4))
thickness_insulation = 0.005  # Thickness of SS Insulation (m)
conductivity_insulation = 0.033  # Conductivity of SS Insulation (W/m.K)
thickness_glass = 0.004  # Thickness of glass m
conductivity_glass = 1.03  # Conductivity of glass W/m.K
gravity = 9.81  # Accelaration due to gravity (m/s^2)

No_months = 1
Hours = 24  # Determine how many hours are being considered
seconds = 3600 * Hours  # Converting hours to seconds

# User input from pyomo.ipynb
initial_salinity = float(Input_data.iloc[-1, 0])  # Water salinity (g/l)
initial_depth = float(Input_data.iloc[-1, 2])  # water depth in SS basin (m)
side_basin = float(Input_data.iloc[-1, 1])  # Side of SS basin (m)


def SS_ZLD_model(TotalYield):
    return TotalYield


area_bottom_basin = side_basin**2  # Area of square basin (m^2)


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

def create_input_arrays():

    global V_wind_hr_list

    T_ambient_hr_list = []
    irradiance_hr_list = []
    V_wind_hr_list = []

    # Collecting hourly weather data
    for day in continuous_day_series:
        if day == 365:
            day = 364
        start_row = day * 24
        end_row = start_row + 24
        for kk in range(24):
            T_ambient_hr = float(Input_data.iloc[start_row + kk, 7])
            irradiance_hr = float(Input_data.iloc[start_row + kk, 4]) + 10
            V_wind_hr = float(Input_data.iloc[start_row + kk, 11])

            if irradiance_hr > 20:
                T_ambient_hr_list.append(T_ambient_hr)
                irradiance_hr_list.append(irradiance_hr)
                V_wind_hr_list.append(V_wind_hr)

    T_ambient_hr_array = pd.Series(T_ambient_hr_list)
    irradiance_hr_array = pd.Series(irradiance_hr_list)
    V_wind_hr_array = pd.Series(V_wind_hr_list)

    return T_ambient_hr_array, irradiance_hr_array, V_wind_hr_array

T_ambient_hr_array, irradiance_hr_array, V_wind_hr_array = create_input_arrays()

# print(V_wind_hr_array.to_list() == V_wind_hr_list)
# assert False

# Allocating the variables
# Second variables
depth = np.zeros(len(irradiance_hr_array) * 3600)
salinity = np.zeros(len(irradiance_hr_array) * 3600)
excess_salinity = np.zeros(len(irradiance_hr_array) * 3600)
volume_scale_formation = np.zeros(len(irradiance_hr_array) * 3600)
thickness_scale_formation = np.zeros(len(irradiance_hr_array) * 3600)
Evap_mass_sw = np.zeros(len(irradiance_hr_array) * 3600)
Evap_mass_sw_sec = np.zeros(len(irradiance_hr_array) * 3600)
irradiance = np.zeros(len(irradiance_hr_array) * 3600)
salt_precipitated = np.zeros(len(irradiance_hr_array) * 3600)
sw_mass = np.zeros(len(irradiance_hr_array) * 3600)
fw_mass = np.zeros(len(irradiance_hr_array) * 3600)
t = np.zeros(len(irradiance_hr_array) * 3600)
V_wind = np.zeros(len(irradiance_hr_array) * 3600)
T_ambient = np.zeros(len(irradiance_hr_array) * 3600)
T_basin = np.zeros(len(irradiance_hr_array) * 3600)
T_glass = np.zeros(len(irradiance_hr_array) * 3600)
T_sky = np.zeros(len(irradiance_hr_array) * 3600)
T_water = np.zeros(len(irradiance_hr_array) * 3600)
Total_water_mass = np.zeros(len(irradiance_hr_array) * 3600)

# Hour variables
Mass_saltwater = len(irradiance_hr_array)
MassWater_Com = len(irradiance_hr_array)
irradiance_hr = len(irradiance_hr_array)
V_wind_hr = len(irradiance_hr_array)
T_ambient_hr = len(irradiance_hr_array)

# Converting hourly data into per second
for hour in range(len(irradiance_hr_array)):
    start_idx = 3600 * hour
    end_idx = start_idx + 3600
    irradiance[start_idx:end_idx] = irradiance_hr_array[hour]
    V_wind[start_idx:end_idx] = V_wind_hr_array.loc[hour]
    T_ambient[start_idx:end_idx] = T_ambient_hr_array[hour]


# Initializing Temperatures
# Initial system is assumed to be in thermal equilibrium with ambient
T_water[1] = T_ambient_hr_array[1]  # Initial water temperature (°C)
T_basin[1] = T_ambient_hr_array[1]  # Initial basin temperature (°C)
T_glass[1] = T_ambient_hr_array[1]  # Initial glass temperature (°C)


B = ((2 * initial_salinity) - 150) / 150
G1 = 0.5
G2 = B
G3 = (2 * (B**2)) - 1
A1 = (4.032 * G1) + (0.115 * G2) + ((3.26e-4) * G3)
A2 = (-0.108 * G1) + ((1.571e-3) * G2) - ((4.23e-4) * G3)
A3 = (-0.012 * G1) + ((1.74e-3) * G2) - ((9e-6) * G3)
A4 = ((6.92e-4) * G1) - ((8.7e-5) * G2) - ((5.3e-5) * G3)
A_rho = ((2 * T_water[1]) - 200) / 160
F1 = 0.5
F2 = A_rho
F3 = (2 * (A_rho**2)) - 1
F4 = (4 * A_rho**3) - 3 * A_rho
rhow_0 = 1e3 * ((A1 * F1) + (A2 * F2) + (A3 * F3) + (A4 * F4))

density_NaCl = 2165  # density of sodium chloride (kg/m^3)
K_NaCl = 3.5  # thermal conductivity of sodium chloride (W/m.K)

# Geometrical properties of squared basin
area_bottom_basin = side_basin**2  # Area of square basin (m^2)
salt_precipitated[1] = 0
salt_precipitated[0] = salt_precipitated[1]
maximum_solubility = 365  # maximum solubility of salt in water (g/l)

salinity[1] = initial_salinity
salinity[0] = initial_salinity
depth[1] = initial_depth
depth[0] = initial_depth
sw_mass[1] = depth[1] * rhow_0 * area_bottom_basin
fw_mass[1] = sw_mass[1] / (1 + salinity[1] / 1000)
Salt_mass = (salinity[1] * fw_mass[1]) / 1000  # Mass of Sodium Chloride (kg)
excess_salinity[1] = salinity[1]  # salinity without maximum solublity (g/l)

T_sky[1] = 0.0552 * (
    (T_ambient[1]) ** 1.5
)  # Initial effective radiation temperature of the sky (deg C), converted into K and then back to Celsius to avoid

t[1] = 1

for i in range(2, len(irradiance), 1):

    if depth[i - 1] <= 0 or fw_mass[i - 1] <= 0:
        depth[i - 1] = initial_depth
        salinity[i - 1] = initial_salinity
        sw_mass[i - 1] = depth[i - 1] * rhow_0 * area_bottom_basin
        fw_mass[i - 1] = sw_mass[i - 1] / (1 + salinity[i] / 1000)

    t[i] = t[i - 1] + 1

    # Avoiding singularities
    Tdif = T_water[i - 1] - T_glass[i - 1]
    if Tdif <= 0:
        Tdif = 0.01

    # Avoiding singularities
    Tdif1 = T_glass[i - 1] - T_ambient[i - 1]
    if Tdif1 <= 0:
        Tdif1 = 0.01

    T_sky[i] = 0.0552 * (
        (T_ambient[i]) ** 1.5
    )  # T_sky Effective radiation temperature of the sky

    area_side_water = (2 * (2 * side_basin)) * depth[
        i - 1
    ]  # Perimeter x depth of water (m^2)

    #        Thermophysical properties calculation as a function of temperature (°C) and salinity (g/l)
    # Calculation of water density as a function of Temperature (°C) and salinity (g/l)
    # Accuracy of correlation is valid for up to 160 g/l; however, intuitive, natural behaviour is reported for up to 350 g/l
    # Atop, according to parametric analysis and available literature, density plays a minimal role in the saltwater evaporation performance
    # Coefficients for density calculation
    B = ((2 * salinity[i - 1]) - 150) / 150
    G1 = 0.5
    G2 = B
    G3 = (2 * (B**2)) - 1
    A1 = (4.032 * G1) + (0.115 * G2) + ((3.26e-4) * G3)
    A2 = (-0.108 * G1) + ((1.571e-3) * G2) - ((4.23e-4) * G3)
    A3 = (-0.012 * G1) + ((1.74e-3) * G2) - ((9e-6) * G3)
    A4 = ((6.92e-4) * G1) - ((8.7e-5) * G2) - ((5.3e-5) * G3)
    A_rho = ((2 * T_water[i - 1]) - 200) / 160
    F1 = 0.5
    F2 = A_rho
    F3 = (2 * (A_rho**2)) - 1
    F4 = (4 * A_rho**3) - 3 * A_rho
    density_sw = 1e3 * (
        (A1 * F1) + (A2 * F2) + (A3 * F3) + (A4 * F4)
    )  # saltwater density (kg/m^3)

    # Calculation of water dynamic viscosity as a function of Temperature (°C) and salinity (g/l)
    # Accuracy of correlation is valid for up to 150 g/l; however, intuitive behaviour is reported for up to 350 g/l
    # Atop, according to parametric analysis and available literature, dynamic viscosity plays a minimal role in the saltwater evaporation performance
    # Coefficients for dynamic viscosity calculation
    dy_visc_sw = (4.2844e-5) + (
        0.157 * ((T_water[i - 1] + 64.993) ** 2) - 91.296
    ) ** -1  # Freshwater dynamic viscosity (Pa.s)

    # Calculation of water specific heat as a function of Temperature (°C) and salinity (g/l)
    # Accuracy of correlation is valid for up to 180 g/l; however, intuitive natural behaviour is reported for up to 240 g/l.
    # Atop, according to parametric analysis and available literature, specific heat plays a minimal role in the saltwater evaporation performance
    # Coefficients for specific heat calculation
    CP_sw = 1e3 * (
        3e-09 * T_water[i - 1] ** 4
        - 7e-07 * T_water[i - 1] ** 3
        + 7e-05 * T_water[i - 1] ** 2
        - 0.0029 * T_water[i - 1]
        + 4.2194
    )  # Specific heat of freshwater  (J / kg.K)

    # Calculation of water thermal conductivity as a function of Temperature (°C) and salinity (g/l)
    # Accuracy of correlation is valid for up to 160 g/l; however, intuitive natural behaviour is reported for up to 350 g/l
    # Atop, according to parametric analysis and available literature, specific heat plays a minimal role in the saltwater evaporation performance
    # Coefficients for thermal conductivity calculation
    coeff_1_Ksw = math.log10(240 + 0.0002 * salinity[i - 1])
    coeff_2_Ksw = 0.434 * (
        2.3 - ((343.5 + (0.037 * salinity[i - 1])) / ((T_water[i - 1] + 273)))
    )
    coeff_3_Ksw = abs(
        1 - ((T_water[i - 1] + 273) / (647 + 0.03 * salinity[i - 1]))
    ) ** (1 / 3)
    if coeff_3_Ksw == "nan":  # Avoiding singularities
        coeff_3_Ksw = 1

    log10_K_sw = coeff_1_Ksw + (coeff_2_Ksw * coeff_3_Ksw)
    K_sw = (10**log10_K_sw) / 1e3  # saltwater thermal conductivity (W / m. K)

    # mass_sw = density_sw * (area_bottom_basin * depth[i-1])  # Water Mass (kg)
    kinem_visc_sw = dy_visc_sw / density_sw
    Prandtl_sw = (CP_sw * dy_visc_sw) / K_sw  # Prandtl number for water (-)

    latent_heat_fw = (
        2501.67 - 2.389 * T_water[i - 1]
    ) * 1000  # Latent heat of vaporization of pure water J/kg

    # Calculation of partial saturated vapor pressure of saltwater
    # According to parametric analysis and available literature, the partial vapor pressure plays a major role in the evaporation of water.
    aw = (
        -0.000566 * salinity[i - 1]
    ) + 0.99853070  # A coefficient obtained from the water molar fraction in salt solutions from 0-350 g/l Paper: Kokya and Kokya
    P_sw = aw * math.exp(
        25.317 - (5144 / (T_water[i - 1] + 273))
    )  # Partial saturated vapor pressure at a saltwater temperature (N/m^2)
    P_glass = math.exp(
        25.317 - (5144 / (T_glass[i - 1] + 273))
    )  # Partial saturated vapor pressure at glass cover temperature (N/m^2)

    # heat transfer coefficients calculation
    C = 0.54
    N = 1 / 4
    beta = 1e-6 * (
        -0.000006 * T_water[i - 1] ** 4
        + 0.001667 * T_water[i - 1] ** 3
        - 0.197796 * T_water[i - 1] ** 2
        + 16.862446 * T_water[i - 1]
        - 64.319951
    )  # Coefficient of volume expansion (1/°C) correlation gotten from Zhutovsky and Kovler (2015)
    Grashof = abs(
        (gravity * beta * (T_basin[i - 1] - T_water[i - 1]) * (depth[i - 1] ** 3))
        / (kinem_visc_sw**2)
    )  # Grashof number
    hw = abs(
        (K_sw / depth[i - 1]) * C * (Grashof * Prandtl_sw) ** N
    )  # Heat transfer coefficient of water layer
    hC_wgi = 0.884 * (
        (
            abs(
                (
                    (T_water[i - 1] - T_glass[i - 1])
                    + (((P_sw - P_glass) * (T_water[i - 1] + 273.15)) / (268900 - P_sw))
                )
            )
        )
        ** (1 / 3)
    )  # Convective heat transfer coefficient (W/m^2.°C) (Dunkel)
    hR_wgi = abs(
        emissivity_water
        * Stephen_Boltzman
        * (
            (((T_water[i - 1] + 273) ** 2) + ((T_glass[i - 1] + 273) ** 2))
            * (T_water[i - 1] + T_glass[i - 1] + 546)
        )
    )  # Radiative heat transfer coefficient  from basin water to glass cover (W/m^2.°C)
    hE_wgi = abs(
        0.01628 * hC_wgi * ((P_sw - P_glass) / (Tdif))
    )  # Evaporative heat transfer coefficient from basin water to glass cover (W/m2.°C)
    hT_wgi = (
        hC_wgi + hR_wgi + hE_wgi
    )  # Total heat transfer coefficient from basin water to glass cover (W/m2°C)
    hRgoa = (
        Stephen_Boltzman
        * emissivity_glass
        * ((((T_glass[i - 1] + 273) ** 4) - ((T_sky[i - 1] + 273) ** 4)) / (Tdif1))
    )  # Radiative heat transfer coefficient from glass cover to ambient (W/m^2.°C)
    if V_wind[i] > 5:
        hcba = 2.8 + (
            3.0 * V_wind[i]
        )  # Convective heat transfer coefficient from basin to ambient (W/m2°C)
        hCgoa = 2.8 + (
            3.0 * V_wind[i]
        )  # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)
    else:
        hcba = 2.8 + (
            3.8 * V_wind[i]
        )  # Convective heat transfer coefficient from basin to ambient (W/m2°C)
        hCgoa = 2.8 + (
            3.8 * V_wind[i]
        )  # Convective heat transfer coefficient from glass cover to ambient (W/m2°C)

    hTgoa = (
        hCgoa + hRgoa
    )  # total heat loss coefficient from the glass cover to the outer atmosphere
    htba = 1 / (
        (thickness_insulation / conductivity_insulation) + (1 / (hcba))
    )  # heat loss coefficient from basin liner to the atmosphere

    alpha_eff = (
        (absorptivity_effective_basin * (hw / (hw + htba + hcba)))
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
    Ub = (hw * htba) / (
        hw + htba
    )  # Overall bottom heat transfer coefficient from bottom to ambient (W/m^2.°C)
    Uss = (area_side_water / area_bottom_basin) * Ub
    Ubs = Ub + Uss
    ULs = UT + Ubs

    # Present [i] temperature calculation
    a = ULs / (sw_mass[i - 1] * CP_sw)
    f = ((alpha_eff * irradiance[i]) + (ULs * T_ambient[i])) / (sw_mass[i - 1] * CP_sw)
    T_glass[i] = (
        (absorptivity_effective_glass * irradiance[i])
        + (hT_wgi * T_water[i - 1])
        + (UTgi_a * T_ambient[i])
    ) / (hT_wgi + UTgi_a)
    T_water[i] = (f / a) * (1 - np.exp(-a * t[i])) + (T_water[i] * np.exp(-a * t[i]))
    T_basin[i] = (
        (absorptivity_effective_basin * irradiance[i])
        + (hw * T_water[i - 1])
        + ((htba + hcba) * T_basin[i - 1])
    ) / (hw + htba + hcba)

    # Evaporation estimation of freshwater and saltwater
    Evap_mass_fw = (
        area_bottom_basin * hE_wgi * (T_water[i] - T_glass[i]) * 3600
    ) / latent_heat_fw  # Distillated (kg)

    Evap_mass_fw_sec = (
        area_bottom_basin * hE_wgi * (Tdif)
    ) / latent_heat_fw  # Distillated (kg)

    Evap_mass_sw[i] = Evap_mass_fw / (
        1 + salinity[i - 1] / 1e3
    )  # Distillated saltwater conversion (Morton) (kg)
    Evap_mass_sw_sec[i] = Evap_mass_fw_sec / (
        1 + salinity[i - 1] / 1e3
    )  # Distillated saltwater conversion (Morton) (kg)

    fw_mass[i] = fw_mass[i - 1] - Evap_mass_sw_sec[i]  # Current Freshwater (kg)
    sw_mass[i] = fw_mass[i] + Salt_mass  # Current saltwater (kg)
    depth[i] = sw_mass[i] / (density_sw * area_bottom_basin)  # Current depth (m)

    if salinity[i - 1] >= maximum_solubility:
        salinity[i] = maximum_solubility

        excess_salinity[i] = (Salt_mass * 1000) / fw_mass[
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
        salinity[i] = (Salt_mass * 1000) / fw_mass[i]
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
