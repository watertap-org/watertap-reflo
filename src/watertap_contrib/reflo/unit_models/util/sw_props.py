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

import math

"""
Utility functions to calculate various saltwater properties for Solar Still Model
"""


def calculate_density(salinity, temperature):
    """
    Calculation of water density as a function of Temperature (°C) and salinity (g/l)
    Accuracy of correlation is valid for up to 160 g/l.
    However, intuitive, natural behaviour is reported for up to 350 g/l.
    """
    dens_coeff_1 = ((2 * salinity) - 150) / 150
    dens_coeff_2 = 0.5
    dens_coeff_3 = dens_coeff_1
    dens_coeff_4 = (2 * (dens_coeff_1**2)) - 1
    dens_coeff_5 = (
        (4.032 * dens_coeff_2) + (0.115 * dens_coeff_3) + ((3.26e-4) * dens_coeff_4)
    )
    dens_coeff_6 = (
        (-0.108 * dens_coeff_2)
        + ((1.571e-3) * dens_coeff_3)
        - ((4.23e-4) * dens_coeff_4)
    )
    dens_coeff_7 = (
        (-0.012 * dens_coeff_2) + ((1.74e-3) * dens_coeff_3) - ((9e-6) * dens_coeff_4)
    )
    dens_coeff_8 = (
        ((6.92e-4) * dens_coeff_2)
        - ((8.7e-5) * dens_coeff_3)
        - ((5.3e-5) * dens_coeff_4)
    )
    dens_coeff_9 = ((2 * temperature) - 200) / 160
    dens_coeff_10 = 0.5
    dens_coeff_11 = dens_coeff_9
    dens_coeff_12 = (2 * (dens_coeff_9**2)) - 1
    dens_coeff_13 = (4 * dens_coeff_9**3) - 3 * dens_coeff_9

    # saltwater density (kg/m^3)
    density = 1e3 * (
        (dens_coeff_5 * dens_coeff_10)
        + (dens_coeff_6 * dens_coeff_11)
        + (dens_coeff_7 * dens_coeff_12)
        + (dens_coeff_8 * dens_coeff_13)
    )

    return density


def calculate_viscosity(salinity, temperature):
    """
    Calculation of water dynamic viscosity as a function of Temperature (°C) and salinity (g/l).
    Accuracy of correlation is valid for up to 150 g/l.
    However, intuitive behaviour is reported for up to 350 g/l.
    """
    freshwater_dyn_visc = (4.2844e-5) + (
        0.157 * ((temperature + 64.993) ** 2) - 91.296
    ) ** -1  # Freshwater dynamic viscosity (Pa.s)
    dyn_visc_coeff_1 = (
        (1.474e-3) + ((1.5e-5) * temperature) - ((3.927e-8) * temperature**2)
    )
    dyn_visc_coeff_2 = (
        (1.073e-5) - ((8.5e-8) * temperature) + ((2.230e-10) * temperature**2)
    )
    dynamic_visc = freshwater_dyn_visc * (
        1 + (dyn_visc_coeff_1 * salinity) + (dyn_visc_coeff_2 * salinity**2)
    )  # Saltwater dynamic viscosity (Pa.s)
    return dynamic_visc


def calculate_specific_heat(salinity, temperature):
    """
    Calculation of water specific heat as a function of Temperature (°C) and salinity (g/l)
    Accuracy of correlation is valid for up to 180 g/l.
    However, intuitive natural behaviour is reported for up to 240 g/l.
    """
    freshwater_specific_heat = (
        3e-09 * temperature**4
        - 7e-07 * temperature**3
        + 7e-05 * temperature**2
        - 0.0029 * temperature
        + 4.2194
    )  # Specific heat of freshwater  (J / kg.K)
    cp_coeff_1 = 5.328 - ((9.76e-2) * salinity) + ((4.04e-4) * salinity**2)
    cp_coeff_2 = (-6.913e-3) + ((7.351e-4) * salinity) - ((3.15e-6) * salinity**2)
    cp_coeff_3 = (9.6e-6) - ((1.927e-6) * salinity) + ((8.23e-9) * salinity**2)
    cp_coeff_4 = (2.5e-9) + ((1.66e-9) * salinity) - ((7.125e-12) * salinity**2)
    specific_heat = 1e3 * (
        cp_coeff_1
        + (cp_coeff_2 * (temperature + 273))
        + ((cp_coeff_3) * (temperature + 273) ** 2)
        + ((cp_coeff_4) * (temperature + 273) ** 3)
    )  # Specific heat of saltwater  (J / kg.K)
    return specific_heat


def calculate_thermal_conductivity(salinity, temperature):
    """
    Calculation of water thermal conductivity as a function of Temperature (°C) and salinity (g/l)
    Accuracy of correlation is valid for up to 160 g/l.
    However, intuitive natural behaviour is reported for up to 350 g/l.
    """
    saltwater_thermal_conductivity_coefficient_1 = math.log10(240 + 0.0002 * salinity)
    saltwater_thermal_conductivity_coefficient_2 = 0.434 * (
        2.3 - ((343.5 + (0.037 * salinity)) / ((temperature + 273)))
    )
    saltwater_thermal_conductivity_coefficient_3 = abs(
        1 - ((temperature + 273) / (647 + 0.03 * salinity))
    ) ** (1 / 3)
    if saltwater_thermal_conductivity_coefficient_3 == "nan":  # Avoiding singularities
        saltwater_thermal_conductivity_coefficient_3 = 1

    log_base_10_thermal_conductivity = saltwater_thermal_conductivity_coefficient_1 + (
        saltwater_thermal_conductivity_coefficient_2
        * saltwater_thermal_conductivity_coefficient_3
    )
    thermal_conductivity = (
        10**log_base_10_thermal_conductivity
    ) / 1e3  # saltwater thermal conductivity (W / m. K)

    return thermal_conductivity
