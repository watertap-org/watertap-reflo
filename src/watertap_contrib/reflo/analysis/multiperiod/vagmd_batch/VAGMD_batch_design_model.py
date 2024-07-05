#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

"""
Calculate the number of periods to reach target recovery rate by solving the system first
"""


def get_n_time_points(
    dt=None,
    feed_flow_rate=600,
    evap_inlet_temp=80,
    cond_inlet_temp=25,
    feed_temp=25,
    feed_salinity=35,
    recovery_ratio=0.5,
    initial_batch_volume=50,
    module_type="AS7C1.5L",
    cooling_system_type="closed",
    cooling_inlet_temp=25,  # not required if cooling system type is "open"
):

    # Identify if the final brine salinity is larger than 175.3 for module "AS7C1.5L"
    final_brine_salinity = feed_salinity / (1 - recovery_ratio)  # g/L
    if module_type == "AS7C1.5L" and final_brine_salinity > 175.3:
        cooling_system_type = "closed"
        feed_flow_rate = 1100
        evap_inlet_temp = 80
        cond_inlet_temp = 25
        high_brine_salinity = True
    else:
        high_brine_salinity = False

    # Identify the module area
    if module_type == "AS7C1.5L":
        module_area = 7.2
    else:  # module_type = "AS26C7.2L"
        module_area = 25.92

    # Identify the time step of the simulation (second)
    if dt is None:
        if module_type == "AS7C1.5L":
            dt = 20352.55 / feed_flow_rate
        else:  # module_type == "AS26C7.2L"
            dt = 73269.19 / feed_flow_rate

    initial_status = _get_membrane_performance(
        evap_inlet_temp,
        feed_flow_rate,
        cond_inlet_temp,
        feed_salinity,
        module_type,
        high_brine_salinity,
    )

    PFlux = [initial_status[0]]  # Permeate flux, kg/h/m2
    TCO = [initial_status[1]]  # Condenser outlet temperature, oC
    TEO = [initial_status[2]]  # Evaporator outlet temperature, oC
    TCI = [cond_inlet_temp]  # Condenser inlet temperature, oC
    Ttank = [feed_temp]  # Feed temperature in thetank, oC
    V = [initial_batch_volume]  # Volume remaining in the tank, L
    Vd = [0]  # Volume of distillate, L
    PFR = [PFlux[0] * module_area]  # Permeate flow rate, L/h
    AccVd = [0]  # Accumulated volume of distillate, L
    RR = [0]  # Accumlated recovery ratio, dimensionless
    S = [feed_salinity]  # Feed salinity in the tank, g/L

    while S[-1] < final_brine_salinity:
        Vd.append(PFR[-1] * dt / 3600)
        V.append(V[-1] - Vd[-1])
        S.append(V[-2] / V[-1] * S[-1])
        Ttank.append(
            (feed_flow_rate * dt / 3600 * TEO[-1] + V[-2] * Ttank[-1])
            / (V[-2] + feed_flow_rate * dt / 3600)
        )

        if cooling_system_type == "closed":
            TCI.append(TCI[-1])

        current_status = _get_membrane_performance(
            evap_inlet_temp,
            feed_flow_rate,
            TCI[-1],
            S[-1],
            module_type,
            high_brine_salinity,
        )

        PFlux.append(current_status[0])
        TCO.append(current_status[1])
        TEO.append(current_status[2])
        PFR.append(PFlux[-1] * module_area)
        AccVd.append(Vd[-1] + AccVd[-1])
        RR.append(
            (
                1
                - (
                    (feed_flow_rate - PFR[-1])
                    / feed_flow_rate
                    * (S[-1] / S[0])
                    * (1 - RR[-1])
                )
            )
        )

    return len(PFlux)


def _get_membrane_performance(TEI, FFR, TCI, SgL, module_type, high_brine_salinity):
    # Model parameters
    PFluxAS26 = [
        0.798993148477908,
        0.314627216640160,
        0.559805181621833,
        -0.146236734128216,
        -0.659197144919924,
        0,
        0.185658514024503,
        -0.107221706014227,
        0,
        0,
        -0.187626469717738,
        0,
        0,
        0,
        0.128420664686447,
    ]
    PFluxAS7_high = [
        9.41014468300000,
        0,
        0,
        0,
        -0.0188989390000000,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    PFluxAS7_low = [
        4.82699491400000,
        1.37479276300000,
        1.91988177200000,
        -0.574212905000000,
        -0.641257664000000,
        0.399259954000000,
        0,
        0,
        0,
        0,
        0,
        0,
        -0.588924321000000,
        0,
        0,
    ]
    TCOAS26 = [
        65.1084685465240,
        9.15474718837607,
        -0.917460918908258,
        0.480070517276181,
        -1.06168979823129,
        0,
        0,
        0,
        0.142552983811052,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    TCOAS7_high = [
        67.0068599900000,
        0,
        0,
        0,
        -0.0145469190000000,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    TCOAS7_low = [
        58.8189910000000,
        8.79691366400000,
        -2.06741401400000,
        1.63187967600000,
        -0.914624645000000,
        -0.536574144000000,
        -0.249657477000000,
        0.398973861000000,
        -0.153760262000000,
        0.102355281000000,
        0.696768080000000,
        -0.300582958000000,
        -0.557410173000000,
    ]
    TEOAS26 = [
        29.2261439896435,
        0.569016152083381,
        0.824636694807529,
        4.62669502530487,
        1.37222105534565,
        0,
        0,
        0.220665657590258,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    TEOAS7_high = [
        36.2497021800000,
        0,
        0,
        0,
        0.0126951860000000,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    TEOAS7_low = [
        34.4309511000000,
        1.55140768500000,
        1.85928314600000,
        4.52887180500000,
        1.10791196800000,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0.808447211000000,
    ]
    coefficients = [
        [-5.68487382500000, 0.0705622560000000, 0.000152146000000000],
        [-1.58460599600000, 0.00102338700000000, 1.20000000000000e-06],
        [-4.27697973100000, 0.175533630000000, -0.000178178000000000],
        [-1.49331349500000, 0.0146627780000000, 5.62000000000000e-06],
        [-7, 0.100000000000000, 0],
        [-2.14285714285714, 0.00285714285714286, 0],
        [-5, 0.200000000000000, 0],
        [-1.33333333333333, 0.00950000000000000, 0],
    ]
    a1 = 0.983930048493388
    a2 = -4.8359231959954e-04
    S_c = a1 * SgL + a2 * SgL**2  # [g/kg]

    # TEI -= 273.15
    # TCI -= 273.15

    surrogate_variables = [
        [1, TEI, TEI**2],
        [1, FFR, FFR**2],
        [1, TCI, TCI**2],
        [1, S_c, S_c**2],
    ]

    # Model calculations
    if module_type == "AS7C1.5L":
        if high_brine_salinity:
            S_r = S_c

            PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high
        else:
            TEI = sum(
                surrogate_variables[0][j] * coefficients[0][j]
                for j in range(len(coefficients[0]))
            )
            FFR = sum(
                surrogate_variables[1][j] * coefficients[1][j]
                for j in range(len(coefficients[1]))
            )
            TCI = sum(
                surrogate_variables[2][j] * coefficients[2][j]
                for j in range(len(coefficients[2]))
            )
            S_r = sum(
                surrogate_variables[3][j] * coefficients[3][j]
                for j in range(len(coefficients[3]))
            )

            PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_low, TCOAS7_low, TEOAS7_low

        VarsAS7 = [
            1,
            TEI,
            FFR,
            TCI,
            S_r,
            FFR * TEI,
            TCI * TEI,
            S_r * TEI,
            FFR * TCI,
            FFR * S_r,
            S_r * TCI,
            TEI**2,
            FFR**2,
            TCI**2,
            S_r**2,
        ]
        VarsAS7_TCO = [
            1,
            TEI,
            FFR,
            TCI,
            S_r,
            FFR * TEI,
            S_r * TEI,
            FFR * TCI,
            FFR * S_r,
            S_r * TCI,
            FFR**2,
            S_r**2,
            FFR**3,
        ]

        PFlux = sum(VarsAS7[j] * PFluxAS7[j] for j in range(len(VarsAS7)))
        TCO = sum(VarsAS7_TCO[j] * TCOAS7[j] for j in range(len(VarsAS7_TCO)))
        TEO = sum(VarsAS7[j] * TEOAS7[j] for j in range(len(VarsAS7)))

    else:  # module_type == "AS26C7.2L":
        TEI = sum(
            surrogate_variables[0][j] * coefficients[4][j]
            for j in range(len(coefficients[0]))
        )
        FFR = sum(
            surrogate_variables[1][j] * coefficients[5][j]
            for j in range(len(coefficients[1]))
        )
        TCI = sum(
            surrogate_variables[2][j] * coefficients[6][j]
            for j in range(len(coefficients[2]))
        )
        S_r = sum(
            surrogate_variables[3][j] * coefficients[7][j]
            for j in range(len(coefficients[3]))
        )

        VarsAS26 = [
            1,
            TEI,
            FFR,
            TCI,
            S_r,
            TCI * TEI,
            FFR * TEI,
            S_r * TEI,
            FFR * TCI,
            S_r * TCI,
            FFR * S_r,
            TEI**2,
            FFR**2,
            TCI**2,
            S_r**2,
        ]

        PFlux = sum(VarsAS26[j] * PFluxAS26[j] for j in range(len(VarsAS26)))
        TCO = sum(VarsAS26[j] * TCOAS26[j] for j in range(len(VarsAS26)))
        TEO = sum(VarsAS26[j] * TEOAS26[j] for j in range(len(VarsAS26)))

    return [PFlux, TCO, TEO, S_c]


#%% Test block

if __name__ == "__main__":

    test_case = get_n_time_points()
    print(test_case)
