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

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
    Set,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.core.solvers import get_solver
from watertap.core.util.initialization import check_dof

from watertap_contrib.reflo.property_models import (
    AirWaterEq,
    AirWaterEqStateBlock,
    MolarVolumeCalculation,
    LiqDiffusivityCalculation,
    VapDiffusivityCalculation,
)

solver = get_solver()


@pytest.fixture(scope="module")
def m1():
    """
    Test NMSU case study 1 results
    """
    props = {
        "solute_list": ["TCA"],
        "mw_data": {"TCA": 0.1334},
        "dynamic_viscosity_data": {"Liq": 0.00115, "Vap": 1.75e-5},
        "henry_constant_data": {"TCA": 0.725},  # salinity adjusted
        "standard_enthalpy_change_data": {"TCA": 28.7e3},
        "temperature_boiling_data": {"TCA": 347},
        "molar_volume_data": {"TCA": 9.81e-5},
        "critical_molar_volume_data": {"TCA": 2.94e-4},
        "density_data": {"Liq": 999.15, "Vap": 1.22},
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block1(m1):
    m = m1
    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 3
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air", "TCA"]
    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    assert isinstance(m.fs.properties.solvent_set, Set)
    assert len(m.fs.properties.solvent_set) == 2
    for j in m.fs.properties.solvent_set:
        assert j in ["H2O", "Air"]
    assert isinstance(m.fs.properties.solute_set, Set)
    assert len(m.fs.properties.solute_set) == 1
    assert len(m.fs.properties.solute_set) == len(m.fs.properties.config.solute_list)
    for j in m.fs.properties.solute_set:
        assert j in ["TCA"]
    assert (len(m.fs.properties.solvent_set) + len(m.fs.properties.solute_set)) == len(
        m.fs.properties.component_list
    )

    assert isinstance(m.fs.properties.liq_solute_set, Set)
    assert len(m.fs.properties.liq_solute_set) == 2
    for (p, j) in m.fs.properties.liq_solute_set:
        assert p == "Liq"
        assert j in m.fs.properties.liq_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.vap_solute_set, Set)
    assert len(m.fs.properties.vap_solute_set) == 2
    for (p, j) in m.fs.properties.vap_solute_set:
        assert p == "Vap"
        assert j in m.fs.properties.vap_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.phase_solute_set, Set)
    assert len(m.fs.properties.phase_solute_set) == len(m.fs.properties.solute_set) * 2
    for (p, j) in m.fs.properties.phase_solute_set:
        assert p in m.fs.properties.phase_list
        assert j in m.fs.properties.solute_set

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["TCA"].value == 133.4e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.critical_molar_volume_comp, Var)
    assert isinstance(m.fs.properties.henry_constant_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert m.fs.properties.config.temp_adjust_henry
    assert isinstance(m.fs.properties.enth_change_dissolution_comp, Var)

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation
        == LiqDiffusivityCalculation.HaydukLaudie
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation
        == VapDiffusivityCalculation.WilkeLee
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation
        == MolarVolumeCalculation.TynCalus
    )


@pytest.mark.component
def test_properties1(m1):
    m = m1
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]

    stream.flow_mass_phase_comp["Liq", "H2O"].fix(157.8657)
    stream.flow_mass_phase_comp["Liq", "TCA"].fix(2.61e-5)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(1.34932)
    stream.flow_mass_phase_comp["Vap", "TCA"].fix(0)
    stream.temperature["Liq"].fix(288)
    stream.temperature["Vap"].fix(288)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mole_phase_comp[...]
    stream.conc_mole_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]

    stream_props = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "dens_mass_phase",
        "mass_frac_phase_comp",
        "flow_mole_phase_comp",
        "conc_mole_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "collision_molecular_separation_comp",
        "collision_function_comp",
        "collision_function_zeta_comp",
        "collision_function_ee_comp",
        "energy_molecular_attraction_phase_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_constant_comp",
    ]

    for prop in stream_props:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325,
        "temperature": {"Liq": 288, "Vap": 288},
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): 157.8657,
            ("Liq", "TCA"): 2.61e-05,
            ("Vap", "Air"): 1.34932,
            ("Vap", "TCA"): 0,
        },
        "conc_mass_phase_comp": {
            ("Liq", "H2O"): 999.1498348101538,
            ("Liq", "TCA"): 0.00016518984610681744,
            ("Vap", "Air"): 1.22,
            ("Vap", "TCA"): 1.1148142932658939e-17,
        },
        "dens_mass_phase": {"Liq": 999.15, "Vap": 1.22},
        "mass_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999998346696081,
            ("Liq", "TCA"): 1.653303769272056e-07,
            ("Vap", "Air"): 0.9999999999999851,
            ("Vap", "TCA"): 1.1148143044160013e-17,
        },
        "flow_mole_phase_comp": {
            ("Liq", "H2O"): 8770.316666666668,
            ("Liq", "TCA"): 0.0001956521739130435,
            ("Vap", "Air"): 46.52827586206896,
            ("Vap", "TCA"): 1.1276171193672585e-16,
        },
        "conc_mole_phase_comp": {
            ("Liq", "H2O"): 55508.32415611966,
            ("Liq", "TCA"): 0.001238304693454404,
            ("Vap", "Air"): 42.06896551724138,
            ("Vap", "TCA"): 8.356928735133821e-17,
        },
        "mole_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999999776915353,
            ("Liq", "TCA"): 2.2308449852439394e-08,
            ("Vap", "Air"): 0.9999999999999851,
            ("Vap", "TCA"): 2.423899789277319e-18,
        },
        "diffus_phase_comp": {
            ("Liq", "TCA"): 7.092554312734387e-10,
            ("Vap", "TCA"): 7.932456576071473e-06,
        },
        "molar_volume_comp": {"TCA": 0.00011007098737879485},
        "collision_molecular_separation_comp": {
            "Air": 0.3711,
            "TCA": 0.5655091393428981,
        },
        "collision_function_comp": {"TCA": 0.5868361666899444},
        "collision_function_zeta_comp": {"TCA": -0.2314831284520491},
        "collision_function_ee_comp": {"TCA": 0.20012379218586487},
        "energy_molecular_attraction_phase_comp": {
            ("Vap", "Air"): 1.085190114e-14,
            ("Vap", "TCA"): 5.7969309563e-14,
        },
        "flow_vol_phase": {"Liq": 0.1580000261222039, "Vap": 1.106},
        "flow_mass_phase": {"Liq": 157.86572610000002, "Vap": 1.34932},
        "henry_constant_comp": {"TCA": 0.4849154431401335},
    }

    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m2():
    """
    Test NMSU case study 2 results
    """
    props = {
        "solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.146},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 0.00011007},
        "critical_molar_volume_data": {"DCP": 2.26e-4},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block2(m2):
    m = m2
    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 3
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air", "DCP"]
    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    assert isinstance(m.fs.properties.solvent_set, Set)
    assert len(m.fs.properties.solvent_set) == 2
    for j in m.fs.properties.solvent_set:
        assert j in ["H2O", "Air"]
    assert isinstance(m.fs.properties.solute_set, Set)
    assert len(m.fs.properties.solute_set) == 1
    assert len(m.fs.properties.solute_set) == len(m.fs.properties.config.solute_list)
    for j in m.fs.properties.solute_set:
        assert j in ["DCP"]
    assert (len(m.fs.properties.solvent_set) + len(m.fs.properties.solute_set)) == len(
        m.fs.properties.component_list
    )

    assert isinstance(m.fs.properties.liq_solute_set, Set)
    assert len(m.fs.properties.liq_solute_set) == 2
    for (p, j) in m.fs.properties.liq_solute_set:
        assert p == "Liq"
        assert j in m.fs.properties.liq_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.vap_solute_set, Set)
    assert len(m.fs.properties.vap_solute_set) == 2
    for (p, j) in m.fs.properties.vap_solute_set:
        assert p == "Vap"
        assert j in m.fs.properties.vap_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.phase_solute_set, Set)
    assert len(m.fs.properties.phase_solute_set) == len(m.fs.properties.solute_set) * 2
    for (p, j) in m.fs.properties.phase_solute_set:
        assert p in m.fs.properties.phase_list
        assert j in m.fs.properties.solute_set

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.critical_molar_volume_comp, Var)
    assert isinstance(m.fs.properties.henry_constant_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert m.fs.properties.config.temp_adjust_henry
    assert isinstance(m.fs.properties.enth_change_dissolution_comp, Var)

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation
        == LiqDiffusivityCalculation.HaydukLaudie
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation
        == VapDiffusivityCalculation.WilkeLee
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation
        == MolarVolumeCalculation.TynCalus
    )


@pytest.mark.component
def test_properties2(m2):
    m = m2
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]

    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mole_phase_comp[...]
    stream.conc_mole_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_constant_comp[...]

    stream_props = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "dens_mass_phase",
        "mass_frac_phase_comp",
        "flow_mole_phase_comp",
        "conc_mole_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "collision_molecular_separation_comp",
        "collision_function_comp",
        "collision_function_zeta_comp",
        "collision_function_ee_comp",
        "energy_molecular_attraction_phase_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_constant_comp",
    ]

    for prop in stream_props:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325,
        "temperature": {"Liq": 283, "Vap": 283},
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): 99.97,
            ("Liq", "DCP"): 0.0001,
            ("Vap", "Air"): 7.482,
            ("Vap", "DCP"): 0,
        },
        "conc_mass_phase_comp": {
            ("Liq", "H2O"): 999.6990000010003,
            ("Liq", "DCP"): 0.0009999989997009106,
            ("Vap", "Air"): 1.247,
            ("Vap", "DCP"): 3.780308869903483e-14,
        },
        "dens_mass_phase": {"Liq": 999.7, "Vap": 1.247},
        "mass_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999989997006304,
            ("Liq", "DCP"): 1.0002990894277389e-06,
            ("Vap", "Air"): 0.99999999999972,
            ("Vap", "DCP"): 3.78030890770692e-14,
        },
        "flow_mole_phase_comp": {
            ("Liq", "H2O"): 5553.888888888889,
            ("Liq", "DCP"): 0.0008851124092759781,
            ("Vap", "Air"): 258.0,
            ("Vap", "DCP"): 2.5034759468456273e-12,
        },
        "conc_mole_phase_comp": {
            ("Liq", "H2O"): 55538.8333333889,
            ("Liq", "DCP"): 0.008851115238988409,
            ("Vap", "Air"): 42.99999999999999,
            ("Vap", "DCP"): 3.345998291647399e-13,
        },
        "mole_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999998406317013,
            ("Liq", "DCP"): 1.593680186847315e-07,
            ("Vap", "Air"): 0.9999999999997105,
            ("Vap", "DCP"): 1.0036286347773076e-14,
        },
        "diffus_phase_comp": {
            ("Liq", "DCP"): 7.210419079256757e-10,
            ("Vap", "DCP"): 8.579253905876393e-06,
        },
        "molar_volume_comp": {"DCP": 8.355077914769306e-05},
        "collision_molecular_separation_comp": {
            "Air": 0.3711,
            "DCP": 0.515860381978062,
        },
        "collision_function_comp": {"DCP": 0.5983130810151872},
        "collision_function_zeta_comp": {"DCP": -0.22307150200399972},
        "collision_function_ee_comp": {"DCP": 0.17911045474853254},
        "energy_molecular_attraction_phase_comp": {
            ("Vap", "Air"): 1.085190114e-14,
            ("Vap", "DCP"): 6.16613030539e-14,
        },
        "flow_vol_phase": {"Liq": 0.100000100030009, "Vap": 6.0},
        "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.482000000000001},
        "henry_constant_comp": {"DCP": 0.07506173060999055},
    }

    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m3():
    """
    Test when molar volume, diffusivity, and Henry constant are all direct user input via CONFIG
    """
    props = {
        "solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "diffusivity_data": {("Liq", "DCP"): 7.21045e-10, ("Vap", "DCP"): 8.57928e-6},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.0525},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 8.355e-5},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
        "temp_adjust_henry": False,
        "molar_volume_calculation": "none",
        "liq_diffus_calculation": "none",
        "vap_diffus_calculation": "none",
    }
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_parameter_block3(m3):
    m = m3
    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 3
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air", "DCP"]
    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    assert isinstance(m.fs.properties.solvent_set, Set)
    assert len(m.fs.properties.solvent_set) == 2
    for j in m.fs.properties.solvent_set:
        assert j in ["H2O", "Air"]
    assert isinstance(m.fs.properties.solute_set, Set)
    assert len(m.fs.properties.solute_set) == 1
    assert len(m.fs.properties.solute_set) == len(m.fs.properties.config.solute_list)
    for j in m.fs.properties.solute_set:
        assert j in ["DCP"]
    assert (len(m.fs.properties.solvent_set) + len(m.fs.properties.solute_set)) == len(
        m.fs.properties.component_list
    )

    assert isinstance(m.fs.properties.liq_solute_set, Set)
    assert len(m.fs.properties.liq_solute_set) == 2
    for (p, j) in m.fs.properties.liq_solute_set:
        assert p == "Liq"
        assert j in m.fs.properties.liq_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.vap_solute_set, Set)
    assert len(m.fs.properties.vap_solute_set) == 2
    for (p, j) in m.fs.properties.vap_solute_set:
        assert p == "Vap"
        assert j in m.fs.properties.vap_comps
        assert j in m.fs.properties.component_list

    assert isinstance(m.fs.properties.phase_solute_set, Set)
    assert len(m.fs.properties.phase_solute_set) == len(m.fs.properties.solute_set) * 2
    for (p, j) in m.fs.properties.phase_solute_set:
        assert p in m.fs.properties.phase_list
        assert j in m.fs.properties.solute_set

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.critical_molar_volume_comp, Var)
    assert isinstance(m.fs.properties.henry_constant_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert not m.fs.properties.config.temp_adjust_henry
    assert not hasattr(m.fs.properties, "enth_change_dissolution_comp")

    assert isinstance(
        m.fs.properties.config.liq_diffus_calculation, LiqDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.liq_diffus_calculation == LiqDiffusivityCalculation.none
    )
    assert isinstance(
        m.fs.properties.config.vap_diffus_calculation, VapDiffusivityCalculation
    )
    assert (
        m.fs.properties.config.vap_diffus_calculation == VapDiffusivityCalculation.none
    )
    assert isinstance(
        m.fs.properties.config.molar_volume_calculation, MolarVolumeCalculation
    )
    assert (
        m.fs.properties.config.molar_volume_calculation == MolarVolumeCalculation.none
    )


@pytest.mark.component
def test_properties3(m3):
    m = m3
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]

    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(1e-3)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mole_phase_comp[...]
    stream.conc_mole_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_constant_comp[...]

    assert not hasattr(stream, "tyn_calus_param")
    assert not hasattr(stream, "tyn_calus_exponent")

    assert (
        value(stream.diffus_phase_comp["Liq", "DCP"])
        == m.fs.properties.config.diffusivity_data["Liq", "DCP"]
    )
    assert (
        value(stream.diffus_phase_comp["Vap", "DCP"])
        == m.fs.properties.config.diffusivity_data["Vap", "DCP"]
    )
    assert (
        value(stream.molar_volume_comp["DCP"])
        == m.fs.properties.config.molar_volume_data["DCP"]
    )
    assert (
        value(stream.henry_constant_comp["DCP"])
        == m.fs.properties.config.henry_constant_data["DCP"]
    )

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325,
        "temperature": {"Liq": 283, "Vap": 283},
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): 99.97,
            ("Liq", "DCP"): 0.0001,
            ("Vap", "Air"): 7.482,
            ("Vap", "DCP"): 0.001,
        },
        "conc_mass_phase_comp": {
            ("Liq", "H2O"): 999.6990000010002,
            ("Liq", "DCP"): 0.0009999989997009106,
            ("Vap", "Air"): 1.2468333556060405,
            ("Vap", "DCP"): 0.00016664439395964187,
        },
        "dens_mass_phase": {"Liq": 999.7, "Vap": 1.247},
        "mass_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999989997009104,
            ("Liq", "DCP"): 1.0002990894277389e-06,
            ("Vap", "Air"): 0.9998663637578511,
            ("Vap", "DCP"): 0.00013363624214887077,
        },
        "flow_mole_phase_comp": {
            ("Liq", "H2O"): 5553.888888888889,
            ("Liq", "DCP"): 0.0008851124092759781,
            ("Vap", "Air"): 258.0,
            ("Vap", "DCP"): 0.008851124092759781,
        },
        "conc_mole_phase_comp": {
            ("Liq", "H2O"): 55538.8333333889,
            ("Liq", "DCP"): 0.008851115238988409,
            ("Vap", "Air"): 42.9942536415876,
            ("Vap", "DCP"): 0.0014749902102995384,
        },
        "mole_frac_phase_comp": {
            ("Liq", "H2O"): 0.9999998406319813,
            ("Liq", "DCP"): 1.593680186847315e-07,
            ("Vap", "Air"): 0.999965694494378,
            ("Vap", "DCP"): 3.4305505621986274e-05,
        },
        "diffus_phase_comp": {("Liq", "DCP"): 7.21045e-10, ("Vap", "DCP"): 8.57928e-06},
        "flow_vol_phase": {"Liq": 0.100000100030009, "Vap": 6.000801924619085},
        "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.483},
    }

    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)
