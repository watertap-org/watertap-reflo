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
from idaes.core.util.exceptions import ConfigurationError

from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.property_models.air_water_equilibrium_properties import (
    AirWaterEq,
    AirWaterEqStateBlock,
    MolarVolumeCalculation,
    LiqDiffusivityCalculation,
    VapDiffusivityCalculation,
    SaturationVaporPressureCalculation,
    VaporPressureCalculation,
    RelativeHumidityCalculation,
    LatentHeatVaporizationCalculation,
    SpecificHeatWaterCalculation,
)

solver = get_solver()


@pytest.fixture(scope="module")
def m1():
    """
    Test NMSU case study 1 results
    """
    props = {
        "volatile_solute_list": ["TCA"],
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
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.ArdenBuck
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
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
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(0.10)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)
    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "TCA"): 2.61e-05,
            ("Liq", "H2O"): 157.8657,
            ("Vap", "TCA"): 0.0,
            ("Vap", "H2O"): 0.1,
            ("Vap", "Air"): 1.34932,
        },
        "conc_mass_phase_comp": {
            ("Liq", "TCA"): 0.00016518,
            ("Liq", "H2O"): 999.1498,
            ("Vap", "TCA"): 1.3e-14,
            ("Vap", "H2O"): 0.084177,
            ("Vap", "Air"): 1.135822,
        },
        "mass_frac_phase_comp": {
            ("Liq", "TCA"): 1.653e-07,
            ("Liq", "H2O"): 0.99999983,
            ("Vap", "TCA"): 1.1e-14,
            ("Vap", "H2O"): 0.068997,
            ("Vap", "Air"): 0.931,
        },
        "flow_mol_phase_comp": {
            ("Liq", "TCA"): 0.0001956,
            ("Liq", "H2O"): 8770.3,
            ("Vap", "TCA"): 1.12e-13,
            ("Vap", "H2O"): 5.5555,
            ("Vap", "Air"): 46.5,
        },
        "conc_mol_phase_comp": {
            ("Liq", "TCA"): 0.001238,
            ("Liq", "H2O"): 55508.3,
            ("Vap", "TCA"): 0,
            ("Vap", "H2O"): 4.67652,
            ("Vap", "Air"): 39.17,
        },
        "mole_frac_phase_comp": {
            ("Liq", "TCA"): 2.2308e-08,
            ("Liq", "H2O"): 0.99999,
            ("Vap", "TCA"): 2e-15,
            ("Vap", "H2O"): 0.10666,
        },
        "diffus_phase_comp": {("Vap", "TCA"): 7.6736e-06, ("Liq", "TCA"): 7.09255e-10},
        "molar_volume_comp": {"TCA": 0.00011007},
        "flow_vol_phase": {"Liq": 0.15800, "Vap": 1.187},
        "flow_mass_phase": {"Liq": 157.8, "Vap": 1.44932},
        "henry_comp": {"TCA": 0.3923},
        "pressure_vap_sat": {"H2O": 1215.6},
        "pressure_vap": {"H2O": 0.0},
        "dh_vap_mass_solvent": 2477.6,
        "cp_mass_solvent": {"Liq": 4197.1, "Vap": 1861.6},
        "collision_molecular_separation": {"TCA": 0.4683},
        "collision_molecular_separation_comp": {"TCA": 0.5655},
        "collision_function_comp": {"TCA": 0.5909},
        "collision_function_zeta_comp": {"TCA": -0.228462},
        "collision_function_ee_comp": {"TCA": 0.19251},
        "energy_molecular_attraction": {"TCA": 2.5e-14},
        "energy_molecular_attraction_air": 1e-14,
        "energy_molecular_attraction_comp": {"TCA": 5.7e-14},
        "arden_buck_exponential_term": 0.6875,
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
        "volatile_solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.146},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 0.00011007},
        "critical_molar_volume_data": {"DCP": 2.26e-4},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
        "relative_humidity_data": 0.5,
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

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 2
    for j in m.fs.properties.liq_comps:
        assert j in ["DCP", "H2O"]

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 2
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j in ["DCP", "H2O"]
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 3
    for j in m.fs.properties.vap_comps:
        assert j in ["DCP", "H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 3
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["DCP", "H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 1
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )
    for j in m.fs.properties.volatile_comps:
        assert j == "DCP"
        assert j == m.fs.properties.config.volatile_solute_list[0]

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) + 1
    )
    for p, j in m.fs.properties.volatile_comp_set:
        assert j == "DCP"
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
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
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.ArdenBuck
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
    )


@pytest.mark.component
def test_properties2(m2):
    m = m2
    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]

    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]
    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)
    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "DCP"): 0.0001,
            ("Liq", "H2O"): 99.97,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.001,
            ("Vap", "Air"): 7.482,
        },
        "conc_mass_phase_comp": {
            ("Liq", "DCP"): 0.000999,
            ("Liq", "H2O"): 999.699,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001666,
            ("Vap", "Air"): 1.24683335,
        },
        "mass_frac_phase_comp": {
            ("Liq", "DCP"): 1.000299e-06,
            ("Liq", "H2O"): 0.99999899,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.00013363,
            ("Vap", "Air"): 0.99986636,
        },
        "flow_mol_phase_comp": {
            ("Liq", "DCP"): 0.0008851,
            ("Liq", "H2O"): 5553.8,
            ("Vap", "DCP"): 2e-12,
            ("Vap", "H2O"): 0.05555,
            ("Vap", "Air"): 258.0,
        },
        "conc_mol_phase_comp": {
            ("Liq", "DCP"): 0.0088511,
            ("Liq", "H2O"): 55538.8,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0092580,
            ("Vap", "Air"): 42.9942,
        },
        "mole_frac_phase_comp": {
            ("Liq", "DCP"): 1.59368e-07,
            ("Liq", "H2O"): 0.9999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0002152,
            ("Vap", "Air"): 0.99978471,
        },
        "diffus_phase_comp": {("Vap", "DCP"): 8.579253e-06, ("Liq", "DCP"): 7.21e-10},
        "molar_volume_comp": {"DCP": 8.3550e-05},
        "flow_vol_phase": {"Liq": 0.1000, "Vap": 6.000},
        "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.483},
        "henry_comp": {"DCP": 0.0750617},
        "pressure_vap_sat": {"H2O": 1215.5768},
        "pressure_vap": {"H2O": 607.7884},
        "dh_vap_mass_solvent": 2477.6833,
        "cp_mass_solvent": {"Liq": 4197.1362, "Vap": 1861.6326},
        "collision_molecular_separation": {"DCP": 0.44348019},
        "collision_molecular_separation_comp": {"DCP": 0.51586},
        "collision_function_comp": {"DCP": 0.59831308},
        "collision_function_zeta_comp": {"DCP": -0.2230715},
        "collision_function_ee_comp": {"DCP": 0.17911045},
        "energy_molecular_attraction": {"DCP": 2.586778e-14},
        "energy_molecular_attraction_comp": {"DCP": 6.1661e-14},
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
        "volatile_solute_list": ["DCP"],
        "mw_data": {"DCP": 0.11298},
        "diffusivity_data": {("Liq", "DCP"): 7.21045e-10, ("Vap", "DCP"): 8.57928e-6},
        "dynamic_viscosity_data": {"Liq": 0.001307, "Vap": 1.79e-5},
        "henry_constant_data": {"DCP": 0.0525},  # salinity adjusted
        "standard_enthalpy_change_data": {"DCP": 31.1e3},
        "temperature_boiling_data": {"DCP": 369.1},
        "molar_volume_data": {"DCP": 8.355e-5},
        "density_data": {"Liq": 999.7, "Vap": 1.247},
        "temp_adjust_henry": False,
        "pressure_vap_sat_data": 1215.5768,
        "pressure_vap_data": 607.7884,
        "latent_heat_of_vaporization_data": 2477.6833,
        "specific_heat_of_water_data": {"Liq": 4197.1362, "Vap": 1861.6326},
        "relative_humidity_data": 0.5,
        "molar_volume_calculation": "none",
        "liq_diffus_calculation": "none",
        "vap_diffus_calculation": "none",
        "saturation_vapor_pressure_calculation": "none",
        "vapor_pressure_calculation": "none",
        "latent_heat_of_vaporization_calculation": "none",
        "specific_heat_of_water_calculation": "none",
        "relative_humidity_calculation": "none",
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

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 2
    for j in m.fs.properties.liq_comps:
        assert j in ["DCP", "H2O"]

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 2
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j in ["DCP", "H2O"]
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 3
    for j in m.fs.properties.vap_comps:
        assert j in ["DCP", "H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 3
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["DCP", "H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 1
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )
    for j in m.fs.properties.volatile_comps:
        assert j == "DCP"
        assert j == m.fs.properties.config.volatile_solute_list[0]

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 2
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) + 1
    )
    for p, j in m.fs.properties.volatile_comp_set:
        assert j == "DCP"
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3
    assert m.fs.properties.mw_comp["DCP"].value == 112.98e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
    assert isinstance(m.fs.properties.temperature_boiling_comp, Var)

    assert not m.fs.properties.config.temp_adjust_henry

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
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.none
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.none
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.none
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.none
    )


@pytest.mark.component
def test_properties3(m3):
    m = m3

    m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
    stream = m.fs.stream[0]
    stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
    stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
    stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
    stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
    stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
    stream.temperature["Liq"].fix(283)
    stream.temperature["Vap"].fix(283)
    stream.pressure.fix(101325)

    stream.conc_mass_phase_comp[...]

    stream.flow_mol_phase_comp[...]
    stream.conc_mol_phase_comp[...]
    stream.mole_frac_phase_comp[...]

    stream.diffus_phase_comp[...]
    stream.flow_vol_phase[...]
    stream.flow_mass_phase[...]
    stream.henry_comp[...]

    stream.pressure_vap_sat[...]
    stream.pressure_vap[...]
    stream.relative_humidity[...]
    stream.dh_vap_mass_solvent[...]
    stream.cp_mass_solvent[...]

    stream_vars = [
        "pressure",
        "temperature",
        "flow_mass_phase_comp",
        "conc_mass_phase_comp",
        "mass_frac_phase_comp",
        "flow_mol_phase_comp",
        "conc_mol_phase_comp",
        "mole_frac_phase_comp",
        "diffus_phase_comp",
        "molar_volume_comp",
        "flow_vol_phase",
        "flow_mass_phase",
        "henry_comp",
        "pressure_vap_sat",
        "pressure_vap",
        "dh_vap_mass_solvent",
        "cp_mass_solvent",
    ]

    for prop in stream_vars:
        assert hasattr(stream, prop)
        assert isinstance(getattr(stream, prop), Var)

    calculate_scaling_factors(m)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)
    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    stream_results = {
        "pressure": 101325.0,
        "temperature": {"Liq": 283.0, "Vap": 283.0},
        "flow_mass_phase_comp": {
            ("Liq", "DCP"): 0.0001,
            ("Liq", "H2O"): 99.97,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.001,
            ("Vap", "Air"): 7.482,
        },
        "conc_mass_phase_comp": {
            ("Liq", "DCP"): 0.000999,
            ("Liq", "H2O"): 999.7,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001666,
            ("Vap", "Air"): 1.2468,
        },
        "mass_frac_phase_comp": {
            ("Liq", "DCP"): 1.0003e-06,
            ("Liq", "H2O"): 0.99999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0001336,
            ("Vap", "Air"): 0.99986,
        },
        "flow_mol_phase_comp": {
            ("Liq", "DCP"): 0.0008851,
            ("Liq", "H2O"): 5553.8,
            ("Vap", "DCP"): 2e-12,
            ("Vap", "H2O"): 0.055555,
            ("Vap", "Air"): 258.0,
        },
        "conc_mol_phase_comp": {
            ("Liq", "DCP"): 0.008851,
            ("Liq", "H2O"): 55538.8,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.009258,
            ("Vap", "Air"): 42.9942,
        },
        "mole_frac_phase_comp": {
            ("Liq", "DCP"): 1.59368e-07,
            ("Liq", "H2O"): 0.9999,
            ("Vap", "DCP"): 0.0,
            ("Vap", "H2O"): 0.0002152,
            ("Vap", "Air"): 0.999784,
        },
        "diffus_phase_comp": {("Vap", "DCP"): 8.5793e-06, ("Liq", "DCP"): 7.21e-10},
        "flow_vol_phase": {"Liq": 0.1000001, "Vap": 6.000},
        "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.483},
    }
    for prop, d in stream_results.items():
        sv = getattr(stream, prop)
        if isinstance(d, dict):
            for i, val in d.items():
                assert value(sv[i]) == pytest.approx(val, rel=1e-3)
        else:
            assert value(sv) == pytest.approx(d, rel=1e-3)


@pytest.fixture(scope="module")
def m_aw():
    """
    Test air-water ONLY system
    """
    props = dict()
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)

    return m


@pytest.mark.unit
def test_air_water_only_config(m_aw):
    pass


@pytest.mark.unit
def test_parameter_block4(m_aw):
    m = m_aw

    assert len(m.fs.properties.config.volatile_solute_list) == 0
    assert len(m.fs.properties.config.non_volatile_solute_list) == 0
    assert isinstance(m.fs.properties.component_list, Set)
    assert len(m.fs.properties.component_list) == 2
    for j in m.fs.properties.component_list:
        assert j in ["H2O", "Air"]

    assert isinstance(m.fs.properties.phase_list, Set)
    assert len(m.fs.properties.phase_list) == 2
    assert m.fs.properties.phase_list == ["Liq", "Vap"]

    ### Liquid components, Sets
    assert isinstance(m.fs.properties.liq_comps, Set)
    assert len(m.fs.properties.liq_comps) == 1
    for j in m.fs.properties.liq_comps:
        assert j == "H2O"

    assert isinstance(m.fs.properties.liq_comp_set, Set)
    assert len(m.fs.properties.liq_comp_set) == 1
    assert len(m.fs.properties.liq_comp_set) == len(m.fs.properties.liq_comps)
    for p, j in m.fs.properties.liq_comp_set:
        assert j == "H2O"
        assert p == "Liq"

    assert isinstance(m.fs.properties.non_volatile_comps, Set)
    assert len(m.fs.properties.non_volatile_comps) == 0
    assert len(m.fs.properties.non_volatile_comps) == len(
        m.fs.properties.config.non_volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 0
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    ### Vapor components, Sets
    assert isinstance(m.fs.properties.vap_comps, Set)
    assert len(m.fs.properties.vap_comps) == 2
    for j in m.fs.properties.vap_comps:
        assert j in ["H2O", "Air"]

    assert isinstance(m.fs.properties.vap_comp_set, Set)
    assert len(m.fs.properties.vap_comp_set) == 2
    assert len(m.fs.properties.vap_comp_set) == len(m.fs.properties.vap_comps)
    for p, j in m.fs.properties.vap_comp_set:
        assert j in ["H2O", "Air"]
        assert p in ["Liq", "Vap"]

    assert isinstance(m.fs.properties.volatile_comps, Set)
    assert len(m.fs.properties.volatile_comps) == 0
    assert len(m.fs.properties.volatile_comps) == len(
        m.fs.properties.config.volatile_solute_list
    )

    assert isinstance(m.fs.properties.volatile_comp_set, Set)
    assert len(m.fs.properties.volatile_comp_set) == 0
    assert (
        len(m.fs.properties.volatile_comp_set)
        == len(m.fs.properties.volatile_comps) * 2
    )

    assert isinstance(m.fs.properties.component_set, Set)
    assert len(m.fs.properties.component_set) == len(
        m.fs.properties.config.non_volatile_solute_list
    ) + len(m.fs.properties.config.volatile_solute_list) + len(["H2O", "Air"])
    assert len(m.fs.properties.component_set) == len(m.fs.properties.solute_set)

    assert m.fs.properties._state_block_class is AirWaterEqStateBlock
    assert isinstance(m.fs.properties.mw_comp, Param)
    assert m.fs.properties.mw_comp["H2O"].value == 18e-3
    assert m.fs.properties.mw_comp["Air"].value == 29e-3

    assert isinstance(m.fs.properties.visc_d_phase, Var)
    assert isinstance(m.fs.properties.molar_volume_comp, Var)
    assert isinstance(m.fs.properties.molar_volume_comp_crit, Var)
    assert isinstance(m.fs.properties.henry_comp, Var)
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
    assert (
        m.fs.properties.config.saturation_vapor_pressure_calculation
        == SaturationVaporPressureCalculation.ArdenBuck
    )
    assert (
        m.fs.properties.config.vapor_pressure_calculation
        == VaporPressureCalculation.FromRelativeHumidity
    )
    assert (
        m.fs.properties.config.relative_humidity_calculation
        == RelativeHumidityCalculation.none
    )
    assert (
        m.fs.properties.config.latent_heat_of_vaporization_calculation
        == LatentHeatVaporizationCalculation.Sharqawy
    )
    assert (
        m.fs.properties.config.specific_heat_of_water_calculation
        == SpecificHeatWaterCalculation.Sharqawy
    )


# @pytest.mark.component
# def test_properties4(m_aw):
#     m = m_aw

#     m.fs.stream = m.fs.properties.build_state_block([0], defined_state=True)
#     stream = m.fs.stream[0]
#     stream.flow_mass_phase_comp["Liq", "H2O"].fix(99.97)
#     stream.flow_mass_phase_comp["Liq", "DCP"].fix(1e-4)
#     stream.flow_mass_phase_comp["Vap", "H2O"].fix(1e-3)
#     stream.flow_mass_phase_comp["Vap", "Air"].fix(7.482)
#     stream.flow_mass_phase_comp["Vap", "DCP"].fix(0)
#     stream.temperature["Liq"].fix(283)
#     stream.temperature["Vap"].fix(283)
#     stream.pressure.fix(101325)

#     stream.conc_mass_phase_comp[...]

#     stream.flow_mol_phase_comp[...]
#     stream.conc_mol_phase_comp[...]
#     stream.mole_frac_phase_comp[...]

#     stream.diffus_phase_comp[...]
#     stream.flow_vol_phase[...]
#     stream.flow_mass_phase[...]
#     stream.henry_comp[...]

#     stream.pressure_vap_sat[...]
#     stream.pressure_vap[...]
#     stream.relative_humidity[...]
#     stream.dh_vap_mass_solvent[...]
#     stream.cp_mass_solvent[...]

#     stream_vars = [
#         "pressure",
#         "temperature",
#         "flow_mass_phase_comp",
#         "conc_mass_phase_comp",
#         "mass_frac_phase_comp",
#         "flow_mol_phase_comp",
#         "conc_mol_phase_comp",
#         "mole_frac_phase_comp",
#         "diffus_phase_comp",
#         "molar_volume_comp",
#         "flow_vol_phase",
#         "flow_mass_phase",
#         "henry_comp",
#         "pressure_vap_sat",
#         "pressure_vap",
#         "dh_vap_mass_solvent",
#         "cp_mass_solvent",
#     ]

#     for prop in stream_vars:
#         assert hasattr(stream, prop)
#         assert isinstance(getattr(stream, prop), Var)

#     calculate_scaling_factors(m)

#     assert_units_consistent(m)

#     check_dof(m, fail_flag=True)
#     m.fs.stream.initialize()

#     results = solver.solve(m)
#     assert_optimal_termination(results)

#     stream_results = {
#         "pressure": 101325.0,
#         "temperature": {"Liq": 283.0, "Vap": 283.0},
#         "flow_mass_phase_comp": {
#             ("Liq", "DCP"): 0.0001,
#             ("Liq", "H2O"): 99.97,
#             ("Vap", "DCP"): 0.0,
#             ("Vap", "H2O"): 0.001,
#             ("Vap", "Air"): 7.482,
#         },
#         "conc_mass_phase_comp": {
#             ("Liq", "DCP"): 0.000999,
#             ("Liq", "H2O"): 999.7,
#             ("Vap", "DCP"): 0.0,
#             ("Vap", "H2O"): 0.0001666,
#             ("Vap", "Air"): 1.2468,
#         },
#         "mass_frac_phase_comp": {
#             ("Liq", "DCP"): 1.0003e-06,
#             ("Liq", "H2O"): 0.99999,
#             ("Vap", "DCP"): 0.0,
#             ("Vap", "H2O"): 0.0001336,
#             ("Vap", "Air"): 0.99986,
#         },
#         "flow_mol_phase_comp": {
#             ("Liq", "DCP"): 0.0008851,
#             ("Liq", "H2O"): 5553.8,
#             ("Vap", "DCP"): 2e-12,
#             ("Vap", "H2O"): 0.055555,
#             ("Vap", "Air"): 258.0,
#         },
#         "conc_mol_phase_comp": {
#             ("Liq", "DCP"): 0.008851,
#             ("Liq", "H2O"): 55538.8,
#             ("Vap", "DCP"): 0.0,
#             ("Vap", "H2O"): 0.009258,
#             ("Vap", "Air"): 42.9942,
#         },
#         "mole_frac_phase_comp": {
#             ("Liq", "DCP"): 1.59368e-07,
#             ("Liq", "H2O"): 0.9999,
#             ("Vap", "DCP"): 0.0,
#             ("Vap", "H2O"): 0.0002152,
#             ("Vap", "Air"): 0.999784,
#         },
#         "diffus_phase_comp": {("Vap", "DCP"): 8.5793e-06, ("Liq", "DCP"): 7.21e-10},
#         "flow_vol_phase": {"Liq": 0.1000001, "Vap": 6.000},
#         "flow_mass_phase": {"Liq": 99.9701, "Vap": 7.483},
#     }
#     for prop, d in stream_results.items():
#         sv = getattr(stream, prop)
#         if isinstance(d, dict):
#             for i, val in d.items():
#                 assert value(sv[i]) == pytest.approx(val, rel=1e-3)
#         else:
#             assert value(sv) == pytest.approx(d, rel=1e-3)
