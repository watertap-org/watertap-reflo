import os
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    value,
    Var,
    Param,
    Constraint,
    Set,
    Var,
    Block,
    SolverFactory,
    TransformationFactory,
    assert_optimal_termination,
    check_optimal_termination,
    log,
    log10,
    units as pyunits,
)
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver

from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *

from idaes.core import UnitModelCostingBlock

from idaes.core.util.testing import initialization_tester

# from watertap.core.zero_order_costing import ZeroOrderCosting, _get_tech_parameters
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
import json
from os.path import join, dirname
from math import floor, ceil

import pytest

from io import StringIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from IPython.display import clear_output
from copy import deepcopy
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

# from watertap_contrib.seto.costing import SETOZeroOrderCosting, SETOWaterTAPCosting
from watertap_contrib.reflo.costing import TreatmentCosting

# from watertap_contrib.seto.solar_models.zero_order import PhotovoltaicZO
# from watertap_contrib.seto.energy import solar_energy
# from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
# from watertap_contrib.seto.unit_models import ChemicalSoftening0D
# from watertap_contrib.seto.property_models.chemical_softening_prop_pack import (
#     ChemSofteningParameterBlock,
# )
from watertap_contrib.reflo.kurby import print_unit_solutions, make_test_dict
from watertap_contrib.reflo.unit_models import SolarStill
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
from watertap.core.solvers import get_solver

# from watertap_contrib.seto.costing.solar import photovoltaic

from pyomo.util.check_units import assert_units_consistent
from watertap_contrib.reflo.unit_models.util.water_yield_calculation import (
    get_solar_still_daily_water_yield,
)
from watertap_contrib.reflo.analysis.case_studies.permian.utils import *

solver = get_solver()

skips = [
    "diffus_phase",
    "diffus_param",
    "dens_mass_param",
    "dh_vap_w_param",
    "cp_phase_param",
    "pressure_sat_param_psatw",
    "enth_mass_param",
    "osm_coeff_param",
    "visc_d_param",
    "therm_cond_phase_param",
    "pressure_sat_param",
    "bpe_",
    "TIC",
    "TPEC",
    "blocks[",
    "yearly_heat_production",
    "yearly_electricity_production",
    "cp_param_NaCl_liq",
    "_translator",
    "permeate_side",
    "properties_interface",
    "material_flow_dx",
    "._flow_terms",
    "pressure_dx",
    "MCAS_properties",
    "cp_param_NaCl_solid",
    "cp_vap_param",
    "temp_sat_solvent",
    "m_mec.fs.properties",
    "cp_mass_phase",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
weather_file = os.path.join(__location__, "data/carlsbad_NM_weather_tmy-2023-full.csv")

rho = 1000 * pyunits.kg / pyunits.m**3
electricity_cost_base = 0.0434618999  # USD_2018/kWh equivalent to 0.0575 USD_2023/kWh
heat_cost_base = 0.00894
save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/SS_results"


def build_ss(initial_salinity=130, daily_water_production=14000):

    daily_water_production = daily_water_production * pyunits.gallons / pyunits.day
    flow_mass_in = pyunits.convert(
        daily_water_production * rho, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_tds = pyunits.convert(
        daily_water_production * initial_salinity * pyunits.g / pyunits.liter,
        to_units=pyunits.kg / pyunits.s,
    )

    inlet_dict = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 31.4038218e-3},
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    water_yield_calc_dict = dict(
        input_weather_file_path=weather_file,
        initial_salinity=initial_salinity,  # initial salinity of influent water; g/L
        initial_water_depth=0.01,  # initial depth of water in solar still basin; m
        length_basin=0.6,  # length of each side of basin (length=width); m
        irradiance_threshold=20,  # irradiance values < threshold assumed to have negligible impact on calculation; W/m2
        temperature_col="Temperature",
        wind_velocity_col="Wind Speed",
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**inlet_dict)

    m.fs.unit = SolarStill(
        property_package=m.fs.properties,
        water_yield_calculation_args=water_yield_calc_dict,
    )

    m.fs.unit.properties_in[0].flow_vol_phase[...]
    m.fs.unit.properties_out[0].flow_vol_phase[...]
    m.fs.unit.properties_in[0].conc_mass_phase_comp[...]
    # m.fs.unit.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_in)
    # m.fs.unit.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix(flow_mass_tds)
    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): daily_water_production,
            ("conc_mass_phase_comp", ("Liq", "TDS")): initial_salinity,
            ("temperature", None): 298.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )
    # m.fs.unit.properties_in[0].pressure.fix(101325)
    # m.fs.unit.properties_in[0].temperature.fix(293.15)

    m.fs.costing = TreatmentCosting()
    m.fs.costing.electricity_cost.fix(electricity_cost_base)
    m.fs.costing.heat_cost.fix(heat_cost_base)
    m.fs.costing.maintenance_labor_chemical_factor.fix(0) # unit model has own relationships/factor for this

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(flow_rate=m.fs.unit.properties_in[0].flow_vol_phase["Liq"])

    return m

def run_permian_ss_water_production_tds_sweep():

    m0 = build_ss()
    rd = build_results_dict(m0, skips=skips)

    water_prod = [7000, 10500, 14000, 17500, 21000]  # gallon per day
    salinity = [100, 120, 130, 165, 200]

    df = pd.DataFrame()

    for wp in water_prod:
        for salt in salinity:
            rd = build_results_dict(m0, skips=skips)
            m = build_ss(initial_salinity=salt, daily_water_production=wp)
            m.fs.unit.initialize()
            results = solver.solve(m)
            assert_optimal_termination(results)
            rd = results_dict_append(m, rd)
            df = pd.concat([df, pd.DataFrame.from_dict(rd)])
            print(f"\nSalinity: {salt} g/L")
            print(f"Water production: {wp} gal/day")
            print(f"LCOW: {m.fs.costing.LCOW()} $/m3\n")

    df.to_csv(f"{save_path}/permian_SS_salinity_water_production_sweep.csv", index=False)

def run_permian_ss_water_production_sweep():

    m = build_ss()
    rd = build_results_dict(m, skips=skips)

    water_prod = [7000, 10500, 14000, 17500, 21000]  # gallon per day

    for wp in water_prod:
        m = build_ss(daily_water_production=wp)
        m.fs.unit.initialize()
        results = solver.solve(m)
        assert_optimal_termination(results)
        rd = results_dict_append(m, rd)
        print(f"Water production: {wp} gal/day")
        print(f"LCOW: {m.fs.costing.LCOW()} $/m3")

    df = pd.DataFrame.from_dict(rd)

    df.to_csv(f"{save_path}/permian_SS_water_production_sweep.csv", index=False)


if __name__ == "__main__":

    run_permian_ss_water_production_tds_sweep()

    # m = build_ss()

    # m.fs.unit.initialize()
    # results = solver.solve(m)
    # assert_optimal_termination(results)
    # m.fs.costing.LCOW.display()
    # m.fs.unit.properties_in[0].conc_mass_phase_comp.display()
