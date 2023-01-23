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
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
import pyomo.util.infeasible as infeas
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import seaborn as sns
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DegeneracyHunter
import idaes.logger as idaeslog
from idaes.core.util.misc import StrEnum
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

# from watertap.property_models.ion_DSPMDE_prop_pack import (
#     DSPMDEParameterBlock,
#     DSPMDEStateBlock,
#     ActivityCoefficientModel,
#     DensityCalculation,
# )
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock, NaClStateBlock
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.unit_models.ion_exchange_0D import (
    IonExchange0D,
    IonExchangeType,
    RegenerantChem,
)

# from watertap.unit_models.zero_order import SolarEnergyZO
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.wt_database import Database
from watertap.core.zero_order_costing import ZeroOrderCosting, _get_tech_parameters
from watertap.core.util.infeasible import *
from watertap.costing import WaterTAPCosting
import watertap.core.zero_order_properties as prop_ZO

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

# import PySAM
# from pysam_seto import *
# print_infeasible_constraints(m)
# print_infeasible_bounds(m)
# print_variables_close_to_bounds(m)
# print_constraints_close_to_bounds(m)
# print_close_to_bounds(m)
# print(infeas.log_infeasible_bounds(m))
# print(infeas.log_infeasible_constraints(m))
# /Users/ksitterl/Documents/Python/watertap-seto/watertap-seto/src/watertap_contrib/seto/costing/seto_zero_order_costing.py
from watertap_contrib.seto.costing import SETOZeroOrderCosting, SETOWaterTAPCosting, TreatmentCosting, EnergyCosting, SETOSystemCosting

from watertap_contrib.seto.solar_models.zero_order import PhotovoltaicZO
from watertap_contrib.seto.energy import solar_energy
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
# from watertap_contrib.seto.costing.solar import photovoltaic
import os

solver = get_solver()

absolute_path = os.path.dirname(__file__)
print(absolute_path)

tech_config_file = "/pysam_data/pvsamv1.json"
tech_config_file = absolute_path + tech_config_file
grid_config_file = "/pysam_data/grid.json"
grid_config_file = absolute_path + grid_config_file
rate_config_file = "/pysam_data/utilityrate5.json"
rate_config_file = absolute_path + rate_config_file
cash_config_file = "/pysam_data/singleowner.json"
cash_config_file = absolute_path + cash_config_file
weather_file = "/pysam_data/phoenix_az_33.450495_-111.983688_psmv3_60_tmy.csv"
weather_file = absolute_path + weather_file


def build_ro_pv():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = SETODatabase()
    m.pysam = PySAMWaterTAP(
        pysam_model="pv",
        tech_config_file=tech_config_file,
        grid_config_file=grid_config_file,
        rate_config_file=rate_config_file,
        cash_config_file=cash_config_file,
        weather_file=weather_file,
    )
    m.fs.properties = NaClParameterBlock()

    treatment = m.fs.treatment = Block()
    energy = m.fs.energy = Block()

    energy.pv = PhotovoltaicZO(property_package=m.fs.properties, database=m.db)
    treatment.feed = Feed(property_package=m.fs.properties)
    treatment.product = Product(property_package=m.fs.properties)
    treatment.disposal = Product(property_package=m.fs.properties)

    treatment.p1 = Pump(property_package=m.fs.properties)

    treatment.ro = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )

    treatment.erd = EnergyRecoveryDevice(property_package=m.fs.properties)

    treatment.a1 = Arc(source=treatment.feed.outlet, destination=treatment.p1.inlet)
    treatment.a2 = Arc(source=treatment.p1.outlet, destination=treatment.ro.inlet)
    treatment.a3 = Arc(
        source=treatment.ro.permeate, destination=treatment.product.inlet
    )
    treatment.a4 = Arc(source=treatment.ro.retentate, destination=treatment.erd.inlet)
    treatment.a5 = Arc(
        source=treatment.erd.outlet, destination=treatment.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(treatment)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    set_scaling_factor(treatment.p1.control_volume.work, 1e-3)
    set_scaling_factor(treatment.ro.area, 1e-2)
    treatment.feed.properties[0].flow_vol_phase["Liq"]
    treatment.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    energy.pv.properties[0].flow_vol_phase["Liq"]
    energy.pv.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    set_scaling_factor(treatment.erd.control_volume.work, 1e-3)
    calculate_scaling_factors(m)

    return m




def set_operating_conditions(
    flow_in=1e-3, tds=35, ro_area_guess=50, water_recovery=0.5
):
    ro = m.fs.treatment.ro
    erd = m.fs.treatment.erd
    pv = m.fs.energy.pv
    p1 = m.fs.treatment.p1
    feed = m.fs.treatment.feed
    feed.properties[0].pressure.fix(101325)
    feed.properties[0].temperature.fix(298.15)
    feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): tds * 1e-3,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )
    pv.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): tds * 1e-3,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )
    p1.efficiency_pump.fix(0.8)
    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.treatment.feed.properties[0],
        solver=solver,
        over_pressure=0.25,
        water_recovery=water_recovery,
        NaCl_passage=0.01,
    )
    operating_pressure_psi = pyunits.convert(operating_pressure * pyunits.Pa, to_units=pyunits.psi)()
    operating_pressure_bar = pyunits.convert(operating_pressure * pyunits.Pa, to_units=pyunits.bar)()
    print(
        f"\nOPERATING PRESSURE ESTIMATE = {round(operating_pressure_bar, 2)} bar = {round(operating_pressure_psi, 2)} psi\n"
    )
    p1.control_volume.properties_out[0].pressure.fix(operating_pressure)
    m.db.get_unit_operation_parameters("solar_energy")
    pv.load_parameters_from_database(use_default_removal=True)

def initialize_treatment(ro_area_guess=50, water_recovery=0.5):
    ro = m.fs.treatment.ro
    p1 = m.fs.treatment.p1
    feed = m.fs.treatment.feed
    erd = m.fs.treatment.erd
    p1 = m.fs.treatment.p1

    ro.A_comp.fix(4.2e-12)
    ro.B_comp.fix(3.5e-8)
    ro.feed_side.channel_height.fix(1e-3)
    ro.feed_side.spacer_porosity.fix(0.97)
    ro.permeate.pressure[0].fix(101325)
    ro.width.fix(5)

    ro.feed_side.properties_in[0].flow_mass_phase_comp[
        "Liq", "H2O"
    ] = p1.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]()
    ro.feed_side.properties_in[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ] = p1.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"]()
    ro.feed_side.properties_in[0].temperature = feed.properties[0].temperature()
    ro.feed_side.properties_in[0].pressure = p1.control_volume.properties_out[
        0
    ].pressure()
    ro.area.fix(ro_area_guess)
    ro.initialize()
    ro.area.unfix()
    ro.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)
    propagate_state(m.fs.treatment.a4)

    erd.efficiency_pump.fix(0.95)
    erd.control_volume.properties_out[0].pressure.fix(101325)
    erd.initialize()

    propagate_state(m.fs.treatment.a1)
    # propagate_state(m.fs.treatment.a6)
    p1.initialize()
    propagate_state(m.fs.treatment.a2)


def initialize_energy():
    m.fs.energy.pv.initialize()


def initialize_sys(solver=None, ro_area_guess=50, water_recovery=0.5):
    if solver is None:
        solver = get_solver()
    optarg = solver.options
    m.fs.treatment.feed.initialize(optarg=optarg)
    initialize_treatment(ro_area_guess=ro_area_guess, water_recovery=water_recovery)
    initialize_energy()

def optimize_setup(
    press_lb=125,
    press_ub=1200,
    area_lb=1,
    area_ub=150,
    prod_salinity=500e-6,
    min_flux=2.5e-4,
):
    ro = m.fs.treatment.ro
    p1 = m.fs.treatment.p1

    press_lb_Pa = pyunits.convert(press_lb * pyunits.psi, to_units=pyunits.Pa)()
    press_ub_Pa = pyunits.convert(press_ub * pyunits.psi, to_units=pyunits.Pa)()
    # m.fs.pv.properties[0].flow_mass_phase_comp['Liq', 'H2O'].fix()
    p1.control_volume.properties_out[0].pressure.unfix()
    p1.control_volume.properties_out[0].pressure.setlb(press_lb_Pa)
    p1.control_volume.properties_out[0].pressure.setub(press_ub_Pa)
    p1.deltaP.setlb(0)

    ro.area.unfix()
    ro.area.setlb(area_lb)
    ro.area.setub(area_ub)

    m.fs.treatment.prod_salinity = Param(initialize=prod_salinity, mutable=True)
    m.fs.treatment.min_flux = Param(initialize=min_flux, mutable=True)

    m.fs.treatment.eq_prod_salinity = Constraint(
        expr=m.fs.treatment.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.treatment.prod_salinity
    )

    constraint_scaling_transform(m.fs.treatment.eq_prod_salinity, 1e4)

    m.fs.treatment.eq_min_flux = Constraint(
        expr=m.fs.treatment.ro.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.treatment.min_flux
    )

    assert_degrees_of_freedom(m.fs.treatment, 1)


def add_costing():
    treatment = m.fs.treatment
    energy = m.fs.energy
    treatment.costing = TreatmentCosting()
    energy.costing = EnergyCosting()
    m.db.get_unit_operation_parameters("photovoltaic")
    energy.pv.load_parameters_from_database(use_default_removal=True)

    energy.pv.costing = UnitModelCostingBlock(flowsheet_costing_block=energy.costing)
    treatment.ro.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )
    treatment.erd.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )
    treatment.p1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )

    treatment.costing.cost_process()
    energy.costing.cost_process()

    m.fs.sys_costing = SETOSystemCosting()
    flow_out = treatment.product.properties[0].flow_vol
    m.fs.sys_costing.add_LCOW(flow_out)
    m.fs.sys_costing.add_specific_electric_energy_consumption(flow_out)

    treatment.costing.initialize()
    energy.costing.initialize()
    # return m

def solve_it(solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")
    return results

def print_ro_results(m, sep="."):
    ro = m.fs.treatment.ro
    erd = m.fs.treatment.erd
    liq = "Liq"
    nacl = "NaCl"
    header = f'{"PARAM":<40s}{"VALUE":<40s}{"UNITS":<40s}'
    prop_in = ro.feed_side.properties_in[0]
    prop_out = ro.feed_side.properties_out[0]
    prop_perm = ro.mixed_permeate[0]
    pv_cost = m.fs.energy.costing.photovoltaic
    line = f'\n{f"{sep*160}":<160s}'
    flux_lmh = pyunits.convert(
        ro.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        / ro.feed_side.properties_in[0].dens_mass_phase["Liq"],
        to_units=pyunits.liter / pyunits.m**2 / pyunits.hr,
    )()
    flow_out_L = pyunits.convert(
        prop_perm.flow_vol, to_units=pyunits.liter / pyunits.hr
    )
    title = f'\n{"=======> INTEGRATED SYSTEM RESULTS <=======":^120}\n'
    print(line)
    print(title)
    print(header)
    if hasattr(m.fs, "sys_costing"):
        print(
            f'{"LCOW":<40s}{f"{m.fs.sys_costing.LCOW():<40.4f}"}{"$/m3":<40s}{f"--":<40}'
        )
        print(
            f'{"Total Capital Cost":<40s}{f"{m.fs.sys_costing.total_capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"Total Operating Cost":<40s}{f"{m.fs.sys_costing.total_operating_cost():<40.2f}"}{"$/yr":<40s}{f"--":<40}'
        )
        print(
            f'{"SEC":<40s}{f"{m.fs.sys_costing.specific_electric_energy_consumption():<40.4f}"}{"kWh/m3":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Avg. Electricity Generation":<40s}{f"{m.fs.energy.pv.electricity():<40.4f}"}{"kW":<40s}{f"--":<40}'
        )
        print(
            f'{"RO Electricity Consumption":<40s}{f"{m.fs.treatment.costing.aggregate_flow_electricity():<40.4f}"}{"kW":<40s}{f"--":<40}'
        )
        print(
            f'{"Overall Electricity Consumption":<40s}{f"{m.fs.sys_costing.aggregate_flow_electricity():<40.4f}"}{"kW":<40s}{f"--":<40}'
        )
        title = f'\n{"=======> PV SYSTEM RESULTS <=======":^120}\n'
        print(title)
        print(header)
        print(
            f'{"PV Capital Cost":<40s}{f"{m.fs.energy.pv.costing.capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Fixed Operating":<40s}{f"{m.fs.energy.pv.costing.fixed_operating_cost():<40.2f}"}{"$/yr":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Fixed Operating by Capacity":<40s}{f"{pv_cost.fixed_operating_by_capacity():<40.2f}"}{"$/yr":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Variable Operating":<40s}{f"{m.fs.energy.pv.costing.variable_operating_cost():<40.2f}"}{"$/yr":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Variable Operating by Annual Gen.":<40s}{f"{pv_cost.variable_operating_by_generation():<40.2f}"}{"$/MWh":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Annual Generation":<40s}{f"{m.fs.energy.pv.costing.annual_generation():<40.4f}"}{"MWh/yr":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Nameplate Capacity":<40s}{f"{m.fs.energy.pv.costing.system_capacity():<40.4f}"}{"W":<40s}{f"--":<40}'
        )
        print(
            f'{"PV Land Required":<40s}{f"{m.fs.energy.pv.costing.land_area():<40.4f}"}{"acres":<40s}{f"--":<40}'
        )
        # if hasattr(m, 'pysam'):
        #     max_gen = max(m.pysam.hourly_energy)
        # print(f'{"PV Max Generation":<40s}{f"{max_gen:<40.4f}"}{"kW":<40s}{f"--":<40}')
        print(
            f'{"PV Avg. Generation":<40s}{f"{-1 * m.fs.energy.pv.electricity():<40.4f}"}{"kW":<40s}{f"--":<40}'
        )
        title = f'\n{"=======> RO SYSTEM RESULTS <=======":^120}\n'
        print(title)
        print(header)
        print(
            f'{"RO Capital Cost":<40s}{f"{m.fs.treatment.ro.costing.capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"Pump Capital Cost":<40s}{f"{m.fs.treatment.p1.costing.capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"ERD Capital Cost":<40s}{f"{m.fs.treatment.erd.costing.capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"Total Treatment Capital Cost":<40s}{f"{m.fs.treatment.costing.total_capital_cost():<40.4f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"RO Fixed Operating Cost":<40s}{f"{m.fs.treatment.ro.costing.fixed_operating_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"Pumping Power":<40s}{f"{m.fs.treatment.p1.control_volume.work[0]():<40.2f}"}{"W":<40s}{f"--":<40}'
        )
        print(
            f'{"RO Capital Cost":<40s}{f"{m.fs.treatment.ro.costing.capital_cost():<40.2f}"}{"$":<40s}{f"--":<40}'
        )
        print(
            f'{"RO Operating Cost":<40s}{f"{m.fs.treatment.costing.total_operating_cost():<40.4f}"}{"$":<40s}{f"--":<40}'
        )

    print(
        f'{"RO Pressure":<40s}{f"{pyunits.convert(ro.inlet.pressure[0], to_units=pyunits.psi)():<40.4f}"}{"psi":<40s}{f"--":<40}'
    )
    print(f'{"Membrane Area":<40s}{f"{ro.area():<40.4f}"}{"m2":<40s}{f"--":<40}')
    print(f'{"Flux":<40s}{f"{flux_lmh:<40.4f}"}{"LMH":<40s}{f"--":<40}')
    # print(f'{"Flux Check":<40s}{f"{flow_out_L() / ro.area():<40.4f}"}{"LMH":<40s}{f"--":<40}')
    print(
        f'{"Vol. Recovery":<40s}{f"{100 * ro.recovery_vol_phase[0, liq]():<40.4f}"}{"%":<40s}{f"--":<40}'
    )
    print(f'{"Flow In":<40s}{f"{prop_in.flow_vol():<40.4f}"}{"m3/s":<40s}{f"--":<40}')
    print(
        f'{"Flow In [MGD]":<40s}{f"{pyunits.convert(prop_in.flow_vol, to_units=pyunits.Mgallons/pyunits.day)():<40.4f}"}{"MGD":<40s}{f"--":<40}'
    )
    print(
        f'{"Flow Out":<40s}{f"{prop_perm.flow_vol():<40.4f}"}{"m3/s":<40s}{f"--":<40}'
    )
    print(
        f'{"Conc. In":<40s}{f"{pyunits.convert(prop_in.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<40.4f}"}{"mg/L":<40s}{f"--":<40}'
    )
    print(
        f'{"Conc. Reject":<40s}{f"{pyunits.convert(prop_out.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<40.4f}"}{"mg/L":<40s}{f"--":<40}'
    )
    print(
        f'{"Conc. Perm":<40s}{f"{pyunits.convert(prop_perm.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<40.4f}"}{"mg/L":<40s}{f"--":<40}'
    )
    print(
        f'{"Pump Pressure":<40s}{f"{pyunits.convert(erd.inlet.pressure[0], to_units=pyunits.psi)():<40.4f}"}{"psi":<40s}{f"--":<40}'
    )
    print(
        f'{"ERD Pressure":<40s}{f"{pyunits.convert(m.fs.treatment.p1.outlet.pressure[0], to_units=pyunits.psi)():<40.4f}"}{"psi":<40s}{f"--":<40}'
    )
    print(
        f'{"ERD Power Recovered":<40s}{f"{-1 * erd.work_mechanical[0]() * 1e-3:<40.4f}"}{"kW":<40s}{f"--":<40}'
    )


# print_ro_results(m)

def fix_pv_costing():
    pvc = m.fs.energy.pv.costing
    pvc.system_capacity.fix(0)
    pvc.annual_generation.fix(0)
    pvc.land_area.fix(0)

def fix_treatment_global_params():
    tc = m.fs.treatment.costing
    tc.factor_total_investment.fix(1)
    # tc.wacc.set_value(0.06)
    tc.factor_maintenance_labor_chemical.fix(0)

flow_in = 4.38e-3  # m3/s
tds = 50  # g/L
ro_area_guess = 50
water_recovery = 0.5
press_lb = 125  # psi
press_ub = 5500  # psi
area_lb = 1  # m2
area_ub = 1200  # m2
prod_salinity = 200e-4
min_flux = 0.5e-6

oversize_factor = 1

m = build_ro_pv()
set_operating_conditions(
    flow_in=flow_in,
    tds=tds,
    ro_area_guess=ro_area_guess,
    water_recovery=water_recovery,
)
initialize_sys()
optimize_setup(
    press_lb=press_lb,
    press_ub=press_ub,
    area_lb=area_lb,
    area_ub=area_ub,
    prod_salinity=prod_salinity,
    min_flux=min_flux,
)
add_costing()

m.fs.obj = Objective(expr=m.fs.sys_costing.LCOW)

m.fs.energy.pv.oversize_factor = oversize_factor

fix_pv_costing()
fix_treatment_global_params()
m.fs.sys_costing.base_currency = pyunits.USD_2018
m.results = solve_it()
print_ro_results(m)
desired_pv_size = m.fs.treatment.costing.aggregate_flow_electricity() * m.fs.energy.pv.oversize_factor

cash_model_kwargs = {"om_fixed": 1e4, "om_production": 20}
tech_model_kwargs = {"subarray1_rear_soiling_loss": 1, "subarray1_rack_shading": 1}
m.pysam.run_pv_single_owner(
    desired_size=desired_pv_size, cash_model_kwargs=cash_model_kwargs
)