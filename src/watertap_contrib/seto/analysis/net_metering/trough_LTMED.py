import os
import sys
from io import StringIO
from os.path import join, dirname, getsize
import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    Var,
    value,
    units as pyunits,
)
from idaes.core.util.scaling import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.seto.costing import (
    EnergyCosting,
    SETOSystemCosting,
    TreatmentCosting,
)
from idaes.core.solvers import get_solver
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap_contrib.seto.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.seto.unit_models.surrogate import LTMEDSurrogate

solver = get_solver()


def build_trough_lt_med(
    salinity=35,
    recovery=0.5,
    feed_temp=25,
    steam_temp=80,
    capacity=80000,
    hours_storage=6
):
    feed_dens = 1000 * pyunits.kg / pyunits.m**3  # g/L = kg/m3

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.water_prop = SeawaterParameterBlock()
    m.fs.steam_prop = WaterParameterBlock()
    m.fs.treatment = Block()
    m.fs.treatment.lt_med = LTMEDSurrogate(
        property_package_water=m.fs.water_prop,
        property_package_steam=m.fs.steam_prop,
    )
    lt_med = m.fs.treatment.lt_med
    feed = lt_med.feed_props[0]
    dist = lt_med.distillate_props[0]
    steam = lt_med.steam_props[0]

    # LT-MED system specification
    feed_salinity = salinity * pyunits.kg / pyunits.m**3  # g/L = kg/m3

    sys_capacity = capacity * pyunits.m**3 / pyunits.day  # m3/day
    recovery_ratio = recovery * pyunits.dimensionless  # dimensionless
    feed_flow = pyunits.convert(
        (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
    )

    feed.flow_mass_phase_comp["Liq", "TDS"].fix(feed_salinity * feed_flow)
    feed.flow_mass_phase_comp["Liq", "H2O"].fix(feed_dens * feed_flow)
    feed.temperature.fix(feed_temp + 273.15)

    steam.temperature.fix(steam_temp + 273.15)
    steam.flow_mass_phase_comp["Liq", "H2O"].fix(0)

    dist.flow_mass_phase_comp["Liq", "TDS"].fix(0)  # salinity in distillate is zero

    lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)

    #  Add LT-MED costing
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.lt_med.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )
    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(dist.flow_vol_phase["Liq"])

    # Build energy block
    m.fs.energy = Block()
    m.fs.energy.trough = trough = TroughSurrogate()
    
    # Add trough costing
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.trough.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.costing.cost_process()

    # Add system costing 
    m.fs.sys_costing = SETOSystemCosting()
    m.fs.sys_costing.add_LCOW(dist.flow_vol_phase["Liq"])
    trough.hours_storage.fix(hours_storage)

    # Linking constraint
    m.fs.energy.lt_med_heat_demand_constr = Constraint(
        expr=trough.heat_load
        == pyunits.convert(lt_med.thermal_power_requirement, to_units=pyunits.MW)
    )

    m.fs.water_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.water_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "TDS")
    )
    m.fs.steam_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.steam_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Vap", "H2O")
    )
    set_scaling_factor(lt_med.thermal_power_requirement, 1e-6)
    constraint_scaling_transform(lt_med.eq_specific_thermal_energy_consumption, 1e6)
    calculate_scaling_factors(m)

    lt_med.initialize()


    return m

for salinity in [30, 45, 60]:
    print(f'SALINITY = {salinity}')
    m = build_trough_lt_med(salinity=salinity)
    trough = m.fs.energy.trough
    lt_med = m.fs.treatment.lt_med
    m.fs.energy.costing.heat_cost.set_value(0)
    m.fs.treatment.costing.heat_cost.set_value(0)

    m.fs.sys_costing.obj = Objective(expr=m.fs.sys_costing.LCOW)

    constraint_scaling_transform(
        m.fs.energy.trough.surrogate_blk.pysmo_constraint["heat_annual"], 1e-3
    )
    constraint_scaling_transform(
        m.fs.energy.trough.surrogate_blk.pysmo_constraint["electricity_annual"], 1e-3
    )

    m.fs.treatment.lt_med.distillate_props[0].flow_vol_phase["Liq"].fix()
    print(f"DOF = {degrees_of_freedom(m)}")

    solver = get_solver()
    solver.options["max_iter"] = 10000
    solver.options["halt_on_ampl_error"] = "yes"
    results = solver.solve(
        m,
        symbolic_solver_labels=True,
        tee=False,
    )
    print(f"SOLVE - {results.solver.termination_condition}")

    print(f"{results.solver.termination_condition}")

    results = solver.solve(m)
    print(m.fs.treatment.lt_med.thermal_power_requirement())
    print(trough.heat_load())
    print(f"DOF = {degrees_of_freedom(m)}")
    print(f"SOLVE - {results.solver.termination_condition}")
    print(m.fs.sys_costing.LCOW())
