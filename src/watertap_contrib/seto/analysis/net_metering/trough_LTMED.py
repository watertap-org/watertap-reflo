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


def main():
    build_solve_trough()
    build_solve_lt_med()
    build_solve_trough_lt_med()


def check_scaling(m, jacob=False):
    if jacob:
        jac, nlp = get_jacobian(m, scaled=False)
        print("Extreme Jacobian entries:")
        for i in extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"derivative of: {i[1]}\nwith respect to: {i[2]}\n\t{i[0]:.2e}\n")
    print("Badly scaled variables:")
    for v, sv in badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-12):
        print(f"var: {v}\n\tscaled value: {sv}\n\tsf: {get_scaling_factor(v)}")


def build_solve_trough():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.energy = Block()
    m.fs.energy.trough = TroughSurrogate()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.trough.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )
    m.fs.energy.costing.cost_process()

    trough = m.fs.energy.trough
    trough.heat_load.fix(216)
    trough.hours_storage.fix(6)

    print(f"DOF = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(f"Trough capital: ${m.fs.energy.costing.total_capital_cost()}")
    print(f"Annual heat generated: {m.fs.energy.trough.heat_annual()} kWh")
    print(f"Annual electricity consumed: {m.fs.energy.trough.electricity_annual()} kWh")


def build_solve_lt_med():
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

    # System specification
    feed_salinity = 35 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    feed_dens = 1000 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    feed_temperature = 25  # degC
    steam_temperature = 80  # degC
    sys_capacity = 80000 * pyunits.m**3 / pyunits.day  # m3/day
    recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless
    feed_flow = pyunits.convert(
        (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
    )

    feed.flow_mass_phase_comp["Liq", "TDS"].fix(feed_salinity * feed_flow)
    feed.flow_mass_phase_comp["Liq", "H2O"].fix(feed_dens * feed_flow)
    feed.temperature.fix(feed_temperature + 273.15)
    steam.temperature.fix(steam_temperature + 273.15)
    # flow rate of liquid steam is zero
    steam.flow_mass_phase_comp["Liq", "H2O"].fix(0)
    dist.flow_mass_phase_comp["Liq", "TDS"].fix(0)  # salinity in distillate is zero

    lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.lt_med.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )
    m.fs.treatment.costing.cost_process()
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
    calculate_scaling_factors(m)
    # check_scaling(m)
    lt_med.initialize()
    check_scaling(m)
    results = solver.solve(m)
    print(f"DOF = {degrees_of_freedom(m)}")
    assert_optimal_termination(results)
    print(f"LT-MED capital: ${m.fs.treatment.costing.total_capital_cost()}")
    print(f"LT-MED thermal power requirement: {lt_med.thermal_power_requirement()} kW")


def build_solve_trough_lt_med():
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

    # System specification
    feed_salinity = 35 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    feed_dens = 1000 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    feed_temperature = 25  # degC
    steam_temperature = 80  # degC
    sys_capacity = 80000 * pyunits.m**3 / pyunits.day  # m3/day
    recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless
    feed_flow = pyunits.convert(
        (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
    )

    feed.flow_mass_phase_comp["Liq", "TDS"].fix(feed_salinity * feed_flow)
    feed.flow_mass_phase_comp["Liq", "H2O"].fix(feed_dens * feed_flow)
    feed.temperature.fix(feed_temperature + 273.15)
    steam.temperature.fix(steam_temperature + 273.15)
    # flow rate of liquid steam is zero
    steam.flow_mass_phase_comp["Liq", "H2O"].fix(0)
    dist.flow_mass_phase_comp["Liq", "TDS"].fix(0)  # salinity in distillate is zero

    lt_med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.lt_med.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )
    m.fs.treatment.costing.cost_process()

    m.fs.energy = Block()
    m.fs.energy.trough = TroughSurrogate()
    trough = m.fs.energy.trough
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.trough.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing
    )

    m.fs.energy.costing.cost_process()
    m.fs.sys_costing = SETOSystemCosting()

    m.fs.energy.lt_med_heat_demand_constr = Constraint(
        expr=trough.heat_load == lt_med.thermal_power_requirement
    )
    trough.hours_storage.fix(6)

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

    calculate_scaling_factors(m)
    lt_med.initialize()

    # check_scaling(m)
    results = solver.solve(m)
    print(f"DOF = {degrees_of_freedom(m)}")
    print(f"SOLVE = {results.solver.termination_condition.swapcase()}")

    print_infeasible_constraints(m)


if __name__ == "__main__":
    main()
