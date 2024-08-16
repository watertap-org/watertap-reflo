import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    TransformationFactory,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import *

# Flowsheet function imports
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_design_model import (
    get_n_time_points,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel


def build_system(model_options,n_time_points):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    
    m.fs.params = SeawaterParameterBlock()

    m.fs.md = FlowsheetBlock(dynamic=False)
    build_streams(m.fs, m.fs.params)
    build_md(m, m.fs.md,model_options,n_time_points)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_streams(blk, prop_package):
    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)


def build_md(m, blk, model_options,n_time_points) -> None:

    print(f'\n{"=======> BUILDING MEMBRANE DISTILLATION SYSTEM <=======":^60}\n')

    # n_time_points = 1

    # Build the multiperiod object for MD
    blk.unit = MultiPeriodModel(
        n_time_points= n_time_points,
        process_model_func=build_vagmd_flowsheet,
        linking_variable_func=get_vagmd_batch_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    feed_flow_rate = pyunits.convert(m.fs.feed.properties[0.0].flow_vol_phase["Liq"],to_units=pyunits.L/pyunits.h)()

    # Because its multiperiod, each instance is assigned an initial value based on model input (kwargs)
    blk.unit.build_multi_period_model(
        model_data_kwargs={t: model_options for t in range(n_time_points)},
        flowsheet_options=model_options,
        unfix_dof_options={"feed_flow_rate": feed_flow_rate},
    )


def init_md(m, blk, model_options, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING MEMBRANE DISTILLATION --------------------\n\n"
    )
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"MD Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    # blk.feed.initialize(optarg=optarg)
    # propagate_state(blk.feed_to_unit)
    feed_flow_rate = pyunits.convert(m.fs.feed.properties[0.0].flow_vol_phase["Liq"],to_units=pyunits.L/pyunits.h)()
    feed_salinity = m.fs.feed.properties[0.0].conc_mass_phase_comp["Liq","TDS"]
    feed_temp = m.fs.feed.properties[0.0].temperature() - 273.15

    add_performance_constraints(m,blk.unit,model_options)

    iscale.calculate_scaling_factors(blk.unit)
    solver = get_solver()
    active_blks = blk.unit.get_active_process_blocks()
    for active in active_blks:
        fix_dof_and_initialize(
            m=active,
            feed_flow_rate=feed_flow_rate,
            feed_salinity=feed_salinity,
            feed_temp=feed_temp,
        )
        result = solver.solve(active)
        unfix_dof(m=active, feed_flow_rate=feed_flow_rate)

    # propagate_state(blk.unit_to_disposal)
    # propagate_state(blk.unit_to_product)
    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)


def set_md_op_conditions(blk):

    active_blks = blk.unit.get_active_process_blocks()

    # Set-up for the first time period
    feed_salinity = m.fs.feed.properties[0.0].conc_mass_phase_comp["Liq","TDS"]
    feed_temp = m.fs.feed.properties[0.0].temperature()

    print("\n--------- MD TIME PERIOD 1 INPUTS ---------\n")
    print('Feed salinity in g/l:',feed_salinity())
    print('Feed temperature in C:', feed_temp-273.15)

    active_blks[0].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(
        feed_salinity
    )
    active_blks[0].fs.vagmd.feed_props[0].temperature.fix(feed_temp)
    active_blks[0].fs.acc_distillate_volume.fix(0)
    active_blks[0].fs.pre_feed_temperature.fix(feed_temp)
    active_blks[0].fs.pre_permeate_flow_rate.fix(0)
    active_blks[0].fs.acc_thermal_energy.fix(0)
    active_blks[0].fs.pre_thermal_power.fix(0)
    active_blks[0].fs.acc_cooling_energy.fix(0)
    active_blks[0].fs.pre_cooling_power.fix(0)
    active_blks[0].fs.acc_electric_energy.fix(0)
    active_blks[0].fs.pre_cooling_pump_power_elec.fix(0)
    active_blks[0].fs.pre_feed_pump_power_elec.fix(0)


def set_system_conditions(blk, model_options):

    feed_flow_rate = model_options['feed_flow_rate']*pyunits.L/pyunits.h
    feed_salinity = model_options['feed_salinity']*pyunits.g/pyunits.L
    feed_temp = model_options['feed_temp']

    flow_mass_in = pyunits.convert(feed_flow_rate, to_units=pyunits.m**3/pyunits.s)*1000*pyunits.kg/pyunits.m**3
    feed_salinity_in = pyunits.convert(feed_salinity*feed_flow_rate, to_units=pyunits.kg/pyunits.s)
    print(feed_salinity(),feed_salinity_in())

    blk.feed.properties[0.0].flow_mass_phase_comp["Liq","H2O"].fix(flow_mass_in)
    blk.feed.properties[0.0].flow_mass_phase_comp["Liq","TDS"].fix(feed_salinity_in)

    blk.feed.properties[0.0].temperature.fix(pyunits.convert_temp_C_to_K(feed_temp))
    blk.feed.properties[0.0].conc_mass_phase_comp["Liq","TDS"]

    solver = get_solver()
    optarg = solver.options
    blk.feed.initialize(optarg=optarg)


def set_system_output(blk,n_time_points,model_options):
    active_blks = blk.md.unit.get_active_process_blocks()

    # Get final variables
    num_modules = active_blks[-1].fs.vagmd.num_modules
    acc_recovery = active_blks[-1].fs.acc_recovery_ratio
    # L/h -> converted to kg/s
    permeate_production_rate = num_modules*value(active_blks[-1].fs.acc_distillate_volume) / (
        value(active_blks[-1].fs.dt) * (n_time_points - 1) / 3600)*pyunits.L/pyunits.h
    permeate_production_mass_rate = pyunits.convert(permeate_production_rate, pyunits.m**3/pyunits.s) *1000*pyunits.kg/pyunits.m**3

    # TODO: Calculate temperature

    blk.product.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(permeate_production_mass_rate)
    blk.product.properties[0].flow_mass_phase_comp["Liq","TDS"].fix(0)
    blk.product.properties[0].temperature.fix()

    # L/h -> converted to kg/s
    brine_production_rate = num_modules*model_options['initial_batch_volume']*(1-acc_recovery)/(
        value(active_blks[-1].fs.dt) * (n_time_points - 1) / 3600)*pyunits.L/pyunits.h

    brine_production_mass_rate = pyunits.convert(brine_production_rate, pyunits.m**3/pyunits.s) *1000*pyunits.kg/pyunits.m**3
    
    brine_salinity = brine_production_rate*pyunits.convert(active_blks[-1].fs.pre_feed_salinity,to_units=pyunits.kg/pyunits.m**3)
    blk.disposal.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(brine_production_mass_rate)
    blk.disposal.properties[0].flow_mass_phase_comp["Liq","TDS"].fix(brine_salinity)
    blk.disposal.properties[0].temperature.fix()

def add_performance_constraints(m,blk,model_options):
    
    # Create accumlative energy terms
    blk.system_capacity = Var(
        initialize=model_options['system_capacity'],
        bounds=(0, None),
        units=pyunits.m**3 / pyunits.day,
        doc="System capacity (m3/day)",
    )

    blk.overall_thermal_power_requirement = Var(
        initialize=2e5,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Thermal power requirement (kW-th)",
    )

    blk.overall_elec_power_requirement = Var(
        initialize=300,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Electric power requirement (kW-e)",
    )

    @blk.Constraint(
        doc="Calculate the overall thermal power requirement through all periods"
    )
    def eq_thermal_power_requirement(b):
        return b.overall_thermal_power_requirement == (
            blk.get_active_process_blocks()[-1].fs.specific_energy_consumption_thermal
            * pyunits.convert(b.system_capacity, to_units=pyunits.m**3 / pyunits.h)
        )

    @blk.Constraint(
        doc="Calculate the overall electric power requirement through all periods"
    )
    def eq_elec_power_requirement(b):
        return b.overall_elec_power_requirement == (
            blk.get_active_process_blocks()[-1].fs.specific_energy_consumption_electric
            * pyunits.convert(b.system_capacity, to_units=pyunits.m**3 / pyunits.h)
        )

    iscale.calculate_scaling_factors(blk)

    if iscale.get_scaling_factor(blk.overall_thermal_power_requirement) is None:
        iscale.set_scaling_factor(blk.overall_thermal_power_requirement, 1e-5)

    if iscale.get_scaling_factor(blk.overall_elec_power_requirement) is None:
        iscale.set_scaling_factor(blk.overall_elec_power_requirement, 1e-3)


def add_md_costing(m, blk):
    """
    This function adds costing module to the target multiperiod module mp,
    by adding the unit costing package to the last time step,
    and overwritting flow values with accumulated ones

    Args:
        mp: A pyomo module that describes the vagmd multiperiod model
        flowsheet_costing_block: Costing block which has been created for the whole flowsheet

    Returns:
        object: A costing module associated to the multiperiod module
    """
    m.fs.costing = REFLOCosting()
    blk.system_capacity.fix()
    # Specify the last time step
    vagmd = blk.get_active_process_blocks()[-1].fs.vagmd
    vagmd.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

    # Overwrite the thermal and electric energy flow with the accumulated values
    vagmd.costing.costing_package.cost_flow(
        blk.overall_thermal_power_requirement - vagmd.thermal_power_requirement,
        "heat",
    )
    vagmd.costing.costing_package.cost_flow(
        blk.overall_elec_power_requirement - vagmd.elec_power_requirement,
        "electricity",
    )

    # Recalculate the number of modules required
    vagmd.eqn_num_modules.deactivate()

    blks = blk.get_active_process_blocks()
    n_time_points = len(blks)

    vagmd.num_modules_constraint = Constraint(
        expr=vagmd.num_modules
        == pyunits.convert(
            vagmd.system_capacity
            / sum(
                blk.get_active_process_blocks()[i].fs.vagmd.permeate_flux
                for i in range(n_time_points)
            )
            * n_time_points
            / blk.get_active_process_blocks()[0].fs.vagmd.module_area,
            to_units=pyunits.dimensionless,
        )
    )



def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results



if __name__ == "__main__":

    model_options = {
        "dt": None,
        "system_capacity": 1000,
        "feed_flow_rate": 600, # L/h
        "evap_inlet_temp": 80,
        "cond_inlet_temp": 25,
        "feed_temp": 50,
        "feed_salinity": 35,
        "initial_batch_volume": 50, #L
        "module_type": "AS7C1.5L",
        "cooling_system_type": "closed",
        "cooling_inlet_temp": 25,
    }

    # Calculate the number of periods to reach target recovery rate by solving the system first
    n_time_points_check = get_n_time_points(
        dt=model_options['dt'],
        feed_flow_rate=model_options['feed_flow_rate'],
        evap_inlet_temp=model_options['evap_inlet_temp'],
        cond_inlet_temp=model_options['cond_inlet_temp'],
        feed_temp=model_options['feed_temp'],
        feed_salinity=model_options['feed_salinity'],
        recovery_ratio= 0.5,
        initial_batch_volume=model_options['initial_batch_volume'],
        module_type=model_options['module_type'],
        cooling_system_type=model_options['cooling_system_type'],
        cooling_inlet_temp=model_options['cooling_inlet_temp'],
        )

    
    n_time_points = 2

    m = build_system(model_options,n_time_points)
    
    set_system_conditions(m.fs,model_options)
    
    init_md(m, m.fs.md, model_options)
    add_md_costing(m, m.fs.md.unit)
    
    set_md_op_conditions(m.fs.md)

    solver = get_solver()
    try:
        results = solver.solve(m.fs.md)
    except ValueError:
        m.fs.md.unit.active_blocks[0].fs.vagmd.display()
        # print_infeasible_constraints(m)
        
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.costing.total_investment_factor.fix(1)
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.capital_recovery_factor.fix(0.08764)
    m.fs.costing.wacc.unfix()

    m.fs.costing.cost_process()

    active_blks = m.fs.md.unit.get_active_process_blocks()
    m.fs.costing.add_annual_water_production(active_blks[-1].fs.vagmd.system_capacity)
    m.fs.costing.add_LCOW(active_blks[-1].fs.vagmd.system_capacity)

    # assert degrees_of_freedom(m) == 0
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    results = solver.solve(m)

    set_system_output(m.fs,n_time_points,model_options)
    print("Overall LCOW ($/m3): ", value(m.fs.costing.LCOW))
    print('Calculate n_time points', n_time_points_check)


# def report_UF(m, blk, stream_table=False):
#     print(f"\n\n-------------------- UF Report --------------------\n")
#     print("\n")
#     print(
#         f'{"Inlet Flow Volume":<30s}{value(m.fs.UF.feed.properties[0.0].flow_vol):<10.3f}{pyunits.get_units(m.fs.UF.feed.properties[0.0].flow_vol)}'
#     )
#     print(f'{"UF Performance:":<30s}')
#     print(
#         f'{"    Recovery":<30s}{100*m.fs.UF.unit.recovery_frac_mass_H2O[0.0].value:<10.1f}{"%"}'
#     )
#     print(
#         f'{"    TDS Removal":<30s}{100*m.fs.UF.unit.removal_frac_mass_comp[0.0,"tds"].value:<10.1f}{"%"}'
#     )
#     print(
#         f'{"    TSS Removal":<30s}{100*m.fs.UF.unit.removal_frac_mass_comp[0.0,"tss"].value:<10.1f}{"%"}'
#     )
#     print(
#         f'{"    Energy Consumption":<30s}{m.fs.UF.unit.electricity[0.0].value:<10.3f}{pyunits.get_units(m.fs.UF.unit.electricity[0.0])}'
#     )
#     print(
#         f'{"    Specific Energy Cons.":<30s}{value(m.fs.UF.unit.energy_electric_flow_vol_inlet):<10.3f}{pyunits.get_units(m.fs.UF.unit.energy_electric_flow_vol_inlet)}'
#     )