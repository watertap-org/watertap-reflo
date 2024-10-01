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
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting

# Flowsheet function imports
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_design_model import (
    get_n_time_points,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet_multiperiod import (
    get_vagmd_batch_variable_pairs,
    unfix_dof,
)

from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel


def build_system(model_options, n_time_points):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Property package
    m.fs.params = SeawaterParameterBlock()

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    build_md(m, m.fs.md, model_options, n_time_points)

    return m


def build_md(m, blk, model_options, n_time_points) -> None:

    print(f'\n{"=======> BUILDING MEMBRANE DISTILLATION SYSTEM <=======":^60}\n')

    # Build the multiperiod object for MD
    blk.unit = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_vagmd_flowsheet,
        linking_variable_func=get_vagmd_batch_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    # Converting to units L/h and supplying value only
    feed_flow_rate = model_options["feed_flow_rate"]

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

    # Converting to units L/h and supplying value only
    feed_flow_rate = model_options["feed_flow_rate"]
    feed_salinity = model_options["feed_salinity"]
    # Converting temperature to C units
    feed_temp = model_options["feed_temp"]

    add_performance_constraints(m, blk.unit, model_options)

    # iscale.calculate_scaling_factors(blk.unit)
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


def set_md_op_conditions(blk):

    active_blks = blk.unit.get_active_process_blocks()

    # Set-up for the first time period
    feed_salinity = model_options["feed_salinity"]
    feed_temp = model_options["feed_temp"]

    print("\n--------- MD TIME PERIOD 1 INPUTS ---------\n")
    print("Feed salinity in g/l:", feed_salinity)
    print("Feed temperature in C:", feed_temp)

    active_blks[0].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(
        feed_salinity
    )
    active_blks[0].fs.vagmd.feed_props[0].temperature.fix(feed_temp + 273.15)
    active_blks[0].fs.acc_distillate_volume.fix(0)
    active_blks[0].fs.pre_feed_temperature.fix(feed_temp + 273.15)
    active_blks[0].fs.pre_permeate_flow_rate.fix(0)
    active_blks[0].fs.acc_thermal_energy.fix(0)
    active_blks[0].fs.pre_thermal_power.fix(0)
    active_blks[0].fs.acc_cooling_energy.fix(0)
    active_blks[0].fs.pre_cooling_power.fix(0)
    active_blks[0].fs.acc_electric_energy.fix(0)
    active_blks[0].fs.pre_cooling_pump_power_elec.fix(0)
    active_blks[0].fs.pre_feed_pump_power_elec.fix(0)


def system_output(blk, n_time_points, model_options):

    active_blks = blk.md.unit.get_active_process_blocks()

    # Get final variables
    num_modules = active_blks[-1].fs.vagmd.num_modules
    acc_recovery = active_blks[-1].fs.acc_recovery_ratio
    # L/h -> converted to kg/s
    permeate_production_rate = (
        num_modules
        * value(active_blks[-1].fs.acc_distillate_volume)
        / (value(active_blks[-1].fs.dt) * (n_time_points - 1) / 3600)
        * pyunits.L
        / pyunits.h
    )
    permeate_flow_rate = pyunits.convert(
        permeate_production_rate, pyunits.m**3 / pyunits.day
    )

    # L/h -> converted to kg/s
    brine_production_rate = (
        num_modules
        * model_options["initial_batch_volume"]
        * (1 - acc_recovery)
        / (value(active_blks[-1].fs.dt) * (n_time_points - 1) / 3600)
        * pyunits.L
        / pyunits.h
    )

    brine_flow_rate = pyunits.convert(brine_production_rate, pyunits.m**3 / pyunits.day)

    brine_salinity = pyunits.convert(
        active_blks[-1].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.kg / pyunits.m**3,
    )

    return permeate_flow_rate, brine_flow_rate, brine_salinity


def add_performance_constraints(m, blk, model_options):

    # Create accumlative energy terms
    blk.system_capacity = Var(
        initialize=model_options["system_capacity"],
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
    vagmd.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

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

    m.fs.costing.total_investment_factor.fix(1)
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.capital_recovery_factor.fix(0.08764)
    m.fs.costing.wacc.unfix()

    m.fs.costing.cost_process()

    active_blks = m.fs.md.unit.get_active_process_blocks()
    m.fs.costing.add_annual_water_production(active_blks[-1].fs.vagmd.system_capacity)
    m.fs.costing.add_LCOW(active_blks[-1].fs.vagmd.system_capacity)


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


def report_MD(m, stream_table=False):
    print(f"\n\n-------------------- MD Report --------------------\n")
    print("\n")

    active_blks = m.fs.md.unit.get_active_process_blocks()

    print("Number of modules:", value(active_blks[-1].fs.vagmd.num_modules))

    print(
        f'{"Inlet Flow Volume":<30s}{value(active_blks[0].fs.vagmd.feed_props[0].flow_vol_phase["Liq"]):<10.5f}{pyunits.get_units(active_blks[-1].fs.vagmd.feed_props[0].flow_vol_phase["Liq"])}'
    )


if __name__ == "__main__":

    Qin = 4 * pyunits.Mgal / pyunits.day
    Qin = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.day)

    overall_recovery = 0.45
    system_capacity = overall_recovery * Qin

    model_options = {
        "dt": None,
        "system_capacity": system_capacity,  # m3/day
        "feed_flow_rate": 750,  # L/h
        "evap_inlet_temp": 80,
        "cond_inlet_temp": 30,
        "feed_temp": 30,
        "feed_salinity": 12,
        "initial_batch_volume": 50,  # L
        "module_type": "AS26C7.2L",
        "cooling_system_type": "closed",
        "cooling_inlet_temp": 25,
    }

    # Calculate the number of periods to reach target recovery rate by solving the system first
    n_time_points_check = get_n_time_points(
        dt=model_options["dt"],
        feed_flow_rate=model_options["feed_flow_rate"],
        evap_inlet_temp=model_options["evap_inlet_temp"],
        cond_inlet_temp=model_options["cond_inlet_temp"],
        feed_temp=model_options["feed_temp"],
        feed_salinity=model_options["feed_salinity"],
        recovery_ratio=overall_recovery,
        initial_batch_volume=model_options["initial_batch_volume"],
        module_type=model_options["module_type"],
        cooling_system_type=model_options["cooling_system_type"],
        cooling_inlet_temp=model_options["cooling_inlet_temp"],
    )

    # n_time_points = 2
    n_time_points = n_time_points_check

    m = build_system(model_options, n_time_points)

    init_md(m, m.fs.md, model_options)
    set_md_op_conditions(m.fs.md)
    add_md_costing(m, m.fs.md.unit)

    solver = get_solver()
    try:
        results = solver.solve(m.fs.md)
    except ValueError:
        m.fs.md.unit.active_blocks[0].fs.vagmd.display()
        # print_infeasible_constraints(m)

    assert degrees_of_freedom(m) == 0
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    results = solver.solve(m)

    active_blks = m.fs.md.unit.get_active_process_blocks()

    permeate_flow_rate, brine_flow_rate, brine_salinity = system_output(
        m.fs, n_time_points, model_options
    )
    print("\nOverall LCOW ($/m3): ", value(m.fs.costing.LCOW))

    print("\nCalculate n_time points", n_time_points)
    print(
        "Time step duration:",
        value(active_blks[0].fs.dt),
        pyunits.get_units(active_blks[0].fs.dt),
    )

    batch_duration = n_time_points * active_blks[0].fs.dt
    print("Batch duration:", value(batch_duration), pyunits.get_units(batch_duration))
    no_of_batches = 24 * pyunits.h / pyunits.convert(batch_duration, to_units=pyunits.h)
    print("Number of batches in 24h:", value(no_of_batches))

    print("Calculated number of modules:", value(active_blks[-1].fs.vagmd.num_modules))
    print(
        "Actual number of modules:",
        value(active_blks[-1].fs.vagmd.num_modules) / value(no_of_batches),
    )

    print(
        "\nSystem Capacity:", value(system_capacity), pyunits.get_units(system_capacity)
    )
    print(
        "Permeate mass flow:",
        value(permeate_flow_rate),
        pyunits.get_units(permeate_flow_rate),
    )
    print(
        "Brine mass flow rate:",
        value(brine_flow_rate),
        pyunits.get_units(brine_flow_rate),
    )
    print("Brine salinity:", value(brine_salinity), pyunits.get_units(brine_salinity))
    print("Accumulate recovery ratio:", value(active_blks[-1].fs.acc_recovery_ratio))
