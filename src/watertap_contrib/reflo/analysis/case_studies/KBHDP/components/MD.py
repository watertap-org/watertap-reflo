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

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_multiperiod_flowsheet import (
    get_vagmd_batch_variable_pairs,
    unfix_dof,
)

from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

import matplotlib.pyplot as plt


def propagate_state(arc):
    _prop_state(arc)


def build_system(overall_recovery=0.5, n_time_points=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Property package
    m.fs.params = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.params)
    m.fs.product = Product(property_package=m.fs.params)
    m.fs.disposal = Product(property_package=m.fs.params)

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    model_options, n_time_points = build_md(m, m.fs.md, overall_recovery, n_time_points)
    add_connections(m)

    return m, model_options, n_time_points


def add_connections(m):

    m.fs.feed_to_md = Arc(source=m.fs.feed.outlet, destination=m.fs.md.feed.inlet)

    m.fs.md_to_product = Arc(
        source=m.fs.md.permeate.outlet, destination=m.fs.product.inlet
    )

    m.fs.md_to_disposal = Arc(
        source=m.fs.md.concentrate.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_md_model_options(
    feed_flow_rate, feed_salinity, overall_recovery, n_time_points=None
):

    system_capacity = overall_recovery * pyunits.convert(
        feed_flow_rate, to_units=pyunits.m**3 / pyunits.day
    )
    feed_salinity = pyunits.convert(feed_salinity, to_units=pyunits.g / pyunits.L)

    model_options = {
        "dt": None,
        "system_capacity": system_capacity,  # m3/day
        "feed_flow_rate": 750,  # L/h
        "evap_inlet_temp": 80,
        "cond_inlet_temp": 30,
        "feed_temp": 30,
        "feed_salinity": feed_salinity(),  # g/L
        "initial_batch_volume": 50,  # L
        "module_type": "AS26C7.2L",
        "cooling_system_type": "closed",
        "cooling_inlet_temp": 25,
    }

    if n_time_points == None:
        # Calculate the number of periods to reach target recovery rate by solving the system first
        n_time_points = get_n_time_points(
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

    return model_options, n_time_points


def build_md(m, blk, overall_recovery=0.5, n_time_points=None) -> None:

    print(f'\n{"=======> BUILDING MEMBRANE DISTILLATION SYSTEM <=======":^60}\n')

    # Build a feed, permeate and brine state function for MD

    blk.feed = StateJunction(property_package=m.fs.params)
    blk.permeate = StateJunction(property_package=m.fs.params)
    blk.concentrate = StateJunction(property_package=m.fs.params)

    model_options, n_time_points = set_md_model_options(
        Qin, feed_salinity, overall_recovery, n_time_points
    )

    # Build the multiperiod object for MD
    blk.unit = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_vagmd_flowsheet,
        linking_variable_func=get_vagmd_batch_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    return model_options, n_time_points


def init_md(blk, model_options, n_time_points, verbose=True, solver=None):
    blk.feed.properties[0]._flow_vol_phase()
    blk.feed.properties[0]._conc_mass_phase_comp()

    blk.feed.initialize()

    feed_flow_rate = pyunits.convert(
        blk.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.h
    )()

    feed_temp = pyunits.convert_temp_K_to_C(blk.feed.properties[0].temperature())

    # Because its multiperiod, each instance is assigned an initial value based on model input (kwargs)
    blk.unit.build_multi_period_model(
        model_data_kwargs={t: model_options for t in range(n_time_points)},
        flowsheet_options=model_options,
        unfix_dof_options={"feed_flow_rate": feed_flow_rate},
    )

    add_performance_constraints(m, blk.unit, model_options)

    blk.unit.system_capacity.fix(model_options["system_capacity"])

    solver = get_solver()
    active_blks = blk.unit.get_active_process_blocks()
    for active in active_blks:
        fix_dof_and_initialize(
            m=active,
            feed_temp=feed_temp,
        )
        result = solver.solve(active)
        unfix_dof(m=active, feed_flow_rate=feed_flow_rate)

    # Build connection to permeate state junction
    blk.permeate.properties[0]._flow_vol_phase()

    @blk.Constraint(
        doc="Assign the permeate flow rate to its respective state junction"
    )
    def get_permeate_flow(b):
        num_modules = blk.unit.get_active_process_blocks()[-1].fs.vagmd.num_modules

        return (
            b.permeate.properties[0].flow_vol_phase["Liq"]
            ==
            # pyunits.convert(
            pyunits.convert(
                (
                    num_modules
                    * b.unit.get_active_process_blocks()[-1].fs.acc_distillate_volume
                    / (
                        b.unit.get_active_process_blocks()[-1].fs.dt
                        * (n_time_points - 1)
                    )
                ),
                to_units=pyunits.m**3 / pyunits.s,
            )
            #     *1000*pyunits.kg/pyunits.m**3
            #    , to_units = pyunits.kg/pyunits.s)
        )

    blk.permeate.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)
    blk.permeate.properties[0].pressure.fix(101325)
    blk.permeate.properties[0].temperature.fix(298.15)

    # Build connection to concentrate state junction

    blk.concentrate.properties[0]._flow_vol_phase()
    blk.concentrate.properties[0]._conc_mass_phase_comp()

    @blk.Constraint(
        doc="Assign the concentrate flow rate to its respective state junction"
    )
    def get_concentrate_flow(b):
        num_modules = blk.unit.get_active_process_blocks()[-1].fs.vagmd.num_modules

        return b.concentrate.properties[0].flow_vol_phase["Liq"] == pyunits.convert(
            (
                num_modules
                * model_options["initial_batch_volume"]
                * pyunits.L
                * (1 - b.unit.get_active_process_blocks()[-1].fs.acc_recovery_ratio)
                / (b.unit.get_active_process_blocks()[-1].fs.dt * (n_time_points - 1))
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    @blk.Constraint(
        doc="Assign the concentrate concentration to its respective state junction"
    )
    def get_concentrate_conc(b):
        return b.concentrate.properties[0].conc_mass_phase_comp[
            "Liq", "TDS"
        ] == pyunits.convert(
            b.unit.get_active_process_blocks()[-1]
            .fs.vagmd.feed_props[0]
            .conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.kg / pyunits.m**3,
        )

    blk.concentrate.properties[0].pressure.fix(101325)
    blk.concentrate.properties[0].temperature.fix(298.15)


def init_system(m, blk, model_options, n_time_points, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING MEMBRANE DISTILLATION --------------------\n\n"
    )
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"MD Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    # # Converting to units L/h and supplying value only
    # feed_flow_rate = model_options["feed_flow_rate"]
    # feed_salinity = model_options["feed_salinity"]
    # # Converting temperature to C units
    # feed_temp = model_options["feed_temp"]

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_md)

    init_md(blk, model_options, n_time_points, verbose=True, solver=None)

    propagate_state(m.fs.md_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.md_to_disposal)
    m.fs.disposal.initialize()


def set_system_op_conditions(blk, model_options):

    feed_flow_rate = model_options["feed_flow_rate"]
    feed_salinity = model_options["feed_salinity"]
    feed_temp = model_options["feed_temp"]

    blk.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): pyunits.convert(
                feed_flow_rate * pyunits.L / pyunits.h,
                to_units=pyunits.m**3 / pyunits.s,
            ),
            ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
            ("temperature", None): feed_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def set_md_op_conditions(blk):

    active_blks = blk.unit.get_active_process_blocks()

    # # Set-up for the first time period
    # feed_salinity = model_options["feed_salinity"]
    # feed_temp = model_options["feed_temp"]
    feed_flow_rate = pyunits.convert(
        m.fs.md.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.h
    )()

    feed_salinity = m.fs.md.feed.properties[0].conc_mass_phase_comp["Liq", "TDS"]()
    feed_temp = pyunits.convert_temp_K_to_C(m.fs.md.feed.properties[0].temperature())

    print("\n--------- MD TIME PERIOD 1 INPUTS ---------\n")
    print("Feed flow rate in L/h:", feed_flow_rate)
    print("Feed salinity in g/l:", feed_salinity)
    print("Feed temperature in C:", feed_temp)
    print("\n")

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


def md_output(blk, n_time_points, model_options):

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
    # blk.system_capacity.fix()
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


def calc_costing(m, blk):

    m.fs.costing.total_investment_factor.fix(1)
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.capital_recovery_factor.fix(0.08764)
    m.fs.costing.wacc.unfix()

    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    # active_blks = m.fs.md.unit.get_active_process_blocks()
    # m.fs.costing.add_annual_water_production(active_blks[-1].fs.vagmd.system_capacity)
    # m.fs.costing.add_LCOW(active_blks[-1].fs.vagmd.system_capacity)

    prod_flow = pyunits.convert(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    m.fs.costing.add_annual_water_production(prod_flow)
    m.fs.costing.add_LCOW(prod_flow)


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


def report_MD(m):

    active_blks = m.fs.md.unit.get_active_process_blocks()

    print(f"\n\n-------------------- MD Operations Report --------------------\n")
    print("\n")

    print(
        f'{"System Capacity":<30s}{value(active_blks[0].fs.vagmd.system_capacity):<10,.2f}{pyunits.get_units(active_blks[0].fs.vagmd.system_capacity)}'
    )

    print(
        f'{"Feed stream salinity":<30s}{value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(m.fs.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )

    print(
        f'{"MD Period 1 Feed salinity":<30s}{value(m.fs.md.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(m.fs.md.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )

    print(
        f'{"Number of modules:":<30s}{value(active_blks[-1].fs.vagmd.num_modules):<10.2f}'
    )

    print(f'{"Membrane type":<30s}{active_blks[-1].fs.vagmd.config.module_type}')

    print(
        f'{"Accumulated recovery":<30s}{value(active_blks[-1].fs.acc_recovery_ratio):<10.2f}{pyunits.get_units(active_blks[-1].fs.acc_recovery_ratio)}'
    )

    perm_flow = pyunits.convert(
        m.fs.md.permeate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Permeate flow rate":<30s}{value(perm_flow):<10,.2f}{pyunits.get_units(perm_flow)}'
    )

    conc_flow = pyunits.convert(
        m.fs.md.concentrate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Concentrate flow rate":<30s}{value(conc_flow):<10,.2f}{pyunits.get_units(conc_flow)}'
    )

    print(
        f'{"Concentrate concentration":<30s}{value(m.fs.md.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(m.fs.md.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )

    prod_flow = pyunits.convert(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Product flow rate":<30s}{value(prod_flow):<10,.2f}{pyunits.get_units(prod_flow)}'
    )

    disp_flow = pyunits.convert(
        m.fs.disposal.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Disposal flow rate":<30s}{value(disp_flow):<10,.2f}{pyunits.get_units(disp_flow)}'
    )

    print(
        f'{"STEC":<30s}{value(active_blks[-1].fs.specific_energy_consumption_thermal):<10.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_thermal)}'
    )

    print(
        f'{"SEC":<30s}{value(active_blks[-1].fs.specific_energy_consumption_electric):<10.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_electric)}'
    )

    print(
        f'{"Overall thermal requirement":<30s}{value(m.fs.md.unit.overall_thermal_power_requirement):<10.2f}{pyunits.get_units(m.fs.md.unit.overall_thermal_power_requirement)}'
    )

    print(
        f'{"Overall elec requirement":<30s}{value(m.fs.md.unit.overall_elec_power_requirement):<10.2f}{pyunits.get_units(m.fs.md.unit.overall_elec_power_requirement)}'
    )

    # Operating parameters to track

    print(f"\n\n-------------------- MD Costing Report --------------------\n")
    print("\n")

    print(
        f'{"LCOW":<30s}{value(m.fs.costing.LCOW):<20.5f}{pyunits.get_units(m.fs.costing.LCOW)}'
    )

    print(
        f'{"Capital Cost":<30s}{value(m.fs.costing.aggregate_capital_cost):<20,.2f}{pyunits.get_units(m.fs.costing.LCOW)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(m.fs.costing.aggregate_fixed_operating_cost):<20,.2f}{pyunits.get_units(m.fs.costing.LCOW)}'
    )

    print(
        f'{"Membrane Cost":<30s}{value(active_blks[-1].fs.vagmd.costing.module_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.module_cost)}'
    )

    print(
        f'{"Aggregated Heat Cost":<30s}{value(m.fs.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_costs["heat"])}'
    )

    print(
        f'{"Aggregated Elec Cost":<30s}{value(m.fs.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_costs["electricity"])}'
    )


if __name__ == "__main__":

    solver = get_solver()

    Qin = 4 * pyunits.Mgal / pyunits.day  # KBHDP flow rate
    Qin = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.day)

    feed_salinity = 12 * pyunits.g / pyunits.L

    overall_recovery = 0.45

    n_time_points = None
    m, model_options, n_time_points = build_system(
        overall_recovery, n_time_points=n_time_points
    )
    print("Number of time points being modeled:", n_time_points)

    set_system_op_conditions(m.fs, model_options)
    init_system(m, m.fs.md, model_options, n_time_points)

    set_md_op_conditions(m.fs.md)

    results = solver.solve(m)

    # Touching variables to solve for volumetric flow rate
    m.fs.product.properties[0].flow_vol_phase
    m.fs.disposal.properties[0].flow_vol_phase

    # results = solver.solve(m)

    add_md_costing(m, m.fs.md.unit)
    calc_costing(m, m.fs)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    report_MD(m)

    # try:
    #     results = solver.solve(m.fs.md)
    # except ValueError:
    #     # m.fs.md.unit.active_blocks[0].fs.vagmd.display()
    #     print_infeasible_constraints(m)

    active_blks = m.fs.md.unit.get_active_process_blocks()

    permeate_flow_rate, brine_flow_rate, brine_salinity = md_output(
        m.fs, n_time_points, model_options
    )

    time_period = [i for i in range(n_time_points)]
    t_minutes = [value(active_blks[i].fs.dt) * i / 60 for i in range(n_time_points)]

    heat_in = [value(active_blks[i].fs.pre_thermal_power) for i in range(n_time_points)]

    # plt.plot(t_minutes,heat_in)
    # plt.show()

    # print("\nCalculate n_time points", n_time_points)
    # print(
    #     "Time step duration:",
    #     value(active_blks[0].fs.dt),
    #     pyunits.get_units(active_blks[0].fs.dt),
    # )

    # batch_duration = n_time_points * active_blks[0].fs.dt
    # print("Batch duration:", value(batch_duration), pyunits.get_units(batch_duration))
    # no_of_batches = 24 * pyunits.h / pyunits.convert(batch_duration, to_units=pyunits.h)
    # print("Number of batches in 24h:", value(no_of_batches))

    # m.fs.md.permeate.properties[0].flow_vol_phase.display()
