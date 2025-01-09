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
    SolverFactory,
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

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_multiperiod_unit_model import (
    VAGMDbatchSurrogate,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_multiperiod_flowsheet import (
    get_vagmd_batch_variable_pairs,
    unfix_dof,
    build_VAGMD_batch_multiperiod_fs,
)

from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

import matplotlib.pyplot as plt

__all__ = [
    "build_md",
    "set_md_model_options",
    "init_md",
    "report_MD",
    "report_md_costing",
]


def propagate_state(arc):
    _prop_state(arc)


def build_system(Qin=4, Cin=12, water_recovery=0.5):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.inlet_flow_rate = pyunits.convert(
        Qin * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = water_recovery

    # Property package
    m.fs.properties = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    build_md(m, m.fs.md, prop_package=m.fs.properties)
    add_connections(m)

    return m


def add_connections(m):

    m.fs.feed_to_md = Arc(source=m.fs.feed.outlet, destination=m.fs.md.feed.inlet)

    m.fs.md_to_product = Arc(
        source=m.fs.md.permeate.outlet, destination=m.fs.product.inlet
    )

    m.fs.md_to_disposal = Arc(
        source=m.fs.md.concentrate.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_md_model_options(m, blk, n_time_points=None):

    m.system_capacity = m.water_recovery * pyunits.convert(
        m.inlet_flow_rate, to_units=pyunits.m**3 / pyunits.day
    )
    m.feed_salinity = pyunits.convert(m.inlet_salinity, to_units=pyunits.g / pyunits.L)

    model_options = {
        "dt": None,
        "system_capacity": value(m.system_capacity),  # m3/day
        "feed_flow_rate": 750,  # L/h
        "evap_inlet_temp": 80,
        "cond_inlet_temp": 20,
        "feed_temp": 25,
        "feed_salinity": value(m.feed_salinity),  # g/L
        "initial_batch_volume": 50,  # L
        "module_type": "AS26C7.2L",
        "cooling_system_type": "closed",
        "cooling_inlet_temp": 25,
        "recovery_ratio": m.water_recovery,
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
            recovery_ratio=model_options["recovery_ratio"],
            initial_batch_volume=model_options["initial_batch_volume"],
            module_type=model_options["module_type"],
            cooling_system_type=model_options["cooling_system_type"],
            cooling_inlet_temp=model_options["cooling_inlet_temp"],
        )

    blk.model_input = model_options
    blk.n_time_points = n_time_points


def build_md(m, blk, prop_package=None) -> None:

    print(f'\n{"=======> BUILDING MEMBRANE DISTILLATION SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties
    # Build a feed, permeate and brine state function for MD

    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = StateJunction(property_package=prop_package)
    blk.concentrate = StateJunction(property_package=prop_package)

    set_md_model_options(m, blk, n_time_points=None)

    blk.unit = VAGMDbatchSurrogate(model_input=blk.model_input)


def init_md(m, blk, verbose=True, solver=None):

    blk.feed.initialize()

    # Build connection to permeate state junction
    blk.permeate.properties[0]._flow_vol_phase

    @blk.Constraint(
        doc="Assign the permeate flow rate to its respective state junction"
    )
    def get_permeate_flow(b):
        # num_modules = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.num_modules

        vagmd = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd
        num_modules = pyunits.convert(
            vagmd.system_capacity
            / sum(
                b.unit.mp.get_active_process_blocks()[i].fs.vagmd.permeate_flux
                for i in range(blk.n_time_points)
            )
            * blk.n_time_points
            / b.unit.mp.get_active_process_blocks()[0].fs.vagmd.module_area,
            to_units=pyunits.dimensionless,
        )

        return b.permeate.properties[0].flow_vol_phase["Liq"] == pyunits.convert(
            (
                num_modules
                * b.unit.mp.get_active_process_blocks()[-1].fs.acc_distillate_volume
                / (
                    b.unit.mp.get_active_process_blocks()[-1].fs.dt
                    * (blk.n_time_points - 1)
                )
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    blk.permeate.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)
    blk.permeate.properties[0].pressure.fix(101325)
    blk.permeate.properties[0].temperature.fix(298.15)

    # # Build connection to concentrate state junction

    blk.concentrate.properties[0]._flow_vol_phase
    blk.concentrate.properties[0]._conc_mass_phase_comp

    @blk.Constraint(
        doc="Assign the concentrate flow rate to its respective state junction"
    )
    def get_concentrate_flow(b):
        num_modules = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.num_modules

        return b.concentrate.properties[0].flow_vol_phase["Liq"] == pyunits.convert(
            (
                b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.system_capacity
                * (1 - m.water_recovery)
                / m.water_recovery
                #     num_modules
                #     * blk.model_input["initial_batch_volume"]
                #     * pyunits.L
                #     * (1 - b.unit.mp.get_active_process_blocks()[-1].fs.acc_recovery_ratio)
                #     / (b.unit.mp.get_active_process_blocks()[-1].fs.dt * (blk.n_time_points - 1))
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
            b.unit.mp.get_active_process_blocks()[-1]
            .fs.vagmd.feed_props[0]
            .conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.kg / pyunits.m**3,
        )

    blk.concentrate.properties[0].pressure.fix(101325)
    blk.concentrate.properties[0].temperature.fix(298.15)


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING MEMBRANE DISTILLATION --------------------\n\n"
    )
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print("\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_md)

    init_md(m, m.fs.md, verbose=True, solver=None)

    propagate_state(m.fs.md_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.md_to_disposal)
    m.fs.disposal.initialize()


def set_system_op_conditions(m):

    feed_flow_rate = m.fs.md.model_input["feed_flow_rate"]
    feed_salinity = m.fs.md.model_input["feed_salinity"]
    feed_temp = m.fs.md.model_input["feed_temp"]

    m.fs.feed.properties.calculate_state(
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

def calc_costing(m, blk):
    m.fs.costing.capital_recovery_factor.fix(0.08764)
    m.fs.costing.wacc.unfix()

    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    # Touching variables to solve for volumetric flow rate
    m.fs.product.properties[0].flow_vol_phase

    prod_flow = pyunits.convert(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    m.fs.costing.add_annual_water_production(prod_flow)
    m.fs.costing.add_LCOW(prod_flow)


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print_infeasible_bounds(m)
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print(msg)
        return results


def report_MD(m, blk):

    active_blks = blk.unit.mp.get_active_process_blocks()

    print(f"\n\n-------------------- MD Operations Report --------------------\n")
    print("\n")

    print(
        f'{"System Capacity":<30s}{value(active_blks[0].fs.vagmd.system_capacity):<10,.2f}{pyunits.get_units(active_blks[0].fs.vagmd.system_capacity)}'
    )

    # print(
    #     f'{"Feed stream salinity":<30s}{value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(m.fs.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    # )

    # print(
    #     f'{"MD Period 1 Feed salinity":<30s}{value(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    # )

    print(
        f'{"Number of modules:":<30s}{value(active_blks[-1].fs.vagmd.num_modules):<10.2f}'
    )

    print(f'{"Membrane type":<30s}{active_blks[-1].fs.vagmd.config.module_type}')

    print(
        f'{"Accumulated recovery":<30s}{value(active_blks[-1].fs.acc_recovery_ratio):<10.2f}{pyunits.get_units(active_blks[-1].fs.acc_recovery_ratio)}'
    )

    perm_flow = pyunits.convert(
        blk.permeate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Permeate flow rate":<30s}{value(perm_flow):<10,.2f}{pyunits.get_units(perm_flow)}'
    )

    conc_flow = pyunits.convert(
        blk.concentrate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Concentrate flow rate":<30s}{value(conc_flow):<10,.2f}{pyunits.get_units(conc_flow)}'
    )

    print(
        f'{"Concentrate concentration":<30s}{value(blk.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"]):<10.2f}{pyunits.get_units(blk.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )

    # prod_flow = pyunits.convert(
    #     m.fs.product.properties[0].flow_vol_phase["Liq"],
    #     to_units=pyunits.m**3 / pyunits.day,
    # )

    # print(
    #     f'{"Product flow rate":<30s}{value(prod_flow):<10,.2f}{pyunits.get_units(prod_flow)}'
    # )

    # disp_flow = pyunits.convert(
    #     m.fs.disposal.properties[0].flow_vol_phase["Liq"],
    #     to_units=pyunits.m**3 / pyunits.day,
    # )

    # print(
    #     f'{"Disposal flow rate":<30s}{value(disp_flow):<10,.2f}{pyunits.get_units(disp_flow)}'
    # )

    print(
        f'{"STEC":<30s}{value(active_blks[-1].fs.specific_energy_consumption_thermal):<10.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_thermal)}'
    )

    print(
        f'{"SEC":<30s}{value(active_blks[-1].fs.specific_energy_consumption_electric):<10.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_electric)}'
    )

    print(
        f'{"Overall thermal requirement":<30s}{value(blk.unit.overall_thermal_power_requirement):<10.2f}{pyunits.get_units(blk.unit.overall_thermal_power_requirement)}'
    )

    print(
        f'{"Overall elec requirement":<30s}{value(blk.unit.overall_elec_power_requirement):<10.2f}{pyunits.get_units(blk.unit.overall_elec_power_requirement)}'
    )


def report_md_costing(m, blk):

    active_blks = blk.md.unit.mp.get_active_process_blocks()

    print(f"\n\n-------------------- MD Costing Report --------------------\n")
    print("\n")

    print(
        f'{"Capital Cost":<30s}{value(active_blks[-1].fs.vagmd.costing.capital_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(active_blks[-1].fs.vagmd.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.fixed_operating_cost)}'
    )

    # print(
    #     f'{"Variable Operating Cost":<30s}{value(blk.costing.variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.variable_operating_cost)}'
    # )

    # print(
    #     f'{"Aggregated Variable Operating Cost":<30s}{value(blk.costing.aggregate_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.aggregate_variable_operating_cost)}'
    # )

    # print(
    #     f'{"Total Operating Cost":<30s}{value(blk.costing.total_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.total_operating_cost)}'
    # )

    print(
        f'{"Membrane Cost":<30s}{value(active_blks[-1].fs.vagmd.costing.module_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.module_cost)}'
    )

    # print(
    # f'{"Heat flow":<30s}{value(blk.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_heat)}'
    # )

    # print(
    #     f'{"Aggregated Heat Cost":<30s}{value(blk.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["heat"])}'
    # )

    # print(
    #     f'{"Elec Flow":<30s}{value(blk.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_electricity)}'
    # )

    # print(
    #     f'{"Aggregated Elec Cost":<30s}{value(blk.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["electricity"])}'
    # )


if __name__ == "__main__":
    m = build_system(Qin=4, Cin=12, water_recovery=0.5)
    set_system_op_conditions(m)
    init_system(m)

    results = solve(m, tee=False)

    m.fs.disposal.properties[0].flow_vol_phase

    m.fs.md.unit.add_costing_module(m.fs.costing)

    calc_costing(m, m.fs)

    print("\nSystem Degrees of Freedom:", degrees_of_freedom(m), "\n")

    assert degrees_of_freedom(m) == 0

    results = solve(m)
    print("\n--------- Cost solve Completed ---------\n")

    print(
        "Inlet flow rate in m3/day:",
        value(
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol_phase["Liq"],
                pyunits.m**3 / pyunits.day,
            )
        ),
    )
    report_MD(m, m.fs.md)
    report_md_costing(m, m.fs)

    print("\n")
    print(
        f'Sys Feed Flow Rate: {value(pyunits.convert(m.fs.feed.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'MD  Feed Flow Rate: {value(pyunits.convert(m.fs.md.feed.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'Sys Perm Flow Rate: {value(pyunits.convert(m.fs.product.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'MD  Perm Flow Rate: {value(pyunits.convert(m.fs.md.permeate.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'Sys Conc Flow Rate: {value(pyunits.convert(m.fs.disposal.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'MD  Conc Flow Rate: {value(pyunits.convert(m.fs.md.concentrate.properties[0].flow_vol_phase["Liq"], pyunits.m ** 3 / pyunits.day)):<10.2f} m3/day'
    )
    print(
        f'Calculated Recovery: {value(m.fs.md.permeate.properties[0].flow_vol_phase["Liq"] / (m.fs.md.permeate.properties[0].flow_vol_phase["Liq"] + m.fs.md.concentrate.properties[0].flow_vol_phase["Liq"])):<10.2f}'
    )