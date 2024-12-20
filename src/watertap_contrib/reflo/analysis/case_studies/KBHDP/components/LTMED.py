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
from watertap.core.util.model_diagnostics import *
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.reflo.unit_models.surrogate import LTMEDSurrogate
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    check_jac,
    calc_scale,
)

__all__ = [
    "build_LTMED",
    "set_LTMED_operating_conditions",
    "init_LTMED",
    "add_LTMED_costing",
    "add_LTMED_scaling",
    "report_LTMED",
]


def _initialize(blk, verbose=False):
    if verbose:
        print("\n")
        print(
            f"{blk.name:<30s}{f'Degrees of Freedom at Initialization = {degrees_of_freedom(blk):<10.0f}'}"
        )
        print("\n")
    try:
        blk.initialize()
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")

        blk.report()
        print_infeasible_bounds(blk)
        print_close_to_bounds(blk)
        assert False


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    print("\n")


# The following 3 functions is a template for building and exporting the LTMED to the overall flowsheet
def build_LTMED(m, blk, liquid_prop, vapor_prop, number_effects=12):
    print(f'\n{"=======> BUILDING LTMED <=======":^60}\n')

    blk.feed = StateJunction(property_package=liquid_prop)
    blk.product = StateJunction(property_package=liquid_prop)
    blk.disposal = StateJunction(property_package=liquid_prop)

    blk.unit = LTMEDSurrogate(
        property_package_liquid=liquid_prop,
        property_package_vapor=vapor_prop,
        number_effects=number_effects,
    )

    # BUG LTMED Surrogate has no inlet port, so can't connect to feed
    blk.feed_to_LTMED = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.feed,
    )

    blk.LTMED_to_product = Arc(
        source=blk.unit.distillate,
        destination=blk.product.inlet,
    )

    blk.LTMED_to_disposal = Arc(
        source=blk.unit.brine,
        destination=blk.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_LTMED_operating_conditions(blk, recovery_ratio=0.45, steam_temperature=80):

    blk.unit.steam_props[0].temperature.fix(steam_temperature + 273.15)
    blk.unit.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    pass


def add_LTMED_scaling(m, blk):
    # set_scaling_factor(blk.product.properties[0.0].flow_mass_phase_comp['Liq','TDS'], 1e-11)
    # set_scaling_factor(blk.unit.thermal_power_requirement, 1e-6)
    # set_scaling_factor(blk.unit.distillate_props[0.0].flow_vol_phase['Liq'], 1)
    set_scaling_factor(blk.unit.specific_energy_consumption_thermal, 1e-2)
    set_scaling_factor(blk.unit.feed_props[0.0].flow_vol_phase["Liq"], 10)
    set_scaling_factor(blk.unit.cooling_out_props[0.0].flow_vol_phase["Liq"], 10)
    set_scaling_factor(blk.unit.distillate_props[0.0].flow_vol_phase["Liq"], 1e2)
    # set_scaling_factor(blk.unit.feed_cool_vol_flow, 1e-4)
    set_scaling_factor(
        blk.unit.distillate_props[0.0].flow_mass_phase_comp["Liq", "H2O"], 1e-3
    )
    set_scaling_factor(blk.unit.brine_props[0.0].flow_vol_phase["Liq"], 100)
    constraint_scaling_transform(blk.unit.eq_specific_area_per_m3_day, 1e-2)
    constraint_scaling_transform(blk.unit.eq_specific_area_kg_s, 1e-2)
    constraint_scaling_transform(blk.unit.eq_steam_mass_flow, 1e-2)
    constraint_scaling_transform(blk.unit.eq_specific_thermal_energy_consumption, 1e-2)
    constraint_scaling_transform(blk.unit.eq_thermal_power_requirement, 1e-4)
    constraint_scaling_transform(blk.unit.eq_feed_cool_vol_flow, 1e-4)
    constraint_scaling_transform(blk.unit.eq_feed_cool_mass_flow, 1e-2)

    for temp_stream in [
        blk.unit.eq_distillate_temp,
        blk.unit.eq_brine_temp,
        blk.unit.eq_cooling_temp,
    ]:
        for e in temp_stream:
            print(temp_stream[e])
            constraint_scaling_transform(temp_stream[e], 1e-2)
        for pressure_stream in [
            blk.unit.eq_feed_to_distillate_isobaric,
            blk.unit.eq_feed_to_brine_isobaric,
            blk.unit.eq_feed_to_cooling_isobaric,
        ]:
            for e in pressure_stream:
                constraint_scaling_transform(pressure_stream[e], 1e-5)

    set_scaling_factor(blk.unit.feed_props[0.0].temperature, 1e-2)


def add_system_costing(m):
    """Add system level costing components"""
    m.fs.costing = ZeroOrderCosting()
    add_LTMED_costing(m, m.fs.treatment.LTMED)
    calc_costing(m, m.fs.treatment.LTMED)

    m.fs.costing.add_annual_water_production(
        m.fs.treatment.LTMED.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.treatment.LTMED.product.properties[0].flow_vol)


def add_LTMED_costing(m, blk, costing_blk=None):
    """Add LTMED model costing components"""
    if costing_blk is None:
        costing_blk = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_blk)


def calc_costing(m, blk):
    """Add system level solve for costing"""
    m.fs.costing.cost_process()
    # m.fs.costing.add_LCOW(blk.ec.properties_treated[0].flow_vol)
    # m.fs.costing.add_electricity_intensity(blk.ec.properties_treated[0].flow_vol)


def init_LTMED(m, blk, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING LTMED --------------------\n\n")
    # assert_no_degrees_of_freedom(m)

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_LTMED)

    _initialize(blk.unit, verbose=True)
    propagate_state(blk.LTMED_to_product)
    propagate_state(blk.LTMED_to_disposal)


def set_system_operating_conditions(m):
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(171.76)
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.6423)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(25 + 273.15)

    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
    )
    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "TDS")
    )

    calculate_scaling_factors(m)


# The following functions are used for testing the component in isolation on thise file
def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    treatment = m.fs.treatment = Block()

    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.liquid_prop)

    treatment.LTMED = FlowsheetBlock(dynamic=False)
    build_LTMED(m, treatment.LTMED, m.fs.liquid_prop, m.fs.vapor_prop)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=treatment.LTMED.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def init_system(m, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    # breakdown_dof(m.fs.treatment.LTMED)
    assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_LTMED(m, m.fs.treatment.LTMED)


def optimize(m, water_recovery=None, objective="LCOW"):

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.treatment.LTMED.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.treatment.LTMED.product.properties[0].flow_vol
    )

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.treatment.LTMED.unit.recovery_vol_phase[0.0, "Liq"].unfix()
        m.fs.water_recovery.fix(water_recovery)


def solve(m, solver=None, tee=True, raise_on_failure=True, debug=False):
    # ---solving---
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            print("\n--------- TREATMENT ---------\n")
            check_jac(m.fs.treatment)

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(m)
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print('\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print('\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print('\n{"=======> CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print("\n--------- FAILED SOLVE!!! ---------\n")
        print(msg)
        assert False


def breakdown_dof(blk):
    equalities = [c for c in activated_equalities_generator(blk)]
    active_vars = variables_in_activated_equalities_set(blk)
    fixed_active_vars = fixed_variables_in_activated_equalities_set(blk)
    unfixed_active_vars = unfixed_variables_in_activated_equalities_set(blk)
    print("\n ===============DOF Breakdown================\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(blk)}")
    print(f"Activated Variables: ({len(active_vars)})")
    for v in active_vars:
        print(f"   {v}")
    print(f"Activated Equalities: ({len(equalities)})")
    for c in equalities:
        print(f"   {c}")

    print(f"Fixed Active Vars: ({len(fixed_active_vars)})")
    for v in fixed_active_vars:
        print(f"   {v}")

    print(f"Unfixed Active Vars: ({len(unfixed_active_vars)})")
    for v in unfixed_active_vars:
        print(f"   {v}")
    print("\n")
    print(f" {f' Active Vars':<30s}{len(active_vars)}")
    print(f"{'-'}{f' Fixed Active Vars':<30s}{len(fixed_active_vars)}")
    print(f"{'-'}{f' Activated Equalities':<30s}{len(equalities)}")
    print(f"{'='}{f' Degrees of Freedom':<30s}{degrees_of_freedom(blk)}")
    print("\nSuggested Variables to Fix:")

    if degrees_of_freedom != 0:
        unfixed_vars_without_constraint = [
            v for v in active_vars if v not in unfixed_active_vars
        ]
        for v in unfixed_vars_without_constraint:
            if v.fixed is False:
                print(f"   {v}")


def report_LTMED(m):
    blk = m.fs.treatment.LTMED
    print(f"\n\n-------------------- LTMED Report --------------------\n")
    print(f'{f"Stream":<20}{f"FLOW RATE H2O":<20}{f"FLOW RATE TDS":<20}')
    print(
        f'{"FEED":<20}{value(blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<20.2f}{value(blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]):<20.2f} kg/s'
    )
    print(
        f'{"PRODUCT":<20}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<20.2f}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "TDS"]):<20.2f} kg/s'
    )
    print(
        f'{"DISPOSAL":<20}{value(blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"]):<20.2f}{value(blk.disposal.properties[0].flow_mass_phase_comp["Liq", "TDS"]):<20.2f} kg/s'
    )
    print("")
    # print(f'Water Recovery: {value(blk.unit.water_recovery[0]):.2f}')
    print(
        f'{"STEC":<20}{value(blk.unit.specific_energy_consumption_thermal):<20.2f}kWh/m3'
    )
    print(
        f'{"GOR":<20}{value(blk.unit.gain_output_ratio):<20.2f}{pyunits.get_units(blk.unit.gain_output_ratio)}'
    )
    print(
        f'{"Thermal Power Req.":<20}{value(blk.unit.thermal_power_requirement):<20,.2f}{pyunits.get_units(blk.unit.thermal_power_requirement)}'
    )
    print(
        f'{"Water Recovery":<20}{value(blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]/blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]) * 100:<20.2f}%'
    )

    # print(
    #     m.fs.treatment.LTMED.product.properties[0.0].flow_mass_phase_comp.display()
    # )  # TDS
    # print(
    #     m.fs.treatment.LTMED.unit.distillate_props[0.0].mass_frac_phase_comp.display()
    # )  # TDS
    # print(
    #     m.fs.treatment.LTMED.unit.distillate_props[0.0].flow_vol_phase.display()
    # )  # TDS
    # print(m.fs.treatment.LTMED.unit.brine_props[0.0].flow_vol_phase.display())


if __name__ == "__main__":

    m = build_system()
    set_system_operating_conditions(m)
    set_LTMED_operating_conditions(m.fs.treatment.LTMED)
    add_system_costing(m)
    add_LTMED_scaling(m, m.fs.treatment.LTMED)
    # add_LTMED_costing(m, m.fs.treatment.LTMED)
    init_system(m)
    solver = get_solver()
    results = solve(m, debug=False)
    optimize(m, water_recovery=0.4)
    results = solve(m, debug=False)
    report_LTMED(m)
    # print(m.fs.treatment.LTMED.unit.costing.display())
