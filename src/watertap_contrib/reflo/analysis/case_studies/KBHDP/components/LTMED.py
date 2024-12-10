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

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.reflo.unit_models.surrogate import LTMEDSurrogate


def propagate_state(arc):
    _prop_state(arc)
    print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    arc.source.display()
    print(arc.destination.name)
    arc.destination.display()
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


def set_LTMED_operating_conditions(m, blk):
    steam_temperature = 80
    recovery_ratio = 0.5

    blk.unit.steam_props[0].temperature.fix(steam_temperature + 273.15)
    blk.unit.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    pass


def init_LTMED(m, blk, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING LTMED --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"LTMED Degrees of Freedom: {degrees_of_freedom(blk)}")
    # assert_no_degrees_of_freedom(m)

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_LTMED)

    blk.unit.initialize(optarg=optarg)
    propagate_state(blk.LTMED_to_product)
    propagate_state(blk.LTMED_to_disposal)


def set_system_operating_conditions(m):
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1000)
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(2)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(25 + 273.15)


# The following functions are used for testing the component in isolation on thise file
def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.liquid_prop)

    m.fs.LTMED = FlowsheetBlock(dynamic=False)
    build_LTMED(m, m.fs.LTMED, m.fs.liquid_prop, m.fs.vapor_prop)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.LTMED.feed.inlet,
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
    print(f"LTMED Degrees of Freedom: {degrees_of_freedom(m.fs.LTMED.unit)}")
    assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_LTMED(m, m.fs.LTMED)


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


if __name__ == "__main__":

    m = build_system()
    set_system_operating_conditions(m)
    set_LTMED_operating_conditions(m, m.fs.LTMED)

    init_system(m)
    solver = get_solver()
    results = solve(m)
