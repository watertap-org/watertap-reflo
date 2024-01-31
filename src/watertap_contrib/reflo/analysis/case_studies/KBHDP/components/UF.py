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
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
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
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

from watertap.unit_models.zero_order.ultra_filtration_zo import UltraFiltrationZO


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_UF(m, blk) -> None:
    print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.RO_properties)
    blk.product = StateJunction(property_package=m.fs.RO_properties)
    blk.disposal = StateJunction(property_package=m.fs.RO_properties)

    # blk.unit = UltraFiltrationZO(property_package=m.fs.RO_properties)

    blk.feed_to_product = Arc(
        source=blk.feed.outlet,
        destination=blk.product.inlet,
    )


def init_UF(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING ULTRAFILTRATION --------------------\n\n"
    )
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Degasifier Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_product)
    blk.product.initialize(optarg=optarg)


def set_UF_op_conditions(blk):
    blk.unit.recovery_frac_mass_H2O.fix(1)
