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

from watertap.unit_models.zero_order.ion_exchange_zo import IonExchangeZO

def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_softener(m, blk) -> None:
    print(f'\n{"=======> BUILDING SOFTENER SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)

    # blk.unit = IonExchangeZO(property_package=m.fs.properties, 
    #                          database=m.db,
    #                          process_subtype="anion_exchange")

    blk.feed_to_product = Arc(
        source=blk.feed.outlet,
        destination=blk.product.inlet,
    )

    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp

def set_softener_operating_conditions(blk):
    # anion exchange
    blk.unit.load_parameters_from_database(use_default_removal=True)
    blk.unit.removal_frac_mass_comp[0, "tds"].fix(0.9)

def init_softener(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING WATER SOFTENER --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    print('\n\n')

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_product)
    blk.product.initialize(optarg=optarg)

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    
    m.fs.softener = FlowsheetBlock(dynamic=False)
    
    build_softener(m,m.fs.softener)
    init_softener(m,m.fs.softener)