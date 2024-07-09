from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
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
    assert_optimal_termination,
    units as pyunits,
)

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver

from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import (
    WaterParameterBlock as WaterParameterBlockZO,
)

from idaes.core.util.initialization import propagate_state as _prop_state
from watertap.unit_models.zero_order import IonExchangeZO
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics import *
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from idaes.core.util.constants import Constants
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from watertap.costing import WaterTAPCosting
import math


def propagate_state(arc):
    _prop_state(arc)


def build_ix(m, prop_package=None):
    '''Function to build IX unit model'''
    
    print(f'\n{"=======> BUILDING EC SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    m.fs.unit.feed = StateJunction(property_package=prop_package)
    m.fs.ix = IonExchangeZO(property_package=m.fs.properties, database=m.db)

    m.fs.feed_to_ix = Arc(
        source=m.fs.unit.feed.outlet,
        destination=m.fs.ix.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def build_system():
    '''Function to create concrete model for individual unit model flowsheet'''
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = WaterParameterBlockZO(
        solute_list=[
        "tds"
        ]
    )
    
    m.fs.feed = Feed(property_package=m.fs.properties)

    m.fs.unit = FlowsheetBlock(dynamic = False)

    build_ix(m)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.unit.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m

def set_system_operating_conditions(m):
    '''This function sets the system operating conditions for individual unit model flowsheet'''
    input = {'q (m3/s)': 0.175,
             'tds (g/l)': 12,
             }
    
    flow_in = input['q (m3/s)'] * pyunits.m**3 / pyunits.s
    flow_in_mass = flow_in * (1000 * pyunits.kg / pyunits.m**3) #kg/s

    tds = input['tds (g/l)']*pyunits.g / pyunits.liter
    tds_in = pyunits.convert(tds, to_units=pyunits.kg / pyunits.m**3)

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(flow_in_mass)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(tds_in * flow_in) #kg/m3 * m3/s = kg/s


def set_ix_operating_conditions(m):
    '''Set EC operating conditions''' 

    m.fs.ix.load_parameters_from_database(use_default_removal=True)
    m.fs.ix.recovery_frac_mass_H2O.fix(1)
    m.fs.ix.removal_frac_mass_comp.fix(0.1)
    m.fs.ix.NaCl_dose.fix()
    m.fs.ix.resin_replacement.fix()

def init_system(m, solver=None):
    '''Initialize system for individual unit process flowsheet'''
    
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"IX Degrees of Freedom: {degrees_of_freedom(m.fs.ix)}")
    assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_ix_model(m)


def init_ix_model(m, solver=None):
    '''Initialize IX model'''

    if solver is None:
        solver = get_solver()

    optarg = solver.options

    m.fs.unit.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_ix)

    m.fs.ix.initialize(optarg=optarg)

def add_system_costing(m):
    '''Add system level costing components'''
    m.fs.costing = ZeroOrderCosting()
    add_ix_costing(m)
    calc_costing(m)


def add_ix_costing(m):
    '''Add IX model costing components'''
    m.fs.ix.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)


def calc_costing(m):
    '''Add system level solve for costing'''
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.ix.properties_treated[0].flow_vol)
    m.fs.costing.add_electricity_intensity(m.fs.ix.properties_treated[0].flow_vol)

if __name__ == "__main__":
    
    m = build_system()
    set_system_operating_conditions(m)
    set_ix_operating_conditions(m)
    # set_scaling(m)

    init_system(m)
    
    solver = get_solver()  
    results = solver.solve(m)

    add_system_costing(m)
    m.fs.objective_lcow = Objective(expr = m.fs.costing.LCOW)
    results = solver.solve(m)

    print(m.fs.objective_lcow())