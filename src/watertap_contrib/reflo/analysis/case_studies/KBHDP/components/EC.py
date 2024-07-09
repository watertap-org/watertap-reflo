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
from watertap.unit_models.zero_order import ElectrocoagulationZO
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


def build_ec(m, prop_package = None):
    '''Function to build EC unit model'''

    print(f'\n{"=======> BUILDING EC SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties
    
    m.fs.unit.feed = StateJunction(property_package=prop_package)

    m.fs.ec =  ElectrocoagulationZO(
    property_package = prop_package,
    database = m.db,
    electrode_material = 'aluminum',
    reactor_material = "pvc",
    overpotential_calculation = 'calculated',
    )
    
    
    m.fs.feed_to_ec = Arc(
        source=m.fs.unit.feed.outlet,
        destination=m.fs.ec.inlet,
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

    build_ec(m)

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



def set_ec_operating_conditions(m):
    '''Set EC operating conditions''' 
    # Check if the set up of the ec inputs is correct
    
    input = {'gap (cm)': 0.32,
             'thickness (cm)': 0.32,
             'ret_time (s)': 23,
             'dose (mg/L)': 2.95,
             'anode_area (cm2)': 184,
             'cd (A/m2)': 218,
             }

    gap = pyunits.convert(input['gap (cm)'] * pyunits.cm, to_units=pyunits.m)()
    e_thick = pyunits.convert(input['thickness (cm)']* pyunits.cm, to_units=pyunits.m)()
    time = pyunits.convert(input['ret_time (s)'] * pyunits.seconds, to_units=pyunits.minutes)()
    
    conv = 5e3 * (pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S)
    tds = m.fs.unit.feed.properties[0].flow_mass_comp["tds"]/(m.fs.unit.feed.properties[0].flow_mass_comp["H2O"]/(1000 * pyunits.kg / pyunits.m**3))
    # kg/s / (kg/s*m3/kg)
    cond = pyunits.convert(pyunits.convert(tds, to_units=pyunits.mg / pyunits.liter)/conv,
                            to_units  = pyunits.S / pyunits.m)
    m.fs.ec.conductivity.fix(cond)

    ec_dose = input['dose (mg/L)']* pyunits.mg / pyunits.liter
    ec_dose = pyunits.convert(
        input['dose (mg/L)']* pyunits.mg / pyunits.liter, to_units=pyunits.kg / pyunits.liter
    )
    
    anode_area = pyunits.convert(input['anode_area (cm2)']*pyunits.cm**2,to_units=pyunits.m**2)

    m.fs.ec.load_parameters_from_database(use_default_removal=True)

    m.fs.ec.electrode_thick.fix(e_thick)
    m.fs.ec.electrode_gap.fix(gap)

    m.fs.ec.current_density.fix(input['cd (A/m2)'])
    m.fs.ec.metal_dose.fix(ec_dose) 
    
    if m.fs.ec.config.electrode_material == 'aluminum':
        m.fs.ec.current_efficiency.fix(1)

    m.fs.ec.overpotential_k1.fix(430)
    m.fs.ec.overpotential_k2.fix(1000)

    
def set_scaling(m):

    calculate_scaling_factors(m)

    set_scaling_factor(m.fs.ec.properties_in[0].flow_vol, 1e7)
    set_scaling_factor(m.fs.ec.properties_in[0].conc_mass_comp['tds'], 1e5)
    set_scaling_factor(m.fs.ec.charge_loading_rate,1e3)
    # set_scaling_factor(m.fs.ec.cell_voltage,1)
    set_scaling_factor(m.fs.ec.anode_area,1e4)
    set_scaling_factor(m.fs.ec.current_density,1e-1)
    # set_scaling_factor(m.fs.ec.applied_current,1e2)
    set_scaling_factor(m.fs.ec.metal_dose,1e4)


def init_system(m, solver=None):
    '''Initialize system for individual unit process flowsheet'''
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"EC Degrees of Freedom: {degrees_of_freedom(m.fs.ec)}")
    assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_ec_model(m)


def init_ec_model(m, solver=None):
    '''Initialize IX model'''

    if solver is None:
        solver = get_solver()

    optarg = solver.options

    m.fs.unit.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_ec)

    m.fs.ec.initialize(optarg=optarg)


def add_system_costing(m):
    '''Add system level costing components'''
    m.fs.costing = ZeroOrderCosting()
    add_ec_costing(m)
    calc_costing(m)


def add_ec_costing(m):
    '''Add EC model costing components'''
    m.fs.ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)


def calc_costing(m):
    '''Add system level solve for costing'''
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.ec.properties_treated[0].flow_vol)
    m.fs.costing.add_electricity_intensity(m.fs.ec.properties_treated[0].flow_vol)


if __name__ == "__main__":
    
    m = build_system()
    set_system_operating_conditions(m)
    set_ec_operating_conditions(m)
    set_scaling(m)

    init_system(m)

    solver = get_solver()  
    results = solver.solve(m)

    add_system_costing(m)

    m.fs.objective_lcow = Objective(expr = m.fs.costing.LCOW)
    results = solver.solve(m)

    print(m.fs.objective_lcow())