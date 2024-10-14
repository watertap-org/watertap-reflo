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

__all__ = [
    "build_ec",
    "set_ec_operating_conditions",
    "init_ec",
    "report_EC",
]

def propagate_state(arc, detailed=True):
    _prop_state(arc)
    if detailed:
        print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print("\n")


def build_ec(m, blk, prop_package=None):
    """Function to build EC unit model"""

    print(f'\n{"=======> BUILDING EC SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.ec = ElectrocoagulationZO(
        property_package=prop_package,
        database=m.db,
        electrode_material="aluminum",
        reactor_material="pvc",
        overpotential_calculation="calculated",
    )

    # print(blk.ec.display())
    # print(blk.ec.report())
    # print(blk.ec.inlet)
    # print(blk.ec.treated)
    # print(blk.ec.byproduct)
    # assert False

    blk.feed_to_ec = Arc(
        source=blk.feed.outlet,
        destination=blk.ec.inlet,
    )

    blk.ec_to_product = Arc(
        source=blk.ec.treated,
        destination=blk.product.inlet,
    )

    blk.ec_to_disposal = Arc(
        source=blk.ec.byproduct,
        destination=blk.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def build_system():
    """Function to create concrete model for individual unit model flowsheet"""
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = WaterParameterBlockZO(solute_list=["tds"])

    m.fs.feed = Feed(property_package=m.fs.properties)

    m.fs.EC = FlowsheetBlock(dynamic=False)

    build_ec(m, m.fs.EC)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.EC.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_system_operating_conditions(m):
    """This function sets the system operating conditions for individual unit model flowsheet"""

    input = {
        "q (m3/s)": 0.175,
        "tds (g/l)": 12.23,
    }

    flow_in = input["q (m3/s)"] * pyunits.m**3 / pyunits.s
    flow_in_mass = flow_in * (1000 * pyunits.kg / pyunits.m**3)  # kg/s

    tds = input["tds (g/l)"] * pyunits.g / pyunits.liter
    tds_in = pyunits.convert(tds, to_units=pyunits.kg / pyunits.m**3)

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(flow_in_mass)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(
        tds_in * flow_in
    )  # kg/m3 * m3/s = kg/s
    # # initialize feed

def set_ec_operating_conditions(m, blk):
    """Set EC operating conditions"""
    # Check if the set up of the ec inputs is correct
    print(f"EC Degrees of Freedom: {degrees_of_freedom(m.fs.EC.ec)}")
    input = {
        "gap (cm)": 0.5,
        "thickness (cm)": 0.1,
        "ret_time (min)": 25,
        "dose (mg/L)": 100,
        "anode_area (cm2)": 184,
        "cd (A/m2)": 500,
    }

    gap = pyunits.convert(input["gap (cm)"] * pyunits.cm, to_units=pyunits.m)()
    e_thick = pyunits.convert(
        input["thickness (cm)"] * pyunits.cm, to_units=pyunits.m
    )()

    time = input["ret_time (min)"] * pyunits.minutes

    conv = 5e3 * (pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S)
    tds = blk.feed.properties[0].flow_mass_comp["tds"] / (
        blk.feed.properties[0].flow_mass_comp["H2O"]
        / (1000 * pyunits.kg / pyunits.m**3)
    )
    # kg/s / (kg/s*m3/kg)
    cond = pyunits.convert(
        pyunits.convert(tds, to_units=pyunits.mg / pyunits.liter) / conv,
        to_units=pyunits.S / pyunits.m,
    )

    blk.ec.conductivity.fix(cond)

    ec_dose = input["dose (mg/L)"] * pyunits.mg / pyunits.liter
    ec_dose = pyunits.convert(
        input["dose (mg/L)"] * pyunits.mg / pyunits.liter,
        to_units=pyunits.kg / pyunits.liter,
    )

    anode_area = pyunits.convert(
        input["anode_area (cm2)"] * pyunits.cm**2, to_units=pyunits.m**2
    )

    blk.ec.load_parameters_from_database(use_default_removal=True)

    blk.ec.electrode_thick.fix(e_thick)
    blk.ec.electrode_gap.fix(gap)

    blk.ec.current_density.fix(input["cd (A/m2)"])
    blk.ec.metal_dose.fix(ec_dose)
    # blk.ec.metal_dose.unfix()

    if blk.ec.config.electrode_material == "aluminum":
        blk.ec.current_efficiency.fix(1.2)

    blk.ec.overpotential.fix(2)
    blk.ec.floc_retention_time.fix(time)
    blk.ec.overpotential_k1.unfix()
    blk.ec.overpotential_k2.unfix()

    fixed_vars = [(v.name, v.value) for v in blk.ec.component_data_objects(ctype=Var, active=True, descend_into=False) if v.fixed]
    unfixed_vars = [(v.name, v.value) for v in blk.ec.component_data_objects(ctype=Var, active=True, descend_into=False) if v.fixed is False]
    print(f"Fixed Vars: ({len(fixed_vars)})")
    for v in fixed_vars:
        print(f"   {v[0]}: {v[1]}")

    print(f"Unfixed Vars: ({len(unfixed_vars)})")
    for v in unfixed_vars:
        print(f"   {v[0]}: {v[1]}")



def set_scaling(m, blk):

    def calc_scale(value):
        return math.floor(math.log(value, 10))

    scale_flow = calc_scale(m.fs.feed.flow_mass_comp[0, "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_comp[0, "tds"].value)

    m.fs.properties.set_default_scaling(
        "flow_mass_comp", 10**-scale_flow, index=("H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_comp", 10**-scale_tds, index=("tds")
    )
    calculate_scaling_factors(m)

    badly_scaled_var_list = list_badly_scaled_variables(m)
    if len(badly_scaled_var_list) > 0:
        [print(i[0].name, i[1]) for i in badly_scaled_var_list]
    else:
        print("Variables are scaled well")


def init_system(m, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"EC Degrees of Freedom: {degrees_of_freedom(m.fs.EC.ec)}")
    # assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_ec(m, m.fs.EC)


def init_ec(m, blk, solver=None):
    """Initialize IX model"""

    if solver is None:
        solver = get_solver()

    optarg = solver.options

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_ec)

    try:
        blk.ec.initialize(optarg=optarg)
    except:
        blk.ec.display()
        assert False
    propagate_state(blk.ec_to_product)
    propagate_state(blk.ec_to_disposal)


def add_system_costing(m):
    """Add system level costing components"""
    m.fs.costing = ZeroOrderCosting()
    add_ec_costing(m, m.fs.EC)
    calc_costing(m, m.fs.EC)


def add_ec_costing(m, blk):
    """Add EC model costing components"""
    blk.ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)


def calc_costing(m, blk):
    """Add system level solve for costing"""
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(blk.ec.properties_treated[0].flow_vol)
    m.fs.costing.add_electricity_intensity(blk.ec.properties_treated[0].flow_vol)


def report_EC(blk):

    print(f"\n\n-------------------- EC Report --------------------\n")
    print(f'{f"Stream":<20}{f"FLOW RATE H2O":<20}{f"FLOW RATE TDS":<20}')
    print(
        f'{"FEED":<20}{value(blk.feed.properties[0].flow_mass_comp["H2O"]):<20.2f}{value(blk.feed.properties[0].flow_mass_comp["tds"]):<20.2f} kg/s'
    )
    print(
        f'{"PRODUCT":<20}{value(blk.product.properties[0].flow_mass_comp["H2O"]):<20.2f}{value(blk.product.properties[0].flow_mass_comp["tds"]):<20.2f} kg/s'
    )
    print(
        f'{"DISPOSAL":<20}{value(blk.disposal.properties[0].flow_mass_comp["H2O"]):<20.2f}{value(blk.disposal.properties[0].flow_mass_comp["tds"]):<20.2f} kg/s'
    )


if __name__ == "__main__":

    m = build_system()
    set_system_operating_conditions(m)
    set_ec_operating_conditions(m, m.fs.EC)
    set_scaling(m, m.fs.EC)
    init_system(m)
    add_system_costing(m)

    solver = get_solver()
    m.fs.objective_lcow = Objective(expr = m.fs.costing.LCOW)
    results = solver.solve(m)

    print(m.fs.objective_lcow())
    report_EC(m.fs.EC)
