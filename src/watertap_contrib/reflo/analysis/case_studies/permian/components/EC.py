import pathlib
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
from watertap.core.solvers import get_solver

from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
import math
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.core import REFLODatabase

rho = 1000 * pyunits.kg / pyunits.m**3
reflo_dir = pathlib.Path(__file__).resolve().parents[4]

case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"

__all__ = [
    "build_ec",
    "set_ec_operating_conditions", 
    "set_ec_scaling",
    "init_ec",
    "add_ec_costing",

]


def propagate_state(arc):
    _prop_state(arc)


def build_ec(m, blk, prop_package=None):
    """Function to build EC unit model"""

    print(f'\n{"=======> BUILDING EC SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = ElectrocoagulationZO(
        property_package=prop_package,
        database=m.db,
        electrode_material="aluminum",
        reactor_material="pvc",
        overpotential_calculation="calculated",
        process_subtype="permian", 
    )

    # print(blk.ec.display())
    # print(blk.ec.report())
    # print(blk.ec.inlet)
    # print(blk.ec.treated)
    # print(blk.ec.byproduct)
    # assert False

    blk.feed_to_ec = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.ec_to_product = Arc(
        source=blk.unit.treated,
        destination=blk.product.inlet,
    )

    blk.ec_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def build_system():
    """Function to create concrete model for individual unit model flowsheet"""
    m = ConcreteModel()
    m.db = REFLODatabase()

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


def set_system_operating_conditions(m, Qin=5, tds=130):
    """This function sets the system operating conditions for individual unit model flowsheet"""


    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3/pyunits.s)
    m.tds = tds

    # TODO: rho should probably be higher for 130 g/L
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    inlet_dict = {"tds": tds * pyunits.kg / pyunits.m**3}
    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)
    m.fs.feed.properties[0].conc_mass_comp

    for solute, solute_conc in inlet_dict.items():
        flow_mass_solute = pyunits.convert(
            Qin * solute_conc, to_units=pyunits.kg / pyunits.s
        )
        sf = 1 / value(flow_mass_solute)
        m.fs.feed.properties[0].flow_mass_comp[solute].fix(flow_mass_solute)
        m.fs.EC.unit.properties_in[0].flow_mass_comp[solute].set_value(flow_mass_solute)
        m.fs.properties.set_default_scaling(
            "flow_mass_comp",
            sf,
            index=(solute),
        )
        m.fs.properties.set_default_scaling(
            "conc_mass_comp",
            1 / solute_conc(),
            index=(solute),
        )

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        1 / value(flow_mass_water),
        index=("H2O"),
    )
    calculate_scaling_factors(m)


def set_ec_operating_conditions(m, blk):
    """Set EC operating conditions"""
    # Check if the set up of the ec inputs is correct

    blk.unit.load_parameters_from_database(use_default_removal=True)

    conv = 5e3 * (pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S)

    cond = pyunits.convert(
        m.tds * pyunits.gram / pyunits.liter / conv,
        to_units=pyunits.S / pyunits.m,
    )
    print(f"cond = {cond()}")
    blk.unit.conductivity.fix(cond)




def set_ec_scaling(m, blk):

    calculate_scaling_factors(m)

    set_scaling_factor(blk.unit.properties_in[0].flow_vol, 1e7)
    set_scaling_factor(blk.unit.properties_in[0].conc_mass_comp["tds"], 1e5)
    set_scaling_factor(blk.unit.charge_loading_rate, 1e3)
    # set_scaling_factor(m.fs.ec.cell_voltage,1)
    set_scaling_factor(blk.unit.anode_area, 1e4)
    set_scaling_factor(blk.unit.current_density, 1e-1)
    # set_scaling_factor(m.fs.ec.applied_current,1e2)
    set_scaling_factor(blk.unit.metal_dose, 1e4)


def init_system(m, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"EC Degrees of Freedom: {degrees_of_freedom(m.fs.EC.unit)}")
    assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_unit)

    init_ec(m, m.fs.EC)


def init_ec(m, blk, solver=None):
    """Initialize EC model"""

    if solver is None:
        solver = get_solver()

    optarg = solver.options

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_ec)

    cvc(blk.unit.overpotential, blk.unit.eq_overpotential)
    cvc(blk.unit.applied_current, blk.unit.eq_applied_current)
    cvc(blk.unit.anode_area, blk.unit.eq_electrode_area_total)
    cvc(blk.unit.ohmic_resistance, blk.unit.eq_ohmic_resistance)

    blk.unit.initialize(optarg=optarg)
    propagate_state(blk.ec_to_product)
    propagate_state(blk.ec_to_disposal)


def add_system_costing(m):
    """Add system level costing components"""
    m.fs.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    # m.fs.costing.electricity_cost.fix(0.07)
    add_ec_costing(m, m.fs.EC)
    calc_costing(m, m.fs.EC)


def add_ec_costing(m, blk):
    """Add EC model costing components"""
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)


def calc_costing(m, blk):
    """Add system level solve for costing"""
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(blk.unit.properties_treated[0].flow_vol)
    m.fs.costing.add_electricity_intensity(blk.unit.properties_treated[0].flow_vol)


def report_EC(blk):
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
    set_system_operating_conditions(m, tds=130)
    set_ec_operating_conditions(m, m.fs.EC)
    set_ec_scaling(m, m.fs.EC)
    init_system(m)
    add_system_costing(m)

    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"LCOW = {m.fs.costing.LCOW()}")
    # m.fs.EC.unit.display()
    # ec = m.fs.EC.unit
    # ec.display()
    # ec.properties_in[0].display()
    # ec.electrolysis_time.display()


