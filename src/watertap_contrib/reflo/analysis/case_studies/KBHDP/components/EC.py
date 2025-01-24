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
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver

from watertap_contrib.reflo.core import REFLODatabase
from watertap.core.zero_order_properties import (
    WaterParameterBlock as WaterParameterBlockZO,
)
from watertap.core.zero_order_properties import WaterParameterBlock
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
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    check_jac,
    calc_scale,
)

__all__ = [
    "build_ec",
    "set_ec_operating_conditions",
    "init_ec",
    "add_ec_costing",
    "add_ec_scaling",
    "report_EC",
    "print_EC_costing_breakdown",
]

__all__ = [
    "build_ec",
    "set_ec_operating_conditions",
    "init_ec",
    "add_ec_costing",
    "add_ec_scaling",
    "report_EC",
    "print_EC_costing_breakdown",
]


def propagate_state(arc, detailed=False):
    _prop_state(arc)
    if detailed:
        print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print("\n")


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


def build_system():
    """Function to create concrete model for individual unit model flowsheet"""
    m = ConcreteModel()
    m.db = REFLODatabase()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.sludge = Product(property_package=m.fs.properties)

    m.fs.EC = FlowsheetBlock(dynamic=False)

    build_ec(m, m.fs.EC)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.EC.feed.inlet,
    )
    m.fs.ec_to_sludge = Arc(
        source=m.fs.EC.disposal.outlet,
        destination=m.fs.sludge.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


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
        process_subtype="kbhdp",
    )

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


def set_system_operating_conditions(m):
    """This function sets the system operating conditions for individual unit model flowsheet"""

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(175.25054)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(2.143156)  # kg/m3 * m3/s = kg/s
    m.fs.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)
    # # initialize feed


def set_ec_operating_conditions(m, blk, conv=5e3):
    """Set EC operating conditions"""
    # Check if the set up of the ec inputs is correct
    print(f"EC Degrees of Freedom: {degrees_of_freedom(blk.ec)}")

    blk.ec.load_parameters_from_database(use_default_removal=True)
    # blk.ec.conductivity.unfix()
    # tds_ec_conversion = conv * (pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S)
    # blk.ec.conductivity_constr = Constraint(
    #     expr=blk.ec.conductivity
    #     == pyunits.convert(
    #         blk.feed.properties[0].conc_mass_comp["tds"] / tds_ec_conversion,
    #         to_units=pyunits.S / pyunits.m,
    #     )
    # )
    # blk.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)
    # blk.ec.overpotential.fix(2)
    print(f"EC Degrees of Freedom: {degrees_of_freedom(blk.ec)}")


def set_scaling(m, blk):

    def calc_scale(value):
        return math.floor(math.log(value, 10))

    scale_flow = calc_scale(m.fs.feed.flow_mass_comp[0, "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_comp[0, "tds"].value)
    scale_tss = calc_scale(m.fs.feed.flow_mass_comp[0, "tss"].value)

    m.fs.properties.set_default_scaling(
        "flow_mass_comp", 10**-scale_flow, index=("H2O")
    )
    m.fs.properties.set_default_scaling("flow_mass_comp", 10**-scale_tds, index=("tds"))
    m.fs.properties.set_default_scaling("flow_mass_comp", 10**-scale_tss, index=("tss"))
    calculate_scaling_factors(m)


def add_ec_scaling(m, blk):
    set_scaling_factor(blk.ec.charge_loading_rate, 1e1)
    set_scaling_factor(blk.ec.electrode_volume, 1)
    # set_scaling_factor(blk.ec.cell_voltage, )
    # set_scaling_factor(blk.ec.applied_current, 1e5)
    # set_scaling_factor(blk.ec.power_required, 1e-6)

    constraint_scaling_transform(blk.ec.eq_power_required, 1e-4)
    set_scaling_factor(blk.ec.properties_in[0.0].flow_mass_comp["H2O"], 1e1)
    set_scaling_factor(blk.ec.properties_in[0.0].flow_mass_comp["tds"], 1)
    set_scaling_factor(blk.ec.properties_treated[0.0].flow_mass_comp["H2O"], 1e1)
    set_scaling_factor(blk.ec.properties_byproduct[0.0].flow_mass_comp["H2O"], 1)


def init_system(m, solver=None):
    """Initialize system for individual unit process flowsheet"""
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"EC Degrees of Freedom: {degrees_of_freedom(m.fs.EC.ec)}")

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

    cvc(blk.ec.overpotential, blk.ec.eq_overpotential)
    cvc(blk.ec.applied_current, blk.ec.eq_applied_current)
    cvc(blk.ec.anode_area, blk.ec.eq_electrode_area_total)
    cvc(blk.ec.ohmic_resistance, blk.ec.eq_ohmic_resistance)

    blk.ec.initialize(optarg=optarg)
    propagate_state(blk.ec_to_product)
    propagate_state(blk.ec_to_disposal)


# def init_ec(m, blk, solver=None):
#     """Initialize IX model"""

#     if solver is None:
#         solver = get_solver()

#     optarg = solver.options
#     # assert_no_degrees_of_freedom(m)
#     _initialize(blk.feed)
#     propagate_state(blk.feed_to_ec)

#     _initialize(blk.ec)

#     propagate_state(blk.ec_to_product)
#     propagate_state(blk.ec_to_disposal)


def add_system_costing(m):
    """Add system level costing components"""
    m.fs.costing = ZeroOrderCosting()
    add_ec_costing(m, m.fs.EC)
    calc_costing(m, m.fs.EC)


def add_ec_costing(m, blk, costing_blk=None):
    """Add EC model costing components"""
    if costing_blk is None:
        costing_blk = m.fs.costing

    blk.ec.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_blk)


def calc_costing(m, blk):
    """Add system level solve for costing"""
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(blk.ec.properties_treated[0].flow_vol)
    m.fs.costing.add_electricity_intensity(blk.ec.properties_treated[0].flow_vol)


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
            check_jac(m)

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


def print_EC_costing_breakdown(blk):
    print(f"\n\n-------------------- EC Costing Breakdown --------------------\n")
    print(f'{"EC Capital Cost":<35s}{f"${blk.ec.costing.capital_cost():<25,.0f}"}')
    print(
        f'{"EC Operating Cost":<35s}{f"${blk.ec.costing.fixed_operating_cost():<25,.0f}"}'
    )


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


if __name__ == "__main__":

    m = build_system()
    set_system_operating_conditions(m)
    set_ec_operating_conditions(m, m.fs.EC)
    add_ec_scaling(m, m.fs.EC)
    set_scaling(m, m.fs.EC)
    init_system(m)
    add_system_costing(m)

    solver = get_solver()
    m.fs.objective_lcow = Objective(expr=m.fs.costing.LCOW)
    results = solve(m, debug=True)

    report_EC(m.fs.EC)
    print_EC_costing_breakdown(m.fs.EC)
    m.fs.EC.ec.conductivity.display()
