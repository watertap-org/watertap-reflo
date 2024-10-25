import pathlib
import os
import math
import numpy as np
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
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

# from idaes.core.solvers import get_solver
from watertap.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.core import Database
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.unit_models.zero_order import CartridgeFiltrationZO
from watertap.core.zero_order_properties import WaterParameterBlock as ZO

# from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

reflo_dir = pathlib.Path(__file__).resolve().parents[4]

case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"

rho = 1000 * pyunits.kg / pyunits.m**3


__all__ = [
    "build_cartridge_filtration", 
    "set_cart_filt_op_conditions",
    "add_cartridge_filtration_costing"

]

def build_system():
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, m.fs.cart_filt, prop_package=m.fs.properties)

    m.fs.feed_to_cart = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.cart_filt.feed.inlet,
    )

    m.fs.cart_to_product = Arc(
        source=m.fs.cart_filt.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.cart_to_disposal = Arc(
        source=m.fs.cart_filt.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_cartridge_filtration(m, blk, prop_package=None) -> None:
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = CartridgeFiltrationZO(
        property_package=prop_package,
        database=m.db,
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.treated,
        destination=blk.product.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_system_operating_conditions(m, Qin=5):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    Qin = Qin * pyunits.Mgal / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    inlet_dict = {"tds": 130 * pyunits.kg / pyunits.m**3}
    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)

    for solute, solute_conc in inlet_dict.items():
        flow_mass_solute = pyunits.convert(
            flow_in * solute_conc, to_units=pyunits.kg / pyunits.s
        )
        sf = 1 / value(flow_mass_solute)
        m.fs.feed.properties[0].flow_mass_comp[solute].fix(flow_mass_solute)
        m.fs.cart_filt.unit.properties_in[0].flow_mass_comp[solute].set_value(
            flow_mass_solute
        )
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


def set_cart_filt_op_conditions(m, cart_filt):

    # data = m.db.get_unit_operation_parameters("chemical_addition")
    cart_filt.load_parameters_from_database()


def add_cartridge_filtration_costing(m, blk):
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)


def init_system(blk, solver=None, flow_out=None):
    if solver is None:
        solver = get_solver()
    if flow_out is None:
        flow_out = blk.fs.feed.properties[0].flow_vol

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Cartridge Filt Degrees of Freedom: {degrees_of_freedom(m.fs.cart_filt)}")

    m.fs.feed.initialize()
    _prop_state(m.fs.feed_to_cart)
    # m.fs.cart_filt.unit.initialize()
    init_cart_filt(m, m.fs.cart_filt)
    _prop_state(m.fs.cart_to_product)
    m.fs.product.initialize()

    m.fs.costing.cost_process()
    m.fs.costing.initialize()
    m.fs.costing.add_LCOW(flow_out)


def init_cart_filt(m, blk):
    blk.unit.initialize()


def print_cart_filt_costing_breakdown(blk):

    print(
        f'{"Chem Addition Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}'
    )
    # print(
    #     f'{"Chem Addition Operating Cost":<35s}{f"${blk.unit.costing.fixed_operating_cost():<25,.0f}"}'
    # )


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
        raise RuntimeError(msg)
    else:
        return results


if __name__ == "__main__":

    m = build_system()
    set_system_operating_conditions(m)
    set_cart_filt_op_conditions(m, m.fs.cart_filt.unit)

    add_cartridge_filtration_costing(m, m.fs.cart_filt)
    init_system(m)
    solve(m)
    print(f"LCOW = {m.fs.costing.LCOW()}")
