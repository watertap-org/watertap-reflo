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
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core import Database
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.unit_models.zero_order import ChemicalAdditionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.util.initialization import *

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

reflo_dir = pathlib.Path(__file__).resolve().parents[4]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3


__all__ = [
    "build_chem_addition",
    "set_chem_addition_op_conditions",
    "set_chem_addition_scaling",
    "add_chem_addition_costing",
    "init_chem_addition",
]


def build_system():
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.chem_addition = FlowsheetBlock(dynamic=False)

    build_chem_addition(m, m.fs.chem_addition, m.fs.properties)

    m.fs.feed_to_chem = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.chem_addition.feed.inlet,
    )

    m.fs.chem_to_product = Arc(
        source=m.fs.chem_addition.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_chem_addition(m, blk, prop_package=None) -> None:

    print(f'\n{"=======> BUILDING CHEMICAL ADDITION SYSTEM <=======":^60}\n')

    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)

    blk.unit = ChemicalAdditionZO(
        property_package=prop_package,
        database=m.db,
        process_subtype="hydrogen_peroxide",
    )
    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_system_operating_conditions(m, Qin=5, tds=130):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    Qin = Qin * pyunits.Mgal / pyunits.day
    tds = tds * pyunits.kg / pyunits.m**3

    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(flow_mass_tds)
    m.fs.chem_addition.unit.properties[0].flow_mass_comp["tds"].set_value(flow_mass_tds)


def set_chem_addition_scaling(m, blk, calc_blk_scaling_factors=False):

    set_scaling_factor(blk.unit.chemical_dosage, 0.1)
    set_scaling_factor(blk.unit.solution_density, 1e-3)
    set_scaling_factor(blk.unit.chemical_flow_vol, 1e6)
    set_scaling_factor(blk.unit.electricity, 1e4)

    # Calculate scaling factors only for chem addition block if in full case study flowsheet
    # so we don't prematurely set scaling factors
    if calc_blk_scaling_factors:
        calculate_scaling_factors(blk)

    # otherwise calculate all scaling factors
    else:
        calculate_scaling_factors(m)


def set_chem_addition_op_conditions(m, blk, **kwargs):

    blk.unit.load_parameters_from_database()
    for k, v in kwargs.items():
        print(blk.name, k)
        try:
            vv = getattr(blk.unit, k)
        except:
            continue
        print(f"{blk.name} {vv.name} {v}")
        if isinstance(vv, Var):
            vv.fix(v)
        if isinstance(vv, Param):
            vv.set_value(v)

    print(f"Chem Addition")
    print(f"\tblock DOF = {degrees_of_freedom(blk)}\n")
    print(f"\tunit DOF = {degrees_of_freedom(blk.unit)}\n")


def add_chem_addition_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )


def init_system(blk, solver=None, flow_out=None):
    if solver is None:
        solver = get_solver()
    if flow_out is None:
        flow_out = blk.fs.feed.properties[0].flow_vol

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Chem Addition Degrees of Freedom: {degrees_of_freedom(m.fs.chem_addition)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_chem)
    init_chem_addition(m, m.fs.chem_addition)
    propagate_state(m.fs.chem_to_product)
    m.fs.product.initialize()

    m.fs.costing.cost_process()
    m.fs.costing.initialize()
    m.fs.costing.add_LCOW(flow_out)


def init_chem_addition(m, blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize()
    propagate_state(blk.unit_to_product)
    print(f"Chem Addition")
    print(f"\tblock DOF after init = {degrees_of_freedom(blk)}\n")
    print(f"\tunit DOF after init = {degrees_of_freedom(blk.unit)}\n")


def print_chem_addition_costing_breakdown(blk):
    print(
        f'{"Hydrogen Peroxide Dose":<35s}{f"{blk.unit.chemical_dosage[0]():<25,.0f} mg/L"}'
    )
    print(
        f'{"Chem Addition Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}'
    )


def solve(m, solver=None, tee=True, raise_on_failure=True):

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
    set_chem_addition_op_conditions(m, m.fs.chem_addition)
    add_chem_addition_costing(m, m.fs.chem_addition)
    set_chem_addition_scaling(m, m.fs.chem_addition, calc_blk_scaling_factors=True)
    init_system(m)
    solve(m)
    print_chem_addition_costing_breakdown(m.fs.chem_addition)
    print(f"LCOW = {m.fs.costing.LCOW()}")
