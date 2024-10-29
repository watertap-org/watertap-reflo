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
    "add_chem_addition_costing",
    "init_chem_addition"
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
        m.fs.chem_addition.unit.properties[0].flow_mass_comp[solute].set_value(
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


def set_chem_addition_op_conditions(m, blk):

    blk.unit.load_parameters_from_database()


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
    set_chem_addition_op_conditions(m, m.fs.chem_addition.unit)

    add_chem_addition_costing(m, m.fs.chem_addition)
    init_system(m)
    solve(m)
    print_chem_addition_costing_breakdown(m.fs.chem_addition)
    print(f"LCOW = {m.fs.costing.LCOW()}")
