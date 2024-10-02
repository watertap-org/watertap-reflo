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
    TransformationFactory,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
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
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.zero_order.ultra_filtration_zo import UltraFiltrationZO
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.wt_database import Database

__all__ = ["build_UF", "init_UF", "set_UF_op_conditions", "add_UF_costing", "report_UF", "print_UF_costing_breakdown"]


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_UF(m, blk, prop_package) -> None:
    print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)
    blk.unit = UltraFiltrationZO(property_package=prop_package, database=m.db)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.treated,
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
    print(f"UF Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize(optarg=optarg)
    propagate_state(blk.unit_to_disposal)
    propagate_state(blk.unit_to_product)
    blk.product.initialize(optarg=optarg)
    blk.disposal.initialize(optarg=optarg)


def set_UF_op_conditions(blk):
    blk.unit.recovery_frac_mass_H2O.fix(0.99)
    blk.unit.removal_frac_mass_comp[0, "tds"].fix(1e-3)
    blk.unit.removal_frac_mass_comp[0, "tss"].fix(0.9)
    blk.unit.energy_electric_flow_vol_inlet.fix(0.05)


def set_system_conditions(blk):
    blk.feed.properties[0.0].flow_mass_comp["H2O"].fix(171.37)
    blk.feed.properties[0.0].flow_mass_comp["tds"].fix(1.96)
    blk.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22E-6)


def add_UF_costing(m, blk):
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    # m.fs.costing.cost_process()


def load_parameters(m, blk):
    m.db.get_unit_operation_parameters("ultra_filtration")
    blk.unit.load_parameters_from_database()


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    # m.fs.costing = ZeroOrderCosting()
    m.db = Database()
    m.fs.RO_properties = NaClParameterBlock()
    m.fs.MCAS_properties = MCASParameterBlock(
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Cl_-",
            "Mg_2+",
            "K_+",
            "SiO2",
            "Na_+",
            "SO2_-4+",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.params = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.UF = FlowsheetBlock(dynamic=False)
    build_UF(m, m.fs.UF, m.fs.params)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

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


def print_stream_table(blk):
    print(f'{"FEED":<20s}')
    print(
        f'{"    H2O":<20s}{m.fs.UF.feed.properties[0.0].flow_mass_comp["H2O"].value:<10.3f}{pyunits.get_units(m.fs.UF.feed.properties[0.0].flow_mass_comp["H2O"])}'
    )
    print(
        f'{"    TDS":<20s}{m.fs.UF.feed.properties[0.0].flow_mass_comp["tds"].value:<10.3f}{pyunits.get_units(m.fs.UF.feed.properties[0.0].flow_mass_comp["tds"])}'
    )
    print(
        f'{"    TSS":<20s}{m.fs.UF.feed.properties[0.0].flow_mass_comp["tss"].value:<10.3f}{pyunits.get_units(m.fs.UF.feed.properties[0.0].flow_mass_comp["tss"])}'
    )
    print(f'{"PRODUCT":<20s}')
    print(
        f'{"    H2O":<20s}{m.fs.UF.product.properties[0.0].flow_mass_comp["H2O"].value:<10.3f}{pyunits.get_units(m.fs.UF.product.properties[0.0].flow_mass_comp["H2O"])}'
    )
    print(
        f'{"    TDS":<20s}{m.fs.UF.product.properties[0.0].flow_mass_comp["tds"].value:<10.3f}{pyunits.get_units(m.fs.UF.product.properties[0.0].flow_mass_comp["tds"])}'
    )
    print(
        f'{"    TSS":<20s}{m.fs.UF.product.properties[0.0].flow_mass_comp["tss"].value:<10.3f}{pyunits.get_units(m.fs.UF.product.properties[0.0].flow_mass_comp["tss"])}'
    )
    print(f'{"DISPOSAL":<20s}')
    print(
        f'{"    H2O":<20s}{m.fs.UF.disposal.properties[0.0].flow_mass_comp["H2O"].value:<10.3f}{pyunits.get_units(m.fs.UF.disposal.properties[0.0].flow_mass_comp["H2O"])}'
    )
    print(
        f'{"    TDS":<20s}{m.fs.UF.disposal.properties[0.0].flow_mass_comp["tds"].value:<10.3f}{pyunits.get_units(m.fs.UF.disposal.properties[0.0].flow_mass_comp["tds"])}'
    )
    print(
        f'{"    TSS":<20s}{m.fs.UF.disposal.properties[0.0].flow_mass_comp["tss"].value:<10.3f}{pyunits.get_units(m.fs.UF.disposal.properties[0.0].flow_mass_comp["tss"])}'
    )


def report_UF(m, blk, stream_table=False):
    print(f"\n\n-------------------- UF Report --------------------\n")
    print("\n")
    print(
        f'{"Inlet Flow Volume":<30s}{value(m.fs.UF.feed.properties[0.0].flow_vol):<10.3f}{pyunits.get_units(m.fs.UF.feed.properties[0.0].flow_vol)}'
    )
    print(f'{"UF Performance:":<30s}')
    print(
        f'{"    Recovery":<30s}{100*m.fs.UF.unit.recovery_frac_mass_H2O[0.0].value:<10.1f}{"%"}'
    )
    print(
        f'{"    TDS Removal":<30s}{100*m.fs.UF.unit.removal_frac_mass_comp[0.0,"tds"].value:<10.1f}{"%"}'
    )
    print(
        f'{"    TSS Removal":<30s}{100*m.fs.UF.unit.removal_frac_mass_comp[0.0,"tss"].value:<10.1f}{"%"}'
    )
    print(
        f'{"    Energy Consumption":<30s}{m.fs.UF.unit.electricity[0.0].value:<10.3f}{pyunits.get_units(m.fs.UF.unit.electricity[0.0])}'
    )
    print(
        f'{"    Specific Energy Cons.":<30s}{value(m.fs.UF.unit.energy_electric_flow_vol_inlet):<10.3f}{pyunits.get_units(m.fs.UF.unit.energy_electric_flow_vol_inlet)}'
    )

def print_UF_costing_breakdown(blk, debug=False):
    print(f'{"UF Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')

    if debug:
        print(blk.unit.costing.display())


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_UF_op_conditions(m.fs.UF)
    set_system_conditions(m.fs.UF)
    # load_parameters(m, m.fs.UF)
    add_UF_costing(m, m.fs.UF)
    m.fs.costing.cost_process()
    # m.fs.costing.initialize()
    init_UF(m, m.fs.UF)
    solve(m)

    report_UF(m, m.fs.UF)
    m.fs.UF.unit.costing.display()
    m.fs.costing.display()
    print_UF_costing_breakdown(m.fs.UF, debug=False)
    # # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
