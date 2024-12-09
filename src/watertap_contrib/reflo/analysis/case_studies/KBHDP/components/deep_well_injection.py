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

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.costing.units.deep_well_injection import (
    blm_costing_params_dict,
)

__all__ = [
    "build_DWI",
    "init_DWI",
    "set_DWI_op_conditions",
    "add_DWI_costing",
    "report_DWI",
    "print_DWI_costing_breakdown",
]


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_DWI(m, blk, prop_package) -> None:
    print(f'\n{"=======> BUILDING DEEP WELL INJECTION SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=prop_package)
    blk.unit = DeepWellInjection(property_package=prop_package)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )


def set_DWI_op_conditions(blk):
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_mgd = 5.08 * pyunits.Mgallons / pyunits.day

    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    prop = blk.unit.properties[0]
    prop.temperature.fix()
    prop.pressure.fix()

    prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
        prop.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / mass_flow_solute),
            index=("Liq", solute),
        )
    prop.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_phase_water),
        index=("Liq", "H2O"),
    )


def init_DWI(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    blk.unit.initialize(optarg=optarg, outlvl=idaeslogger.INFO)


def add_DWI_costing(m, blk):
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={
            "cost_method": "as_opex"
        },  # could be "as_capex" or "blm"
    )


def report_DWI(blk):
    print(f"\n\n-------------------- DWI Report --------------------\n")
    print("\n")
    print(
        f'{"Injection Well Depth":<30s}{value(blk.unit.config.injection_well_depth):<10.3f}{pyunits.get_units(blk.unit.config.injection_well_depth)}'
    )


def print_DWI_costing_breakdown(blk):
    print(f'{"DWI Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')
    print(
        f'{"DWI Operating Cost":<35s}{f"${blk.unit.costing.variable_operating_cost():<25,.0f}"}'
    )
    print("\n")


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )

    m.fs.DWI = FlowsheetBlock(dynamic=False)
    build_DWI(m, m.fs.DWI, m.fs.properties)

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


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_DWI_op_conditions(m.fs.DWI)

    init_DWI(m, m.fs.DWI)
    add_DWI_costing(m, m.fs.DWI)
    m.fs.costing.cost_process()
    solve(m)

    report_DWI(m.fs.DWI)
    print_DWI_costing_breakdown(m.fs.DWI)
