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
    assert_optimal_termination,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
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
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.costing.units.deep_well_injection import (
    blm_costing_params_dict,
)


# reflo_dir = pathlib.Path(__file__).resolve().parents[4]
# case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3


__all__ = [
    "build_dwi",
    "build_and_run_dwi",
    "init_dwi",
    "add_dwi_costing",
    "report_DWI",
    "print_DWI_costing_breakdown",
]


def build_and_run_dwi(Qin=5, tds=130, **kwargs):

    m = build_system()
    m.fs.optimal_solve = Var(initialize=1)
    m.fs.optimal_solve.fix()
    add_dwi_costing(m, m.fs.DWI)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol_phase["Liq"])
    set_system_operating_conditions(m, Qin=Qin, tds=tds, **kwargs)
    print(f"dof = {degrees_of_freedom(m)}")
    init_system(m, m.fs.DWI)

    solver = get_solver()
    results = solver.solve(m)
    # assert_optimal_termination(results)
    if not check_optimal_termination(results):
        m.fs.optimal_solve.fix(0)

    print(f"LCOW = {m.fs.costing.LCOW()}")

    return m


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)

    m.fs.DWI = FlowsheetBlock(dynamic=False)
    build_dwi(m, m.fs.DWI, m.fs.properties)

    m.fs.feed_to_dwi = Arc(source=m.fs.feed.outlet, destination=m.fs.DWI.feed.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_dwi(m, blk, prop_package, injection_well_depth=5000):

    print(f'\n{"=======> BUILDING DEEP WELL INJECTION SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=prop_package)
    blk.unit = DeepWellInjection(
        property_package=prop_package, injection_well_depth=injection_well_depth
    )

    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)


def set_system_operating_conditions(m, Qin=5, tds=130):

    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_water)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(flow_mass_tds)
    m.fs.feed.properties[0].temperature.fix(293)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_water),
        index=("Liq", "H2O"),
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_tds),
        index=("Liq", "TDS"),
    )


def init_system(m, blk):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_dwi)
    init_dwi(m, blk)
    print("DWI")
    print(f"\tblock DOF after init = {degrees_of_freedom(blk)}\n")
    print(f"\tunit DOF after init = {degrees_of_freedom(blk.unit)}\n")


def init_dwi(m, blk, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize(optarg=optarg, outlvl=idaeslogger.INFO)


def add_dwi_costing(m, blk, flowsheet_costing_block=None):
    """
    Add DWI costing using BLM approach
    """
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block,
        costing_method_arguments={
            "cost_method": "as_opex"
        },  # could be "as_capex" or "blm"
    )
    flowsheet_costing_block.deep_well_injection.dwi_lcow.set_value(0.1)


def report_DWI(m, blk):
    print(f"\n\n-------------------- DWI Report --------------------\n")
    print("\n")
    print(
        f'{"Injection Well Depth":<30s}{value(blk.unit.config.injection_well_depth):<10.3f}{pyunits.get_units(blk.unit.config.injection_well_depth)}'
    )


def print_DWI_costing_breakdown(m, blk):
    print(f"\n\n-------------------- DWI Costing Breakdown --------------------\n")
    print("\n")
    print(f'{"Capital Cost":<30s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')
    print(
        f'{"Capital Cost":<30s}{f"${blk.unit.costing.variable_operating_cost():<25,.0f}"}'
    )


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_and_run_dwi(Qin=5, tds=130)

    # m = build_system()
    # add_dwi_costing(m, m.fs.DWI)
    # m.fs.costing.cost_process()
    # m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol_phase["Liq"])
    # set_system_operating_conditions(m)
    print(f"dof = {degrees_of_freedom(m)}")
    # init_system(m, m.fs.DWI)

    # solver = get_solver()
    # results = solver.solve(m)
    # assert_optimal_termination(results)

    # print(f"LCOW = {m.fs.costing.LCOW()}")

    # # init_DWI(m, m.fs.DWI)
    # # add_DWI_costing(m, m.fs.DWI)
    # # m.fs.costing.cost_process()
    # # solve(m)

    # # report_DWI(m, m.fs.DWI)
    # # print_DWI_costing_breakdown(m, m.fs.DWI)
