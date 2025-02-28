import os
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
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
from watertap.core.util.initialization import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap_contrib.reflo.core import REFLODatabase
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
    "add_DWI_scaling",
    "report_DWI",
    "print_DWI_costing_breakdown",
]


def propagate_state(arc):
    _prop_state(arc)


def build_DWI(m, blk, prop_package) -> None:
    print(f'\n{"=======> BUILDING DEEP WELL INJECTION SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=prop_package)
    blk.unit = DeepWellInjection(property_package=prop_package)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )


def set_system_op_conditions(blk):
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_mgd = 2.08 * pyunits.Mgallons / pyunits.day

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


def set_DWI_op_conditions(blk):
    pass


def init_DWI(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    # assert_no_degrees_of_freedom(m)
    blk.feed.initialize(optarg=optarg, outlvl=idaeslogger.INFO)
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize(optarg=optarg, outlvl=idaeslogger.INFO)


def add_DWI_costing(m, blk, costing_blk=None):
    if costing_blk is None:
        costing_blk = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_blk,
        costing_method_arguments={
            "cost_method": "as_opex"
        },  # could be "as_capex" or "blm"
    )


def add_DWI_scaling(m, blk):
    set_scaling_factor(blk.feed.properties[0.0].mass_frac_phase_comp["Liq", "H2O"], 0.1)
    # set_scaling_factor(blk.feed.properties[0.0].mass_frac_phase_comp['Liq','TDS'], 1e-2)

    set_scaling_factor(blk.unit.properties[0.0].mass_frac_phase_comp["Liq", "H2O"], 0.1)
    set_scaling_factor(blk.unit.properties[0.0].flow_vol_phase["Liq"], 1e-2)
    # set_scaling_factor(blk.unit.properties[0.0].mass_frac_phase_comp['Liq','TDS'], 1e-2)
    # set_scaling_factor(blk.unit.properties[0.0].dens_mass_phase['Liq'], 1e3)

    # constraint_scaling_transform(blk.feed.properties[0].eq_conc_mass_phase_comp, 1e-1)
    # constraint_scaling_transform(blk.feed.properties[0.0].eq_conc_mass_phase_comp, 1)
    # set_scaling_factor(blk.feed.properties[0.0].dens_mass_phase['Liq'], 1e3)


def report_DWI(blk):
    print(f"\n\n-------------------- DWI Report --------------------\n")
    print("\n")
    print(
        f'{"Injection Well Depth":<30s}{value(blk.unit.config.injection_well_depth):<10.3f}{pyunits.get_units(blk.unit.config.injection_well_depth)}'
    )


def print_DWI_costing_breakdown(blk):
    print(f"\n\n-------------------- DWI Costing Breakdown --------------------\n")
    print(f'{"DWI Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')
    print(
        f'{"DWI Operating Cost":<35s}{f"${blk.unit.costing.variable_operating_cost():<25,.0f}"}'
    )
    print("\n")


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    m.db = REFLODatabase()
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }
    # BUG MCAS for RO flowsheets WaterParameterBlock for LTMED causing issues
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
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_system_op_conditions(m.fs.DWI)

    init_DWI(m, m.fs.DWI)
    # add_DWI_costing(m, m.fs.DWI)
    # solve(m)
    # print_DWI_costing_breakdown(m.fs.DWI)

    # print(m.fs.DWI.display())