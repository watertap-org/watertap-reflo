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
from idaes.core.solvers import get_solver
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
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap_contrib.reflo.unit_models.zero_order.chemical_softening_zo import (
    ChemicalSofteningZO,
)


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_softener(m, blk, prop_package=None) -> None:
    print(f'\n{"=======> BUILDING SOFTENER SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = ChemicalSofteningZO(
        property_package=prop_package,
        silica_removal=True,
        softening_procedure_type="excess_lime_soda",
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.waste,
        destination=blk.disposal.inlet,
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_system_operating_conditions(m):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )

    Q_basis = 1 * pyunits.L / pyunits.s  # m3/d
    ca_in = 0.143 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    mg_in = 0.0314 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    sio2_in = 0.031 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    alk_in = 0.0821 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
    rho = 1000 * pyunits.kg / pyunits.m**3

    prop_in = m.fs.feed.properties[0.0]
    prop_out = m.fs.disposal.properties[0.0]

    flow_mass_phase_water = pyunits.convert(
        Q_basis * rho, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_ca = pyunits.convert(
        Q_basis * ca_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_mg = pyunits.convert(
        Q_basis * mg_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_si = pyunits.convert(
        Q_basis * sio2_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_alk = pyunits.convert(
        Q_basis * alk_in, to_units=pyunits.kg / pyunits.s
    )

    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(flow_mass_phase_ca)
    prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(flow_mass_phase_mg)
    prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(flow_mass_phase_si)
    prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(flow_mass_phase_alk)

    prop_in.temperature.fix(298)
    prop_in.pressure.fix(101356)

    m.fs.feed.properties[0].conc_mass_phase_comp

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_softener_operating_conditions(blk):
    print(
        "\n\n-------------------- SETTING SOFTENER OPERATING CONDITIONS --------------------\n\n"
    )
    Q_basis = 1 * pyunits.m**3 / pyunits.day
    ca_in = 0.13 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    mg_in = 0.03 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    sio2_in = 0.031 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    alk_in = 0.08 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    rho = 1000 * pyunits.kg / pyunits.m**3

    prop_in = blk.feed.properties[0.0]

    flow_mass_phase_water = pyunits.convert(
        Q_basis * rho, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_ca = pyunits.convert(
        Q_basis * ca_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_mg = pyunits.convert(
        Q_basis * mg_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_si = pyunits.convert(
        Q_basis * sio2_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_alk = pyunits.convert(
        Q_basis * alk_in, to_units=pyunits.kg / pyunits.s
    )

    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    # prop_in.conc_mass_phase_comp["Liq", "Ca_2+"].fix(ca_in)
    # prop_in.conc_mass_phase_comp["Liq", "Mg_2+"].fix(mg_in)
    # prop_in.conc_mass_phase_comp["Liq", "SiO2"].fix(sio2_in)
    # prop_in.conc_mass_phase_comp["Liq", "Alkalinity_2-"].fix(alk_in)
    prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(flow_mass_phase_ca)
    prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(flow_mass_phase_mg)
    prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(flow_mass_phase_si)
    prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(flow_mass_phase_alk)

    prop_in.temperature.fix(298)
    prop_in.pressure.fix(101325)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_removal_eff(blk):
    print(
        "\n\n-------------------- SETTING SOFTENER REMOVAL EFFICIENCY --------------------\n\n"
    )
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3

    blk.unit.ca_eff_target.fix(0.3)
    blk.unit.mg_eff_target.fix(0.2)
    blk.unit.no_of_mixer.fix(1)
    blk.unit.no_of_floc.fix(2)
    blk.unit.retention_time_mixer.fix(0.4)
    blk.unit.retention_time_floc.fix(25)
    blk.unit.retention_time_sed.fix(130)
    blk.unit.retention_time_recarb.fix(20)
    blk.unit.frac_vol_recovery.fix()
    blk.unit.removal_efficiency.fix({
        "Cl_-": 0.5,
    })  # Dict for components (non-hardness)
    blk.unit.CO2_CaCO3.fix(CO2_in)
    # blk.unit.sedimentation_overflow.fix(90)
    # blk.unit.vel_gradient_mix.fix(500)
    # blk.unit.vel_gradient_floc.fix(30)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_scaling(m):
    print("\n\n-------------------- SCALING WATER SOFTENER --------------------\n\n")
    m.fs.properties.set_default_scaling(
        "conc_mass_phase_comp",
        1,
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "conc_mass_phase_comp",
        1,
        index=("Liq", "Ca_2+"),
    )
    m.fs.properties.set_default_scaling(
        "conc_mass_phase_comp",
        1,
        index=("Liq", "Mg_2+"),
    )
    m.fs.properties.set_default_scaling(
        "conc_mass_phase_comp", 1, index=("Liq", "SiO2")
    )
    m.fs.properties.set_default_scaling(
        "conc_mass_phase_comp",
        1,
        index=("Liq", "Alkalinity_2-"),
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def init_system(blk, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_softener)

    try:
        init_softener(m, m.fs.softener)
    except:
        print("========== HERE ==========")
        print_infeasible_bounds(m.fs.softener)
        print_infeasible_constraints(m.fs.softener)
        print_close_to_bounds(m.fs.softener)
        print("========== HERE2 ==========")

    print("========== HERE3 ==========")
    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)


def init_softener(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING WATER SOFTENER --------------------\n\n"
    )

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_unit)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    
    try:
        blk.unit.initialize(optarg=optarg)
    except:
        print("========== HERE ==========")
        print_infeasible_bounds(m.fs.softener.unit)
        print_infeasible_constraints(m.fs.softener.unit)
        print_close_to_bounds(m.fs.softener.unit)
        print("========== HERE2 ==========")
        assert False

    print("========== HERE3 ==========")
    propagate_state(blk.unit_to_product)
    blk.product.initialize(optarg=optarg)
    blk.disposal.initialize(optarg=optarg)
    propagate_state(blk.unit_to_disposal)
    propagate_state(blk.unit_to_product)
    print(blk.unit.report())



def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["Alkalinity_2-", "Ca_2+", "Mg_2+", "SiO2"],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.softener = FlowsheetBlock(dynamic=False)
    build_softener(m, m.fs.softener)

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.feed.inlet,
    )

    m.fs.softener_to_product = Arc(
        source=m.fs.softener.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.unit_to_disposal = Arc(
        source=m.fs.softener.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")

    return m


def set_softener_op_conditions(m, blk, ca_eff=3, mg_eff=0.2):
    print(
        "\n\n-------------------- SETTING SOFTENER REMOVAL EFFICIENCY --------------------\n\n"
    )
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
    blk.unit.ca_eff_target.fix(ca_eff)
    blk.unit.mg_eff_target.fix(mg_eff)

    blk.unit.no_of_mixer.fix(1)
    blk.unit.no_of_floc.fix(2)
    blk.unit.retention_time_mixer.fix(0.4)
    blk.unit.retention_time_floc.fix(25)
    blk.unit.retention_time_sed.fix(130)
    blk.unit.retention_time_recarb.fix(20)
    blk.unit.frac_vol_recovery.fix()
    blk.unit.removal_efficiency.fix()
    # blk.unit.removal_efficiency = ({
    #     "Cl_-": 0.5,
    # })  # Dict for components (non-hardness)
    blk.unit.CO2_CaCO3.fix(CO2_in)
    blk.unit.vel_gradient_mix.fix(300)
    blk.unit.vel_gradient_floc.fix(50)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


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


def display_dof_breakdown(blk, decend=False, report=False):
    print(
        "\n\n-------------------- DEGREE OF FREEDOM BREAKDOWN --------------------\n\n"
    )
    print(f'{"BLOCK":<40s}{"DEGREES OF FREEDOM":<30s}')
    print(f'{"m.fs":<40s}{degrees_of_freedom(m)}')
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(f"{v.name:<40s}{degrees_of_freedom(v)}")


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f'{"m.fs":<40s}{number_unused_variables(m)}')
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")

def report_MCAS_stream_conc(m):
    solute_set = m.fs.properties.solute_set
    print("\n\n-------------------- FEED CONCENTRATIONS --------------------\n\n")
    print(f'{"Component":<15s}{"Conc.":<10s}{"Units":10s}')
    for i in solute_set:
        print(f"{i:<15s}: {m.fs.feed.properties[0].conc_mass_phase_comp['Liq', i].value:<10.3f}{pyunits.get_units(m.fs.feed.properties[0].conc_mass_phase_comp['Liq', i])}")
    print(f'{"Overall TDS":<15s}: {sum(value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq", i]) for i in solute_set):<10.3f}')
    print(f"{'Vol. Flow Rate':<15s}: {m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value:<10.3f}{pyunits.get_units(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'])}")


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_system_operating_conditions(m)
    set_removal_eff(m.fs.softener)
    set_scaling(m)
    init_system(m)
    results = solve(m)
    assert_optimal_termination(results)
    print(m.fs.softener.report())
    report_MCAS_stream_conc(m)