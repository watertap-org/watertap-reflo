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
from watertap.core.util.initialization import *
from watertap_contrib.reflo.unit_models.zero_order.chemical_softening_zo import (
    ChemicalSofteningZO,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    m.fs.properties = MCASParameterBlock(
        solute_list=["Alkalinity_2-", "Ca_2+", "Mg_2+", "SiO2"],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    # m.fs.product = Product(property_package=m.fs.properties)
    # m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.softener = ChemicalSofteningZO(
        property_package=m.fs.properties,
        silica_removal=True,
        softening_procedure_type="excess_lime",
    )

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.inlet,
    )

    # m.fs.softener_to_product = Arc(
    #     source=m.fs.softener.outlet,
    #     destination=m.fs.product.inlet,
    # )

    # m.fs.unit_to_disposal = Arc(
    #     source=m.fs.softener.waste,
    #     destination=m.fs.disposal.inlet,
    # )

    TransformationFactory("network.expand_arcs").apply_to(m)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")

    return m


def build_softener(m, blk, prop_package=None) -> None:
    print(f'\n{"=======> BUILDING SOFTENER SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    blk.unit = ChemicalSofteningZO(
        property_package=prop_package,
        silica_removal=True,
        softening_procedure_type="excess_lime",
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(blk.unit)}")


def set_system_operating_conditions(m):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    soft = m.fs.softener
    ca_in = 0.61 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    mg_in = 0.161 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    sio2_in = 0.13 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    alk_in = 0.0821 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
    # q_in = 3785 * pyunits.m**3 / pyunits.day  # m3/d
    q_in = 10 * pyunits.L / pyunits.s
    rho = 1000 * pyunits.kg / pyunits.m**3

    prop_in = soft.properties_in[0]
    prop_out = soft.properties_out[0.0]

    flow_mass_phase_water = pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_phase_ca = pyunits.convert(q_in * ca_in, to_units=pyunits.kg / pyunits.s)
    flow_mass_phase_mg = pyunits.convert(q_in * mg_in, to_units=pyunits.kg / pyunits.s)
    flow_mass_phase_si = pyunits.convert(
        q_in * sio2_in, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_phase_alk = pyunits.convert(
        q_in * alk_in, to_units=pyunits.kg / pyunits.s
    )
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)

    prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(flow_mass_phase_ca)
    prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(flow_mass_phase_mg)
    prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(flow_mass_phase_si)
    prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(flow_mass_phase_alk)

    prop_in.temperature.fix(298)
    prop_in.pressure.fix(101356)
    prop_out.temperature.fix(298)
    prop_out.pressure.fix(101356)

    # m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.softener.properties_waste[0].conc_mass_phase_comp["Liq", "Mg_2+"]
    m.fs.softener.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_softener_op_conditions(m, blk, ca_eff=30, mg_eff=2):
    print(
        "\n\n-------------------- SETTING SOFTENER REMOVAL EFFICIENCY --------------------\n\n"
    )
    soft = blk
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
    # blk.ca_eff_target.fix(ca_eff)
    # blk.mg_eff_target.fix(mg_eff)

    # #Ca to CaCO3 conv = 2.5
    # ca_eff_target = ca_in()*10
    # #Mg to CaCO3 conv = 4.12
    # mg_eff_target = mg_in()*0.1

    blk.ca_eff_target.fix(1.312938301706484)
    blk.mg_eff_target.fix(4.2551166976399655)

    soft.no_of_mixer.fix(1)
    soft.no_of_floc.fix(2)
    soft.retention_time_mixer.fix(0.4)
    soft.retention_time_floc.fix(25)
    soft.retention_time_sed.fix(130)
    soft.retention_time_recarb.fix(20)
    soft.frac_vol_recovery.fix()
    soft.removal_efficiency.fix()
    soft.CO2_CaCO3.fix(CO2_in)
    soft.vel_gradient_mix.fix(300)
    soft.vel_gradient_floc.fix(50)
    # soft.excess_CaO.fix(0)
    soft.CO2_second_basin.fix(0)
    soft.Na2CO3_dosing.fix(0)

    prop_in = soft.properties_in[0.0]
    prop_out = soft.properties_out[0.0]
    prop_out.temperature.fix(298)
    prop_out.pressure.fix(value(prop_in.pressure))

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    # assert False


def set_scaling(m):
    print("\n\n-------------------- SCALING WATER SOFTENER --------------------\n\n")
    prop_in = m.fs.feed.properties[0.0]
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / value(prop_in.flow_mass_phase_comp["Liq", "H2O"])),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / value(prop_in.flow_mass_phase_comp["Liq", "Ca_2+"])),
        index=("Liq", "Ca_2+"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / value(prop_in.flow_mass_phase_comp["Liq", "Mg_2+"])),
        index=("Liq", "Mg_2+"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / value(prop_in.flow_mass_phase_comp["Liq", "SiO2"])),
        index=("Liq", "SiO2"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / value(prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"])),
        index=("Liq", "Alkalinity_2-"),
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def add_softener_costing(m, blk):
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    # m.fs.costing.cost_process()
    # m.fs.costing.add_LCOW(blk.properties_in[0].flow_vol)


def init_system(blk, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    # assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_softener)

    try:
        init_softener(m, m.fs.softener)
    except:
        print_infeasible_bounds(m.fs.softener)
        print_infeasible_constraints(m.fs.softener)
        print_close_to_bounds(m.fs.softener)

        assert False

    # m.fs.product.initialize(optarg=optarg)
    # m.fs.disposal.initialize(optarg=optarg)


def init_softener(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING WATER SOFTENER --------------------\n\n"
    )
    # assert_units_consistent(m)
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    blk.initialize(optarg=optarg)
    # try:
    #     blk.initialize(optarg=optarg)
    # except:
    #     print_infeasible_bounds(m.fs.softener)
    #     print_infeasible_constraints(m.fs.softener)
    #     print_close_to_bounds(m.fs.softener)

    # propagate_state(blk.unit_to_product)
    # blk.product.initialize(optarg=optarg)
    # blk.disposal.initialize(optarg=optarg)
    # propagate_state(blk.unit_to_disposal)
    # propagate_state(blk.unit_to_product)

    blk.report()


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
        print(
            f"{i:<15s}: {m.fs.feed.properties[0].conc_mass_phase_comp['Liq', i].value:<10.3f}{pyunits.get_units(m.fs.feed.properties[0].conc_mass_phase_comp['Liq', i])}"
        )
    print(
        f'{"Overall TDS":<15s}: {sum(value(m.fs.feed.properties[0].conc_mass_phase_comp["Liq", i]) for i in solute_set):<10.3f}'
    )
    print(
        f"{'Vol. Flow Rate':<15s}: {m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value:<10.3f}{pyunits.get_units(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'])}"
    )


def report_softener(m):
    print(f"\n\n-------------------- Softener Report --------------------\n")
    print("Stream Table:")
    print(
        f'{"Component":<20s}{"Flow In":<20s}{"Flow Product":<20s}{"Flow Waste":20s}{"Removal %":20s}'
    )
    for idx, val in m.fs.softener.unit.properties_in[0.0].flow_mass_phase_comp.items():
        flow_in = val.value
        flow_out = (
            m.fs.softener.unit.properties_out[0.0].flow_mass_phase_comp[idx].value
        )
        flow_waste = (
            m.fs.softener.unit.properties_waste[0.0].flow_mass_phase_comp[idx].value
        )
        removal = (flow_in - flow_out) / flow_in * 100
        print(
            f"{idx[1]:20s}{flow_in:<20.3f}{flow_out:<20.3f}{flow_waste:<20.3f}{removal:<20.1f}"
        )

    print("\nDosing Details:")
    print(
        f'{"CaO Dose":<20s}{m.fs.softener.unit.CaO_dosing.value:<20.3f}{pyunits.get_units(m.fs.softener.unit.CaO_dosing)}'
    )
    print(
        f'{"MgCl Dose":<20s}{value(m.fs.softener.unit.MgCl2_dosing):<20.3f}{pyunits.get_units(m.fs.softener.unit.MgCl2_dosing)}'
    )
    print(
        f'{"Na2CO3 Dose":<20s}{m.fs.softener.unit.Na2CO3_dosing.value:<20.3f}{pyunits.get_units(m.fs.softener.unit.Na2CO3_dosing)}'
    )
    print(
        f'{"Excess CaO":<20s}{m.fs.softener.unit.excess_CaO.value:<20.3f}{pyunits.get_units(m.fs.softener.unit.excess_CaO)}'
    )
    print(
        f'{"CO2 CaCO3":<20s}{m.fs.softener.unit.CO2_CaCO3.value:<20.3f}{pyunits.get_units(m.fs.softener.unit.CO2_CaCO3)}'
    )
    print(
        f'{"Sludge Produced":<20s}{m.fs.softener.unit.sludge_prod.value:<20.3f}{pyunits.get_units(m.fs.softener.unit.sludge_prod)}'
    )

    print("\nCosting Report:")
    print(
        f'{"Capital Cost":<19s}{f"${m.fs.softener.unit.costing.capital_cost.value:<15,.0f}"}{pyunits.get_units(m.fs.softener.unit.costing.capital_cost)}'
    )
    print(
        f'{"Operating Cost":<19s}{f"${m.fs.softener.unit.costing.fixed_operating_cost.value:<15,.0f}"}{pyunits.get_units(m.fs.softener.unit.costing.fixed_operating_cost)}'
    )
    print(
        f'{"Electricity Flow":<20s}{m.fs.softener.unit.costing.electricity_flow.value:<15.2f}{pyunits.get_units(m.fs.softener.unit.costing.electricity_flow)}'
    )


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_system_operating_conditions(m)
    set_softener_op_conditions(m, m.fs.softener)
    # add_softener_costing(m, m.fs.softener)
    # set_scaling(m)
    # m.fs.softener.unit.display()
    init_system(m)
    results = solve(m)
    # assert_optimal_termination(results)
