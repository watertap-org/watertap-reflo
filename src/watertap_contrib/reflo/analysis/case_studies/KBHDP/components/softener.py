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

__all__ = [
    "build_softener",
    "set_softener_op_conditions",
    "add_softener_costing",
    "report_softener",
    "init_softener",
]


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
        solute_list=["Alkalinity_2-", "Ca_2+", "Mg_2+", "SiO2", "Na_+", "Cl_-"],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    # m.fs.product = Product(property_package=m.fs.properties)
    # m.fs.disposal = Product(property_package=m.fs.properties)
    # m.fs.softener = ChemicalSofteningZO(
    #     property_package=m.fs.properties,
    #     silica_removal=True,
    #     softening_procedure_type="excess_lime",
    # )

    m.fs.softener = FlowsheetBlock(dynamic=False)
    build_softener(m, m.fs.softener, m.fs.properties)

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.unit.inlet,
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
        silica_removal=False,
        softening_procedure_type="excess_lime",
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(blk.unit)}")


def set_system_operating_conditions(m):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )

    soft = m.fs.softener.unit
    ca_in = 0.61 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    mg_in = 0.161 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    sio2_in = 0.13 * pyunits.kg / pyunits.m**3  # g/L = kg/m3
    alk_in = 0.0821 * pyunits.kg / pyunits.m**3  # g/L = kg/m3

    q_in = 4 * pyunits.Mgal / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3

    nacl_in = 5.6
    na_in = nacl_in * pyunits.kg / pyunits.m**3
    cl_in = nacl_in * pyunits.kg / pyunits.m**3

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

    flow_mass_phase_na = pyunits.convert(q_in * na_in, to_units=pyunits.kg / pyunits.s)
    flow_mass_phase_cl = pyunits.convert(q_in * cl_in, to_units=pyunits.kg / pyunits.s)

    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)

    prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].fix(flow_mass_phase_ca)
    prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].fix(flow_mass_phase_mg)
    prop_in.flow_mass_phase_comp["Liq", "SiO2"].fix(flow_mass_phase_si)
    prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].fix(flow_mass_phase_alk)

    prop_in.flow_mass_phase_comp["Liq", "Na_+"].fix(flow_mass_phase_na)
    prop_in.flow_mass_phase_comp["Liq", "Cl_-"].fix(flow_mass_phase_cl)

    print(prop_in.flow_mass_phase_comp["Liq", "Ca_2+"].value)
    print(prop_in.flow_mass_phase_comp["Liq", "Mg_2+"].value)
    print(prop_in.flow_mass_phase_comp["Liq", "SiO2"].value)
    print(prop_in.flow_mass_phase_comp["Liq", "Alkalinity_2-"].value)

    prop_in.temperature.fix(298)
    prop_in.pressure.fix(101356)
    prop_out.temperature.fix(298)
    prop_out.pressure.fix(101356)

    # m.fs.feed.properties[0].conc_mass_phase_comp
    # soft.properties_waste[0].conc_mass_phase_comp["Liq", "Mg_2+"]
    # soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
    print(degrees_of_freedom(m))

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_softener_op_conditions(m, soft, ca_effluent=0.03, mg_effluent=0.02, non_important_removals=0.01):
    print(
        "\n\n-------------------- SETTING SOFTENER REMOVAL EFFICIENCY --------------------\n\n"
    )
    CO2_in = 0.10844915 * pyunits.kg / pyunits.m**3
    
    non_important_comps = ["Na_+", "Cl_-", "K_+", "SO2_-4+"]
    # fix removal efficiency for all comps...
    soft.removal_efficiency.fix()
    # ...then refix non important ones
    for comp in non_important_comps:
        m.fs.softener.unit.removal_efficiency[comp].fix(non_important_removals)

    prop_out = soft.properties_out[0.0]
    prop_out.temperature.fix(293)
    prop_out.pressure.fix(101325)

    soft.ca_eff_target.fix(ca_effluent)
    soft.mg_eff_target.fix(mg_effluent)

    soft.no_of_mixer.fix(2)
    soft.no_of_floc.fix(4)
    soft.retention_time_mixer.fix(0.4)
    soft.retention_time_floc.fix(25)
    soft.retention_time_sed.fix(120)
    soft.retention_time_recarb.fix(20)
    soft.frac_vol_recovery.fix(0.99)
    soft.vel_gradient_mix.fix(300)
    soft.vel_gradient_floc.fix(50)
    
    # soft.removal_efficiency["SiO2"].fix(0)
    # soft.CO2_CaCO3.fix(CO2_in)

    # soft.excess_CaO.fix(0)
    soft.CO2_second_basin.fix(0)
    soft.Na2CO3_dosing.fix(0)
    soft.MgCl2_dosing.fix(0)

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
    # m.fs.costing.add_LCOW(blk.unit.properties_in[0].flow_vol)


def init_system(blk, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")
    # assert_no_degrees_of_freedom(m)
    print("\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_softener)

    try:
        init_softener(m, m.fs.softener.unit)
    except:
        print_infeasible_bounds(m.fs.softener)
        print_infeasible_constraints(m.fs.softener)
        print_close_to_bounds(m.fs.softener)

        # assert False

    # m.fs.product.initialize()
    # m.fs.disposal.initialize()


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
    blk.initialize()
    # try:
    #     blk.initialize()
    # except:
    #     print_infeasible_bounds(m.fs.softener)
    #     print_infeasible_constraints(m.fs.softener)
    #     print_close_to_bounds(m.fs.softener)

    # propagate_state(blk.unit_to_product)
    # blk.product.initialize()
    # blk.disposal.initialize()
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


def report_softener(m, blk=None):
    if blk is None:
        blk = m.fs.softener.unit

    comps = blk.config.property_package.solute_set
    print(f"\n\n-------------------- Softener Report --------------------\n")

    print("Stream Table:\n")

    print(
        f'{"Flow Basis":<20s}{value(pyunits.convert(blk.properties_in[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day)):<10.3f}{pyunits.get_units(pyunits.convert(blk.properties_in[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day))}'
    )
    print(
        f'{"Overall TDS:":<20s}{sum(value(blk.properties_in[0].conc_mass_phase_comp["Liq", i]) for i in comps):<10.3f}{"g/L"}\n'
    )

    print(
        f'{"Component":<20s}{"Flow In (g/L)":<20s}{"Flow Product (g/L)":<20s}{"Flow Waste (g/L)":20s}{"Removal %":20s}'
    )

    for solute in comps:
        flow_in = blk.properties_in[0.0].conc_mass_phase_comp["Liq", solute].value
        flow_out = blk.properties_out[0.0].conc_mass_phase_comp["Liq", solute].value
        flow_waste = blk.properties_waste[0.0].conc_mass_phase_comp["Liq", solute].value
        removal = (flow_in - flow_out) / flow_in * 100
        print(
            f"{solute:20s}{flow_in:<20.3f}{flow_out:<20.3f}{flow_waste:<20.3f}{removal:<20.1f}"
        )

    print("\nDosing Details:")
    print(
        f'{"CaO Dose":<20s}{blk.CaO_dosing.value:<10.3f}{pyunits.get_units(blk.CaO_dosing)}'
    )
    print(
        f'{"MgCl Dose":<20s}{value(blk.MgCl2_dosing):<10.3f}{pyunits.get_units(blk.MgCl2_dosing)}'
    )
    print(
        f'{"Na2CO3 Dose":<20s}{blk.Na2CO3_dosing.value:<10.3f}{pyunits.get_units(blk.Na2CO3_dosing)}'
    )
    print(
        f'{"Excess CaO":<20s}{blk.excess_CaO.value:<10.3f}{pyunits.get_units(blk.excess_CaO)}'
    )
    print(
        f'{"CO2 CaCO3":<20s}{blk.CO2_CaCO3.value:<10.3f}{pyunits.get_units(blk.CO2_CaCO3)}'
    )
    print(
        f'{"Sludge Produced":<20s}{blk.sludge_prod.value:<10.3f}{pyunits.get_units(blk.sludge_prod)}'
    )

    print("\nCosting Report:")
    print(
        f'{"Capital Cost":<19s}{f"${blk.costing.capital_cost.value:<15,.0f}"}{pyunits.get_units(blk.costing.capital_cost)}'
    )
    print(
        f'{"Operating Cost":<19s}{f"${blk.costing.fixed_operating_cost.value:<15,.0f}"}{pyunits.get_units(blk.costing.fixed_operating_cost)}'
    )
    print(
        f'{"Electricity Flow":<20s}{blk.costing.electricity_flow.value:<15.2f}{pyunits.get_units(blk.costing.electricity_flow)}'
    )


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_system_operating_conditions(m)
    set_softener_op_conditions(m, m.fs.softener.unit)
    add_softener_costing(m, m.fs.softener)
    # set_scaling(m)
    # report_MCAS_stream_conc(m)
    init_system(m)