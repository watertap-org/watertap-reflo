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
from watertap_contrib.reflo.unit_models.chemical_softening import (
    ChemicalSoftening,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

__all__ = [
    "build_softener",
    "set_softener_op_conditions",
    "add_softener_costing",
    "report_softener",
    "init_softener",
    "print_softening_costing_breakdown",
]
rho = 1000 * pyunits.kg / pyunits.m**3


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
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Mg_2+",
            "SiO2",
            "Na_+",
            "Cl_-",
            "K_+",
            "SO2_-4+",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.softener = FlowsheetBlock(dynamic=False)
    build_softener(m, m.fs.softener, m.fs.properties)

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.unit.inlet,
    )

    m.fs.softener_to_product = Arc(
        source=m.fs.softener.unit.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.softener_to_disposal = Arc(
        source=m.fs.softener.unit.waste,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")

    return m


def build_softener(m, blk, prop_package=None) -> None:
    print(f'\n{"=======> BUILDING SOFTENER SYSTEM <=======":^60}\n')
    if prop_package is None:
        prop_package = m.fs.properties

    blk.unit = ChemicalSoftening(
        property_package=prop_package,
        silica_removal=False,
        softening_procedure_type="excess_lime_soda",
    )

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(blk.unit)}")


def set_system_operating_conditions(
    m,
    Qin=4,
):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )

    soft = m.fs.softener.unit
    prop_in = soft.properties_in[0]
    prop_out = soft.properties_out[0]
    Qin = Qin * pyunits.Mgal / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)

    flow_mass_phase_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)

    inlet_dict = {
        "Ca_2+": 0.61 * pyunits.kg / pyunits.m**3,
        "Mg_2+": 0.161 * pyunits.kg / pyunits.m**3,
        "Alkalinity_2-": 0.0821 * pyunits.kg / pyunits.m**3,
        "SiO2": 0.13 * pyunits.kg / pyunits.m**3,
        "Cl_-": 5.5 * pyunits.kg / pyunits.m**3,
        "Na_+": 5.5 * pyunits.kg / pyunits.m**3,
        "K_+": 0.016 * pyunits.kg / pyunits.m**3,
        "SO2_-4+": 0.23 * pyunits.kg / pyunits.m**3,
    }
    calc_state_dict = {
        ("flow_vol_phase", "Liq"): value(flow_in),
        ("pressure", None): 101325,
        ("temperature", None): 298,
    }

    for solute, solute_conc in inlet_dict.items():
        calc_state_dict[("conc_mass_phase_comp", ("Liq", solute))] = solute_conc
        flow_mass_solute = pyunits.convert(
            flow_in * solute_conc, to_units=pyunits.kg / pyunits.s
        )
        sf = 1 / value(flow_mass_solute)
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].set_value(
            flow_mass_solute
        )
        m.fs.softener.unit.properties_in[0].flow_mass_phase_comp[
            "Liq", solute
        ].set_value(flow_mass_solute)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            sf,
            index=("Liq", solute),
        )
        m.fs.properties.set_default_scaling(
            "conc_mass_phase_comp",
            1 / solute_conc(),
            index=("Liq", solute),
        )
        # m.fs.properties.set_default_scaling(
        #     "mass_frac_phase_comp",
        #     1 / value(flow_mass_solute / flow_mass_phase_water),
        #     index=("Liq", solute),
        # )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )
    calculate_scaling_factors(m)

    m.fs.feed.properties.calculate_state(var_args=calc_state_dict, hold_state=True)
    calculate_scaling_factors(m)

    print(f"DOF = {degrees_of_freedom(m)}")

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener)}")


def set_softener_op_conditions(
    m, soft, ca_effluent=0.03, mg_effluent=0.02, non_important_removals=0.01
):
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
    # prop_out.temperature.fix(293)
    # prop_out.pressure.fix(101325)

    soft.ca_eff_target.fix(ca_effluent)
    soft.mg_eff_target.fix(mg_effluent)
    # soft.ca_eff_target.set_value(ca_effluent)
    # soft.mg_eff_target.set_value(mg_effluent)

    soft.number_mixers.value = 2
    soft.number_floc.value = 4
    soft.retention_time_mixer.fix(0.4)
    soft.retention_time_floc.fix(25)
    soft.retention_time_sed.fix(120)
    soft.retention_time_recarb.fix(20)
    soft.frac_mass_water_recovery.fix(0.99)
    soft.vel_gradient_mix.fix(300)
    soft.vel_gradient_floc.fix(50)

    # soft.removal_efficiency["SiO2"].fix(0)
    # soft.CO2_CaCO3.fix(0.10)

    # # soft.excess_CaO.fix(0)
    # soft.CO2_second_basin.fix(0)
    # soft.Na2CO3_dosing.fix(0)
    soft.MgCl2_dosing.fix(0)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Softener Degrees of Freedom: {degrees_of_freedom(m.fs.softener.unit)}")
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


def add_softener_costing(m, blk, costing_blk=None):
    if costing_blk is None:
        costing_blk = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_blk,
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

    # try:
    init_softener(m, m.fs.softener.unit)
    # except:
    #     print_infeasible_bounds(m.fs.softener)
    #     print_infeasible_constraints(m.fs.softener)
    #     print_close_to_bounds(m.fs.softener)

    # assert False
    propagate_state(m.fs.softener_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.softener_to_disposal)
    m.fs.disposal.initialize()


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

    # blk.report()


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


def print_softening_costing_breakdown(blk):
    print(
        f'{"Softening Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}'
    )
    print(
        f'{"Softening Operating Cost":<35s}{f"${blk.unit.costing.fixed_operating_cost():<25,.0f}"}'
    )


# calculate_variable_from_constraint
if __name__ == "__main__":
    from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.utils import check_jac

    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()

    soft = m.fs.softener.unit
    soft.properties_out[0].conc_mass_phase_comp[...]
    soft.sedimentation_overflow.fix()
    soft.CO2_CaCO3.fix(0.063)
    soft.ca_eff_target.fix(0.01)
    soft.mg_eff_target.fix(0.01)

    set_system_operating_conditions(m)
    set_softener_op_conditions(m, m.fs.softener.unit)
    print(f"\n\n\n\n\nDOF = {degrees_of_freedom(m)}\n\n\n\n")

    add_softener_costing(m, m.fs.softener)
    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    m.fs.costing.add_LCOW(0.2 * m.fs.softener.unit.properties_in[0].flow_vol)
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    init_system(m)
    solve(m)

    print_softening_costing_breakdown(m.fs.softener)
    print(f'{"Softening LCOW":<35s}{f"${m.fs.costing.LCOW():<25,.2f}"}')
    report_softener(m, m.fs.softener.unit)
    # print(m.fs.costing.display())
    print(
        f'{"Product Flow":<35s}{f"{value(pyunits.convert(m.fs.feed.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
    )
    print(
        f'{"Total Capital Cost":<35s}{f"${m.fs.costing.total_capital_cost():<25,.0f}"}'
    )
    print(
        f'{"Total Operating Cost":<35s}{f"${m.fs.costing.total_operating_cost():<25,.0f}"}'
    )
    for flow in m.fs.costing.aggregate_flow_costs:
        print(
            f'{f"    Flow Cost [{flow}]":<35s}{f"${m.fs.costing.aggregate_flow_costs[flow]():<25,.3f}"}'
        )
    # assert_units_consistent(soft)
    # print(pyunits.get_units(soft.Mg_CaCO3))
    # assert False
    # soft = m.fs.softener.unit
    # # soft.display()
    # # soft.excess_CaO_coeff.display()
    # soft.CO2_CaCO3.display()
    # soft.CaO_dosing.display()
    # soft.Ca_CaCO3.display()
    # soft.Mg_CaCO3.display()
    # soft.excess_CaO.display()
    # soft.Mg_hardness_CaCO3.display()
    # soft.Ca_hardness_CaCO3.display()
    # # soft.ca_eff_target.display()
    # # soft.mg_eff_target.display()
    # # soft.CO2_first_basin.display()
    # soft.MgCl2_dosing.display()
    # soft.MgCl2_SiO2_ratio.display()
    # soft.CO2_second_basin.display()
    # soft.properties_in[0].conc_mass_phase_comp.display()
    # soft.properties_out[0].conc_mass_phase_comp.display()
    # soft.properties_out[0].flow_mass_phase_comp.display()
    # soft.excess_CaO.fix()
    # soft.CaO_dosing.fix()
    # print(f"DOF = {degrees_of_freedom(m)}")
    # print(f"LCOW = {m.fs.costing.LCOW()}")
    # print(soft.config.softening_procedure_type)
    # e = (
    #     soft.CO2_CaCO3
    #     + soft.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
    #     + soft.Mg_CaCO3
    # ) * soft.excess_CaO_coeff
    # e = pyunits.convert(
    #     (
    #         soft.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
    #         - (soft.Ca_hardness_CaCO3 + soft.Mg_hardness_CaCO3)
    #         + soft.excess_CaO
    #         # + soft.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
    #         + soft.ca_eff_target * soft.Ca_CaCO3_conv
    #         # + soft.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
    #         + soft.mg_eff_target * soft.Mg_CaCO3_conv
    #     )
    #     * soft.properties_in[0].flow_vol_phase["Liq"]
    #     * soft.CO2_mw
    #     / soft.CaCO3_mw,
    #     to_units=pyunits.kg / pyunits.d,
    # )
    # e = soft.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"] - (
    #     soft.Ca_hardness_CaCO3 + soft.Mg_hardness_CaCO3
    # )
    # e = soft.properties_in[0].conc_mass_phase_comp["Liq", "SiO2"] * 2.35
    # print(value(e))
