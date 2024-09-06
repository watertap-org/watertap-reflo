import os
import math

from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import Product, Feed, StateJunction, Separator

from watertap.core.wt_database import Database
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
import watertap.core.zero_order_properties as prop_ZO
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.pressure_changer import Pump

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.analysis.case_studies.KBHDP import *

__all__ = [
    "build_system",
    "add_connections",
    "add_constraints",
    "add_costing",
    "relax_constraints",
    "set_inlet_conditions",
    "set_operating_conditions",
    "report_MCAS_stream_conc",
    "display_system_stream_table",
    "display_system_build",
    "init_system",
]


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print("\n")


def main():
    file_dir = os.path.dirname(os.path.abspath(__file__))

    m = build_system()
    display_system_build(m)
    add_connections(m)
    add_constraints(m)
    relax_constraints(m, m.fs.RO)
    set_operating_conditions(m)
    init_system(m)
    add_costing(m)
    solve(m)
    display_system_stream_table(m)
    report_softener(m)
    report_UF(m, m.fs.UF)
    report_RO(m, m.fs.RO)
    display_costing_breakdown(m)


def build_system():
    m = ConcreteModel()
    m.db = Database()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REFLOCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

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

    m.fs.RO_properties = NaClParameterBlock()
    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.feed = Feed(property_package=m.fs.MCAS_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.disposal = Product(property_package=m.fs.RO_properties)

    # Define the Unit Models
    m.fs.softener = FlowsheetBlock(dynamic=False)
    m.fs.UF = FlowsheetBlock(dynamic=False)
    # m.fs.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.RO = FlowsheetBlock(dynamic=False)

    # Define the Translator Blocks
    m.fs.MCAS_to_NaCl_translator = Translator_MCAS_to_NACL(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.MCAS_to_TDS_translator = Translator_MCAS_to_TDS(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.TDS_to_NaCl_translator = Translator_TDS_to_NACL(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_softener(m, m.fs.softener, prop_package=m.fs.MCAS_properties)
    build_UF(m, m.fs.UF, prop_package=m.fs.UF_properties)
    build_ro(m, m.fs.RO, prop_package=m.fs.RO_properties)

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def add_connections(m):

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.unit.inlet,
    )

    m.fs.softener_to_translator = Arc(
        source=m.fs.softener.unit.outlet,
        destination=m.fs.MCAS_to_TDS_translator.inlet,
    )

    m.fs.translator_to_UF = Arc(
        source=m.fs.MCAS_to_TDS_translator.outlet,
        destination=m.fs.UF.feed.inlet,
    )

    m.fs.UF_to_translator3 = Arc(
        source=m.fs.UF.product.outlet,
        destination=m.fs.TDS_to_NaCl_translator.inlet,
    )

    # m.fs.translator_to_pump = Arc(
    #     source=m.fs.TDS_to_NaCl_translator.outlet,
    #     destination=m.fs.pump.inlet,
    # )

    # m.fs.pump_to_ro = Arc(
    #     source=m.fs.pump.outlet,
    #     destination=m.fs.RO.feed.inlet,
    # )

    m.fs.translator_to_ro = Arc(
        source=m.fs.TDS_to_NaCl_translator.outlet,
        destination=m.fs.RO.feed.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.RO.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.ro_to_disposal = Arc(
        source=m.fs.RO.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.feed_flow_vol = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.L / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.product.properties[0].conc_mass_phase_comp
    m.fs.disposal.properties[0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_in[0.0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_out[0.0].conc_mass_phase_comp


def add_costing(m):
    # m.fs.pump.costing = UnitModelCostingBlock(
    #     flowsheet_costing_block=m.fs.costing,
    # )

    add_softener_costing(m, m.fs.softener)
    add_UF_costing(m, m.fs.UF)
    add_ro_costing(m, m.fs.RO)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)


def relax_constraints(m, blk):
    # Release constraints related to low concentration
    for idx, stage in blk.stage.items():
        stage.module.width.setub(10000)
        # for item in [stage.module.permeate_side, stage.module.feed_side.properties_interface]:
        #     for idx, param in item.items():
        #         if idx[1] > 0:
        #             param.molality_phase_comp["Liq", "NaCl"].setlb(0)
        #             param.pressure_osm_phase["Liq"].setlb(0)
        #             param.conc_mass_phase_comp["Liq", "NaCl"].setlb(0)

    # for idx, param in blk.module.feed_side.friction_factor_darcy.items():
    #     # if idx[1] > 0:
    #     param.setub(100)

    # # Release constraints related to low velocity and low flux
    # for idx1, item in enumerate([blk.module.feed_side.K, blk.module.feed_side.cp_modulus]):
    #     for idx2, param in item.items():
    #         if idx1 > 0:
    #             if idx2[1] > 0:
    #                 param.setub(4)
    #         else:
    #             if idx2[1] > 0:
    #                 param.setlb(0)


def define_inlet_composition(m):

    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=[
            "cod",
            "nonbiodegradable_cod",
            "ammonium_as_nitrogen",
            "phosphates",
        ]
    )


def set_inlet_conditions(
    m,
    Qin=None,
    Cin=None,
    water_recovery=None,
    supply_pressure=1e5,
    primary_pump_pressure=15e5,
):
    """Sets operating condition for the PV-RO system

    Args:
        m (obj): Pyomo model
        flow_in (float, optional): feed volumetric flow rate [m3/s]. Defaults to 1e-2.
        conc_in (int, optional): solute concentration [g/L]. Defaults to 30.
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    q_in = Qin * pyunits.Mgal / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_mass_phase_water = pyunits.convert(q_in * rho, to_units=pyunits.kg / pyunits.s)

    if Qin is None:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            flow_mass_phase_water
        )

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

    # inlet_dict = {
    #     "Ca_2+": 0.13 * pyunits.kg / pyunits.m**3,
    #     "Mg_2+": 0.03 * pyunits.kg / pyunits.m**3,
    #     "Alkalinity_2-": 0.08 * pyunits.kg / pyunits.m**3,
    #     "SiO2": 0.031 * pyunits.kg / pyunits.m**3,
    #     "Cl_-": 1.18 * pyunits.kg / pyunits.m**3,
    #     "Na_+": 0.77 * pyunits.kg / pyunits.m**3,
    #     "K_+": 0.016 * pyunits.kg / pyunits.m**3,
    #     "SO2_-4+": 0.23 * pyunits.kg / pyunits.m**3,
    # }

    for solute, solute_conc in inlet_dict.items():
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            pyunits.convert(
                (
                    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
                    / (1000 * pyunits.kg / pyunits.m**3)
                )
                * solute_conc,
                to_units=pyunits.kg / pyunits.s,
            )
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute]),
            index=("Liq", solute),
        )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )

    # if Cin is None:
    #     m.fs.feed_salinity.fix(10)
    # else:
    #     m.fs.feed_salinity.fix(Cin)

    # if water_recovery is not None:
    #     m.fs.water_recovery.fix(water_recovery)
    #     m.fs.primary_pump.control_volume.properties_out[0].pressure.unfix()
    # else:
    #     m.fs.water_recovery.unfix()
    #     m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(primary_pump_pressure)

    # m.fs.pump.control_volume.properties_out[0].pressure.fix(primary_pump_pressure)

    # # iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)
    # iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    # iscale.set_scaling_factor(m.fs.feed_salinity, 1)

    # m.fs.feed_flow_constraint = Constraint(
    #         expr=m.fs.feed_flow_mass == m.fs.perm_flow_mass / m.fs.water_recovery
    #     )
    # iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)

    feed_temperature = 273.15 + 20
    pressure_atm = 101325

    # # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)
    # m.fs.disposal.pressure[0].fix(101356)
    # m.fs.disposal.temperature[0].fix(feed_temperature)

    # m.fs.primary_pump.efficiency_pump.fix(0.85)
    # iscale.set_scaling_factor(m.fs.primary_pump.control_volume.work, 1e-3)

    # m.fs.feed.properties[0].flow_vol_phase["Liq"]
    # m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    # m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
    #     m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    # )
    # m.fs.feed.flow_mass_phase_comp[
    #     0, "Liq", "H2O"
    # ].value = m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)

    # scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    # scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    # )
    # m.fs.properties.set_default_scaling(
    #     "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    # )

    # assert_units_consistent(m)
    # m.fs.feed.properties[0].display()
    # report_MCAS_stream_conc(m)


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def set_operating_conditions(m):
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4, supply_pressure=1e5, primary_pump_pressure=10e5)
    set_softener_op_conditions(m, m.fs.softener.unit, ca_eff=0.3, mg_eff=0.2)
    # # inlet_dict = {
    # #     "Ca_2+": 0.13 * pyunits.kg / pyunits.m**3,
    # #     "Mg_2+": 0.03 * pyunits.kg / pyunits.m**3,
    # #     "Alkalinity_2-": 0.08 * pyunits.kg / pyunits.m**3,
    # #     "SiO2": 0.031 * pyunits.kg / pyunits.m**3,
    # #     "Cl_-": 1.18 * pyunits.kg / pyunits.m**3,
    # #     "Na_+": 0.77 * pyunits.kg / pyunits.m**3,
    # #     "K_+": 0.016 * pyunits.kg / pyunits.m**3,
    # #     "SO2_-4+": 0.23 * pyunits.kg / pyunits.m**3,}
    non_important_removals = 0.01
    m.fs.softener.unit.removal_efficiency["Cl_-"].fix(non_important_removals)
    m.fs.softener.unit.removal_efficiency["Na_+"].fix(non_important_removals)
    m.fs.softener.unit.removal_efficiency["K_+"].fix(non_important_removals)
    m.fs.softener.unit.removal_efficiency["SO2_-4+"].fix(non_important_removals)

    set_UF_op_conditions(m.fs.UF)
    set_ro_system_operating_conditions(
        m, m.fs.RO, mem_area=10000, RO_pump_pressure=20e5
    )


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize(optarg=optarg)
    # propagate_state(m.fs.feed_to_primary_pump)
    propagate_state(m.fs.feed_to_softener)
    report_MCAS_stream_conc(m, m.fs.feed.properties[0.0])
    init_softener(m, m.fs.softener.unit)
    propagate_state(m.fs.softener_to_translator)
    m.fs.MCAS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(m.fs.translator_to_UF)
    init_UF(m, m.fs.UF)
    propagate_state(m.fs.UF_to_translator3)
    # BUG Need to add UF to disposal

    m.fs.TDS_to_NaCl_translator.initialize(optarg=optarg)

    # propagate_state(m.fs.translator_to_pump)
    # m.fs.pump.initialize(optarg=optarg)

    # propagate_state(m.fs.pump_to_ro)
    propagate_state(m.fs.translator_to_ro)
    init_ro_system(m, m.fs.RO)

    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_disposal)

    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)
    display_system_stream_table(m)
    m.fs.softener.unit.properties_waste[0].conc_mass_phase_comp


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
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
        print_infeasible_bounds(m)
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print(msg)
        return results


def report_MCAS_stream_conc(m, stream):
    solute_set = m.fs.MCAS_properties.solute_set
    print(f"\n\n-------------------- {stream} CONCENTRATIONS --------------------\n\n")
    print(f'{"Component":<15s}{"Conc.":<10s}{"Units":10s}')
    for i in solute_set:
        print(
            f"{i:<15s}: {stream.conc_mass_phase_comp['Liq', i].value:<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp['Liq', i])}"
        )
    print(
        f'{"Overall TDS":<15s}: {sum(value(stream.conc_mass_phase_comp["Liq", i]) for i in solute_set):<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp["Liq", "Ca_2+"])}'
    )
    print(
        f"{'Vol. Flow Rate':<15s}: {stream.flow_mass_phase_comp['Liq', 'H2O'].value:<10.3f}{pyunits.get_units(stream.flow_mass_phase_comp['Liq', 'H2O'])}"
    )


def report_stream_ion_conc(m, stream):
    print(f"\n\n-------------------- {stream} CONCENTRATIONS --------------------\n\n")
    for ion_conc in stream.conc_mass_phase_comp:
        print(
            f"{ion_conc[1]:<15s}: {stream.conc_mass_phase_comp[ion_conc].value:<10.3f}{str(pyunits.get_units(stream.conc_mass_phase_comp[ion_conc]))}"
        )


def display_system_stream_table(m):
    print("\n\n-------------------- SYSTEM STREAM TABLE --------------------\n\n")
    print(
        f'{"NODE":<20s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<20s}{m.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}'
    )
    print(
        f'{"Softener Inlet":<20s}{m.fs.softener.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.softener.unit.properties_in[0.0].pressure, to_units=pyunits.bar)():<30.1f}'
    )
    print(
        f'{"Softener Outlet":<20s}{m.fs.softener.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.softener.unit.properties_out[0.0].pressure, to_units=pyunits.bar)():<30.1f}'
    )
    display_flow_table(m.fs.RO)
    print("\n\n")


def display_system_build(m):
    blocks = []
    for v in m.fs.component_data_objects(ctype=Block, active=True, descend_into=False):
        print(v)


def display_costing_breakdown(m):
    print("\n\n-------------------- SYSTEM COSTING BREAKDOWN --------------------\n\n")
    header = f'{"PARAM":<25s}{"VALUE":<25s}{"UNITS":<25s}'
    print(header)
    print(
        f'{"Product Flow":<25s}{f"{value(pyunits.convert(m.fs.product.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
    )
    print(f'{"LCOW":<24s}{f"${m.fs.costing.LCOW():<25.3f}"}{"$/m3":<25s}')


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()

# TODO Add costing to the system
# TODO Use case study input values
