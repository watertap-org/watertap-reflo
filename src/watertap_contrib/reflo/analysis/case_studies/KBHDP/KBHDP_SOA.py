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
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc

from idaes.core.util.scaling import *
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
    "set_inlet_conditions",
    "set_operating_conditions",
    "report_MCAS_stream_conc",
    "display_system_stream_table",
    "display_system_build",
    "init_system",
]


def propagate_state(arc, detailed=True):
    _prop_state(arc)
    if detailed:
        print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print("\n")


def main():
    file_dir = os.path.dirname(os.path.abspath(__file__))

    m = build_system()
    display_system_build(m)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    init_system(m)
    add_costing(m)
    display_system_build(m)
    optimize(m, ro_mem_area=None, water_recovery=0.5)
    solve(m, debug=True)
    display_system_stream_table(m)
    report_softener(m)
    report_UF(m, m.fs.UF)
    report_RO(m, m.fs.RO)
    display_costing_breakdown(m)
    print_all_results(m)


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
    m.fs.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.RO = FlowsheetBlock(dynamic=False)
    m.fs.DWI = FlowsheetBlock(dynamic=False)

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
    build_ro(m, m.fs.RO, prop_package=m.fs.RO_properties, number_of_stages=2)
    build_DWI(m, m.fs.DWI, prop_package=m.fs.RO_properties)

    m.fs.units = [m.fs.softener, m.fs.UF, m.fs.pump, m.fs.RO, m.fs.DWI]

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

    m.fs.translator_to_pump = Arc(
        source=m.fs.TDS_to_NaCl_translator.outlet,
        destination=m.fs.pump.inlet,
    )

    m.fs.pump_to_ro = Arc(
        source=m.fs.pump.outlet,
        destination=m.fs.RO.feed.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.RO.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.ro_to_disposal = Arc(
        source=m.fs.RO.disposal.outlet,
        destination=m.fs.DWI.feed.inlet,
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
        expr=m.fs.feed.properties[0].flow_vol_phase["Liq"] * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.product.properties[0].conc_mass_phase_comp
    m.fs.disposal.properties[0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_in[0.0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_out[0.0].conc_mass_phase_comp


def add_costing(m):
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    add_softener_costing(m, m.fs.softener)
    add_UF_costing(m, m.fs.UF)
    add_ro_costing(m, m.fs.RO)
    add_DWI_costing(m, m.fs.DWI)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    m.fs.costing.initialize()

    m.fs.costing.total_annualized_cost = pyo.Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        ),
        doc="Total annualized cost of operation",
    )


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
    Qin=4,
    supply_pressure=1.1e5,
    primary_pump_pressure=15e5,
):

    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')
    feed_temperature = 293

    Qin = Qin * pyunits.Mgal / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    rho = 1000 * pyunits.kg / pyunits.m**3
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
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            sf,
            index=("Liq", solute),
        )
        m.fs.MCAS_properties.set_default_scaling(
            "conc_mass_phase_comp",
            1 / solute_conc(),
            index=("Liq", solute),
        )
        m.fs.MCAS_properties.set_default_scaling(
            "mass_frac_phase_comp",
            1 / value(flow_mass_solute / flow_mass_phase_water),
            index=("Liq", solute),
        )

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )
    calculate_scaling_factors(m)

    m.fs.feed.properties.calculate_state(var_args=calc_state_dict, hold_state=True)
    calculate_scaling_factors(m)

    # # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def set_operating_conditions(m, RO_pressure=20e5, supply_pressure=1.1e5):
    pump_efi = 0.8  # pump efficiency [-]

    set_inlet_conditions(m, Qin=4, supply_pressure=1.1e5)
    set_softener_op_conditions(m, m.fs.softener.unit)
    set_UF_op_conditions(m.fs.UF)
    m.fs.pump.efficiency_pump.fix(pump_efi)
    m.fs.pump.control_volume.properties_in[0].pressure.fix(supply_pressure)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(
        m,
        m.fs.RO,
        mem_area=10000,
    )


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_softener)
    report_MCAS_stream_conc(m, m.fs.feed.properties[0.0])
    init_softener(m, m.fs.softener.unit)
    propagate_state(m.fs.softener_to_translator)
    m.fs.MCAS_to_TDS_translator.initialize()
    propagate_state(m.fs.translator_to_UF)
    init_UF(m, m.fs.UF)
    propagate_state(m.fs.UF_to_translator3)

    m.fs.TDS_to_NaCl_translator.initialize()

    propagate_state(m.fs.translator_to_pump)
    m.fs.pump.initialize()

    propagate_state(m.fs.pump_to_ro)

    init_ro_system(m, m.fs.RO)
    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_disposal)

    m.fs.product.initialize()
    m.fs.disposal.initialize()
    display_system_stream_table(m)
    m.fs.softener.unit.properties_waste[0].conc_mass_phase_comp


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        m.fs.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        print(f"------- Unfixed RO Pump Pressure -------")
        m.fs.pump.control_volume.properties_out[0].pressure.unfix()

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for idx, stage in m.fs.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        print(f"\n------- Unfixed RO Membrane Area -------\n")
        for idx, stage in m.fs.RO.stage.items():
            stage.module.area.unfix()


def solve(model, solver=None, tee=False, raise_on_failure=False, debug=False):
    print(f"DEGREES OF FREEDOM BEFORE SOLVING: {degrees_of_freedom(model)}")

    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    solver.options["max_iter"] = 3000
    results = solver.solve(
        model,
        tee=tee,
    )

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        if debug:
            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(model)

            print("\n--------- INFEASIBLE BOUNDS ---------\n")
            print_infeasible_bounds(model)

            print("\n--------- CHECKING JACOBIAN ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CHECKING JACOBIAN ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)
        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)
        raise RuntimeError(msg)
    else:
        print(msg)
        assert False


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
    print(
        f'{"UF Inlet":<20s}{m.fs.UF.feed.properties[0.0].flow_mass_comp["H2O"].value:<30.3f}'
    )
    print(
        f'{"UF Product":<20s}{m.fs.UF.product.properties[0.0].flow_mass_comp["H2O"].value:<30.3f}'
    )
    print(
        f'{"UF Disposal":<20s}{m.fs.UF.disposal.properties[0.0].flow_mass_comp["H2O"].value:<30.3f}'
    )

    print(
        f'{"Pump Inlet":<20s}{m.fs.pump.control_volume.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.pump.control_volume.properties_in[0.0].pressure, to_units=pyunits.bar)():<30.1f}'
    )
    print(
        f'{"Pump Outlet":<20s}{m.fs.pump.control_volume.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.pump.control_volume.properties_out[0.0].pressure, to_units=pyunits.bar)():<30.1f}'
    )
    print(m.fs.pump.report())
    display_flow_table(m.fs.RO)
    print("\n\n")


def display_system_build(m):
    blocks = []
    for unit in m.fs.component_data_objects(
        ctype=Block, active=True, descend_into=False
    ):
        print(unit)
        for component in unit.component_data_objects(
            ctype=Block, active=True, descend_into=False
        ):
            print("   ", component)


def display_costing_breakdown(m):
    print("\n\n-------------------- SYSTEM COSTING BREAKDOWN --------------------\n\n")
    header = f'{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<25s}'
    print(header)
    print(
        f'{"Product Flow":<35s}{f"{value(pyunits.convert(m.fs.product.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
    )
    print(f'{"LCOW":<34s}{f"${m.fs.costing.LCOW():<25.3f}"}{"$/m3":<25s}')
    print("\n")
    print_RO_costing_breakdown(m.fs.RO)
    print_softening_costing_breakdown(m.fs.softener)
    print_UF_costing_breakdown(m.fs.UF)
    print_DWI_costing_breakdown(m.fs.DWI)
    print(
        f'{"Pump Capital Cost":<35s}{f"${value(m.fs.pump.costing.capital_cost):<25,.0f}"}'
    )
    print("\n")
    print(
        f'{"Total Capital Cost":<35s}{f"${m.fs.costing.total_capital_cost():<25,.3f}"}'
    )
    print(
        f'{"Total Operating Cost":<35s}{f"${m.fs.costing.total_operating_cost():<25,.3f}"}'
    )
    print(
        f'{"Total Annualized Cost":<35s}{f"${m.fs.costing.total_annualized_cost():<25,.3f}"}'
    )
    print("\nOperating Costs:")
    print(
        f'{f"Total Fixed Operating Costs":<35s}{f"${m.fs.costing.total_fixed_operating_cost():<25,.3f}"}'
    )
    print(
        f'{f"    Agg Fixed Operating Costs":<35s}{f"${m.fs.costing.aggregate_fixed_operating_cost():<25,.3f}"}'
    )
    print(
        f'{f"    MLC Operating Costs":<35s}{f"${m.fs.costing.maintenance_labor_chemical_operating_cost():<25,.3f}"}'
    )
    print(
        f'{f"Total Variable Operating Costs":<35s}{f"${m.fs.costing.total_variable_operating_cost():<25,.3f}"}'
    )
    print(
        f'{f"    Agg Variable Operating Costs":<35s}{f"${m.fs.costing.aggregate_variable_operating_cost():<25,.3f}"}'
    )
    for flow in m.fs.costing.aggregate_flow_costs:
        print(
            f'{f"    Flow Cost [{flow}]":<35s}{f"${m.fs.costing.aggregate_flow_costs[flow]():>25,.0f}"}'
        )


def print_all_results(m):
    display_system_stream_table(m)
    report_softener(m)
    report_UF(m, m.fs.UF)
    report_RO(m, m.fs.RO)
    report_DWI(m.fs.DWI)
    display_costing_breakdown(m)


if __name__ == "__main__":
    main()
