import os

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

from watertap_contrib.reflo.core import REFLODatabase
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
    REFLOSystemCosting,
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

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    optimize(m, ro_mem_area=None, water_recovery=0.75)
    solve(m, debug=True)
    report_UF(m, m.fs.UF)

    # print(m.fs.treatment.product.display())
    # print(m.fs.treatment.product.properties[0].flow_vol_phase.display())
    # print(m.fs.treatment.costing.display())
    # print(m.fs.treatment.costing.LCOW.display())


def build_system():
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)

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

    build_treatment(m)

    return m


def build_treatment(m):
    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.RO_properties)
    treatment.sludge = Product(property_package=m.fs.MCAS_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    # Define the Unit Models
    treatment.softener = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.RO = FlowsheetBlock(dynamic=False)
    treatment.DWI = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_TDS_translator = Translator_MCAS_to_TDS(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    treatment.TDS_to_NaCl_translator = Translator_TDS_to_NACL(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_softener(m, treatment.softener, prop_package=m.fs.MCAS_properties)
    build_UF(m, treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(m, treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=2)
    build_DWI(m, treatment.DWI, prop_package=m.fs.RO_properties)

    m.fs.units = [
        treatment.softener,
        treatment.UF,
        treatment.pump,
        treatment.RO,
        treatment.DWI,
    ]

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e-2, index=("H2O"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1, index=("tds"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e5, index=("tss"))

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def build_sweep():

    m = build_system()
    display_system_build(m)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    display_system_build(m)
    optimize(m, ro_mem_area=None, water_recovery=0.75)

    return m


def add_connections(m):
    treatment = m.fs.treatment

    treatment.feed_to_softener = Arc(
        source=treatment.feed.outlet,
        destination=treatment.softener.unit.inlet,
    )

    treatment.softener_to_translator = Arc(
        source=treatment.softener.unit.outlet,
        destination=treatment.MCAS_to_TDS_translator.inlet,
    )

    treatment.softener_to_sludge = Arc(
        source=treatment.softener.unit.waste,
        destination=treatment.sludge.inlet,
    )

    treatment.translator_to_UF = Arc(
        source=treatment.MCAS_to_TDS_translator.outlet,
        destination=treatment.UF.feed.inlet,
    )

    treatment.UF_to_translator3 = Arc(
        source=treatment.UF.product.outlet,
        destination=treatment.TDS_to_NaCl_translator.inlet,
    )

    treatment.UF_to_waste = Arc(
        source=treatment.UF.disposal.outlet,
        destination=treatment.UF_waste.inlet,
    )

    treatment.translator_to_pump = Arc(
        source=treatment.TDS_to_NaCl_translator.outlet,
        destination=treatment.pump.inlet,
    )

    treatment.pump_to_ro = Arc(
        source=treatment.pump.outlet,
        destination=treatment.RO.feed.inlet,
    )

    treatment.ro_to_product = Arc(
        source=treatment.RO.product.outlet,
        destination=treatment.product.inlet,
    )

    treatment.ro_to_disposal = Arc(
        source=treatment.RO.disposal.outlet,
        destination=treatment.DWI.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_constraints(m):
    treatment = m.fs.treatment

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=treatment.feed.properties[0].flow_vol * m.fs.water_recovery
        == treatment.product.properties[0].flow_vol
    )


def apply_scaling(m):
    treatment = m.fs.treatment
    add_UF_scaling(treatment.UF)
    add_ro_scaling(m, treatment.RO)
    calculate_scaling_factors(m)


def add_costing(m):
    treatment = m.fs.treatment
    # treatment.costing = TreatmentCosting()
    treatment.costing = TreatmentCosting()
    elec_cost = pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)()
    treatment.costing.electricity_cost.fix(elec_cost)

    # energy = m.fs.energy = Block()
    # energy.costing = EnergyCosting()

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )

    add_softener_costing(m, treatment.softener, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_ro_costing(m, treatment.RO, treatment.costing)
    add_DWI_costing(m, treatment.DWI, treatment.costing)

    # treatment.costing.cost_process()
    # treatment.costing.initialize()

    # energy.costing.cost_process()
    # energy.costing.initialize()

    # m.fs.costing = REFLOCosting()
    # m.fs.costing.cost_process()

    treatment.costing.cost_process()

    treatment.costing.add_annual_water_production(
        treatment.product.properties[0].flow_vol
    )
    treatment.costing.add_LCOW(treatment.product.properties[0].flow_vol_phase["Liq"])

    # m.fs.costing.initialize()
    treatment.costing.initialize()


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
    water_recovery=None,
    supply_pressure=101325,
):

    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    treatment = m.fs.treatment

    # Convert Q_in from MGD to kg/s
    Qin = pyunits.convert(
        Qin * pyunits.Mgallon * pyunits.day**-1, to_units=pyunits.L / pyunits.s
    )
    feed_density = 1000 * pyunits.kg / pyunits.m**3
    print('\n=======> SETTING FEED CONDITIONS <======="\n')
    print(f"Flow Rate: {value(Qin):<10.2f}{pyunits.get_units(Qin)}")

    if Qin is None:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            Qin * feed_density
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

    for solute, solute_conc in inlet_dict.items():
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            pyunits.convert(
                (
                    treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
                    / (1000 * pyunits.kg / pyunits.m**3)
                )
                * solute_conc,
                to_units=pyunits.kg / pyunits.s,
            )
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute]),
            index=("Liq", solute),
        )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )

    feed_temperature = 273.15 + 20

    # # initialize feed
    treatment.feed.pressure[0].fix(supply_pressure)
    treatment.feed.temperature[0].fix(feed_temperature)


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def set_operating_conditions(m, RO_pressure=20e5, supply_pressure=1.1e5):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]

    set_inlet_conditions(m, Qin=4)
    set_softener_op_conditions(m, treatment.softener.unit)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(pump_efi)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m, treatment.RO, mem_area=10000)


def init_system(m, verbose=True, solver=None):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    if degrees_of_freedom(m) != 0:
        breakdown_dof(m, detailed=True)
    assert_no_degrees_of_freedom(m)
    init_treatment(m)


def init_treatment(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    treatment.feed.initialize()
    propagate_state(treatment.feed_to_softener)
    report_MCAS_stream_conc(m, treatment.feed.properties[0.0])

    init_softener(m, treatment.softener.unit)
    propagate_state(treatment.softener_to_translator)
    propagate_state(treatment.softener_to_sludge)
    treatment.sludge.initialize()

    treatment.MCAS_to_TDS_translator.initialize()
    propagate_state(treatment.translator_to_UF)
    init_UF(m, treatment.UF)
    propagate_state(treatment.UF_to_translator3)
    propagate_state(treatment.UF_to_waste)
    treatment.UF_waste.initialize()

    treatment.TDS_to_NaCl_translator.initialize()

    propagate_state(treatment.translator_to_pump)
    treatment.pump.initialize()

    propagate_state(treatment.pump_to_ro)

    init_ro_system(m, treatment.RO)
    propagate_state(treatment.ro_to_product)
    propagate_state(treatment.ro_to_disposal)

    treatment.product.initialize()
    init_DWI(m, treatment.DWI)


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))
    treatment = m.fs.treatment
    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=treatment.costing.LCOW)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        treatment.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        treatment.pump.control_volume.properties_out[0].pressure.unfix()
        treatment.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        treatment.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for idx, stage in treatment.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for idx, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)


def solve(model, solver=None, tee=False, raise_on_failure=True, debug=False):
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

            # print("\n--------- CHECKING JACOBIAN ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

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
    for v in m.fs.component_data_objects(ctype=Block, active=True, descend_into=False):
        print(v)


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
