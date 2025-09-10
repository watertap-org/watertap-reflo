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
from idaes.core.util.initialization import propagate_state
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


def main():

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    optimize(m, ro_mem_area=None, water_recovery=0.8)
    solve(m, debug=True)
    print_all_results(m)
    m.fs.costing.SEC.display()


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
            "SO4_2-",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )

    m.fs.RO_properties = NaClParameterBlock()
    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])

    build_treatment(m)

    return m


def build_treatment(m):

    m.fs.feed = Feed(property_package=m.fs.MCAS_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.sludge = Product(property_package=m.fs.MCAS_properties)
    m.fs.UF_waste = Product(property_package=m.fs.UF_properties)

    # Define the Unit Models
    m.fs.softener = FlowsheetBlock(dynamic=False)
    m.fs.UF = FlowsheetBlock(dynamic=False)
    m.fs.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.RO = FlowsheetBlock(dynamic=False)
    m.fs.DWI = FlowsheetBlock(dynamic=False)

    m.fs.MCAS_to_ZO_TDS_translator = TranslatorMCAStoZO(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.ZO_TDS_to_NaCl_translator = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_softener(m.fs.softener, prop_package=m.fs.MCAS_properties)
    build_UF(m, m.fs.UF, prop_package=m.fs.UF_properties)
    build_ro(m, m.fs.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(m, m.fs.DWI, prop_package=m.fs.RO_properties)

    m.fs.units = [
        m.fs.softener,
        m.fs.UF,
        m.fs.pump,
        m.fs.RO,
        m.fs.DWI,
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


# def build_sweep():

#     m = build_system()
#     display_system_build(m)
#     add_connections(m)
#     add_constraints(m)
#     set_operating_conditions(m)
#     apply_scaling(m)
#     init_system(m)
#     add_costing(m)
#     display_system_build(m)
#     optimize(m, ro_mem_area=None, water_recovery=0.8)

#     return m


def add_connections(m):

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.unit.inlet,
    )

    m.fs.softener_to_translator = Arc(
        source=m.fs.softener.unit.outlet,
        destination=m.fs.MCAS_to_ZO_TDS_translator.inlet,
    )

    m.fs.softener_to_sludge = Arc(
        source=m.fs.softener.unit.waste,
        destination=m.fs.sludge.inlet,
    )

    m.fs.translator_to_UF = Arc(
        source=m.fs.MCAS_to_ZO_TDS_translator.outlet,
        destination=m.fs.UF.feed.inlet,
    )

    m.fs.UF_to_translator3 = Arc(
        source=m.fs.UF.product.outlet,
        destination=m.fs.ZO_TDS_to_NaCl_translator.inlet,
    )

    m.fs.UF_to_waste = Arc(
        source=m.fs.UF.disposal.outlet,
        destination=m.fs.UF_waste.inlet,
    )

    m.fs.translator_to_pump = Arc(
        source=m.fs.ZO_TDS_to_NaCl_translator.outlet,
        destination=m.fs.pump.inlet,
    )

    m.fs.pump_to_ro = Arc(
        source=m.fs.pump.outlet,
        destination=m.fs.RO.feed.inlet,
    )

    m.fs.RO_to_product = Arc(
        source=m.fs.RO.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.RO_to_disposal = Arc(
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

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )


def apply_scaling(m):
    add_UF_scaling(m.fs.UF)
    add_ro_scaling(m.fs.RO)
    calculate_scaling_factors(m)


def add_costing(m):
    # treatment = m.fs.treatment
    # treatment.costing = TreatmentCosting()
    m.fs.costing = TreatmentCosting()
    elec_cost = value(pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018))
    m.fs.costing.electricity_cost.fix(elec_cost)

    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    add_softener_costing(m.fs.softener, costing_blk=m.fs.costing)
    add_UF_costing(m.fs.UF, costing_blk=m.fs.costing)
    add_ro_costing(m.fs.RO, costing_blk=m.fs.costing)
    add_DWI_costing(m.fs.DWI, costing_blk=m.fs.costing)

    # treatment.costing.cost_process()
    # treatment.costing.initialize()

    # energy.costing.cost_process()
    # energy.costing.initialize()

    # m.fs.costing = REFLOCosting()
    # m.fs.costing.cost_process()

    m.fs.costing.cost_process()

    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol, name="SEC")


    # m.fs.costing.initialize()
    m.fs.costing.initialize()


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


    # Convert Q_in from MGD to kg/s
    Qin = pyunits.convert(
        Qin * pyunits.Mgallon * pyunits.day**-1, to_units=pyunits.L / pyunits.s
    )
    feed_density = 1000 * pyunits.kg / pyunits.m**3
    print('\n=======> SETTING FEED CONDITIONS <======="\n')
    print(f"Flow Rate: {value(Qin):<10.2f}{pyunits.get_units(Qin)}")

    if Qin is None:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
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
        "SO4_2-": 0.23 * pyunits.kg / pyunits.m**3,
    }

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

    feed_temperature = 273.15 + 20

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

    set_inlet_conditions(m, Qin=4)
    set_softener_op_conditions(m.fs.softener)
    set_UF_op_conditions(m.fs.UF)
    m.fs.pump.efficiency_pump.fix(pump_efi)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m.fs.RO, mem_area=10000)


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


    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_softener)
    report_MCAS_stream_conc(m, m.fs.feed.properties[0.0])

    init_softener(m.fs.softener)
    propagate_state(m.fs.softener_to_translator)
    propagate_state(m.fs.softener_to_sludge)
    m.fs.sludge.initialize()

    m.fs.MCAS_to_ZO_TDS_translator.initialize()
    propagate_state(m.fs.translator_to_UF)
    init_UF(m, m.fs.UF)
    propagate_state(m.fs.UF_to_translator3)
    propagate_state(m.fs.UF_to_waste)
    m.fs.UF_waste.initialize()

    m.fs.ZO_TDS_to_NaCl_translator.initialize()

    propagate_state(m.fs.translator_to_pump)
    m.fs.pump.initialize()

    propagate_state(m.fs.pump_to_ro)

    init_ro_system(m.fs.RO)
    propagate_state(m.fs.RO_to_product)
    propagate_state(m.fs.RO_to_disposal)

    m.fs.product.initialize()
    init_DWI(m, m.fs.DWI)


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
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        m.fs.pump.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        m.fs.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for idx, stage in m.fs.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for idx, stage in m.fs.RO.stage.items():
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
    # print_softening_costing_breakdown(m.fs.treatment.softener)
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
    # display_system_stream_table(m)
    # report_softener(m)
    report_UF(m.fs.UF)
    report_RO(m.fs.RO)
    report_DWI(m.fs.DWI)
    print(m.fs.pump.display())
    print(m.fs.product.display())
    # print(
    #     value(
    #         pyunits.convert(
    #             m.fs.treatment.pump.control_volume.work[0], to_units=pyunits.kW
    #             )
    #         ) /
    #     value(
    #         pyunits.convert(
    #             m.fs.treatment.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.m**3 / pyunits.hr
    #             )
    #         )
    # )
    print(
        value(
            pyunits.convert(
                m.fs.costing.aggregate_flow_electricity, to_units=pyunits.kW
            )
        )
        / value(
            pyunits.convert(
                m.fs.product.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
    )
    display_costing_breakdown(m)


if __name__ == "__main__":
    main()
    # build_sweep()
