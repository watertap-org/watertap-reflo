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
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap.core.wt_database import Database

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
    REFLOSystemCosting,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.analysis.case_studies.KBHDP import *
from watertap_contrib.reflo.solar_models.surrogate.pv import PVSurrogate


def propagate_state(arc):
    _prop_state(arc)
    print(f"\nPropogation of {arc.source.name} to {arc.destination.name} successful.")
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
    apply_scaling(m)
    init_system(m)
    solve(m, debug=True)

    add_costing(m)
    # # get_scaling_factors(m)
    # # set_pv_constraints(m, focus="Energy")
    scale_costing(m)

    optimize(m, ro_mem_area=None, water_recovery=0.8, objective=None)
    solve(m, debug=True)



    report_RO(m, m.fs.treatment.RO)
    # # report_pump(m, m.fs.treatment.pump)
    # report_PV(m)
    # # m.fs.treatment.costing.display()
    # # m.fs.energy.costing.display()""
    # # m.fs.costing.display()
    # display_costing_breakdown(m)
    # # print(m.fs.energy.pv.display())


def build_system():
    m = ConcreteModel()
    m.db = Database()
    m.fs = FlowsheetBlock(dynamic=False)

    treatment = m.fs.treatment = Block()
    energy = m.fs.energy = Block()

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

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.RO_properties)
    treatment.disposal = Product(property_package=m.fs.RO_properties)

    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
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

    build_ec(m, treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(m, treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(m, treatment.RO, prop_package=m.fs.RO_properties)
    build_DWI(m, treatment.DWI, prop_package=m.fs.RO_properties)
    build_pv(m)

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
    treatment = m.fs.treatment

    treatment.feed_to_translator = Arc(
        source=treatment.feed.outlet,
        destination=treatment.MCAS_to_TDS_translator.inlet,
    )

    treatment.translator_to_EC = Arc(
        source=treatment.MCAS_to_TDS_translator.outlet,
        destination=treatment.EC.feed.inlet,
    )

    treatment.EC_to_UF = Arc(
        source=treatment.EC.product.outlet,
        destination=treatment.UF.feed.inlet,
    )

    treatment.UF_to_translator3 = Arc(
        source=treatment.UF.product.outlet,
        destination=treatment.TDS_to_NaCl_translator.inlet,
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
    energy = m.fs.energy

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

    m.fs.annual_treatment_energy = Var(
        initialize=10000000,
        domain=NonNegativeReals,
        units=pyunits.kWh / pyunits.year,
        doc="Annual Energy Consumption of Treatment System",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=treatment.feed.properties[0].flow_vol * m.fs.water_recovery
        == treatment.product.properties[0].flow_vol
    )


def add_treatment_costing(m):
    treatment = m.fs.treatment
    treatment.costing = TreatmentCosting()

    # print_fixed_and_unfixed_vars(treatment.costing)

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )
    add_ec_costing(m, treatment.EC, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_ro_costing(m, treatment.RO, treatment.costing)
    add_DWI_costing(m, treatment.DWI, treatment.costing)

    treatment.costing.ultra_filtration.capital_a_parameter.fix(500000)
    treatment.costing.total_investment_factor.fix(1)
    treatment.costing.maintenance_labor_chemical_factor.fix(0)

    treatment.costing.cost_process()
    treatment.costing.initialize()

    m.fs.annual_treatment_energy = Expression(
        expr=pyunits.convert(
            m.fs.treatment.costing.aggregate_flow_electricity,
            to_units=pyunits.kWh / pyunits.year,
        )
    )

    # print_fixed_and_unfixed_vars(treatment.costing)


def add_energy_costing(m):
    energy = m.fs.energy
    energy.costing = EnergyCosting()

    energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
    )

    # energy.costing.total_investment_factor.fix(1)
    # energy.costing.maintenance_labor_chemical_factor.fix(0)

    # set_pv_constraints(m, focus="Energy")

    energy.costing.cost_process()
    energy.costing.initialize()

    # set_pv_constraints(m)


def add_costing(m):
    treatment = m.fs.treatment
    energy = m.fs.energy

    add_treatment_costing(m)
    add_energy_costing(m)

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.costing.cost_process()

    m.fs.costing.add_annual_water_production(treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(treatment.product.properties[0].flow_vol)

    m.fs.costing.initialize()


def get_scaling_factors(m):
    for var in [m.fs.treatment.costing.aggregate_flow_electricity]:
        val = value(var)
        scale = calc_scale(val)
        sf = iscale.get_scaling_factor(var)
        if sf is None:
            sf = scale
            iscale.set_scaling_factor(var, sf)

        print(f"{var.name:<50s}{val:<20.3f}{scale:<20.3f}{sf:<20.3f}")


def scale_costing(m):
    treatment = m.fs.treatment
    energy = m.fs.energy

    iscale.set_scaling_factor(m.fs.energy.pv.electricity, 1e-10)
    iscale.set_scaling_factor(m.fs.energy.pv.annual_energy, 1e-10)
    iscale.set_scaling_factor(m.fs.energy.pv.costing.annual_generation, 1e-10)
    iscale.set_scaling_factor(m.fs.energy.pv.costing.system_capacity, 1e-6)

    # iscale.constraint_scaling_transform(
    #     m.fs.energy.pv.costing.annual_generation_constraint, 1e-8
    # )

    iscale.constraint_scaling_transform(m.fs.energy.pv.electricity_constraint, 1e-3)
    # iscale.constraint_scaling_transform(m.fs.energy.pv.costing.system_capacity_constraint, 1e-3)

    # agg_elec_scale = calc_scale(pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)())
    # iscale.set_scaling_factor(m.fs.energy.pv.annual_energy, 1e-8)
    # iscale.constraint_scaling_transform(energy.pv_design_constraint, 1e-2)
    # for e in m.fs.treatment.RO.stage[1].module.feed_side.eq_K:
    #     iscale.constraint_scaling_transform(m.fs.treatment.RO.stage[1].module.feed_side.eq_K[e], 1e6)

def apply_system_scaling(m):
    iscale.set_scaling_factor(m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3)
    iscale.constraint_scaling_transform(m.fs.treatment.product.properties[0.0].eq_flow_vol_phase["Liq"], 1e-6)
    iscale.set_scaling_factor(m.fs.treatment.pump.control_volume.work, 1e-6)

def apply_scaling(m):

    add_ec_scaling(m, m.fs.treatment.EC)
    add_ro_scaling(m, m.fs.treatment.RO)
    add_pv_scaling(m, m.fs.energy.pv)
    apply_system_scaling(m)
    iscale.calculate_scaling_factors(m)


def define_inlet_composition(m):
    import watertap.core.zero_order_properties as prop_ZO

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
):
    """Sets operating condition for the PV-RO system

    Args:
        m (obj): Pyomo model
        flow_in (float, optional): feed volumetric flow rate [m3/s]. Defaults to 1e-2.
        conc_in (int, optional): solute concentration [g/L]. Defaults to 30.
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    treatment = m.fs.treatment

    # Convert Q_in from MGD to kg/s
    Qin = pyunits.convert(Qin*pyunits.Mgallon * pyunits.day ** -1, to_units=pyunits.L / pyunits.s)
    feed_density = 1000 * pyunits.kg / pyunits.m**3
    print('\n=======> SETTING FEED CONDITIONS <======="\n')
    print(f"Flow Rate: {value(Qin):<10.2f}{pyunits.get_units(Qin)}")

    if Qin is None:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(Qin * feed_density)

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
    treatment.feed.pressure[0].fix(supply_pressure)
    treatment.feed.temperature[0].fix(feed_temperature)
    # m.fs.disposal.pressure[0].fix(101356)
    # m.fs.disposal.temperature[0].fix(feed_temperature)

    # m.fs.primary_pump.efficiency_pump.fix(0.85)
    

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


def set_operating_conditions(m, RO_pressure=30e5, supply_pressure=1.1e5):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4, supply_pressure=1e5)
    set_ec_operating_conditions(m, treatment.EC)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(pump_efi)
    treatment.pump.control_volume.properties_in[0].pressure.fix(supply_pressure)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m, treatment.RO, mem_area=10000)
    # # set__ED_op_conditions


def initialize_energy(m, train=False):
    if train:
        train_pv_surrogate(m)

    m.fs.energy.pv.load_surrogate()


def init_treatment(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    treatment.feed.initialize(optarg=optarg)
    propagate_state(treatment.feed_to_translator)

    treatment.MCAS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_EC)

    init_ec(m, treatment.EC)
    propagate_state(treatment.EC_to_UF)

    init_UF(m, treatment.UF)
    propagate_state(treatment.UF_to_translator3)

    treatment.TDS_to_NaCl_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_pump)

    treatment.pump.initialize(optarg=optarg)
    propagate_state(treatment.pump_to_ro)

    init_ro_system(m, treatment.RO)
    propagate_state(treatment.ro_to_product)
    propagate_state(treatment.ro_to_disposal)

    treatment.product.initialize(optarg=optarg)
    treatment.disposal.initialize(optarg=optarg)
    display_system_stream_table(m)


def init_system(m, verbose=True, solver=None):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    # assert_no_degrees_of_freedom(m)
    initialize_energy(m)
    init_treatment(m)


def solve(m, solver=None, tee=True, raise_on_failure=True, debug=False):
    # ---solving---
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000
        optarg = {
            # "constr_viol_tol": 1e-8,
            # "nlp_scaling_method": "user-scaling",
            # "linear_solver": "ma57",
            # "OF_ma57_automatic_scaling": "yes",
            # "max_iter": 350,
            # "tol": 1e-8,
            "halt_on_ampl_error": "yes",
        }

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            print("\n--------- TREATMENT ---------\n")
            check_jac(m.fs.treatment)
            print("\n--------- ENERGY ---------\n")
            check_jac(m.fs.energy)

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(m)
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print('\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print('\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print('\n{"=======> CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print("\n--------- FAILED SOLVE!!! ---------\n")
        print(msg)
        assert False


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    treatment = m.fs.treatment
    energy = m.fs.energy
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    else:
        m.fs.membrane_area_objective = Objective(expr=treatment.RO.stage[1].module.area)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.01
        upper_bound = 0.99
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
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


def report_pump(m, pump):
    print(f"\n\n-------------------- SYSTEM PUMP --------------------\n\n")
    print(
        f'{"Pump Pressure":<30s}{value(pyunits.convert(pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}'
    )
    print(
        f'{"Pump Work":<30s}{value(pyunits.convert(pump.control_volume.work[0], to_units=pyunits.kW)):<10.3f}{"kW"}'
    )


def display_system_stream_table(m):
    treatment = m.fs.treatment
    print("\n\n-------------------- SYSTEM STREAM TABLE --------------------\n\n")
    print(
        f'{"NODE":<20s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<20s}{treatment.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(treatment.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}'
    )
    display_flow_table(treatment.RO)
    print("\n\n")


def display_system_build(m):
    blocks = []
    for v in m.fs.component_data_objects(ctype=Block, active=True, descend_into=False):
        print(v)


def display_costing_breakdown(m):
    header = f'\n{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<25s}'
    print(header)
    print(
        f'{"Product Flow":<35s}{f"{value(pyunits.convert(m.fs.treatment.product.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
    )
    print(f'{"LCOW":<34s}{f"${m.fs.costing.LCOW():<25.3f}"}{"$/m3":<25s}\n')

    print_RO_costing_breakdown(m.fs.treatment.RO)
    print_EC_costing_breakdown(m.fs.treatment.EC)
    print_UF_costing_breakdown(m.fs.treatment.UF)
    print_DWI_costing_breakdown(m.fs.treatment.DWI)
    print_PV_costing_breakdown(m.fs.energy.pv)


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()

# BUG: RO pressure and membrane area are changing based on the pv energy!!!
