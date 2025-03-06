import os
import math
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap_contrib.reflo.core import REFLODatabase
import idaes.logger as idaeslog
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import *
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
from watertap.core.zero_order_properties import WaterParameterBlock

_log = idaeslog.getLogger(__name__)


def propagate_state(arc, detailed=True):
    _prop_state(arc)

    if detailed:
        print(
            f"\nPropogation of {arc.source.name} to {arc.destination.name} successful."
        )
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print("\n")


def main():
    file_dir = os.path.dirname(os.path.abspath(__file__))

    m = build_system(RE=True)
    display_system_build(m)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    scale_costing(m)
    # box_solve_problem(m)
    # solve(m, debug=True)
    # scale_costing(m)
    optimize(m, ro_mem_area=20000, water_recovery=0.8, grid_frac=0.5, objective="LCOW")
    solve(m, debug=True)
    # # display_flow_table(m)
    # display_system_stream_table(m)
    # report_RO(m, m.fs.treatment.RO)
    # # # # # report_pump(m, m.fs.treatment.pump)
    # # report_PV(m)
    # # # # # m.fs.treatment.costing.display()
    # # # # # m.fs.energy.costing.display()""
    # # # # # m.fs.costing.display()
    # display_costing_breakdown(m)
    # # # # # print(m.fs.energy.pv.display())
    # # # print_system_scaling_report(m)
    report_PV(m)
    report_pump(m, m.fs.treatment.pump)
    print(m.fs.costing.frac_elec_from_grid.display())
    print(m.fs.costing.aggregate_flow_electricity_purchased.display())

    return m


def build_sweep(
    grid_frac=None,
    elec_price=None,
    water_recovery=None,
    ro_mem_area=None,
    objective="LCOT",
):

    m = build_system(RE=True)
    display_system_build(m)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    scale_costing(m)
    optimize(
        m,
        ro_mem_area=ro_mem_area,
        water_recovery=water_recovery,
        grid_frac=grid_frac,
        elec_price=elec_price,
        objective=objective,
    )

    return m


def build_system(RE=True):
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
    # if RE:
    #     print('Building System with Renewable Energy')
    #     m.fs.RE = RE
    build_energy(m)
    # else:
    #     print('Building System without Renewable Energy')
    #     m.fs.RE = RE

    m.fs.RE = RE

    return m


def build_treatment(m):
    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.RO_properties)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.RO = FlowsheetBlock(dynamic=False)
    treatment.DWI = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_TDS_translator = Translator_MCAS_to_TDS(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    treatment.TDS_to_NaCl_translator = Translator_TDS_to_NACL(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    build_ec(m, treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(m, treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(m, treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(m, treatment.DWI, prop_package=m.fs.RO_properties)

    m.fs.units = [
        treatment.feed,
        treatment.EC,
        treatment.UF,
        treatment.pump,
        treatment.RO,
        treatment.DWI,
        # treatment.product,
        treatment.sludge,
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


def build_energy(m):
    energy = m.fs.energy = Block()
    build_pv(m)


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

    treatment.EC_to_sludge = Arc(
        source=treatment.EC.disposal.outlet,
        destination=treatment.sludge.inlet,
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

    treatment.ro_to_dwi = Arc(
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


def add_treatment_costing(m):
    treatment = m.fs.treatment
    treatment.costing = TreatmentCosting()

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )
    add_ec_costing(m, treatment.EC, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_ro_costing(m, treatment.RO, treatment.costing)
    add_DWI_costing(m, treatment.DWI, treatment.costing)

    treatment.costing.cost_process()
    treatment.costing.initialize()


def add_energy_costing(m):
    energy = m.fs.energy
    energy.costing = EnergyCosting()
    energy.costing.has_electricity_generation = m.fs.RE

    energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
        costing_method_arguments={
            "cost_method": "simple"
        }
    )

    energy.costing.cost_process()
    energy.costing.add_LCOE()
    energy.costing.initialize()


def add_costing(m):
    treatment = m.fs.treatment
    # energy = m.fs.energy

    add_treatment_costing(m)
    # if m.fs.RE:
    add_energy_costing(m)

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()

    m.fs.costing.add_annual_water_production(treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOT(treatment.product.properties[0].flow_vol)

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

    add_pv_costing_scaling(m, energy.costing)
    # iscale.set_scaling_factor(m.fs.energy.pv.electricity, 1e-5)
    # iscale.set_scaling_factor(m.fs.energy.pv.annual_energy, 1/10000000)
    # # iscale.set_scaling_factor(m.fs.energy.pv.costing.annual_generation, 1e-10)
    # # iscale.set_scaling_factor(m.fs.energy.pv.costing.system_capacity, 1e-6)

    # # iscale.constraint_scaling_transform(
    # #     m.fs.energy.pv.costing.annual_generation_constraint, 1e-8
    # # )

    # iscale.constraint_scaling_transform(m.fs.energy.pv.electricity_constraint, 1e-1)
    # iscale.constraint_scaling_transform(m.fs.energy.pv.costing.system_capacity_constraint, 1e-3)

    # agg_elec_scale = calc_scale(pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)())
    # iscale.set_scaling_factor(m.fs.energy.pv.annual_energy, 1e-8)
    # iscale.constraint_scaling_transform(energy.pv_design_constraint, 1e-2)
    # for e in m.fs.treatment.RO.stage[1].module.feed_side.eq_K:
    #     iscale.constraint_scaling_transform(m.fs.treatment.RO.stage[1].module.feed_side.eq_K[e], 1e6)


def apply_system_scaling(m):
    # iscale.set_scaling_factor(m.fs.treatment.sludge.properties[0.0].flow_mass_comp["H2O"], 1)
    # iscale.set_scaling_factor(m.fs.treatment.sludge.properties[0.0].flow_mass_comp["tds"], 1)
    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )

    # iscale.set_scaling_factor(
    #     m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3
    # )
    # iscale.constraint_scaling_transform(
    #     m.fs.treatment.product.properties[0.0].eq_flow_vol_phase["Liq"], 1e-6
    # )
    # iscale.set_scaling_factor(m.fs.treatment.pump.control_volume.work, 1e-6)


def apply_scaling(m):

    add_ec_scaling(m, m.fs.treatment.EC)
    add_UF_scaling(m.fs.treatment.UF)
    add_ro_scaling(m, m.fs.treatment.RO)
    # if m.fs.RE:
    add_pv_scaling(m, m.fs.energy.pv)
    apply_system_scaling(m)
    iscale.calculate_scaling_factors(m)


def set_inlet_conditions(
    m,
    Qin=None,
    Cin=None,
    water_recovery=None,
    supply_pressure=101325,
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

    # treatment.product.properties[0.0].pressure.fix(101356)
    # treatment.product.properties[0.0].temperature.fix(feed_temperature)

    # assert_units_consistent(m)


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def set_operating_conditions(m, RO_pressure=20e5):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4)
    set_ec_operating_conditions(m, treatment.EC)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(pump_efi)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m, treatment.RO, mem_area=10000)
    set_pv_constraints(m, focus="Energy")


def init_treatment(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    assert_no_degrees_of_freedom(m)
    treatment.feed.initialize(optarg=optarg)
    propagate_state(treatment.feed_to_translator)
    report_MCAS_stream_conc(m, treatment.feed.properties[0.0])
    treatment.MCAS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_EC)

    init_ec(m, treatment.EC)
    propagate_state(treatment.EC_to_UF)
    propagate_state(treatment.EC_to_sludge)
    treatment.sludge.initialize(optarg=optarg)

    init_UF(m, treatment.UF)
    propagate_state(treatment.UF_to_translator3)
    propagate_state(treatment.UF_to_waste)
    treatment.UF_waste.initialize(optarg=optarg)

    treatment.TDS_to_NaCl_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_pump)

    treatment.pump.initialize(optarg=optarg)

    propagate_state(treatment.pump_to_ro)

    init_ro_system(m, treatment.RO)
    propagate_state(treatment.ro_to_product)
    # propagate_state(treatment.ro_to_dwi)

    treatment.product.initialize(optarg=optarg)
    init_DWI(m, treatment.DWI)
    display_system_stream_table(m)


def init_system(m, verbose=True, solver=None):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    if degrees_of_freedom(m) != 0:
        breakdown_dof(m, detailed=True)
    assert_no_degrees_of_freedom(m)
    init_treatment(m)


def solve(m, solver=None, tee=True, raise_on_failure=True, debug=False):
    # ---solving---
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000

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


def set_prob_for_box_solve(m):
    treatment = m.fs.treatment
    # m.fs.water_recovery.unfix()
    # m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(75e5)
    treatment.pump.control_volume.properties_out[0].pressure.unfix()
    for idx, stage in treatment.RO.stage.items():
        stage.module.recovery_vol_phase[0.0, "Liq"].fix(0.5)


def box_solve_problem(m):
    set_prob_for_box_solve(m)
    breakdown_dof(m, detailed=False)
    fsTools.standard_solve(
        m,
        tee=True,
        check_close_to_bounds=True,
        check_var_scailing=True,
        check_dofs=True,
        expected_DOFs=0,
    )


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    grid_frac=None,
    elec_price=None,
    objective="LCOW",
):
    treatment = m.fs.treatment
    # energy = m.fs.energy
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    elif objective == "LCOT":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOT)
    else:
        m.fs.membrane_area_objective = Objective(expr=treatment.RO.stage[1].module.area)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.5
        upper_bound = 0.8
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(lower_bound)
        m.fs.water_recovery.setub(upper_bound)

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
        # for idx, stage in treatment.RO.stage.items():
        #     stage.module.area.fix(ro_mem_area)
        treatment.RO.total_membrane_area.fix(ro_mem_area)
        for idx, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)
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

    if grid_frac is not None:
        m.fs.costing.frac_elec_from_grid.fix(grid_frac)
        m.fs.energy.pv.design_size.unfix()
        m.fs.energy.pv.annual_energy.unfix()

    if elec_price is not None:
        m.fs.costing.frac_elec_from_grid.unfix()
        m.fs.energy.pv.design_size.unfix()
        m.fs.energy.pv.annual_energy.unfix()

    for idx, stage in treatment.RO.stage.items():
        stage.module.width.setub(5000)
        stage.module.feed_side.velocity[0, 0].unfix()
        stage.module.feed_side.velocity[0, 1].setlb(0.0)
        stage.module.feed_side.K.setlb(1e-6)
        stage.module.feed_side.friction_factor_darcy.setub(50)
        stage.module.flux_mass_phase_comp.setub(1)
        # stage.module.flux_mass_phase_comp.setlb(1e-5)
        stage.module.feed_side.cp_modulus.setub(10)
        stage.module.rejection_phase_comp.setlb(1e-4)
        stage.module.feed_side.N_Re.setlb(1)
        stage.module.recovery_mass_phase_comp.setlb(1e-7)


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


def print_system_scaling_report(m):
    badly_scaled_var_list = iscale.list_badly_scaled_variables(m)
    if len(badly_scaled_var_list) > 0:
        print("Variables are not scaled well")
        print(
            f'{"Variable":<83s}{"Val":<15s}{"Val Scale":<10s}{"SF":<10s}{"Diff":<10s}'
        )
        print("Treatment:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "treatment"
        ]
        print("Energy:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "energy"
        ]
        print("Costing:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "costing"
        ]

        # for var in badly_scaled_var_list:
        #     keys = var[0].name.split(".")
        #     val_scale = -1 * calc_scale(var[1])
        #     sf_scale = -1 * calc_scale(iscale.get_scaling_factor(var[0]))
        #     scale_diff = val_scale - sf_scale
        #     print(f"{var[0].name:<80s}{val_scale:<10.1f}{sf_scale:<10.1f}{scale_diff:<10.1f}")
        #     print(f"{var[0].name:<80s}{var[1]}")
        # [print(i[0].name, i[1]) for i in badly_scaled_var_list]
        # [print(f'Variable: {var.name:<80s} Value Scale:{-1*calc_scale(val):<5.1f} Scale Factor:{-1*calc_scale(iscale.get_scaling_factor(var))}') for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)]
    else:
        print("Variables are scaled well")


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

    print_EC_costing_breakdown(m.fs.treatment.EC)
    print_UF_costing_breakdown(m.fs.treatment.UF)
    print_RO_costing_breakdown(m.fs.treatment.RO)
    print_DWI_costing_breakdown(m.fs.treatment.DWI)
    print_PV_costing_breakdown(m.fs.energy.pv)


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = main()