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
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock as vapor_prop

_log = idaeslog.getLogger(__name__)


def propagate_state(arc, detailed=False):
    _prop_state(arc)

    if detailed:
        print(
            f"\nPropogation of {arc.source.name} to {arc.destination.name} successful."
        )
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print("\n")


def main(frac_heat_from_grid=0.75):

    m = build_system(RE=True)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    solve(m.fs.treatment, debug=False)
    solve(m.fs.treatment, debug=False)

    heat_load = value(
        pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat, to_units=pyunits.MW)
    )
    pysam_results = run_pysam_fpc_model(m, heat_load=heat_load)
    m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
    m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
    m.fs.energy.FPC.heat_load.fix(heat_load)
    print(f"dof = {degrees_of_freedom(m)}")
    results = solve(m)
    display_costing_breakdown(m)

    heat_annual_required = value(
        pyunits.convert(
            m.fs.treatment.costing.aggregate_flow_heat,
            to_units=pyunits.kWh * pyunits.year**-1,
        )
    )  * (1 - frac_heat_from_grid)

    heat_load_required, pysam_results = get_fpc_heat_load(
        m, heat_annual_required, increment_heat_load=1
    )
    print(heat_annual_required, pysam_results["heat_annual"], heat_load_required)
    m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
    m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
    m.fs.energy.FPC.heat_load.fix(heat_load_required)
    assert degrees_of_freedom(m) == 0
    results = solve(m)

    display_costing_breakdown(m)

    return m


def build_sweep(
    frac_heat_from_grid=None,
    heat_load=None,
):
    m = build_system(RE=True)
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    apply_scaling(m)
    init_system(m)
    add_costing(m)
    solve(m.fs.treatment, debug=False)
    solve(m.fs.treatment, debug=False)
    if frac_heat_from_grid is not None:
        m.fs.costing.frac_heat_from_grid.set_value(frac_heat_from_grid)
        m.fs.costing.frac_heat_from_grid.setub(None)
    
        heat_annual_required = value(
            pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat, to_units=pyunits.kWh * pyunits.year**-1)
        ) * (1 - value(m.fs.costing.frac_heat_from_grid))
        # pysam_results = run_pysam_fpc_model(m, heat_load=heat_load)
        heat_load_required, pysam_results = get_fpc_heat_load(
            m, heat_annual_required, increment_heat_load=1
        )
        m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
        m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
        m.fs.energy.FPC.heat_load.fix(heat_load_required)

    elif heat_load is not None:
        pysam_results = run_pysam_fpc(m, heat_load=heat_load)
        m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
        m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
        m.fs.energy.FPC.heat_load.fix(heat_load)
    else:

        heat_load = value(
            pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat, to_units=pyunits.MW)
        )
        pysam_results = run_pysam_fpc(m, heat_load=heat_load)
        m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
        m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
        m.fs.energy.FPC.heat_load.fix(heat_load)
    
    results = solve(m)

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

    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])
    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = vapor_prop()

    build_treatment(m)
    build_energy(m)

    return m


def build_treatment(m):
    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.product = Product(property_package=m.fs.liquid_prop)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.liquid_prop)
    treatment.LTMED = FlowsheetBlock(dynamic=False)
    treatment.DWI = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_TDS_translator = Translator_MCAS_to_TDS(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    treatment.TDS_to_TDS_translator = Translator_TDS_to_TDS(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.liquid_prop,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    build_ec(m, treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(m, treatment.UF, prop_package=m.fs.UF_properties)
    build_LTMED(m, treatment.LTMED, m.fs.liquid_prop, m.fs.vapor_prop)
    build_DWI(m, treatment.DWI, prop_package=m.fs.liquid_prop)

    m.fs.units = [
        treatment.feed,
        treatment.EC,
        treatment.UF,
        treatment.pump,
        treatment.LTMED,
        treatment.DWI,
        treatment.product,
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

    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def build_energy(m):
    energy = m.fs.energy = Block()
    build_fpc_pysam(m)


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
        destination=treatment.TDS_to_TDS_translator.inlet,
    )

    treatment.UF_to_waste = Arc(
        source=treatment.UF.disposal.outlet,
        destination=treatment.UF_waste.inlet,
    )

    treatment.translator_to_pump = Arc(
        source=treatment.TDS_to_TDS_translator.outlet,
        destination=treatment.pump.inlet,
    )

    treatment.pump_to_LTMED = Arc(
        source=treatment.pump.outlet,
        destination=treatment.LTMED.feed.inlet,
    )

    treatment.LTMED_to_product = Arc(
        source=treatment.LTMED.product.outlet,
        destination=treatment.product.inlet,
    )

    treatment.LTMED_to_dwi = Arc(
        source=treatment.LTMED.disposal.outlet,
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
    elec_cost = pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)()
    treatment.costing.electricity_cost.fix(elec_cost)

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )
    add_ec_costing(m, treatment.EC, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_LTMED_costing(m, treatment.LTMED, treatment.costing)
    add_DWI_costing(m, treatment.DWI, treatment.costing)

    treatment.costing.cost_process()
    treatment.costing.initialize()


def add_energy_costing(m):
    energy = m.fs.energy
    energy.costing = EnergyCosting()

    elec_cost = pyunits.convert(0.066 * pyunits.USD_2023, to_units=pyunits.USD_2018)()
    m.fs.energy.costing.electricity_cost.fix(elec_cost)

    add_fpc_pysam_costing(m, costing_block=energy.costing)

    energy.costing.cost_process()
    energy.costing.add_LCOH()
    energy.costing.initialize()


def add_costing(m):
    treatment = m.fs.treatment
    energy = m.fs.energy

    add_treatment_costing(m)
    add_energy_costing(m)

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.electricity_cost_buy.set_value(0.066)
    m.fs.costing.heat_cost_buy.set_value(0.00894)
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
    # TBD

    print(m.fs.costing.display())

    # iscale.set_scaling_factor(m.fs.costing.aggregate_flow_electricity_sold, 0)
    # iscale.set_scaling_factor(m.fs.costing.aggregate_flow_heat_sold, 0)
    # iscale.set_scaling_factor(m.fs.costing.aggregate_flow_heat_purchased, 1e-6)

    # iscale.constraint_scaling_transform(m.fs.costing.aggregate_electricity_balance, 0)
    # iscale.constraint_scaling_transform(m.fs.costing.aggregate_heat_balance, 0)
    # iscale.constraint_scaling_transform(m.fs.costing.aggregate_heat_complement, 0)

    # iscale.calculate_scaling_factors(m.fs.treatment.costing)
    # iscale.calculate_scaling_factors(m.fs.energy.costing)
    # iscale.calculate_scaling_factors(m.fs.costing)


def apply_system_scaling(m):
    iscale.set_scaling_factor(
        m.fs.treatment.sludge.properties[0.0].flow_mass_comp["tss"], 1e1
    )
    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )

    iscale.set_scaling_factor(
        m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3
    )
    iscale.set_scaling_factor(
        m.fs.treatment.product.properties[0.0].dens_mass_phase["Liq"], 1e-3
    )


def apply_scaling(m):

    add_ec_scaling(m, m.fs.treatment.EC)
    add_UF_scaling(m.fs.treatment.UF)
    add_LTMED_scaling(m, m.fs.treatment.LTMED)
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


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def set_operating_conditions(m, RO_pressure=101325):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=4)
    set_ec_operating_conditions(m, treatment.EC)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(pump_efi)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_LTMED_operating_conditions(treatment.LTMED)
    set_fpc_pysam_op_conditions(m)


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

    treatment.MCAS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_EC)

    init_ec(m, treatment.EC)
    propagate_state(treatment.EC_to_UF)

    init_UF(m, treatment.UF)
    propagate_state(treatment.UF_to_translator3)
    propagate_state(treatment.UF_to_waste)

    treatment.TDS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(treatment.translator_to_pump)

    treatment.pump.initialize(optarg=optarg)

    propagate_state(treatment.pump_to_LTMED)

    init_LTMED(m, treatment.LTMED)
    propagate_state(treatment.LTMED_to_product)

    propagate_state(treatment.LTMED_to_dwi)

    treatment.product.initialize(optarg=optarg)
    init_DWI(m, treatment.DWI)
    # display_system_stream_table(m)


def init_system(m, verbose=True, solver=None):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    if degrees_of_freedom(m) != 0:
        breakdown_dof(m, detailed=True)
    assert_no_degrees_of_freedom(m)
    init_treatment(m)


def solve(m, solver=None, tee=False, raise_on_failure=False, debug=False):
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

            print("\n--------- COSTING ---------\n")
            check_jac(m.fs.costing)

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
        # assert False
        return results


def set_prob_for_box_solve(m):
    treatment = m.fs.treatment
    treatment.pump.control_volume.properties_out[0].pressure.unfix()
    for idx, stage in treatment.RO.stage.items():
        stage.module.recovery_vol_phase[0.0, "Liq"].fix(0.5)


def box_solve_problem(m):
    # set_prob_for_box_solve(m)
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
    water_recovery=None,
    fixed_pressure=None,
    grid_frac_heat=None,
    heat_price=None,
    objective="LCOT",
):
    treatment = m.fs.treatment
    energy = m.fs.energy
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    if objective == "LCOT":
        m.fs.lcot_objective = Objective(expr=m.fs.costing.LCOT)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.treatment.LTMED.unit.recovery_vol_phase[0.0, "Liq"].unfix()
        m.fs.water_recovery.fix(water_recovery)

    if grid_frac_heat is not None:
        energy.FPC.heat_load.unfix()
        energy.FPC.hours_storage.unfix()
        m.fs.costing.frac_heat_from_grid.fix(grid_frac_heat)

    if heat_price is not None:
        energy.FPC.heat_load.unfix()
        m.fs.costing.frac_heat_from_grid.unfix()
        m.fs.costing.heat_cost_buy.fix(heat_price)

    print(f"Degrees of Feedom: {degrees_of_freedom(m)}")
    assert degrees_of_freedom(m) >= 0


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
    else:
        print("Variables are scaled well")


def display_system_stream_table(m):
    treatment = m.fs.treatment
    print("\n\n-------------------- SYSTEM STREAM TABLE --------------------\n\n")
    print(
        f'{"NODE":<20s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<20s}{treatment.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.1f}{value(pyunits.convert(treatment.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}'
    )
    print(
        f'{"Product":<20s}{treatment.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.1f}{value(pyunits.convert(treatment.product.properties[0.0].pressure, to_units=pyunits.bar)):<20.1f}{treatment.product.properties[0.0].flow_mass_phase_comp["Liq", "TDS"].value:<30.3f}{treatment.product.properties[0.0].conc_mass_phase_comp["Liq", "TDS"].value:<30.3f}'
    )
    print(
        f'{"Disposal":<20s}{treatment.DWI.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.1f}{value(pyunits.convert(treatment.DWI.feed.properties[0.0].pressure, to_units=pyunits.bar)):<20.1f}{treatment.DWI.feed.properties[0.0].flow_mass_phase_comp["Liq", "TDS"].value:<30.3f}{treatment.DWI.feed.properties[0.0].conc_mass_phase_comp["Liq", "TDS"].value:<30.3f}'
    )
    print("\n\n")

    report_EC(treatment.EC)
    report_UF(m, treatment.UF)
    report_pump(m, treatment.pump)
    report_LTMED(m)
    report_DWI(treatment.DWI)
    report_fpc(m)


def display_system_build(blk):
    blocks = []
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=False):
        print(v)


def display_costing_breakdown(m):
    energy = m.fs.energy
    treatment = m.fs.treatment
    header = f'\n{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<25s}'
    print(header)
    # print(
    #     f'{"Product Flow":<35s}{f"{value(pyunits.convert(m.fs.treatment.product.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
    # )
    print(f'{"LCOW":<34s}{f"${m.fs.costing.LCOW():<25.3f}"}{"$/m3":<25s}\n')

    print_EC_costing_breakdown(m.fs.treatment.EC)
    print_UF_costing_breakdown(m.fs.treatment.UF)
    print_DWI_costing_breakdown(m.fs.treatment.DWI)
    print_FPC_costing_breakdown(m, m.fs.energy.FPC)

    # print(
    #     f'{"Grid Heat Frac.":<30s}{value(m.fs.costing.frac_heat_from_grid):<20,.2f}{pyunits.get_units(energy.costing.frac_heat_from_grid)}'
    # )
    print(
        f'{"Treatment Agg. Flow Heat":<30s}{value(treatment.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(treatment.costing.aggregate_flow_heat)}'
    )
    print(
        f'{"Energy Agg. Flow Heat":<30s}{value(energy.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(energy.costing.aggregate_flow_heat)}'
    )
    print(
        f'{"Agg. Flow Heat Purchased":<30s}{value(m.fs.costing.aggregate_flow_heat_purchased):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_heat_purchased)}'
    )
    print("")
    print(
        f'{"Treatment Agg. Flow Elec.":<30s}{value(treatment.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(treatment.costing.aggregate_flow_electricity)}'
    )
    print(
        f'{"Energy Agg. Flow Elec.":<30s}{value(energy.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(energy.costing.aggregate_flow_electricity)}'
    )
    print(
        f'{"Agg. Flow Elec. Purchased":<30s}{value(m.fs.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_electricity)}'
    )
    print("")
    print(
        f'{"Grid Frac Elec.":<30s}{value(m.fs.costing.frac_elec_from_grid):<20,.2f}{pyunits.get_units(m.fs.costing.frac_elec_from_grid)}'
    )
    print(
        f'{"Grid Frac Heat":<30s}{value(m.fs.costing.frac_heat_from_grid):<20,.2f}{pyunits.get_units(m.fs.costing.frac_heat_from_grid)}'
    )


if __name__ == "__main__":

    # m = main()
    # m.fs.costing.display()
    m = build_sweep(frac_heat_from_grid=None, heat_load=33.5)
    results = solve(m, raise_on_failure=True)
    display_costing_breakdown(m)
    m.fs.costing.LCOT.display()
    m.fs.energy.FPC.storage_volume.display()
    # m.fs.energy.FPC.display()
    # m.fs.energy.costing.display()
    # m.fs.energy.costing.LCOH.display()

