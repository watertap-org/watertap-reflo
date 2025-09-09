import os
import math
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
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
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap_contrib.reflo.core import REFLODatabase
import idaes.logger as idaeslog
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.unit_specific.cryst_prop_pack import NaClParameterBlock as cryst_prop_pack
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
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

from idaes.models.unit_models import (
    Product,
    Feed,
    StateJunction,
    Separator,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core import MaterialBalanceType

def kbhdp_zld_ro(ro_recovery=0.5,Qin=4):

    m = build_zld_ro_treatment()
    add_zld_ro_connections(m)
    set_zld_ro_operating_conditions(m,Qin)
    apply_zld_ro_scaling(m)
    init_zld_ro_treatment(m, verbose=False)
    add_ro_recovery_constraint(m, m.fs.treatment.RO,ro_recovery)
    optimize_zld_ro(m, water_recovery=ro_recovery)
    solve(m, debug=False)
    
    return m

def add_ro_recovery_constraint(m, blk,ro_recovery):
    m.fs.treatment.ro_water_recovery = Var(
        initialize=ro_recovery,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="RO Water Recovery",
    )

    blk.eq_water_recovery = Constraint(
        expr=blk.feed.properties[0].flow_vol * m.fs.treatment.ro_water_recovery
        == blk.product.properties[0].flow_vol
    )


def build_zld_ro_treatment():
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
    m.fs.properties_NaCl = cryst_prop_pack()

    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.RO_waste = Product(property_package=m.fs.RO_properties)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.RO = FlowsheetBlock(dynamic=False)


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

    treatment.NaCl_cryst_product_translator = Translator_NaCl_to_NaCl(
        inlet_property_package = m.fs.RO_properties,
        outlet_property_package= m.fs.properties_NaCl,
    )

    build_ec(m, treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(m, treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(m, treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)

    m.fs.units = [
        treatment.feed,
        treatment.EC,
        treatment.UF,
        treatment.pump,
        treatment.RO,

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

    return m


def add_zld_ro_connections(m):
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
        destination=treatment.NaCl_cryst_product_translator.inlet,
    )

    treatment.ro_to_md = Arc(
        source=treatment.RO.disposal.outlet,
        destination=treatment.RO_waste.inlet,
    )


    TransformationFactory("network.expand_arcs").apply_to(m)


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

    # assert_units_consistent(m)


def set_zld_ro_operating_conditions(m, Qin,RO_pressure=20e5):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin)
    set_ec_operating_conditions(m, treatment.EC)
    set_UF_op_conditions(treatment.UF)
    treatment.pump.efficiency_pump.fix(pump_efi)
    treatment.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m, treatment.RO, mem_area=10000)


def set_ec_scaling(m, blk, calc_blk_scaling_factors=False):

    set_scaling_factor(blk.ec.properties_in[0].flow_vol, 1e-2)
    set_scaling_factor(blk.ec.properties_in[0].conc_mass_comp["tds"], 1)
    set_scaling_factor(blk.ec.charge_loading_rate, 1e-2)
    set_scaling_factor(blk.ec.cell_voltage, 1)
    set_scaling_factor(blk.ec.anode_area, 1e-2)
    set_scaling_factor(blk.ec.cathode_area, 1e-2)
    set_scaling_factor(blk.ec.current_density, 1e-2)
    set_scaling_factor(blk.ec.applied_current, 1e-6)
    set_scaling_factor(blk.ec.metal_dose, 1e4)
    set_scaling_factor(blk.ec.electrode_thick, 1e3)
    set_scaling_factor(blk.ec.electrode_mass, 1e-2)
    set_scaling_factor(blk.ec.electrode_volume, 1e1)
    set_scaling_factor(blk.ec.electrode_gap, 1e3)
    set_scaling_factor(blk.ec.floc_basin_vol, 1e-2)
    set_scaling_factor(blk.ec.ohmic_resistance, 1e7)
    set_scaling_factor(blk.ec.power_required, 1e-6)

    # Calculate scaling factors only for EC block if in full case study flowsheet
    # so we don't prematurely set scaling factors
    if calc_blk_scaling_factors:
        calculate_scaling_factors(blk)

    # otherwise calculate all scaling factors
    else:
        calculate_scaling_factors(m)


def apply_zld_ro_scaling(m):

    add_ec_scaling(m, m.fs.treatment.EC)
    set_ec_scaling(m, m.fs.treatment.EC, calc_blk_scaling_factors=False)
    add_UF_scaling(m.fs.treatment.UF)
    add_ro_scaling(m, m.fs.treatment.RO)

    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )
    iscale.calculate_scaling_factors(m)


def init_zld_ro_treatment(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    assert_no_degrees_of_freedom(m)

    # treatment.feed.properties[0].flow_vol
    # treatment.feed.properties[0].flow_vol_phase

    treatment.feed.initialize(optarg=optarg)
    propagate_state(treatment.feed_to_translator)

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
    propagate_state(treatment.ro_to_md)

    treatment.RO.disposal.properties[0].flow_vol
    treatment.RO.disposal.properties[0].flow_vol_phase
    treatment.RO.disposal.initialize()

    treatment.NaCl_cryst_product_translator.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treatment.NaCl_cryst_product_translator.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)


def solve(m, solver=None, tee=False, raise_on_failure=True, debug=False):
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
            # print("\n--------- ENERGY ---------\n")
            # check_jac(m.fs.energy)

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


def optimize_zld_ro(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    treatment = m.fs.treatment
    # energy = m.fs.energy
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    # if objective == "LCOW":
    #     m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    # elif objective == "LCOT":
    #     m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOT)
    # else:
    #     m.fs.membrane_area_objective = Objective(expr=treatment.RO.stage[1].module.area)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.treatment.ro_water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.5
        upper_bound = 0.8
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.treatment.ro_water_recovery.unfix()
        m.fs.treatment.ro_water_recovery.setlb(lower_bound)
        m.fs.treatment.ro_water_recovery.setub(upper_bound)

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


def add_zld_md(m=None, Q_md=0.22478, Cin=118, md_water_recovery=0.5):
    if m == None:
        m = ConcreteModel()
        # m.db = REFLODatabase()
        m.fs = FlowsheetBlock(dynamic=False)

        treat = m.fs.treatment = Block()
    
    else:
        treat = m.fs.treatment 

    m.fs.properties_md = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.inlet_flow_rate = pyunits.convert(
        Q_md * pyunits.m**3 / pyunits.s, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    
    m.water_recovery = md_water_recovery


    m.fs.treatment.sw_to_nacl_product = Translator_SW_to_NaCl(
        inlet_property_package = m.fs.properties_md,
        outlet_property_package= m.fs.properties_NaCl,
    )

    m.fs.treatment.sw_to_nacl_disposal = Translator_SW_to_NaCl(
        inlet_property_package = m.fs.properties_md,
        outlet_property_package= m.fs.properties_NaCl,
    )

    m.fs.treatment.md = FlowsheetBlock(dynamic=False)
    build_md(m, m.fs.treatment.md, m.fs.properties_md)

    treat.md_to_product = Arc(
        source=treat.md.permeate.outlet, destination=treat.sw_to_nacl_product.inlet
    )

    treat.md_disposal_to_nacl_translator = Arc(
        source=treat.md.concentrate.outlet,
        destination=treat.sw_to_nacl_disposal.inlet,
    )


    TransformationFactory("network.expand_arcs").apply_to(m)

    init_md(m, treat.md)
    
    propagate_state(treat.md_to_product)
    propagate_state(treat.md_disposal_to_nacl_translator)

    treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)
    treat.sw_to_nacl_disposal.initialize()

    results = solve(m.fs.treatment.md)
    results = solve(m)

    report_MD(m, treat.md)

    return m


def add_zld_mec(m,permian_cryst_config):

    treat = m.fs.treatment

    print("\n--------- Adding MEC---------\n")
    treat.mec = FlowsheetBlock(dynamic=False)
    build_mec(m, treat.mec, 
            #   prop_package = m.fs.properties_NaCl,
            #   prop_package_vapor = m.fs.properties_vapor
            )
    
    total_feed_H2O_mass = treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ].value
    total_feed_NaCl_mass = treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[
        0, "Liq", "NaCl"
    ].value
    
    set_mec_op_conditions(
                            m, 
                            treat.mec,
                            operating_pressures=permian_cryst_config["operating_pressures"],
                            nacl_yield=permian_cryst_config["nacl_yield"],
                            heat_transfer_coeff = permian_cryst_config["heat_transfer_coefficient"],
                        )
    
    init_mec(treat.mec)

    unfix_mec(treat.mec)

    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())
    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]())

    treat.mec.unit.inlet.temperature[0].fix(treat.sw_to_nacl_disposal.outlet.temperature[0]())
    treat.mec.unit.inlet.pressure[0].fix(treat.sw_to_nacl_disposal.outlet.pressure[0]())

    print("Water",treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())
    print("Nacl",treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]())
    print("Temp",treat.sw_to_nacl_disposal.outlet.temperature[0]())
    print("Pressure",treat.sw_to_nacl_disposal.outlet.pressure[0]())

    mec_rescaling(
        m.fs.treatment.mec,
        flow_mass_phase_water_total=total_feed_H2O_mass,
        flow_mass_phase_salt_total=total_feed_NaCl_mass,
    )

    treat.product_NaCl_mixer = Mixer(
        property_package=m.fs.properties_NaCl,
        num_inlets=3,
        inlet_list=["ro_product","md_product", "cryst_product"],
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.product = Product(property_package=m.fs.properties_NaCl)


    treat.md_translator_to_product_NaCl_mixer = Arc(
        source=treat.sw_to_nacl_product.outlet, destination=treat.product_NaCl_mixer.md_product,
    ) 

    treat.cryst_to_product_NaCl_mixer = Arc(
        source=treat.mec.unit.outlet, destination=treat.product_NaCl_mixer.cryst_product,
    )

    treat.ro_to_product_NaCl_mixer = Arc(
        source=treat.NaCl_cryst_product_translator.outlet, destination=treat.product_NaCl_mixer.ro_product,
    )  

    treat.product_NaCl_mixer_to_product = Arc(
        source=treat.product_NaCl_mixer.outlet, destination=treat.product.inlet,
    ) 

    TransformationFactory("network.expand_arcs").apply_to(m)

    treat.sw_to_nacl_product.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treat.sw_to_nacl_product.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)

    treat.sw_to_nacl_product.initialize()
    
    propagate_state(treat.md_translator_to_product_NaCl_mixer)
    propagate_state(treat.ro_to_product_NaCl_mixer)
    propagate_state(treat.cryst_to_product_NaCl_mixer)

    treat.product_NaCl_mixer.initialize()
    treat.product_NaCl_mixer.outlet.pressure[0].fix()
    
    propagate_state(treat.product_NaCl_mixer_to_product)

    treat.product.properties[0].flow_vol
    treat.product.properties[0].flow_vol_phase
    treat.product.initialize()


def add_zld_cst(m):
    m.fs.energy.cst = FlowsheetBlock()
    build_cst(m.fs.energy.cst)
    set_cst_op_conditions(m.fs.energy.cst, hours_storage=24)
    init_cst(m.fs.energy.cst)


def add_zld_pv(m):

    build_pv(m)
    add_pv_scaling(m, m.fs.energy.pv)
    set_pv_constraints(m, focus="Size")


def add_zld_treatment_costing(m,heat_price,electricity_price):

    # Add treatment costing
    treatment = m.fs.treatment
    treatment.costing = TreatmentCosting()

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )
    m.fs.treatment.costing.heat_cost.fix(heat_price)
    m.fs.treatment.costing.electricity_cost.fix(electricity_price)

    add_ec_costing(m, treatment.EC, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_ro_costing(m, treatment.RO, treatment.costing)

    constraint_scaling_transform(m.fs.treatment.EC.ec.costing.capital_cost_floc_constraint, 1e-6)
    constraint_scaling_transform(m.fs.treatment.EC.ec.costing.capital_cost_reactor_constraint, 1e-3)
    constraint_scaling_transform(m.fs.treatment.EC.ec.costing.capital_cost_power_supply_constraint, 1e-6)
    constraint_scaling_transform(m.fs.treatment.EC.ec.costing.capital_cost_electrodes_constraint, 1e-3)

    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)

    add_mec_costing(m, m.fs.treatment.mec, flowsheet_costing_block=m.fs.treatment.costing)

    constraint_scaling_transform(m.fs.treatment.mec.unit.costing.total_capital_cost_effect_1_constraint,1e-7)
    constraint_scaling_transform(m.fs.treatment.mec.unit.costing.total_capital_cost_effect_2_constraint,1e-7)
    constraint_scaling_transform(m.fs.treatment.mec.unit.costing.total_capital_cost_effect_3_constraint,1e-7)
    constraint_scaling_transform(m.fs.treatment.mec.unit.costing.total_capital_cost_effect_4_constraint,1e-7)
    
    treatment.costing.cost_process()
    treatment.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)

    treatment.costing.initialize()

    print("\n--------- Treatment Costing Initialization Complete ---------\n")
    

def add_zld_heat_energy_costing(m,cost_per_total_aperture_area,cost_per_storage_capital,heat_price,electricity_price):
    # Add energy costing
    energy = m.fs.energy

    add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)
    add_cst_costing_scaling(m,m.fs.energy.cst.unit)


def add_zld_electricity_energy_costing(m,cost_per_watt_installed,heat_price,electricity_price):
    energy = m.fs.energy

    energy.costing.has_electricity_generation = True
    m.fs.energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing,
        costing_method_arguments={
            "cost_method": "simple"
        }
    )


def add_zld_system_energy_costing(m,heat_price,electricity_price):
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)

    print("\n--------- System Costing Initialization Complete ---------\n")


def kbhdp_zld_md_reporting_variables(m):
    # For reporting purposes
    m.fs.treatment.md.unit.capital_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.capital_cost), mutable=True
    )
    m.fs.treatment.md.unit.fixed_operating_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.fixed_operating_cost), mutable=True
    )
    m.fs.treatment.md.unit.module_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.module_cost), mutable=True
    )
    m.fs.treatment.md.unit.other_capital_cost = Param(
        initialize=value(m.fs.treatment.md.unit.costing.other_capital_cost), mutable=True
    )    


def zld_main(
        Qin=4,ro_recovery=0.5, md_water_recovery = 0.7, nacl_recovery_price=0,
        heat_price=0.0166,electricity_price=0.04989, grid_frac_heat=0.5,
        cost_per_total_aperture_area=373,cost_per_storage_capital=62,
        cost_per_watt_installed = 1.6,
        ):
    
    if grid_frac_heat==1:
        treatment_only = True
    else:
        treatment_only = False
    
    m = kbhdp_zld_ro(ro_recovery,Qin)

    print(
        f'RO Recovery: {100 * (value(m.fs.treatment.RO.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]) / value(m.fs.treatment.RO.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"])):<5.2f}%'
    )
    print(f"\nDOF after RO = {degrees_of_freedom(m)}")
    print("\n")

    ro_flow = pyunits.convert(m.fs.treatment.RO.feed.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)
    md_flow = pyunits.convert(m.fs.treatment.RO.disposal.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)
    md_conc = pyunits.convert(m.fs.treatment.RO.disposal.properties[0].conc_mass_phase_comp['Liq', 'NaCl'],to_units=pyunits.g / pyunits.L,)

    
    print("\n---------RO outputs after RO solve ---------\n")
    for idx, stage in m.fs.treatment.RO.stage.items():
        print('RO feed side velocity:',stage.module.feed_side.velocity[0, 0]())
    
    print('RO pump pressure:',m.fs.treatment.pump.control_volume.properties_out[0].pressure())

    print("RO flow:",ro_flow())
    print('MD flow:', md_flow())
    print('MD Conc:', md_conc())

    # Get RO waste stream flow rate and TDS
    # Add MD
    m = add_zld_md(m, Q_md=md_flow(), Cin = md_conc(), md_water_recovery=md_water_recovery)
    print(f"\nDOF after MD = {degrees_of_freedom(m)}")
    print("\n")
    
    # Add MEC
    permian_cryst_config = {
    "operating_pressures": [0.45, 0.25, 0.208, 0.095], # Operating pressure of each effect (bar)
    "nacl_yield": 0.9, # Yield
    "heat_transfer_coefficient": 1300
    }

    add_zld_mec(m,permian_cryst_config)
    results = solve(m.fs.treatment.mec)

    print(f"\nDOF after MEC = {degrees_of_freedom(m)}")
    print("\n")


    if treatment_only == True:
        add_zld_treatment_costing(m,heat_price,electricity_price)
    
        try:
            results = solve(m)
            print(f"\nDOF after Costing = {degrees_of_freedom(m)}")
            print("\n")
        except:
            print_infeasible_constraints(m)

        results = solve(m)

        m.fs.treatment.costing.nacl_recovered.cost.set_value(nacl_recovery_price)

        # feed_density = 1000 * pyunits.kg / pyunits.m**3 
        feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
        )

        product_m3h = pyunits.convert(
            m.fs.treatment.product.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
        )

        print("\nFeed flow in m3/h:",feed_m3h())
        print("Product flow in m3/h:",product_m3h())
        print("\nAggregate flow electricity:",m.fs.treatment.costing.aggregate_flow_electricity())
        print("Aggregate flow heat:",m.fs.treatment.costing.aggregate_flow_heat())
        
        print("LCOW:",m.fs.treatment.costing.LCOW(),pyunits.get_units(m.fs.treatment.costing.LCOW))
        print("SEC (electricity) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_electricity()/product_m3h())
        print("SEC (heat) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_heat()/product_m3h())
        print("System recovery (%):", product_m3h()/feed_m3h()*100)
        print("Capex ($M):",m.fs.treatment.costing.total_capital_cost()/1e6, 
            pyunits.get_units(m.fs.treatment.costing.total_capital_cost))
        print("Opex ($M/yr):",m.fs.treatment.costing.total_operating_cost()/1e6, 
            pyunits.get_units(m.fs.treatment.costing.total_operating_cost))
        # print('Electricity demand (MWh/year):',pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity,to_units=pyunits.MW*pyunits.h/pyunits.year)())
        # print('Heat demand (MWh/year):',pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat,to_units=pyunits.MW*pyunits.h/pyunits.year)())

    else:

        add_zld_treatment_costing(m,heat_price=0,electricity_price=0)
        try:
            results = solve(m)
            print(f"\nDOF after Costing = {degrees_of_freedom(m)}")
            print("\n")
        except:
            print_infeasible_constraints(m)

        results = solve(m)
        
        # Calculate required heat load and CST sizing
        m.fs.energy = Block()
        m.fs.energy.costing = EnergyCosting()

        add_zld_cst(m)
        add_zld_heat_energy_costing(m,cost_per_total_aperture_area,cost_per_storage_capital,heat_price,electricity_price)
        m.fs.energy.cst.unit.heat_load.unfix()
        results = solve(m)

        print(f"\nDOF after Energy = {degrees_of_freedom(m)}")
        print("\n")

        print("\n--------- CST Initialized ---------\n")

        print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
        print('CST Heat:', value(m.fs.energy.cst.unit.heat))
        print("\n")
        
        add_zld_pv(m)

        add_zld_electricity_energy_costing(m,cost_per_watt_installed,heat_price,electricity_price)

        m.fs.energy.cst.unit.heat_load.unfix()
        m.fs.energy.pv.design_size.unfix()

        m.fs.energy.costing.cost_process()
        m.fs.energy.costing.initialize()

        print("\n--------- Energy Costing Initialization Complete ---------\n")

        add_zld_system_energy_costing(m,heat_price,electricity_price)
        
        # CST heat load calculated
        m.fs.energy.cst.unit.heat_load.unfix()
        m.fs.costing.frac_heat_from_grid.fix(grid_frac_heat)

        # m.fs.energy.pv.annual_energy.unfix()
        m.fs.energy.pv.design_size.unfix()
        m.fs.costing.frac_elec_from_grid.fix(0.5)
        
        try:
            results = solve(m)
            print(f"\nDOF after Costing = {degrees_of_freedom(m)}")
            print("\n")
        except:
            print_infeasible_constraints(m)

        print("\n--------- CST Calculated ---------\n")

        print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
        print('CST Heat:', value(m.fs.energy.cst.unit.heat))
        print("\n")
        print(f"\nDOF after CST sizing = {degrees_of_freedom(m)}")
        print("\n")
        
        m.fs.treatment.costing.nacl_recovered.cost.set_value(nacl_recovery_price)

        m.fs.energy.costing.trough_surrogate.cost_per_total_aperture_area.fix(cost_per_total_aperture_area)
        m.fs.energy.costing.trough_surrogate.cost_per_storage_capital.fix(cost_per_storage_capital)
        m.fs.energy.costing.pv_surrogate.cost_per_watt_installed.fix(cost_per_watt_installed)

        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOT)
        
        results = solve(m)
        print(f"\nDOF after PV sizing = {degrees_of_freedom(m)}")

        print("LCOT:", m.fs.costing.LCOT())

        kbhdp_zld_md_reporting_variables(m)

        print("CST Heat load:", m.fs.energy.cst.unit.heat_load())
        print("Heat grid fraction:",m.fs.costing.frac_heat_from_grid())
        print("Heat purchased:",m.fs.costing.aggregate_flow_heat_purchased())
        print("Heat generated:",m.fs.energy.costing.aggregate_flow_heat())
        print("Heat required:",m.fs.treatment.costing.aggregate_flow_heat())
        print("MD heat:", m.fs.treatment.md.unit.overall_thermal_power_requirement() )
        print("MEC heat:",pyunits.convert(m.fs.treatment.mec.unit.effects[1].effect.work_mechanical[0],
                                        to_units=pyunits.kW)())    
        
        print("\nPV design size:",m.fs.energy.pv.design_size())
        print("PV annual energy:",m.fs.energy.pv.annual_energy())
        print("Electricity grid fraction:",m.fs.costing.frac_elec_from_grid())

        print("Electricity purchased:",m.fs.costing.aggregate_flow_electricity_purchased())
        print("Electricity generated (PV-CST):",m.fs.energy.costing.aggregate_flow_electricity())
        print("Electricity required:",m.fs.treatment.costing.aggregate_flow_electricity())
        print('CST electricity required:',m.fs.energy.cst.unit.electricity())

        print("\n---------RO outputs after complete solve ---------\n")
        for idx, stage in m.fs.treatment.RO.stage.items():
            print('RO feed side velocity:',stage.module.feed_side.velocity[0, 0]())
        
        print('RO pump pressure:',m.fs.treatment.pump.control_volume.properties_out[0].pressure(),pyunits.get_units(m.fs.treatment.pump.control_volume.properties_out[0].pressure))
        
        print("Cost per aperature area:", m.fs.energy.costing.trough_surrogate.cost_per_total_aperture_area())

        product_m3h = pyunits.convert(
        m.fs.treatment.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.m**3 / pyunits.h
        )
        
        print('Electricity demand (MWh/year):',pyunits.convert(m.fs.costing.aggregate_flow_electricity,to_units=pyunits.MW*pyunits.h/pyunits.year)())
        print('Heat demand (MWh/year):',pyunits.convert(m.fs.costing.aggregate_flow_heat,to_units=pyunits.MW*pyunits.h/pyunits.year)())
        print("SEC (electricity) in kWh/m3:",m.fs.costing.aggregate_flow_electricity()/product_m3h())
        print("SEC (heat) in kWh/m3:",m.fs.costing.aggregate_flow_heat()/product_m3h())
        print("Capex ($M):",m.fs.costing.total_capital_cost()/1e6, 
            pyunits.get_units(m.fs.costing.total_capital_cost))
        print("Opex ($M/yr):",m.fs.costing.total_operating_cost()/1e6, 
            pyunits.get_units(m.fs.costing.total_operating_cost))
        
    # Adding SEC
    feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
        )

    m.fs.treatment.costing._add_flow_component_breakdowns(
        "heat",
        "SEC_th",
        feed_m3h,
        period=pyunits.hr 
        )

    m.fs.treatment.costing._add_flow_component_breakdowns(
        "electricity",
        "SEC_elec",
        feed_m3h,
        period=pyunits.hr 
        )
    results = solve(m)
    
    return m

def recovery_check(m):
    # RO recovery
    ro_product = m.fs.treatment.RO.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    ro_feed = m.fs.treatment.RO.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]

    print("\nRO")
    print("RO feed:",ro_feed())
    print("RO product:", ro_product())
    print(
        f'RO Recovery: {100 * (value(ro_product) / value(ro_feed)):<5.2f}%'
    )

    # MD recovery
    md_feed = m.fs.treatment.RO.disposal.properties[0].flow_vol
    md_product = m.fs.treatment.md.permeate.properties[0.0].flow_vol_phase["Liq"]

    print("\nMD")
    print("MD feed:",md_feed())
    print("MD product:", md_product())
    print(
        f'MD Recovery: {100 * (value(md_product) / value(md_feed)):<5.2f}%'
    )

    # MD to MEC transition
    print("\nMD to MEC transition")
    print("MD SW to NaCl:",m.fs.treatment.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())
    # print("Normalized flow:",m.fs.treatment.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())

    # MEC recovery
    mec_feed = m.fs.treatment.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
    mec_product = m.fs.treatment.mec.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]

    print("\nMEC")
    print("MEC feed:",mec_feed())
    print("MEC product:", mec_product())

    print(
        f'MEC Recovery: {100 * (value(mec_product) / value(mec_feed)):<5.2f}%'
    )

    # System recovery
    feed_density = 1000 * pyunits.kg / pyunits.m**3 
    system_feed = pyunits.convert(
            m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]/feed_density, to_units=pyunits.m**3 / pyunits.s
        )

    system_product = m.fs.treatment.product.properties[0].flow_vol
    
    print("\nSystem")
    print("System feed:",system_feed())
    print("System product:", system_product())
    print(
        f'System Recovery: {100 * (value(system_product) / value(system_feed)):<5.2f}%'
    )
    

if __name__ == "__main__":

    m = zld_main(
        ro_recovery = 0.8,
        md_water_recovery = 0.78,
        nacl_recovery_price = 0,
        heat_price = 0.00894,
        electricity_price = 0.04989,
        grid_frac_heat = 1,
        cost_per_total_aperture_area = 297,
        cost_per_storage_capital = 62,
        cost_per_watt_installed = 1.6
        )
    
    recovery_check(m)

    total_salt = 0
    for effect_number, eff in m.fs.treatment.mec.unit.effects.items():
        print(effect_number, value(eff.effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"]), 
              pyunits.get_units(eff.effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"]))

        total_salt+=value(eff.effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"])

    print(total_salt)

    feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
    )

    print("\nFeed Basis Treatment SEC (heat) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_heat()/feed_m3h())
    print("Feed Basis Treatment SEC (electricity) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_electricity()/feed_m3h())
    
    print('\nThermal SEC Breakdown')
    print('VAGMD SEC (heat):',m.fs.treatment.md.unit.overall_thermal_power_requirement()/feed_m3h())
    print('MEC SEC (heat):',m.fs.treatment.costing.SEC_th_component["fs.treatment.mec.unit.effects[1].effect"]())
    
    print('\nElectrical SEC Breakdown')
    print('VAGMD SEC (electric):',m.fs.treatment.md.unit.overall_elec_power_requirement()/feed_m3h())
    print('EC SEC (electric):',m.fs.treatment.costing.SEC_elec_component["fs.treatment.EC.ec"]())

    mec_total_sec =  0
    for key in m.fs.treatment.costing.SEC_elec_component.keys():
        if "fs.treatment.mec.unit" in key:
            mec_total_sec += m.fs.treatment.costing.SEC_elec_component[key]()

    print('MEC SEC (electric):',mec_total_sec)
    # print(m.fs.treatment.costing.SEC_elec_component.display())



