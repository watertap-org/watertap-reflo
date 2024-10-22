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
    relax_constaints(m, m.fs.treatment.RO)
    set_operating_conditions(m)
    init_system(m)
    add_costing(m)
    optimize(m, ro_mem_area=None, water_recovery=0.5)
    solve(m)
    display_system_stream_table(m)
    # display_costing_breakdown(m)
    report_EC(m.fs.treatment.EC)
    report_UF(m, m.fs.treatment.UF)
    report_RO(m, m.fs.treatment.RO)
    report_pump(m, m.fs.treatment.pump)
    report_PV(m)
    # m.fs.treatment.costing.display()
    # m.fs.energy.costing.display()
    # m.fs.costing.display()
    # display_costing_breakdown(m)
    # print(m.fs.energy.pv.display())

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

    energy.pv = PVSurrogate(
        surrogate_model_file="/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/solar_models/surrogate/pv/pv_surrogate.json",
        dataset_filename="/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/solar_models/surrogate/pv/data/dataset.pkl",
        input_variables={
            "labels": ["design_size"],
            "bounds": {"design_size": [1, 200000]},
            "units": {"design_size": "kW"},
        },
        output_variables={
            "labels": ["annual_energy", "land_req"],
            "units": {"annual_energy": "kWh", "land_req": "acre"},
        },
        scale_training_data=False,
    )

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

def train_pv_surrogate(m):
    energy = m.fs.energy
    
    energy.pv.create_rbf_surrogate()

    assert False

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
        destination=treatment.disposal.inlet,
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

    m.fs.eq_water_recovery = Constraint(
        expr=treatment.feed.properties[0].flow_vol * m.fs.water_recovery
        == treatment.product.properties[0].flow_vol
    )

    # m.fs.feed.properties[0].conc_mass_phase_comp
    # m.fs.product.properties[0].conc_mass_phase_comp
    # m.fs.disposal.properties[0].conc_mass_phase_comp


def add_treatment_costing(m):
    treatment = m.fs.treatment
    treatment.costing = TreatmentCosting()

    treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing,
    )
    add_ec_costing(m, treatment.EC, treatment.costing)
    add_UF_costing(m, treatment.UF, treatment.costing)
    add_ro_costing(m, treatment.RO, treatment.costing)

    treatment.costing.ultra_filtration.capital_a_parameter.fix(500000)

    treatment.costing.cost_process()
    treatment.costing.initialize()


def add_energy_costing(m):
    energy = m.fs.energy
    energy.costing = EnergyCosting()

    energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
    )

    energy.pv_design_constraint = Constraint(
        expr=m.fs.energy.pv.design_size
        == m.fs.treatment.costing.aggregate_flow_electricity
    )

    energy.costing.cost_process()
    energy.costing.initialize()


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


def relax_constaints(m, blk):
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

    if Qin is None:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(Qin)

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


def set_operating_conditions(m, RO_pressure=20e5, supply_pressure=1.1e5):
    treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=1000, supply_pressure=1e5)
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
    initialize_energy(m)
    init_treatment(m)


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 1000

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
        print(f"------- Unfixed RO Pump Pressure -------")
        treatment.pump.control_volume.properties_out[0].pressure.unfix()

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for idx, stage in treatment.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        print(f"\n------- Unfixed RO Membrane Area -------\n")
        for idx, stage in treatment.RO.stage.items():
            stage.module.area.unfix()

    # energy.pv_design_constraint = Constraint(
    #     expr=energy.pv.design_size
    #     == treatment.costing.aggregate_flow_electricity
    # )


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
    print(f'{"Pump Pressure":<30s}{value(pyunits.convert(pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}')
    print(f'{"Pump Work":<30s}{value(pyunits.convert(pump.control_volume.work[0], to_units=pyunits.kW)):<10.3f}{"kW"}')


def report_PV(m):
    elec = "electricity"
    print(f"\n\n-------------------- PHOTOVOLTAIC SYSTEM --------------------\n\n")
    print(f'{"System Agg. Flow Electricity":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}')
    print(f'{"PV Agg. Flow Elec.":<30s}{value(m.fs.energy.pv.design_size):<10.1f}{pyunits.get_units(m.fs.energy.pv.design_size)}')
    print(f'{"Treatment Agg. Flow Elec.":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}')
    print(f'{"Land Requirement":<30s}{value(m.fs.energy.pv.land_req):<10.1f}{pyunits.get_units(m.fs.energy.pv.land_req)}')
    print(f'{"PV Annual Energy":<30s}{value(m.fs.energy.pv.annual_energy):<10,.1f}{pyunits.get_units(m.fs.energy.pv.annual_energy)}')
    print(f'{"Treatment Annual Energy":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}')
    print('\n')
    print(f'{"PV Annual Generation":<25s}{f"{pyunits.convert(-1*m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}')
    print(f'{"Treatment Annual Demand":<25s}{f"{pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}')

    print(f'{"Treatment Elec Cost":<25s}{f"${value(m.fs.treatment.costing.aggregate_flow_costs[elec]):<25,.0f}"}{"$/yr":<10s}')
    print(f'{"Energy Elec Cost":<25s}{f"${value(m.fs.energy.costing.aggregate_flow_costs[elec]):<25,.0f}"}{"$/yr":<10s}')


def print_PV_costing_breakdown(pv):
    print(
        f'{"PV Capital Cost":<35s}{f"${value(pv.costing.capital_cost):<25,.0f}"}'
    )
    print(
        f'{"PV Operating Cost":<35s}{f"${value(pv.costing.fixed_operating_cost):<25,.0f}"}'
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
    print_PV_costing_breakdown(m.fs.energy.pv)


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()
