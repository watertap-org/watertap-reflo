import pathlib
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
import pandas as pd
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import (
    Product,
    Feed,
    StateJunction,
    Separator,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *
from idaes.core.util.initialization import propagate_state

from watertap.core.solvers import get_solver
# from idaes.core.solvers import get_solver
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import TreatmentCosting, EnergyCosting, REFLOSystemCosting
from watertap_contrib.reflo.analysis.case_studies.permian.components import *
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.MD import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.CST import *

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
# rho = 1125 * pyunits.kg / pyunits.m**3
# rho_water = 997 * pyunits.kg / pyunits.m**3

solver = get_solver()

__all__ = [
    "build_treatment_permian_st1_md",
    "build_energy_permian_st1_md",
    "set_treatment_operating_conditions_st1_md",
    "set_energy_operating_conditions_st1_md",
    "add_treatment_costing_st1_md",
    "add_system_costing_st1_md",
    "set_permian_pretreatment_scaling_st1_md",
    "init_treatment_system_st1_md",
    "init_energy_system_st1_md",
    "run_permian_st1_md",
]


def get_stream_density(Qin=5, tds=130, **kwargs):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.feed_sw = Feed(property_package=m.fs.properties_feed)
    m.fs.feed_sw.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,
            ("conc_mass_phase_comp", ("Liq", "TDS")): tds * pyunits.g / pyunits.liter,
            ("temperature", None): 300,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    m.fs.feed_sw.initialize()
    rho = (
        value(m.fs.feed_sw.properties[0].dens_mass_phase["Liq"])
        * pyunits.kg
        / pyunits.m**3
    )

    return rho


def build_treatment_permian_st1_md(Qin=5, Q_md=0.22478, Cin=118, water_recovery=0.2, rho=None):
    """
    Build Permian pretreatment flowsheet
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.properties.dens_mass_default = rho

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    # Begin building Treatment Block
    m.fs.treatment = treat = Block()

    treat.feed = Feed(property_package=m.fs.properties)
    treat.product = Product(property_package=m.fs.properties_feed)

    m.inlet_flow_rate = pyunits.convert(
        Q_md * pyunits.m**3 / pyunits.s, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )

    # MD water recovery
    m.water_recovery = water_recovery
    m.fs.water_recovery =  Param(
        initialize=water_recovery, mutable=True
    )

    # Add translator blocks
    treat.zo_to_sw_feed = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )

    treat.zo_to_sw_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )

    treat.disposal_ZO_mixer = Mixer(
        property_package=m.fs.properties,
        num_inlets=2,
        inlet_list=["ec_disposal", "cart_filt_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.disposal_SW_mixer = Mixer(
        property_package=m.fs.properties_feed,
        num_inlets=2,
        inlet_list=["zo_mixer", "md_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.EC = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.EC)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.md = FlowsheetBlock(dynamic=False)
    build_md(m, treat.md, m.fs.properties_feed)

    treat.DWI = FlowsheetBlock(dynamic=False)
    build_dwi(m, treat.DWI, m.fs.properties_feed)

    # BUILD PRODUCT STREAM
    # feed > chem_addition > EC > cart_filt > ZO_to_SW_translator > desal unit > product
    treat.feed_to_chem_addition = Arc(
        source=treat.feed.outlet, destination=treat.chem_addition.feed.inlet
    )

    treat.chem_addition_to_ec = Arc(
        source=treat.chem_addition.product.outlet, destination=treat.EC.feed.inlet
    )

    treat.ec_to_cart_filt = Arc(
        source=treat.EC.product.outlet, destination=treat.cart_filt.feed.inlet
    )

    # from ZO to SW - feed
    treat.cart_filt_to_translator = Arc(
        source=treat.cart_filt.product.outlet, destination=treat.zo_to_sw_feed.inlet
    )

    treat.cart_filt_translated_to_md = Arc(
        source=treat.zo_to_sw_feed.outlet, destination=treat.md.feed.inlet
    )

    treat.md_to_product = Arc(
        source=treat.md.permeate.outlet, destination=treat.product.inlet
    )

    # BUILD DISPOSAL STREAM
    #        EC > ZO_mixer > ZO_to_SW_translator > disposal_mixer > disposal_mixer > DWI
    # cart_filt > ZO_mixer
    #                                   MD unit  > disposal_mixer

    treat.ec_to_disposal_mix = Arc(
        source=treat.EC.disposal.outlet, destination=treat.disposal_ZO_mixer.ec_disposal
    )

    treat.cart_filt_to_disposal_mix = Arc(
        source=treat.cart_filt.disposal.outlet,
        destination=treat.disposal_ZO_mixer.cart_filt_disposal,
    )

    treat.disposal_ZO_mix_to_translator = Arc(
        source=treat.disposal_ZO_mixer.outlet, destination=treat.zo_to_sw_disposal.inlet
    )

    treat.disposal_ZO_mix_translated_to_disposal_SW_mixer = Arc(
        source=treat.zo_to_sw_disposal.outlet,
        destination=treat.disposal_SW_mixer.zo_mixer,
    )

    treat.md_disposal_to_translator = Arc(
        source=treat.md.concentrate.outlet,
        destination=treat.disposal_SW_mixer.md_disposal,
    )

    treat.disposal_SW_mixer_to_dwi = Arc(
        source=treat.disposal_SW_mixer.outlet, destination=treat.DWI.feed.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Add treatment costing 
    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)

    return m

def build_energy_permian_st1_md(m):
    # Build energy block
    m.fs.energy = energy = Block()
    m.fs.energy.cst = FlowsheetBlock()
    # Energy block selection based on treatment train energy range
    build_cst(m.fs.energy.cst)

    # Add energy costing
    m.fs.energy.costing = EnergyCosting()

    return m

def set_energy_operating_conditions_st1_md(m,heat_load=10,hours_storage=24):
    set_cst_op_conditions(m.fs.energy.cst,heat_load, hours_storage)


def set_treatment_operating_conditions_st1_md(m, rho, Qin=5, tds=130, **kwargs):

    global flow_mass_water, flow_mass_tds, flow_in

    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )

    m.fs.treatment.feed.properties[0].flow_mass_comp["H2O"].fix(
        flow_mass_water - flow_mass_tds
    )
    m.fs.treatment.feed.properties[0].flow_mass_comp["tds"].fix(flow_mass_tds)
    m.fs.treatment.feed.properties[0].conc_mass_comp[...]

    set_chem_addition_op_conditions(m, m.fs.treatment.chem_addition, **kwargs)
    set_ec_operating_conditions(m, m.fs.treatment.EC, **kwargs)
    set_cart_filt_op_conditions(m, m.fs.treatment.cart_filt)


def set_permian_pretreatment_scaling_st1_md(
    m, calculate_m_scaling_factors=False, **kwargs
):

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        # 1 / value(flow_mass_water),
        1e-2,
        index=("H2O"),
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        # 1 / value(flow_mass_tds),
        0.1,
        index=("tds"),
    )

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        0.1,
        index=("Liq", "TDS"),
    )

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1e-2,
        index=("Liq", "H2O"),
    )

    set_chem_addition_scaling(
        m, m.fs.treatment.chem_addition, calc_blk_scaling_factors=True
    )

    set_cart_filt_scaling(m, m.fs.treatment.cart_filt, calc_blk_scaling_factors=True)

    set_ec_scaling(m, m.fs.treatment.EC, calc_blk_scaling_factors=True)

    set_scaling_factor(
        m.fs.treatment.product.properties[0].flow_mass_phase_comp["Liq", "H2O"], 1e-2
    )
    set_scaling_factor(
        m.fs.treatment.product.properties[0].flow_mass_phase_comp["Liq", "TDS"], 1e5
    )

    # ZO to SW feed translator
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_feed.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-2,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_feed.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     0.1,
    # )

    # # ZO to SW disposal translator
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_disposal.properties_in[0].flow_mass_comp["H2O"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_disposal.properties_in[0].flow_mass_comp["tds"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_disposal.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_to_sw_disposal.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     1,
    # )

    # # ZO DISPOSAL MIXER
    # # CF inlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.cart_filt_disposal_state[0].flow_mass_comp[
    #         "H2O"
    #     ],
    #     100,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.cart_filt_disposal_state[0].flow_mass_comp[
    #         "tds"
    #     ],
    #     1e8,
    # )

    # # EC inlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.ec_disposal_state[0].flow_mass_comp["H2O"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.ec_disposal_state[0].flow_mass_comp["tds"],
    #     1,
    # )

    # # mixed state
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.mixed_state[0].flow_mass_comp["H2O"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_ZO_mixer.mixed_state[0].flow_mass_comp["tds"],
    #     1,
    # )
    # # SW DISPOSAL MIXER
    # # ZO mixer inlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.zo_mixer_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     100,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.zo_mixer_state[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     10,
    # )

    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.md_disposal_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-3,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.md_disposal_state[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     1e-2,
    # )

    # # mixed state outlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.mixed_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_SW_mixer.mixed_state[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     1e-2,
    # )

    # DWI
    set_scaling_factor(
        m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-5,
    )
    set_scaling_factor(
        m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        1e-1,
    )

    if calculate_m_scaling_factors:
        print("calculate_m_scaling_factors\n\n\n")
        calculate_scaling_factors(m)


def init_treatment_system_st1_md(m, **kwargs):

    treat = m.fs.treatment

    # Touch flow_vol variable
    treat.feed.properties[0].flow_vol

    treat.feed.initialize()
    propagate_state(treat.feed_to_chem_addition)

    init_chem_addition(m, treat.chem_addition)

    propagate_state(treat.chem_addition_to_ec)

    init_ec(m, treat.EC)
    propagate_state(treat.ec_to_cart_filt)
    propagate_state(treat.ec_to_disposal_mix)

    init_cart_filt(m, treat.cart_filt)
    propagate_state(treat.cart_filt_to_translator)
    propagate_state(treat.cart_filt_to_disposal_mix)

    treat.disposal_ZO_mixer.initialize()
    propagate_state(treat.disposal_ZO_mix_to_translator)

    # NOTE: If the selected temperature is similar to md_disposal temperature, the mixer has issues in the solve
    treat.zo_to_sw_disposal.outlet.temperature[0].fix(300)
    treat.zo_to_sw_disposal.outlet.pressure[0].fix(101325)
    treat.zo_to_sw_disposal.initialize()

    treat.zo_to_sw_feed.properties_out[0].temperature.fix(300)
    treat.zo_to_sw_feed.properties_out[0].pressure.fix(101325)
    treat.zo_to_sw_feed.initialize()

    propagate_state(treat.cart_filt_translated_to_md)

    init_md(m, treat.md)
    propagate_state(treat.md_to_product)
    propagate_state(treat.md_disposal_to_translator)

    propagate_state(treat.disposal_ZO_mix_translated_to_disposal_SW_mixer)
    # NOTE: variable that affects DOF in unclear way

    treat.disposal_SW_mixer.initialize()
    treat.disposal_SW_mixer.mixed_state[0].temperature.fix(300)
    treat.disposal_SW_mixer.mixed_state[0].pressure.fix()

    propagate_state(treat.disposal_SW_mixer_to_dwi)
    # NOTE: variables that affect DOF in unclear way
    # treat.DWI.feed.properties[0].temperature.fix()
    # treat.DWI.feed.properties[0].pressure.fix()
    treat.DWI.unit.properties[0].flow_vol_phase
    treat.DWI.unit.properties[0].conc_mass_phase_comp
    init_dwi(m, treat.DWI)

    treat.product.properties[0].flow_vol_phase
    treat.product.properties[0].conc_mass_phase_comp
    treat.product.initialize()


def init_energy_system_st1_md(m):
    init_cst(m.fs.energy.cst)


def add_treatment_costing_st1_md(m, heat_price, electricity_price):
   # Treatment costing only. To be used when grid_frac = 1
    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.EC, flowsheet_costing_block=m.fs.treatment.costing)

    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )

    add_dwi_costing(
        m, m.fs.treatment.DWI, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)
    m.fs.treatment.costing.heat_cost.fix(heat_price)
    m.fs.treatment.costing.electricity_cost.fix(electricity_price)

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.product.properties[0].flow_vol
    )

    # Add energy costing
    add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)

    m.fs.energy.costing.heat_cost.set_value(0)
    m.fs.energy.costing.cost_process()
    # m.fs.energy.costing.maintenance_labor_chemical_factor.fix(0)
    # # m.fs.energy.costing.add_LCOH()

    # Add system costing
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")
    m.fs.treatment.costing.initialize()
    m.fs.energy.costing.initialize()
    m.fs.costing.initialize()

    print("\n--------- INITIALIZING SYSTEM COSTING COMPLETE---------\n")

    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)



def add_system_costing_st1_md(m, heat_price, electricity_price):

    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.EC, flowsheet_costing_block=m.fs.treatment.costing)

    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )

    add_dwi_costing(
        m, m.fs.treatment.DWI, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.product.properties[0].flow_vol
    )

    # Add energy costing
    add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)

    m.fs.energy.costing.heat_cost.set_value(0)
    m.fs.energy.costing.electricity_cost.fix(electricity_price)
    m.fs.energy.costing.cost_process()
    # m.fs.energy.costing.maintenance_labor_chemical_factor.fix(0)

    # Add system costing
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")
    m.fs.treatment.costing.initialize()
    m.fs.energy.costing.initialize()
    m.fs.costing.initialize()

    print("\n--------- INITIALIZING SYSTEM COSTING COMPLETE---------\n")

    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)


def set_md_cost_scaling(m):
    iscale.calculate_scaling_factors(m.fs.treatment.md.unit.mp)
    if (
        iscale.get_scaling_factor(
            m.fs.treatment.md.unit.overall_thermal_power_requirement
        )
        is None
    ):
        iscale.set_scaling_factor(
            m.fs.treatment.md.unit.overall_thermal_power_requirement, 1e-6
        )

    if (
        iscale.get_scaling_factor(m.fs.treatment.md.unit.overall_elec_power_requirement)
        is None
    ):
        iscale.set_scaling_factor(
            m.fs.treatment.md.unit.overall_elec_power_requirement, 1e-4
        )

    if (
        iscale.get_scaling_factor(
            m.fs.treatment.md.unit.mp.get_active_process_blocks()[
                -1
            ].fs.vagmd.system_capacity
        )
        is None
    ):
        iscale.set_scaling_factor(
            m.fs.treatment.md.unit.mp.get_active_process_blocks()[
                -1
            ].fs.vagmd.system_capacity,
            1e-4,
        )

def solve(
    m, solver=None, tee=False, raise_on_failure=True, symbolic_solver_labels=True
):
    # ---solving---
    if solver is None:
        solver = get_solver()

    solver.options["max_iter"] = 1000
    solver.options["halt_on_ampl_error"] = "yes"

    print(f"\n--------- SOLVING {m.name} ---------\n")

    results = solver.solve(m, tee=tee, symbolic_solver_labels=True)

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


def run_permian_st1_md(Qin=5, tds=130, grid_frac_heat = 0.5, water_recovery = 0.3, 
                       heat_price=0.00894, electricity_price=0.0575, dwi_lcow = 8.4, **kwargs):
    
    """
    Permian pretreatment flowsheet to get input for MD
    """

    rho = get_stream_density(Qin, tds)
    m_pretreatment = build_and_run_permian_pretreatment(Qin,tds)

    _ = solve(m_pretreatment)

    print(
        f"Pretreatment Product Flow: {pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)():.4f} m3/s"
    )
    print(
        f"Pretreatment Product Flow: {pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].conc_mass_phase_comp['Liq', 'TDS'],to_units=pyunits.g / pyunits.L,)():.4f} g/L"
    )

    pretreatment_system_recovery = (
        m_pretreatment.fs.treatment.product.properties[0].flow_vol() / m_pretreatment.fs.treatment.feed.properties[0].flow_vol()
    )

    print("Pretreatment Recovery:", pretreatment_system_recovery)

    md_flow = pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)
    md_conc = pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].conc_mass_phase_comp['Liq', 'TDS'],to_units=pyunits.g / pyunits.L,)

    print("md_flow:",md_flow())
    print("md_conc:",md_conc())

    '''
    Complete flowsheet
    '''

    m = build_treatment_permian_st1_md(Q_md=md_flow(), Cin=md_conc(), water_recovery = water_recovery, rho=rho)
    treat = m.fs.treatment

    set_treatment_operating_conditions_st1_md(m, rho, Qin, tds)
    set_permian_pretreatment_scaling_st1_md(
        m, calculate_m_scaling_factors=True
    )  # Doesn't solve without this even before costing

    # Initialize system
    init_treatment_system_st1_md(m)
    print(f"DOF = {degrees_of_freedom(m)}")

    results = solve(m)

    print("\n--------- Before costing solve Completed ---------\n")
    report_MD(m, treat.md)
    print("\n")

    set_md_cost_scaling(m)

    # Add energy block
    m = build_energy_permian_st1_md(m)
    set_energy_operating_conditions_st1_md(m)
    init_energy_system_st1_md(m)

    results = solve(m)

    print("\n--------- Before costing solve Completed ---------\n")

    print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
    print('CST Heat:', value(m.fs.energy.cst.unit.heat))
    print("\n")

    # Add costing  

    if grid_frac_heat==1:
        add_treatment_costing_st1_md(m, heat_price, electricity_price)

        if dwi_lcow!= None:
            m.fs.treatment.costing.deep_well_injection.dwi_lcow.set_value(dwi_lcow)

        print(f"DOF = {degrees_of_freedom(m)}")

        results = solve(m)

        # Update fs.costing block results for only treatment costs
        # Update total heat/electric operating to be treatment aggregate_flow_costs
        m.fs.costing.total_heat_operating_cost = m.fs.treatment.costing.aggregate_flow_costs['heat']
        m.fs.costing.total_electric_operating_cost = m.fs.treatment.costing.aggregate_flow_costs['electricity']
        
        # Update the total_capital, total_operating
        m.fs.energy.cst.unit.costing.capital_cost.fix(1e-20)
        m.fs.energy.cst.unit.costing.fixed_operating_cost.fix(1e-20)
        m.fs.energy.cst.unit.costing.variable_operating_cost.fix(1e-20)

        # Update LCOT to be LCOW
        m.fs.costing.LCOT.fix(m.fs.treatment.costing.LCOW())
        # Update grid fraction reported
        m.fs.costing.frac_heat_from_grid.fix(1)

    else:
        
        add_system_costing_st1_md(m, heat_price, electricity_price)
        add_cst_costing_scaling(m,m.fs.energy.cst.unit)

        results = solve(m)

        m.fs.energy.cst.unit.heat_load.unfix()
        m.fs.costing.frac_heat_from_grid.fix(grid_frac_heat)

        if dwi_lcow!= None:
            m.fs.treatment.costing.deep_well_injection.dwi_lcow.set_value(dwi_lcow)

        print(f"DOF = {degrees_of_freedom(m)}")

        results = solve(m)

        print("\n--------- Costing solve with fixed grid fraction ---------\n")

        print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
        print('CST Heat:', value(m.fs.energy.cst.unit.heat))
        print('Grid fraction:',value(m.fs.costing.frac_heat_from_grid))

    permian_md_reporting_variables(m)

    return m


def permian_md_reporting_variables(m):
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

def report_costing(blk):

    print(f"\n\n-------------------- System Costing Report --------------------\n")
    print("\n")

    print(f'{"LCOT":<30s}{value(blk.LCOT):<20,.2f}{pyunits.get_units(blk.LCOT)}')

    print(
        f'{"Capital Cost":<30s}{value(blk.total_capital_cost):<20,.2f}{pyunits.get_units(blk.total_capital_cost)}'
    )

    print(
        f'{"Total Operating Cost":<30s}{value(blk.total_operating_cost):<20,.2f}{pyunits.get_units(blk.total_operating_cost)}'
    )

    print(
        f'{"Agg Fixed Operating Cost":<30s}{value(blk.aggregate_fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.aggregate_fixed_operating_cost)}'
    )

    print(
        f'{"Agg Variable Operating Cost":<30s}{value(blk.aggregate_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.aggregate_variable_operating_cost)}'
    )

    print(
        f'{"Heat flow":<30s}{value(blk.aggregate_flow_heat):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat)}'
    )

    # print(
    #     f'{"Total heat cost":<30s}{value(blk.total_heat_operating_cost):<20,.2f}{pyunits.get_units(blk.total_heat_operating_cost)}'
    # )

    print(
        f'{"Heat purchased":<30s}{value(blk.aggregate_flow_heat_purchased):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat_purchased)}'
    )

    # print(
    #     f'{"Heat sold":<30s}{value(blk.aggregate_flow_heat_sold):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat_sold)}'
    # )

    print(
        f'{"Elec Flow":<30s}{value(blk.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity)}'
    )

    # print(
    #     f'{"Total elec cost":<30s}{value(blk.total_electric_operating_cost):<20,.2f}{pyunits.get_units(blk.total_electric_operating_cost)}'
    # )

    print(
        f'{"Elec purchased":<30s}{value(blk.aggregate_flow_electricity_purchased):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity_purchased)}'
    )

    # print(
    #     f'{"Elec sold":<30s}{value(blk.aggregate_flow_electricity_sold):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity_sold)}'
    # )

def print_results(m):
    treat = m.fs.treatment
    report_MD(m, treat.md)

    print(f"\nDOF = {degrees_of_freedom(m)}")

    report_cst(m, m.fs.energy.cst.unit)

    system_recovery = (
        treat.product.properties[0].flow_vol() / treat.feed.properties[0].flow_vol()
    )

    print(f"\n\n-------------------- System Cost Report --------------------\n")
    print("\n")

    print(
        f'{"Treatment LCOW":<30s}{value(m.fs.treatment.costing.LCOW):<10.2f}{pyunits.get_units(m.fs.treatment.costing.LCOW)}'
    )

    # print("\n")
    # print(
    #     f'{"Energy LCOH":<30s}{value(m.fs.energy.costing.LCOH):<10.2f}{pyunits.get_units(m.fs.energy.costing.LCOH)}'
    # )

    print("\n")
    print(
        f'{"System LCOT":<30s}{value(m.fs.costing.LCOT) :<10.2f}{pyunits.get_units(m.fs.costing.LCOT)}'
    )

    print("\n")
    print(
        f'{"Percent from the grid":<30s}{value(m.fs.costing.frac_heat_from_grid):<10.2f}{pyunits.get_units(m.fs.costing.frac_heat_from_grid)}'
    )


    print(f"\n\n-------------------- Pretreatment Report --------------------\n")

    print("\n")
    print(
        f'{"System Recovery":<30s}{system_recovery:.2f}'
    )

    print(
        f'{"Inlet flow_vol":<30s} {treat.feed.properties[0].flow_vol():<10.2f} {pyunits.get_units(treat.feed.properties[0].flow_vol)}'
    )
    print(
        f'{"Inlet TDS conc":<30s} {treat.feed.properties[0].conc_mass_comp["tds"]():<10.2f} {pyunits.get_units(treat.feed.properties[0].conc_mass_comp["tds"])}'
    )

    print(
        f'{"EC feed TDS conc":<30s} {treat.EC.feed.properties[0].conc_mass_comp["tds"]():.<10.2f} {pyunits.get_units(treat.EC.feed.properties[0].conc_mass_comp["tds"])}'
    )

    print(
        f'{"EC product TDS conc":<30s} {treat.EC.product.properties[0].conc_mass_comp["tds"]():<10.2f} { pyunits.get_units(treat.EC.product.properties[0].conc_mass_comp["tds"])}'
    )

    print(
        f'{"EC disposal TDS conc":<30s} {treat.EC.disposal.properties[0].conc_mass_comp["tds"]():<10.2f} {pyunits.get_units(treat.EC.disposal.properties[0].conc_mass_comp["tds"])}'
    )

    print(
        f'{"CF feed TDS conc":<30s} {treat.cart_filt.product.properties[0].conc_mass_comp["tds"]():<10.2f} {pyunits.get_units(treat.cart_filt.product.properties[0].conc_mass_comp["tds"])}'
    )

    print(
        f'{"Product TDS conc":<30s} {treat.product.properties[0].conc_mass_phase_comp["Liq", "TDS"]():.<10.2f} {pyunits.get_units(treat.product.properties[0].conc_mass_phase_comp["Liq", "TDS"]())}'
    )

    print(
        f'{"Product flow_vol":<30s} {treat.product.properties[0].flow_vol_phase["Liq"]():<10.2f} {pyunits.get_units(treat.product.properties[0].flow_vol_phase["Liq"])}'
    )

    print(
        f'{"DWI flow_vol":<30s} {treat.DWI.unit.properties[0].flow_vol_phase["Liq"]():<10.2f} {pyunits.get_units(treat.DWI.unit.properties[0].flow_vol_phase["Liq"])}'
    )

    print(
        f'{"DWI TDS conc":<30s} {treat.DWI.unit.properties[0].conc_mass_phase_comp["Liq", "TDS"]():<10.2f} {pyunits.get_units(treat.DWI.unit.properties[0].conc_mass_phase_comp["Liq", "TDS"])}'
    )
    print(
        f'{"DWI pressure":<30s} {treat.DWI.feed.properties[0].pressure()} Pa'
        )

    print(
        f'{"Translator pressure":<30s} {treat.disposal_SW_mixer.zo_mixer_state[0].pressure()} Pa'
    )

    print(
        f'{"Aggregated Heat Cost":<30s}{value(m.fs.treatment.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(m.fs.treatment.costing.aggregate_flow_costs["heat"])}'
    )

    report_costing(m.fs.costing)

def sweep_feed_flow_salinity():
    heat_price = 0.00894
    electricity_price = 0.04346

    sweep_dict = {
        'Qin':[1,5,9],
        'tds': [100,130,200],
        'recovery': [0.6,0.48,0.23]
    }

    for flow in sweep_dict["Qin"]:
        for i in range(0,len(sweep_dict['tds'])):
            m = run_permian_st1_md(
                            Qin=flow, 
                            tds=sweep_dict['tds'][i], 
                            water_recovery = sweep_dict['recovery'][i],
                            grid_frac_heat = 0.5,
                            heat_price=heat_price, 
                            electricity_price=electricity_price, 
                            dwi_lcow  = 8.4
                            )


def main():
    heat_price = 0.00894
    electricity_price = 0.04346  # Updated 0.0575 in USD 2018 to USD 2023

    m = run_permian_st1_md(
        Qin=5, 
        tds=200, 
        water_recovery = 0.23,
        grid_frac_heat = 1,
        heat_price=heat_price, 
        electricity_price=electricity_price, 
        dwi_lcow  = 8.4
        )
    
    print_results(m)

    print('Treatment heat:', m.fs.treatment.costing.aggregate_flow_heat())
    print('Energy heat:', m.fs.energy.costing.aggregate_flow_heat())

if __name__ == "__main__":
    main()