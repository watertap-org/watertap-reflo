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
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.property_models.unit_specific.cryst_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from idaes.core import MaterialBalanceType
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *

from watertap_contrib.reflo.costing import TreatmentCosting,  EnergyCosting, REFLOSystemCosting
from watertap_contrib.reflo.analysis.case_studies.permian.components import *
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.MD import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.translator_zo_to_nacl import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.translator_sw_to_nacl import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.normalizer_cryst import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.denormalizer_cryst import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.multi_effect_crystallizer import *
from watertap_contrib.reflo.analysis.case_studies.permian.components.CST import *

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"


solver = get_solver()

__all__ = [
    "build_permian_st2_md",
    "set_operating_conditions_st2_md",
    "add_costing_st2_md",
    "set_permian_pretreatment_scaling_st2_md",
    "init_system_st2_md",
    "run_permian_st2_md",
]


def get_stream_density(Qin=5, tds=130, **kwargs):
    # global rho

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

    # rho = m.fs.feed_sw.properties[0].dens_mass_phase["Liq"]

    return rho


def build_permian_st2_md(Qin=5, Q_md=0.22478, Cin=118, water_recovery=0.2, rho=None):
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
    m.fs.properties_NaCl = NaClParameterBlock()

    # Begin building Treatment Block
    m.fs.treatment = treat = Block()

    treat.feed = Feed(property_package=m.fs.properties)
    # treat.product = Product(property_package=m.fs.properties_NaCl)

    m.inlet_flow_rate = pyunits.convert(
        Q_md * pyunits.m**3 / pyunits.s, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )

    m.water_recovery = water_recovery
    m.fs.water_recovery =  Param(
        initialize=water_recovery, mutable=True
    )

    # Add translator blocks
    treat.zo_to_sw_feed = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )

    treat.zo_mix_to_nacl_disposal = Translator_ZO_to_NaCl(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_NaCl,
    )

    treat.disposal_ZO_mixer = Mixer(
        property_package=m.fs.properties,
        num_inlets=2,
        inlet_list=["ec_disposal", "cart_filt_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.sw_to_nacl_product = Translator_SW_to_NaCl(
        inlet_property_package = m.fs.properties_feed,
        outlet_property_package= m.fs.properties_NaCl,
    )

    treat.sw_to_nacl_disposal = Translator_SW_to_NaCl(
        inlet_property_package = m.fs.properties_feed,
        outlet_property_package= m.fs.properties_NaCl,
    )

    treat.disposal_NaCl_mixer = Mixer(
        property_package=m.fs.properties_NaCl,
        num_inlets=2,
        inlet_list=["zo_mixer", "md_disposal"],
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.norm_feed = Normalizer_Cryst(
        inlet_property_package = m.fs.properties_NaCl,
        outlet_property_package= m.fs.properties_NaCl,
    )

    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.EC = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.EC)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.md = FlowsheetBlock(dynamic=False)
    build_md(m, treat.md, m.fs.properties_feed)


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
        source=treat.md.permeate.outlet, destination=treat.sw_to_nacl_product.inlet
    )

    # BUILD DISPOSAL STREAM
    #        EC > ZO_mixer > ZO_to_NaCl_translator > disposal_mixer > disposal_mixer > DWI
    # cart_filt > ZO_mixer
    #                                   MD unit  > disposal_mixer

    treat.ec_to_disposal_mix = Arc(
        source=treat.EC.disposal.outlet, destination=treat.disposal_ZO_mixer.ec_disposal
    )

    treat.cart_filt_to_disposal_mix = Arc(
        source=treat.cart_filt.disposal.outlet,
        destination=treat.disposal_ZO_mixer.cart_filt_disposal,
    )

    treat.disposal_ZO_mix_to_nacl_translator = Arc(
        source=treat.disposal_ZO_mixer.outlet, 
        destination=treat.zo_mix_to_nacl_disposal.inlet
    )

    treat.md_disposal_to_nacl_translator = Arc(
        source=treat.md.concentrate.outlet,
        destination=treat.sw_to_nacl_disposal.inlet,
    )

    treat.disposal_ZO_mix_translated_to_disposal_NaCl_mixer = Arc(
        source=treat.zo_mix_to_nacl_disposal.outlet,
        destination=treat.disposal_NaCl_mixer.zo_mixer,
    )

    treat.md_disposal_translated_to_disposal_NaCl_mixer = Arc(
        source=treat.sw_to_nacl_disposal.outlet,
        destination=treat.disposal_NaCl_mixer.md_disposal,
    )

    treat.mixer_to_normalized_feed = Arc(
        source=treat.disposal_NaCl_mixer.outlet, destination=treat.norm_feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Add treatment costing 
    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)

    return m

def build_energy_permian_st2_md(m):
    # Build energy block
    m.fs.energy = energy = Block()
    m.fs.energy.cst = FlowsheetBlock()
    # Energy block selection based on treatment train energy range
    build_cst(m.fs.energy.cst)

    # Add energy costing
    m.fs.energy.costing = EnergyCosting()

    return m


def set_operating_conditions_st2_md(m, rho, Qin=5, tds=130, **kwargs):

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


def set_energy_operating_conditions_st2_md(m,heat_load=10,hours_storage=24):
    set_cst_op_conditions(m.fs.energy.cst, heat_load, hours_storage)


def add_treatment_costing_st2_md(m,heat_price=0.018, electricity_price=0.0626):

    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)

    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.EC, flowsheet_costing_block=m.fs.treatment.costing)
    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )
    
    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)

    add_mec_costing(m, m.fs.treatment.mec, flowsheet_costing_block=m.fs.treatment.costing)

    m.fs.treatment.costing.heat_cost.fix(heat_price)
    m.fs.treatment.costing.electricity_cost.fix(electricity_price)

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)

    m.fs.treatment.costing.initialize()

    # # Add energy costing

    # add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)

    # m.fs.energy.costing.heat_cost.set_value(0)
    # # m.fs.energy.costing.electricity_cost.fix(electricity_price)
    # m.fs.energy.costing.cost_process()

    # # Add system costing
    # m.fs.costing = REFLOSystemCosting()
    # m.fs.costing.heat_cost_buy.fix(heat_price)
    # m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    # m.fs.costing.cost_process()

    # print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")
    
    # m.fs.energy.costing.initialize()
    # m.fs.treatment.costing.initialize()
    # m.fs.costing.initialize()

    # print("\n--------- INITIALIZING SYSTEM COSTING COMPLETE---------\n")

    # m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)


def add_system_costing_st2_md(m,heat_price=0.018, electricity_price=0.0626):

    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.EC, flowsheet_costing_block=m.fs.treatment.costing)
    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )
    
    m.fs.treatment.md.unit.add_costing_module(m.fs.treatment.costing)

    add_mec_costing(m, m.fs.treatment.mec, flowsheet_costing_block=m.fs.treatment.costing)

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)

    # Add energy costing

    add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)

    m.fs.energy.costing.heat_cost.set_value(0)
    m.fs.energy.costing.cost_process()

    # Add system costing
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")
    
    m.fs.energy.costing.initialize()
    m.fs.treatment.costing.initialize()
    m.fs.costing.initialize()

    print("\n--------- INITIALIZING SYSTEM COSTING COMPLETE---------\n")

    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    

def set_permian_pretreatment_scaling_st2_md(
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

    # SW to NaCl product tranlator 
    set_scaling_factor(
        m.fs.treatment.sw_to_nacl_product.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-2,
    )

    # ZO to SW feed translator
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_feed.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-2,
    )
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_feed.properties_out[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        0.1,
    )


    # ZO to NaCl disposal translator
    # set_scaling_factor(
    #     m.fs.treatment.zo_mix_to_nacl_disposal.properties_in[0].flow_mass_comp["H2O"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_mix_to_nacl_disposal.properties_in[0].flow_mass_comp["tds"],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_mix_to_nacl_disposal.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.zo_mix_to_nacl_disposal.properties_out[0].flow_mass_phase_comp[
    #         "Liq", "NaCl"
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

    # EC inlet
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
    # NaCl DISPOSAL MIXER
    # ZO mixer inlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.zo_mixer_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     100,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.zo_mixer_state[0].flow_mass_phase_comp[
    #         "Liq", "NaCl"
    #     ],
    #     10,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.zo_mixer_state[0].flow_mass_phase_comp[
    #         "Sol", "NaCl"
    #     ],
    #     10,
    # )

    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.md_disposal_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-3,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.md_disposal_state[0].flow_mass_phase_comp[
    #         "Liq", "NaCl"
    #     ],
    #     1e-2,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.md_disposal_state[0].flow_mass_phase_comp[
    #         "Sol", "NaCl"
    #     ],
    #     1e-2,
    # )

    # # mixed state outlet
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.mixed_state[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-1,
    # )
    set_scaling_factor(
        m.fs.treatment.disposal_NaCl_mixer.mixed_state[0].flow_mass_phase_comp[
            "Vap", "H2O"
        ],
        1e-2,
    )
    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.mixed_state[0].flow_mass_phase_comp[
    #         "Liq", "NaCl"
    #     ],
    #     1e-2,
    # )

    # set_scaling_factor(
    #     m.fs.treatment.disposal_NaCl_mixer.mixed_state[0].flow_mass_phase_comp[
    #         "Sol", "NaCl"
    #     ],
    #     10,
    # )

    if calculate_m_scaling_factors:
        print("calclate_m_scaling_factors\n\n\n")
        calculate_scaling_factors(m)


def init_system_st2_md(m, **kwargs):

    treat = m.fs.treatment

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
    propagate_state(treat.disposal_ZO_mix_to_nacl_translator)

    # NOTE: If the selected temperature is similar to md_disposal temperature, the mixer has issues in the solve
    treat.zo_mix_to_nacl_disposal.outlet.temperature[0].fix(300)
    treat.zo_mix_to_nacl_disposal.outlet.pressure[0].fix(101325)
    treat.zo_mix_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treat.zo_mix_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)
    treat.zo_mix_to_nacl_disposal.initialize()

    treat.zo_to_sw_feed.properties_out[0].temperature.fix(300)
    treat.zo_to_sw_feed.properties_out[0].pressure.fix(101325)
    treat.zo_to_sw_feed.initialize()

    propagate_state(treat.cart_filt_translated_to_md)

    init_md(m, treat.md)
    propagate_state(treat.md_to_product)
    propagate_state(treat.md_disposal_to_nacl_translator)

    treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)

    propagate_state(treat.disposal_ZO_mix_translated_to_disposal_NaCl_mixer)
    propagate_state(treat.md_disposal_translated_to_disposal_NaCl_mixer)
    # # NOTE: variable that affects DOF in unclear way

    treat.disposal_NaCl_mixer.initialize()
    treat.disposal_NaCl_mixer.mixed_state[0].temperature.fix(300)
    treat.disposal_NaCl_mixer.mixed_state[0].pressure.fix(101325)

    propagate_state(treat.mixer_to_normalized_feed)


def init_energy_system_st2_md(m):
    init_cst(m.fs.energy.cst)

def add_mec_st2_md(m,permian_cryst_config):

    treat = m.fs.treatment

    print("\n--------- Adding MEC---------\n")
    treat.mec = FlowsheetBlock(dynamic=False)
    build_mec(m, treat.mec, 
              prop_package = m.fs.properties_NaCl,
              prop_package_vapor = m.fs.properties_vapor
            )
    
    set_mec_op_conditions(
                            m, 
                            treat.mec,
                            operating_pressures=permian_cryst_config["operating_pressures"],
                            nacl_yield=permian_cryst_config["nacl_yield"],
                            heat_transfer_coefficient = permian_cryst_config["heat_transfer_coefficient"],
                        )
    
    init_mec(m, treat.mec)

    unfix_mec(treat.mec)

    treat.cryst_feed_H2O_constraint = Constraint(
    expr = treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        == treat.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    treat.cryst_feed_NaCl_constraint = Constraint(
    expr = treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        == treat.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
    )
    treat.cryst_feed_temp_constraint = Constraint(
    expr = treat.mec.unit.inlet.temperature[0]
        == treat.norm_feed.outlet.temperature[0]
    )
    treat.cryst_feed_pressure_constraint = Constraint(
    expr = treat.mec.unit.inlet.pressure[0]
        == treat.norm_feed.outlet.pressure[0]
    )

    print("Water",treat.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())
    print("Nacl",treat.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]())
    print("Temp",treat.norm_feed.outlet.temperature[0]())
    print("Pressure",treat.norm_feed.outlet.pressure[0]())

    mec_rescaling(m.fs,
                  flow_mass_phase_water_total = treat.norm_feed.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]())

    treat.denorm_cryst_product = Denormalizer_Cryst(
        inlet_property_package = m.fs.properties_NaCl,
        outlet_property_package= m.fs.properties_NaCl,
    )

    treat.product_NaCl_mixer = Mixer(
        property_package=m.fs.properties_NaCl,
        num_inlets=2,
        inlet_list=["md_product", "cryst_product"],
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.product = Product(property_package=m.fs.properties_NaCl)

    treat.cryst_product_to_denomalizer = Arc(
        source=treat.mec.unit.outlet, destination=treat.denorm_cryst_product.inlet,
    ) # (8)
    treat.md_translator_to_product_NaCl_mixer = Arc(
        source=treat.sw_to_nacl_product.outlet, destination=treat.product_NaCl_mixer.md_product,
    ) # (7)
    treat.cryst_denomalizer_to_product_NaCl_mixer = Arc(
        source=treat.denorm_cryst_product.outlet, destination=treat.product_NaCl_mixer.cryst_product,
    ) # (9)
    treat.product_NaCl_mixer_to_product = Arc(
        source=treat.product_NaCl_mixer.outlet, destination=treat.product.inlet,
    ) # (10)

    TransformationFactory("network.expand_arcs").apply_to(m)

    treat.sw_to_nacl_product.outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    treat.sw_to_nacl_product.outlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(0)
    treat.sw_to_nacl_product.initialize()
    
    propagate_state(treat.cryst_product_to_denomalizer)
    treat.denorm_cryst_product.initialize()
    
    propagate_state(treat.md_translator_to_product_NaCl_mixer)
    propagate_state(treat.cryst_denomalizer_to_product_NaCl_mixer)

    treat.product_NaCl_mixer.outlet.pressure[0].fix()
    treat.product_NaCl_mixer.initialize()
    
    propagate_state(treat.product_NaCl_mixer_to_product)

    treat.product.properties[0].flow_vol
    treat.product.properties[0].flow_vol_phase
    treat.product.initialize()


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

def run_permian_st2_md(permian_cryst_config, Qin=5, tds=130, grid_frac_heat =0.5, 
                       water_recovery = 0.3, heat_price = 0.00894, electricity_price = 0.04346,
                        treatment_only=False, **kwargs):
    """
    Run Permian pretreatment flowsheet
    """
    rho = get_stream_density(Qin, tds)

    m_pretreatment = build_and_run_permian_pretreatment(Qin, tds)

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

    '''
    Pretreatment + MD flowsheet
    '''

    m = build_permian_st2_md(Q_md=md_flow(), Cin=md_conc(), water_recovery=water_recovery*pretreatment_system_recovery, rho=rho)
    treat = m.fs.treatment

    set_operating_conditions_st2_md(m, rho, Qin, tds)
    set_permian_pretreatment_scaling_st2_md(
        m, calculate_m_scaling_factors=True
    )  # Doesn't solve without this even before costing

    
    init_system_st2_md(m)
    print(f"DOF = {degrees_of_freedom(m)}")

    results = solve(m)
    # print_infeasible_constraints(m)
    assert_optimal_termination(results)

    # print('MD product:', m.fs.treatment.sw_to_nacl_product.display())
    # assert False

    print("\n--------- Before costing solve Completed ---------\n")
    report_MD(m, treat.md)

    '''
    Pretreatment + MD + MEC
    '''

    add_mec_st2_md(m,permian_cryst_config)

    print(f"DOF = {degrees_of_freedom(m)}")

    results = solve(m.fs.treatment.mec)
    results = solve(m)

    print_infeasible_constraints(m)
    assert_optimal_termination(results)
    # assert False
    
    if treatment_only == True:

    # Add costing
        add_treatment_costing_st2_md(m, heat_price, electricity_price)

        print(f"DOF = {degrees_of_freedom(m)}")

        results = solve(m)

        # Update fs.costing block results for only treatment costs
        # Update total heat/electric operating to be treatment aggregate_flow_costs
        # m.fs.costing.total_heat_operating_cost = m.fs.treatment.costing.aggregate_flow_costs['heat']
        # m.fs.costing.total_electric_operating_cost = m.fs.treatment.costing.aggregate_flow_costs['electricity']
        
        # Update the total_capital, total_operating
        # m.fs.energy.cst.unit.costing.capital_cost.fix(1e-25)
        # m.fs.energy.cst.unit.costing.fixed_operating_cost.fix(1e-25)
        # m.fs.energy.cst.unit.costing.variable_operating_cost.fix(1e-25)

        # # Update LCOT to be LCOW
        # m.fs.costing.LCOT.fix(m.fs.treatment.costing.LCOW())
        # # Update grid fraction reported
        # m.fs.costing.frac_heat_from_grid.fix(1)

    else:
        build_energy_permian_st2_md(m)
        set_energy_operating_conditions_st2_md(m)
        init_energy_system_st2_md(m)

        print("\n--------- CST Inputs Completed ---------\n")

        print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
        print('CST Heat:', value(m.fs.energy.cst.unit.heat))
        print("\n")
            
        add_system_costing_st2_md(m, heat_price, electricity_price)
        add_cst_costing_scaling(m,m.fs.energy.cst.unit)

        results = solve(m)

        m.fs.energy.cst.unit.heat_load.unfix()
        m.fs.costing.frac_heat_from_grid.fix(grid_frac_heat)

        print(f"DOF = {degrees_of_freedom(m)}")

        results = solve(m)

        print("\n--------- Costing solve with fixed grid fraction ---------\n")

        print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
        print('CST Heat:', value(m.fs.energy.cst.unit.heat))
        print('Grid fraction:',value(m.fs.costing.frac_heat_from_grid))

    permian_md_reporting_variables(m)

    return m

def print_results(m):
    treat = m.fs.treatment
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

    print("\n--------- CST Inputs Completed ---------\n")

    print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
    print('CST Heat:', value(m.fs.energy.cst.unit.heat))
    print("\n")


    print(f"\n\n-------------------- Pretreatment Report --------------------\n")

    print("\n")
    print(
        f'{"Pretreatment Recovery":<30s}{system_recovery:.2f}'
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
        f'{"Product TDS conc":<30s} {treat.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"]():.<10.2f} {pyunits.get_units(treat.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"]())}'
    )

    print(
        f'{"Product flow_vol":<30s} {treat.product.properties[0].flow_vol_phase["Liq"]():<10.2f} {pyunits.get_units(treat.product.properties[0].flow_vol_phase["Liq"])}'
    )

    print(
        f'{"Aggregated Heat Cost":<30s}{value(m.fs.treatment.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(m.fs.treatment.costing.aggregate_flow_costs["heat"])}'
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

def sweep_feed_flow_salinity():

    permian_cryst_config = {
    "operating_pressures": [0.45, 0.25, 0.208, 0.095], # Operating pressure of each effect (bar)
    "nacl_yield": 0.9, # Yield
    "heat_transfer_coefficient": 0.13
    }

    heat_price = 0.00894
    electricity_price = 0.04346

    sweep_dict = {
        'Qin':[1,5,9],
        'tds': [100,130,200],
        'recovery': [0.6,0.48,0.23]
    }

    for flow in sweep_dict["Qin"]:
        for i in range(0,len(sweep_dict['tds'])):
            m = run_permian_st2_md(
                        Qin=flow, 
                        tds=sweep_dict['tds'][i], 
                        water_recovery = sweep_dict['recovery'][i],
                        grid_frac_heat = 0.5,
                        heat_price=heat_price, 
                        electricity_price=electricity_price, 
                        permian_cryst_config=permian_cryst_config
                        )

def main():
    permian_cryst_config = {
    "operating_pressures": [0.45, 0.25, 0.208, 0.095], # Operating pressure of each effect (bar)
    "nacl_yield": 0.9, # Yield
    "heat_transfer_coefficient": 0.13
    }

    heat_price = 0.00894
    electricity_price = 0.04346  # Updated 0.0575 in USD 2018 to USD 2023

    m = run_permian_st2_md(
                        Qin=5, 
                        tds=130, 
                        water_recovery = 0.5,
                        grid_frac_heat = 1,
                        heat_price=heat_price, 
                        electricity_price=electricity_price, 
                        permian_cryst_config=permian_cryst_config,
                        treatment_only=True
                        )
    
   
    print(f"DOF = {degrees_of_freedom(m)}")

    # report_costing(m.fs.costing)

    feed_m3h = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.h
    )

    product_m3h = pyunits.convert(
        m.fs.treatment.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.m**3 / pyunits.h
    )
    print("\nProduct flow in m3/h:",product_m3h())
    
    print("LCOW:",m.fs.treatment.costing.LCOW(),pyunits.get_units(m.fs.treatment.costing.LCOW))
    print("SEC (electricity) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_electricity()/product_m3h())
    print("SEC (heat) in kWh/m3:",m.fs.treatment.costing.aggregate_flow_heat()/product_m3h())
    print("System recovery (%):", product_m3h()/feed_m3h()*100)
    print("Capex ($M):",m.fs.treatment.costing.total_capital_cost()/1e6, 
          pyunits.get_units(m.fs.treatment.costing.total_capital_cost))
    print("Opex ($M/yr):",m.fs.treatment.costing.total_operating_cost()/1e6, 
          pyunits.get_units(m.fs.treatment.costing.total_operating_cost))
    print('Electricity demand (MWh/year):',pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity,to_units=pyunits.MW*pyunits.h/pyunits.year)())
    print('Heat demand (MWh/year):',pyunits.convert(m.fs.treatment.costing.aggregate_flow_heat,to_units=pyunits.MW*pyunits.h/pyunits.year)())


            

if __name__ == "__main__":

    main()

    
    