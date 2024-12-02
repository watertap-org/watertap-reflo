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
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state as _prop_state
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
from watertap.core import Database
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock as MCAS,
)

# from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.unit_specific.cryst_prop_pack import NaClParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

from watertap_contrib.reflo.analysis.example_flowsheets.fo_trevi_flowsheet import (
    build_fo_trevi_flowsheet,
    fix_dof_and_initialize,
    get_flowsheet_performance,
)
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import FODrawSolutionParameterBlock

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3

solver = get_solver()

def build_permian_FO():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.treatment = treat = Block()

    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_draw = FODrawSolutionParameterBlock()
    m.fs.properties_nacl = NaClParameterBlock()

    treat.feed = Feed(property_package=m.fs.properties)
    treat.product = Product(property_package=m.fs.properties_feed)

    # Add translator blocks
    treat.zo_to_sw_feed = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )
    treat.zo_to_sw_ec_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )
    treat.zo_to_sw_cart_filt_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )
    treat.draw_to_sw = Translator_Draw_to_SW(
        inlet_property_package = m.fs.properties_draw,
        outlet_property_package= m.fs.properties_feed,
    )
    treat.sw_to_nacl = Translator_SW_to_NaCl(
        inlet_property_package = m.fs.properties_feed,
        outlet_property_package= m.fs.properties_nacl,
    )

    # Add components
    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.ec = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.ec)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.FO = build_fo_trevi_flowsheet(feed_vol_flow=0.22, # m3/s
                                        feed_TDS_mass=0.119, # mass fraction, 0.119 is about 130 g/L
                                        recovery_ratio=0.45,
                                        RO_recovery_ratio=0.9,  # RO recovery ratio
                                        NF_recovery_ratio=0.8,  # Nanofiltration recovery ratio
                                        feed_temperature=25,
                                        strong_draw_temp=25,  # Strong draw solution inlet temperature (C)
                                        strong_draw_mass=0.95,  # Strong draw solution mass fraction
                                        )

    treat.mec = FlowsheetBlock(dynamic=False)
    build_mec(m, treat.mec)

    treat.disposal_SW_mixer = Mixer(
        property_package=m.fs.properties_feed,
        num_inlets=3,
        inlet_list=["ec_disposal", "cart_filt_disposal", "fo_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    # BUILD PRODUCT STREAM
    # feed (1)> chem_addition (2)> EC (3)> cart_filt (4)> ZO_to_SW_translator (5)> FO (6)> Draw_to_SW_translator (7)> product
    
    treat.feed_to_chem_addition = Arc(
        source=treat.feed.outlet, destination=treat.chem_addition.feed.inlet
    ) # (1)
    treat.chem_addition_to_ec = Arc(
        source=treat.chem_addition.product.outlet, destination=treat.ec.feed.inlet
    ) # (2)
    treat.ec_to_cart_filt = Arc(
        source=treat.ec.product.outlet, destination=treat.cart_filt.feed.inlet
    ) # (3)
    treat.cart_filt_to_translator = Arc(
        source=treat.cart_filt.product.outlet, destination=treat.zo_to_sw_feed.inlet
    ) # (4)
    treat.cart_filt_translated_to_fo = Arc(
        source=treat.zo_to_sw_feed.outlet, destination=treat.FO.fs.fo.feed
    ) # (5)
    treat.fo_to_translator = Arc(
        source=treat.FO.fs.S2.fresh_water, destination=treat.draw_to_sw.inlet
    ) # (6)
    treat.fo_translator_to_product = Arc(
        source=treat.draw_to_sw.outlet, destination=treat.product.inlet
    ) # (7)

    # BUILD DISPOSAL STREAM
    #        EC (1)> ZO_to_SW_translator (3)> disposal_mixer (6)> SW_to_NaCl_translator (7)> crystallizer
    # cart_filt (2)> ZO_to_SW_translator (4)> disposal_mixer
    #                                FO  (5)> disposal_mixer

    treat.ec_disposal_to_translator = Arc(
        source=treat.ec.disposal.outlet, destination=treat.zo_to_sw_ec_disposal.inlet
    ) # (1)
    treat.cart_filt_disposal_to_translator = Arc(
        source=treat.cart_filt.disposal.outlet, destination=treat.zo_to_sw_cart_filt_disposal.inlet,
    ) # (2)
    treat.ec_disposal_translator_to_SW_mixer = Arc(
        source=treat.zo_to_sw_ec_disposal.outlet, destination=treat.disposal_SW_mixer.ec_disposal,
    ) # (3)
    treat.cart_filt_disposal_translator_to_SW_mixer = Arc(
        source=treat.zo_to_sw_cart_filt_disposal.outlet, destination=treat.disposal_SW_mixer.cart_filt_disposal,
    ) # (4)
    treat.fo_disposal_translator_to_SW_mixer = Arc(
        source=treat.FO.fs.fo.brine, destination=treat.disposal_SW_mixer.fo_disposal,
    ) # (5)
    treat.SW_mixer_to_translator = Arc(
        source=treat.disposal_SW_mixer.fo_disposal, destination=treat.sw_to_nacl.inlet,
    ) # (6)
    treat.SW_mixer_translator_to_mec = Arc(
        source=treat.sw_to_nacl.outlet, destination=treat.mec.unit.inlet,
    ) # (7)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m

def set_operating_conditions(m, Qin=5, tds=130, **kwargs):

    global flow_mass_water, flow_mass_tds, flow_in

    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_water = 211.62 * pyunits.kg / pyunits.s
    flow_mass_tds = 28.58 * pyunits.kg / pyunits.s

    m.fs.treatment.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)
    m.fs.treatment.feed.properties[0].flow_mass_comp["tds"].fix(flow_mass_tds)

    m.fs.treatment.feed.properties[0].conc_mass_comp[...]

    set_chem_addition_op_conditions(m, m.fs.treatment.chem_addition, **kwargs)
    set_ec_operating_conditions(m, m.fs.treatment.ec, **kwargs)
    set_cart_filt_op_conditions(m, m.fs.treatment.cart_filt)

    operating_pressures = [0.45, 0.25, 0.208, 0.095]
    set_mec_op_conditions(m, m.fs.treatment.mec, operating_pressures=operating_pressures)

    # m.fs.treatment.mec.unit.inlet.temperature[0].fix(273.15 + 20)

def set_permian_scaling(m, **kwargs):

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        1 / value(flow_mass_water),
        index=("H2O"),
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        1 / value(flow_mass_tds),
        index=("tds"),
    )

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_water),
        index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_tds),
        index=("Liq", "TDS")
    )

    calculate_scaling_factors(m)

def init_system(m):
    treat = m.fs.treatment

    treat.feed.initialize()
    propagate_state(treat.feed_to_chem_addition)   

    init_chem_addition(m, treat.chem_addition)
    propagate_state(treat.chem_addition_to_ec)

    init_ec(m, treat.ec)
    propagate_state(treat.ec_to_cart_filt)
    propagate_state(treat.ec_disposal_to_translator)

    init_cart_filt(m, treat.cart_filt)
    propagate_state(arc = treat.cart_filt_to_translator)
    propagate_state(arc = treat.cart_filt_translated_to_fo)
    propagate_state(treat.cart_filt_disposal_to_translator)

    treat.zo_to_sw_feed.initialize()

    fix_dof_and_initialize(treat.FO,                
                           strong_draw_mass_frac=0.90,
                           product_draw_mass_frac=0.01,
                           RO_recovery_ratio=0.9,
                           NF_recovery_ratio=0.8,)
    
    # unfix FO fs and set operating point
    treat.FO.fs.HX1.area.unfix()
    treat.FO.fs.HX2.area.unfix()
    treat.FO.fs.HX1.weak_draw_outlet.temperature.fix(78 + 273.15)
    treat.FO.fs.HX1.product_water_outlet.temperature.fix(32 + 273.15)

    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    propagate_state(arc = treat.fo_to_translator)
    propagate_state(arc = treat.fo_translator_to_product)

    treat.product.initialize()

    treat.zo_to_sw_ec_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_ec_disposal.outlet.pressure[0].fix(101325)
    treat.zo_to_sw_cart_filt_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_cart_filt_disposal.outlet.pressure[0].fix(101325)
    treat.FO.fs.fo.brine.pressure[0].fix(101325)

    treat.zo_to_sw_ec_disposal.initialize()
    treat.zo_to_sw_cart_filt_disposal.initialize()
    
    treat.disposal_SW_mixer.initialize()
    propagate_state(arc=treat.SW_mixer_to_translator)
    propagate_state(arc=treat.SW_mixer_translator_to_mec)

    print('')
    print('')
    print('mec init started')
    init_mec_and_unfix(treat.mec)
    print('')
    print('')
    print('mec init finished')
    mec_rescaling(treat.mec,
                  flow_mass_phase_water_total = 112,
                  flow_mass_phase_salt_total = 28)

def add_treatment_costing(m):

    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.ec, flowsheet_costing_block=m.fs.treatment.costing)
    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.FO.fs.fo.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.treatment.costing)
    
    # add_dwi_costing(
    #     m, m.fs.treatment.DWI, flowsheet_costing_block=m.fs.treatment.costing
    # )

    m.fs.treatment.costing.cost_process()

def run_permian_FO(recovery=0.5):
    m = build_permian_FO()
    treat = m.fs.treatment

    set_operating_conditions(m)

if __name__ == "__main__":

    m = build_permian_FO()
    treat = m.fs.treatment

    set_operating_conditions(m)
    set_permian_scaling(m)

    treat.feed.properties[0].flow_vol

    init_system(m)

    print('DOF before solving: ', degrees_of_freedom(m))
    results = solver.solve(m)
    assert_optimal_termination(results)

    add_treatment_costing(m)
    flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]
    treat.costing.electricity_cost.fix(0.07)
    treat.costing.add_LCOW(flow_vol)
    treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
    treat.costing.initialize()

    print(f"DOF after add costing: {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)