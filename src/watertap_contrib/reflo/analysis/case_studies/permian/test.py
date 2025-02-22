import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Block,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import (
    Product,
    Feed,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *
from idaes.core.util.initialization import propagate_state

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
# from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

from watertap_contrib.reflo.analysis.example_flowsheets.fo_trevi_flowsheet import (
    build_fo_trevi_flowsheet,
    fix_dof_and_initialize,
)
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import FODrawSolutionParameterBlock

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3

solver = get_solver()

__all__ = [
    "build_permian_FO",
    "set_operating_conditions",
    "add_treatment_costing",
    "add_energy_costing",
    "set_permian_scaling",
    "init_system",
    "run_permian_FO",
]

def build_permian_FO(permian_fo_config):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.treatment = treat = Block()
    m.fs.energy = energy = Block()

    # m.fs.energy.cst = FlowsheetBlock()
    # build_cst(m.fs.energy.cst)

    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_draw = FODrawSolutionParameterBlock()

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

    # # Add components
    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.ec = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.ec)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.FO = build_fo_trevi_flowsheet(feed_vol_flow     =permian_fo_config["feed_vol_flow"], 
                                        feed_TDS_mass     =permian_fo_config["feed_TDS_mass"], 
                                        recovery_ratio    =permian_fo_config["recovery_ratio"],
                                        RO_recovery_ratio =permian_fo_config["RO_recovery_ratio"], 
                                        NF_recovery_ratio =permian_fo_config["NF_recovery_ratio"], 
                                        feed_temperature  =permian_fo_config["feed_temperature"],
                                        strong_draw_temp  =permian_fo_config["strong_draw_temp"],  
                                        strong_draw_mass  =permian_fo_config["strong_draw_mass_frac"],  
                                        )

    treat.DWI = DeepWellInjection(
        property_package=m.fs.properties_feed, injection_well_depth=5000
    )

    treat.disposal_SW_mixer = Mixer(
        property_package=m.fs.properties_feed,
        num_inlets=3,
        inlet_list=["ec_disposal", "cart_filt_disposal", "fo_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    # BUILD PRODUCT STREAM
    # feed (1)> chem_addition (2)> EC (3)> cart_filt 
    #      (4)> ZO_to_SW_translator (5)> FO (6)> Draw_to_SW_translator (7)> product
    
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
    #        EC (1)> ZO_to_SW_translator (3)> disposal_mixer (6)> DWI
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
    treat.SW_mixer_to_DWI = Arc(
        source=treat.disposal_SW_mixer.outlet, destination=treat.DWI.inlet,
    ) # (6)

    TransformationFactory("network.expand_arcs").apply_to(m)
    
    return m

def set_operating_conditions(m, operating_condition, **kwargs):
    Qin, tds = operating_condition["feed_vol_flow"], operating_condition["feed_tds"]

    global flow_mass_water, flow_mass_tds, flow_in

    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )

    m.fs.treatment.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)
    m.fs.treatment.feed.properties[0].flow_mass_comp["tds"].fix(flow_mass_tds)

    m.fs.treatment.feed.properties[0].conc_mass_comp[...]

    set_chem_addition_op_conditions(m, m.fs.treatment.chem_addition, **kwargs)
    set_ec_operating_conditions(m, m.fs.treatment.ec, **kwargs)
    set_cart_filt_op_conditions(m, m.fs.treatment.cart_filt)

    # Set energy system condition
    # set_cst_op_conditions(m.fs.energy.cst)

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

    m.fs.properties_draw.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_water),
        index=("Liq", "H2O")
    )

    m.fs.properties_draw.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_water),
        index=("Liq", "DrawSolution")
    )

    calculate_scaling_factors(m)

def init_system(m, permian_fo_config, CST_config):
    treat = m.fs.treatment

    treat.feed.initialize()
    propagate_state(arc = treat.feed_to_chem_addition)   

    init_chem_addition(m, treat.chem_addition)
    propagate_state(arc = treat.chem_addition_to_ec)

    init_ec(m, treat.ec)
    propagate_state(arc = treat.ec_to_cart_filt)
    propagate_state(arc = treat.ec_disposal_to_translator)
    propagate_state(arc = treat.ec_disposal_translator_to_SW_mixer)

    init_cart_filt(m, treat.cart_filt)
    propagate_state(arc = treat.cart_filt_to_translator)
    treat.zo_to_sw_feed.initialize()

    # Unfix FO feed properties so that it's propogated from filtration outlet
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
    propagate_state(arc = treat.cart_filt_translated_to_fo)
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].fix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].fix()


    propagate_state(arc = treat.cart_filt_disposal_to_translator)
    treat.zo_to_sw_cart_filt_disposal.initialize()
    propagate_state(arc = treat.cart_filt_disposal_translator_to_SW_mixer)

    fix_dof_and_initialize(treat.FO,                
                           strong_draw_mass_frac =permian_fo_config["strong_draw_mass_frac"],
                           product_draw_mass_frac=permian_fo_config["product_draw_mass_frac"],
                           RO_recovery_ratio     =permian_fo_config["RO_recovery_ratio"],
                           NF_recovery_ratio     =permian_fo_config["NF_recovery_ratio"],)
    
    # unfix FO fs and set operating point
    treat.FO.fs.HX1.area.unfix()
    treat.FO.fs.HX2.area.unfix()
    treat.FO.fs.HX1.weak_draw_outlet.temperature.fix(permian_fo_config["HX1_cold_out_temp"])
    treat.FO.fs.HX1.product_water_outlet.temperature.fix(permian_fo_config["HX1_hot_out_temp"])

    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    propagate_state(arc = treat.fo_to_translator)
    treat.draw_to_sw.initialize()

    # @treat.Constraint(doc="FO_product")
    # def FO_product_temp(b):
    #     return b.draw_to_sw.inlet.temperature[0] == b.FO.fs.S2.fresh_water.temperature[0] 
    
    # @treat.Constraint(doc="FO_product")
    # def FO_product_pres(b):
    #     return b.draw_to_sw.inlet.pressure[0] == b.FO.fs.S2.fresh_water.pressure[0] 
    
    # treat.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'H2O'].fix(100)
    # treat.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'DrawSolution'].fix(0)
    
    # @treat.Constraint(doc="FO_product")
    # def FO_product_H2O(b):
    #     return b.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'H2O'] == b.FO.fs.S2.fresh_water.flow_mass_phase_comp[0,"Liq",'H2O'] 

    # @treat.Constraint(doc="FO_product")
    # def FO_product_TDS(b):
    #     return b.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'DrawSolution'] == b.FO.fs.S2.fresh_water.flow_mass_phase_comp[0,"Liq",'DrawSolution'] 
    
    # treat.draw_to_sw.inlet.temperature[0] = treat.FO.fs.fo.brine.temperature[0]
    # treat.draw_to_sw.inlet.pressure[0] = treat.FO.fs.fo.brine.pressure[0]
    # treat.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'DrawSolution'] = 0
    # treat.draw_to_sw.inlet.flow_mass_phase_comp[0,"Liq",'H2O'] = treat.FO.fs.fo.brine.flow_mass_phase_comp[0,"Liq",'H2O']

    propagate_state(arc = treat.fo_translator_to_product)
    treat.product.initialize()

    # @treat.Constraint(doc="FO_product")
    # def FO_product_H2O(b):
    #     return b.product.flow_mass_phase_comp[0,"Liq",'H2O'] == 100#b.FO.fs.S2.fresh_water.flow_mass_phase_comp[0,"Liq",'H2O'] 

    # @treat.Constraint(doc="FO_product")
    # def FO_product_TDS(b):
    #     return b.product.flow_mass_phase_comp[0,"Liq",'TDS'] == 10#b.FO.fs.S2.fresh_water.flow_mass_phase_comp[0,"Liq",'DrawSolution'] 
    
    # @treat.Constraint(doc="FO_product")
    # def FO_product_temp(b):
    #     return b.product.temperature[0] == b.FO.fs.S2.fresh_water.temperature[0] 
    
    propagate_state(arc = treat.fo_disposal_translator_to_SW_mixer)

    treat.zo_to_sw_ec_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_ec_disposal.outlet.pressure[0].fix(101325)
    treat.zo_to_sw_cart_filt_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_cart_filt_disposal.outlet.pressure[0].fix(101325)
    treat.FO.fs.fo.brine.pressure[0].fix(101325)

    treat.zo_to_sw_ec_disposal.initialize()
    treat.zo_to_sw_cart_filt_disposal.initialize()
    
    treat.disposal_SW_mixer.initialize()
    propagate_state(arc = treat.SW_mixer_to_DWI)

    treat.DWI.properties[0].temperature.fix()
    treat.DWI.properties[0].pressure.fix()
    treat.DWI.initialize()

    # init_cst(m.fs.energy.cst, 
    #          storage=CST_config['storage'], 
    #          heat_load=CST_config['heat_load'])

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
    
    m.fs.treatment.DWI.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.treatment.costing)

    m.fs.treatment.costing.cost_process()

def add_energy_costing(m, CST_config):
    energy = m.fs.energy
    energy.costing = EnergyCosting()
    energy.costing.base_currency = pyunits.USD_2023

    add_cst_costing(m.fs.energy.cst, costing_block=m.fs.energy.costing)


    energy.costing.heat_cost.set_value(0)
    energy.costing.cost_process()
    energy.costing.initialize()


    energy.cst.unit.heat_load.unfix()
    energy.costing.aggregate_flow_heat.fix(CST_config["heat_flow"])


def run_permian_FO(operating_condition,
                   permian_fo_config,
                   CST_config):
    m = build_permian_FO(permian_fo_config)
    treat = m.fs.treatment

    set_operating_conditions(m, operating_condition)
    set_permian_scaling(m)

    treat.feed.properties[0].flow_vol

    init_system(m, permian_fo_config, CST_config)

    print('DOF after init: ', degrees_of_freedom(m))
    results = solver.solve(m)
    assert_optimal_termination(results)

    add_treatment_costing(m)
    # add_energy_costing(m,CST_config)

    flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]
    # flow_vol = treat.FO.fs.S2.fresh_water_state[0].flow_vol_phase["Liq"]

    treat.costing.base_currency = pyunits.USD_2023
    treat.costing.electricity_cost.fix(0.0575)
    treat.costing.heat_cost.fix(0.00894)
    treat.costing.add_LCOW(flow_vol)
    treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
    treat.costing.initialize()

    # scaling (based on grid participation), setup order
    # deactivate constraints, 
    heat_price = 0.01    # $/kWh-th
    electricity_price = 0.07 # $/kWh-e
    # m.fs.costing = REFLOSystemCosting()
    # m.fs.costing.heat_cost_buy.fix(heat_price)
    # m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    # m.fs.costing.cost_process()

    # m.fs.energy.cst.unit.heat_load.unfix()
    # m.fs.energy.costing.aggregate_flow_heat.unfix()
    # m.fs.costing.frac_heat_from_grid.fix(0.05)

    # m.fs.costing.initialize()
    # m.fs.costing.add_LCOH()
    # m.fs.costing.add_LCOW(flow_vol)
    # m.fs.costing.add_LCOT(flow_vol)


    print(f"DOF after add costing: {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


if __name__ == "__main__":
    fail=[]
    heat=[]
    fo_feed=[]
    brine=[]
    grid_frac =[]
    # for rr in [0.52 + 0.005*i for i in range(8)]:
    rr = 0.4
    permian_fo_config = {
    "feed_vol_flow": 0.22, # initial value for fo model setup
    "feed_TDS_mass": 0.039, # mass fraction, 0.119 is about 130 g/L, 0.092 for 100 g/L, 0.19 for 200 g/L
    "recovery_ratio": rr,
    "RO_recovery_ratio":1,  # RO recovery ratio
    "NF_recovery_ratio":0.8,  # Nanofiltration recovery ratio
    "feed_temperature":25,
    "strong_draw_temp":25,  # Strong draw solution inlet temperature (C)
    "strong_draw_mass_frac":0.9,  # Strong draw solution mass fraction
    "product_draw_mass_frac": 0.01,   # FO product draw solution mass fraction
    "HX1_cold_out_temp": 78 + 273.15, # HX1 coldside outlet temperature
    "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
    }

    CST_config = {
        "storage":12, # hr
        "heat_load":25, # MW
        "heat_flow": -5000, # kW
    }

    operating_condition = {
    "feed_vol_flow": 5, # MGD
    "feed_tds": 130 # g/L
    }
    # try:
    m = run_permian_FO(operating_condition,
                    permian_fo_config,
                    CST_config,
                    )
    # heat.append((rr,value(m.fs.treatment.FO.fs.fo.costing.thermal_energy_flow)))
    fo_feed.append((rr, value(m.fs.treatment.FO.fs.fo.feed_props[0].conc_mass_phase_comp["Liq","TDS"])))
    brine.append((rr, value(m.fs.treatment.FO.fs.fo.brine_props[0].conc_mass_phase_comp["Liq","TDS"])))
        # grid_frac.append((rr,m.fs.costing.frac_heat_from_grid.value))
    # except:
    #     fo_feed.append((rr,'fail'))
    #     brine.append((rr,'fail'))
            # heat.append((rr,'fail'))
            # grid_frac.append((rr,'fail'))
        