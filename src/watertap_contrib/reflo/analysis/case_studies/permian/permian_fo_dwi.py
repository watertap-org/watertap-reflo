#%%
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
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

from watertap_contrib.reflo.analysis.example_flowsheets.fo_trevi_flowsheet import (
    build_fo_trevi_flowsheet,
    fix_dof_and_initialize,
    get_flowsheet_performance,
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

    m.fs.energy.cst = FlowsheetBlock()
    build_cst(m.fs.energy.cst)

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

    # Add components
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

    calculate_scaling_factors(m)

def init_system(m, permian_fo_config):
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
    propagate_state(arc = treat.SW_mixer_to_DWI)

    treat.DWI.properties[0].temperature.fix()
    treat.DWI.properties[0].pressure.fix()
    treat.DWI.initialize()

    init_cst(m.fs.energy.cst)

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

def add_energy_costing(m):
    energy = m.fs.energy
    energy.costing = EnergyCosting()

    energy.cst.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
    )

    energy.costing.heat_cost.set_value(0)
    energy.costing.cost_process()
    energy.costing.initialize()


    energy.cst.unit.heat_load.unfix()
    energy.costing.aggregate_flow_heat.fix(0)

def run_permian_FO(operating_condition,
                   permian_fo_config):
    m = build_permian_FO(permian_fo_config)
    treat = m.fs.treatment

    set_operating_conditions(m, operating_condition)
    set_permian_scaling(m)

    treat.feed.properties[0].flow_vol

    init_system(m, permian_fo_config)

    print('DOF after init: ', degrees_of_freedom(m))
    results = solver.solve(m)
    assert_optimal_termination(results)

    add_treatment_costing(m)
    # add_energy_costing(m)

    flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]
    treat.costing.base_currency = pyunits.USD_2021
    treat.costing.electricity_cost.fix(0.07)
    treat.costing.heat_cost.fix(0.02)
    treat.costing.add_LCOW(flow_vol)
    treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
    treat.costing.initialize()

    # scaling (based on grid participation), setup order
    # deactivate constraints, 
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.base_currency = pyunits.USD_2021
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(flow_vol)
    m.fs.costing.add_LCOW(flow_vol)
    m.fs.costing.initialize()
    

    print(f"DOF after add costing: {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m

if __name__ == "__main__":
    permian_fo_config = {
    "feed_vol_flow": 0.22, # initial value for fo model setup
    "feed_TDS_mass": 0.119, # mass fraction, 0.119 is about 130 g/L
    "recovery_ratio": 0.45,
    "RO_recovery_ratio":1,  # RO recovery ratio
    "NF_recovery_ratio":0.8,  # Nanofiltration recovery ratio
    "feed_temperature":25,
    "strong_draw_temp":25,  # Strong draw solution inlet temperature (C)
    "strong_draw_mass_frac":0.9,  # Strong draw solution mass fraction
    "product_draw_mass_frac": 0.01,   # FO product draw solution mass fraction
    "HX1_cold_out_temp": 78 + 273.15, # HX1 coldside outlet temperature
    "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
    }

    operating_condition = {
    "feed_vol_flow": 5, # MGD
    "feed_tds": 130 # g/L
    }

    m = run_permian_FO(operating_condition,
                       permian_fo_config,
                       )
    
    #%%
    lcow = value(m.fs.treatment.costing.LCOW)
    capex_total = value(m.fs.treatment.costing.total_capital_cost)
    chem_capex = m.fs.treatment.chem_addition.unit.costing.capital_cost()
    filt_capex = m.fs.treatment.cart_filt.unit.costing.capital_cost()
    ec_capex = m.fs.treatment.ec.unit.costing.capital_cost()
    fo_capex = m.fs.treatment.FO.fs.fo.costing.capital_cost()
    dwi_capex = m.fs.treatment.DWI.costing.capital_cost()

    opex_total = value(m.fs.treatment.costing.total_operating_cost)
    fix_opex = value(m.fs.treatment.costing.maintenance_labor_chemical_operating_cost)
    chem_opex = chem_capex * 0.03
    filt_opex = filt_capex * 0.03
    ec_opex = ec_capex * 0.03
    fo_opex = fo_capex * 0.03
    dwi_opex = dwi_capex * 0.03

    # fo_opex = m.fs.treatment.FO.fs.fo.costing.fixed_operating_cost()
    # dwi_opex = m.fs.treatment.DWI.costing.fixed_operating_cost()


    var_opex_total = value(m.fs.treatment.costing.total_variable_operating_cost)
    elec_cost = value(m.fs.treatment.costing.aggregate_flow_costs["electricity"])
    heat_cost = value(m.fs.treatment.costing.aggregate_flow_costs["heat"])
    alum_cost = value(m.fs.treatment.costing.aggregate_flow_costs["aluminum"])
    h2o2_cost = value(m.fs.treatment.costing.aggregate_flow_costs["hydrogen_peroxide"])


    chem_elec_cost = 0.07 * value(m.fs.treatment.chem_addition.unit.electricity[0]) * 8766
    ec_elec_cost = 0.07 * value(m.fs.treatment.ec.unit.costing.electricity_flow) * 8766
    filt_elec_cost = 0.07 * value(m.fs.treatment.cart_filt.unit.electricity[0]) * 8766
    fo_elec_cost = 0.07 * value(m.fs.treatment.FO.fs.fo.costing.electricity_flow) * 8766
    fo_heat_cost = 0.02 * value(m.fs.treatment.FO.fs.fo.costing.thermal_energy_flow) * 8766
    dwi_elec_cost = 0.07 * value(m.fs.treatment.DWI.costing.pumping_power_required) * 8766

    print('LCOW                     ', lcow)

    print('TOTAL CAPEX              ', capex_total)
    print('    chem_addition capex     ', chem_capex)
    print('    elec_coag capex       ', ec_capex)
    print('    cart_filtration capex ', filt_capex)
    print('    FO capex              ', fo_capex)
    print('    DWI capex             ', dwi_capex)
    print('')

    print('TOTAL OPEX           ', opex_total)
    print('     Total fixed OPEX ', fix_opex)
    print('           chem_opex    ', chem_opex)
    print('           filt_opex   ', filt_opex)
    print('           ec_opex     ', ec_opex)
    print('           fo_opex    ', fo_opex)
    print('           dwi_opex   ', dwi_opex)
    print('')
    print('     TOTAL VAR OPEX      ', var_opex_total)
    print('           Elec cost      ', elec_cost)
    print('                chem elec      ', chem_elec_cost)
    print('                ec elec   ', ec_elec_cost)
    print('                filt elec    ', filt_elec_cost)
    print('                FO elec   ', fo_elec_cost)
    print('                DWI elec   ', dwi_elec_cost)
    print('')
    print('           Heat cost     ', heat_cost)
    print('                FO heat  ', fo_heat_cost)
    print('')
    print('           Al cost       ', alum_cost)
    print('           H2O2 cost      ', h2o2_cost)

    #%%
    chem_capexs = []
    ec_capexs = []
    filt_capexs = []
    fo_capexs = []
    dwi_capexs = []

    chem_opexs =[]
    ec_opexs = []
    filt_opexs = []
    fo_opexs = []
    dwi_opexs = []

    elecs = []
    heats = []
    alums = []
    h2o2s = []

    LCOWs =[]
    failed= []

    rr = [ 0.2,0.24,0.28,0.32,0.36,0.4,0.44, 0.45]
    strong_draw_mass = [i*0.03 + 0.80 for i in range(6)]

    for v in strong_draw_mass:
        permian_fo_config = {
    "feed_vol_flow": 0.22, # initial value for fo model setup
    "feed_TDS_mass": 0.119, # mass fraction, 0.119 is about 130 g/L
    "recovery_ratio": 0.4,
    "RO_recovery_ratio":1,  # RO recovery ratio
    "NF_recovery_ratio":0.8,  # Nanofiltration recovery ratio
    "feed_temperature":25,
    "strong_draw_temp":25,  # Strong draw solution inlet temperature (C)
    "strong_draw_mass_frac":v,  # Strong draw solution mass fraction
    "product_draw_mass_frac": 0.01,   # FO product draw solution mass fraction
    "HX1_cold_out_temp": 78 + 273.15, # HX1 coldside outlet temperature
    "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
    }
        operating_condition = {
    "feed_vol_flow": 5, # MGD
    "feed_tds": 130 # g/L
    }
        try:
            m = run_permian_FO(operating_condition,
                        permian_fo_config,
                        )
        except:
            failed.append(v)

        lcow = value(m.fs.treatment.costing.LCOW)
        capex_total = value(m.fs.treatment.costing.total_capital_cost)
        chem_capex = m.fs.treatment.chem_addition.unit.costing.capital_cost()
        filt_capex = m.fs.treatment.cart_filt.unit.costing.capital_cost()
        ec_capex = m.fs.treatment.ec.unit.costing.capital_cost()
        fo_capex = m.fs.treatment.FO.fs.fo.costing.capital_cost()
        dwi_capex = m.fs.treatment.DWI.costing.capital_cost()

        opex_total = value(m.fs.treatment.costing.total_operating_cost)
        fix_opex = value(m.fs.treatment.costing.maintenance_labor_chemical_operating_cost)
        chem_opex = chem_capex * 0.03
        filt_opex = filt_capex * 0.03
        ec_opex = ec_capex * 0.03
        fo_opex = fo_capex * 0.03
        dwi_opex = dwi_capex * 0.03

        fo_elec_cost = 0.07 * value(m.fs.treatment.FO.fs.fo.costing.electricity_flow) * 8766
        fo_heat_cost = 0.02 * value(m.fs.treatment.FO.fs.fo.costing.thermal_energy_flow) * 8766


        var_opex_total = value(m.fs.treatment.costing.total_variable_operating_cost)
        elec_cost = value(m.fs.treatment.costing.aggregate_flow_costs["electricity"])
        heat_cost = value(m.fs.treatment.costing.aggregate_flow_costs["heat"])
        alum_cost = value(m.fs.treatment.costing.aggregate_flow_costs["aluminum"])
        h2o2_cost = value(m.fs.treatment.costing.aggregate_flow_costs["hydrogen_peroxide"])

        capital_recovery_rate = value(m.fs.treatment.costing.capital_recovery_factor)
        flow_vol = value(pyunits.convert(m.fs.treatment.product.properties[0].flow_vol_phase["Liq"],
                                         to_units=pyunits.m**3/pyunits.year))

        LCOWs.append(lcow)
        chem_capexs.append(chem_capex*capital_recovery_rate/flow_vol)
        ec_capexs.append(ec_capex*capital_recovery_rate/flow_vol)
        filt_capexs.append(filt_capex*capital_recovery_rate/flow_vol)
        fo_capexs.append(fo_capex*capital_recovery_rate/flow_vol)
        dwi_capexs.append(dwi_capex*capital_recovery_rate/flow_vol)

        chem_opexs.append(chem_opex /flow_vol)
        ec_opexs.append(ec_opex /flow_vol)
        filt_opexs.append(filt_opex/flow_vol)
        fo_opexs.append(fo_opex /flow_vol)
        dwi_opexs.append(dwi_opex /flow_vol)

        elecs.append(elec_cost/flow_vol)
        heats.append(heat_cost/flow_vol)
        alums.append(alum_cost/flow_vol)
        h2o2s.append(h2o2_cost/flow_vol)
    
    print('failed',failed)

#%%
    import matplotlib.pyplot as plt

    plt.stackplot(strong_draw_mass,
                chem_capexs, chem_opexs,
                ec_capexs, ec_opexs,
                filt_capexs, filt_opexs,
                fo_capexs, fo_opexs,
                dwi_capexs, dwi_opexs,
                elecs, heats, alums, h2o2s,
                labels=['Chem add CAPEX', 'Chem add OPEX',
                        'EC CAPEX', 'EC OPEX',
                        'Cart filt CAPEX', 'Cart filt OPEX',
                        'FO CAPEX', 'FO OPEX',
                        'DWI CAPEX', 'DWI OPEX',
                        'Elec', 'Heat','Aluminum','H2O2'
                        ],
                hatch =['', '\\\\',
                        '', '\\\\',
                        '', '\\\\',
                        '', '\\\\',
                        '', '\\\\',
                        '','','','',
                        ],
                colors=['gray','gray',
                        'tomato', 'tomato',
                        'sandybrown','sandybrown',
                        'khaki','khaki',
                        'lightgreen','lightgreen',
                        'gold','indianred','royalblue','darkviolet'
                        ],
                edgecolor='black',
                      )

    plt.rcParams['figure.dpi']=300
    # Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Stacked Line Chart')

    # Show the legend
    plt.legend(loc='upper right', ncol =2,prop={'size': 6})

    plt.ylabel('LCOW ($/m3)')
    plt.xlabel('FO strong draw solution mass fraction')
    plt.title('')
    # Display the chart
    plt.show()

# %%
