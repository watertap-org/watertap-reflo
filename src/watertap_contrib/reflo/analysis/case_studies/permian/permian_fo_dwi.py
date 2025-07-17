import pathlib
import pandas as pd

# Import Pyomo packages
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

# Import IDAES packages
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialBalanceType
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
    get_scaling_factor,
    extreme_jacobian_columns,
    extreme_jacobian_rows,
    badly_scaled_var_generator,
    unscaled_variables_generator,
    extreme_jacobian_entries,
)
from idaes.models.unit_models import (
    Product,
    Feed,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *
from idaes.core.util.initialization import propagate_state

# Import WaterTAP packages
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)

# Import WaterTAP-REFLO packages
from watertap_contrib.reflo.analysis.example_flowsheets.fo_trevi_flowsheet import (
    build_fo_trevi_flowsheet,
    fix_dof_and_initialize,
)
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import (
    FODrawSolutionParameterBlock,
)

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
save_dir = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/sweep_results"
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


skips = [
    "diffus_phase",
    "diffus_param",
    "dens_mass_param",
    "dh_vap_w_param",
    "cp_phase_param",
    "pressure_sat_param_psatw",
    "enth_mass_param",
    "osm_coeff_param",
    "visc_d_param",
    "therm_cond_phase_param",
    "pressure_sat_param",
    "bpe_",
    "TIC",
    "TPEC",
    "blocks[",
    "yearly_heat_production",
    "yearly_electricity_production",
    "cp_param_NaCl_liq",
    "_translator",
    "permeate_side",
    "properties_interface",
    "material_flow_dx",
    "._flow_terms",
    "pressure_dx",
    "MCAS_properties",
    "cp_param_NaCl_solid",
    "cp_vap_param",
    "temp_sat_solvent",
    "cp_mass_phase",
    ".properties_NaCl",
    ".properties_draw",
    ".properties_sw",
    ".vapor_properties",
    ".properties.",
    ".seawater_properties."
]


def build_permian_FO_DWI(permian_fo_config):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.treatment = treat = Block()

    m.fs.properties = ZO(solute_list=["tds"])
    m.fs.properties_sw = SeawaterParameterBlock()
    m.fs.properties_draw = FODrawSolutionParameterBlock()

    treat.feed = Feed(property_package=m.fs.properties)
    treat.product = Product(property_package=m.fs.properties_sw)

    # Add translator blocks
    treat.zo_to_sw_feed = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_sw,
    )
    treat.zo_to_sw_ec_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_sw,
    )
    treat.zo_to_sw_cart_filt_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_sw,
    )
    treat.draw_to_sw = Translator_Draw_to_SW(
        inlet_property_package=m.fs.properties_draw,
        outlet_property_package=m.fs.properties_sw,
    )

    # Add components
    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.ec = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.ec)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.FO = build_fo_trevi_flowsheet(
        feed_vol_flow=permian_fo_config["feed_vol_flow"],
        feed_TDS_mass=permian_fo_config["feed_TDS_mass"],
        recovery_ratio=permian_fo_config["recovery_ratio"],
        RO_recovery_ratio=permian_fo_config["RO_recovery_ratio"],
        NF_recovery_ratio=permian_fo_config["NF_recovery_ratio"],
        feed_temperature=permian_fo_config["feed_temperature"],
        strong_draw_temp=permian_fo_config["strong_draw_temp"],
        strong_draw_mass=permian_fo_config["strong_draw_mass_frac"],
    )

    treat.disposal_SW_mixer = Mixer(
        property_package=m.fs.properties_sw,
        num_inlets=3,
        inlet_list=["ec_disposal", "cart_filt_disposal", "fo_disposal"],
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.DWI = FlowsheetBlock(dynamic=False)
    build_dwi(m, treat.DWI, prop_package=m.fs.properties_sw)

    # BUILD PRODUCT STREAM
    # feed (1)> chem_addition (2)> EC (3)> cart_filt
    #      (4)> ZO_to_SW_translator (5)> FO (6)> Draw_to_SW_translator (7)> product

    treat.feed_to_chem_addition = Arc(
        source=treat.feed.outlet, destination=treat.chem_addition.feed.inlet
    )  # (1)
    treat.chem_addition_to_ec = Arc(
        source=treat.chem_addition.product.outlet, destination=treat.ec.feed.inlet
    )  # (2)
    treat.ec_to_cart_filt = Arc(
        source=treat.ec.product.outlet, destination=treat.cart_filt.feed.inlet
    )  # (3)
    treat.cart_filt_to_translator = Arc(
        source=treat.cart_filt.product.outlet, destination=treat.zo_to_sw_feed.inlet
    )  # (4)
    treat.cart_filt_translated_to_fo = Arc(
        source=treat.zo_to_sw_feed.outlet, destination=treat.FO.fs.fo.feed
    )  # (5)
    treat.fo_to_translator = Arc(
        source=treat.FO.fs.S2.fresh_water, destination=treat.draw_to_sw.inlet
    )  # (6)
    treat.fo_translator_to_product = Arc(
        source=treat.draw_to_sw.outlet, destination=treat.product.inlet
    )  # (7)

    # BUILD DISPOSAL STREAM
    #        EC (1)> ZO_to_SW_translator (3)> disposal_mixer (6)> DWI
    # cart_filt (2)> ZO_to_SW_translator (4)> disposal_mixer
    #                                FO  (5)> disposal_mixer

    treat.ec_disposal_to_translator = Arc(
        source=treat.ec.disposal.outlet, destination=treat.zo_to_sw_ec_disposal.inlet
    )  # (1)
    treat.cart_filt_disposal_to_translator = Arc(
        source=treat.cart_filt.disposal.outlet,
        destination=treat.zo_to_sw_cart_filt_disposal.inlet,
    )  # (2)
    treat.ec_disposal_translator_to_SW_mixer = Arc(
        source=treat.zo_to_sw_ec_disposal.outlet,
        destination=treat.disposal_SW_mixer.ec_disposal,
    )  # (3)
    treat.cart_filt_disposal_translator_to_SW_mixer = Arc(
        source=treat.zo_to_sw_cart_filt_disposal.outlet,
        destination=treat.disposal_SW_mixer.cart_filt_disposal,
    )  # (4)
    treat.fo_disposal_translator_to_SW_mixer = Arc(
        source=treat.FO.fs.fo.brine,
        destination=treat.disposal_SW_mixer.fo_disposal,
    )  # (5)
    treat.SW_mixer_to_DWI = Arc(
        source=treat.disposal_SW_mixer.outlet,
        destination=treat.DWI.feed.inlet,
    )  # (6)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_operating_conditions(m, operating_condition, **kwargs):
    Qin, tds = operating_condition["feed_vol_flow"], operating_condition["feed_tds"]

    global flow_mass_water, flow_mass_tds, flow_in

    # Calculate feed flow density
    rho = get_stream_density(tds=tds)
    m.fs.properties.dens_mass_default = rho

    # Calculate feed mass flow rate
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

    # Set up components
    set_chem_addition_op_conditions(m, m.fs.treatment.chem_addition, **kwargs)
    set_ec_operating_conditions(m, m.fs.treatment.ec, **kwargs)
    set_cart_filt_op_conditions(m, m.fs.treatment.cart_filt)


def set_permian_scaling(m, **kwargs):
    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        # 1 / value(flow_mass_water),
        1e-2,
        index=("H2O"),
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_comp",
        # 1 / value(flow_mass_tds),
        1e-1,
        index=("tds"),
    )

    m.fs.properties_sw.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(flow_mass_water),
        1e-1,
        index=("Liq", "H2O"),
    )
    m.fs.properties_sw.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(flow_mass_tds),
        1e-2,
        index=("Liq", "TDS"),
    )

    m.fs.properties_draw.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(flow_mass_water),
        1e-1,
        index=("Liq", "H2O"),
    )

    m.fs.properties_draw.set_default_scaling(
        "flow_mass_phase_comp",
        # 1 / value(flow_mass_water),
        1e-1,
        index=("Liq", "DrawSolution"),
    )

    # set_scaling_factor(
    #     m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
    #         "Liq", "H2O"
    #     ],
    #     1e-5,
    # )
    # set_scaling_factor(
    #     m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
    #         "Liq", "TDS"
    #     ],
    #     1e-1,
    # )

    set_scaling_factor(
        m.fs.treatment.zo_to_sw_cart_filt_disposal.properties_in[0].flow_mass_comp[
            "H2O"
        ],
        1e3,
    )
    set_scaling_factor(
        m.fs.treatment.cart_filt.disposal.properties[0].flow_mass_comp["H2O"], 1e3
    )
    set_scaling_factor(
        m.fs.treatment.cart_filt.unit.properties_byproduct[0].flow_mass_comp["H2O"], 1e3
    )

    set_chem_addition_scaling(
        m, m.fs.treatment.chem_addition, calc_blk_scaling_factors=True
    )

    set_cart_filt_scaling(m, m.fs.treatment.cart_filt, calc_blk_scaling_factors=True)

    set_ec_scaling(m, m.fs.treatment.ec, calc_blk_scaling_factors=True)

    calculate_scaling_factors(m)

    # set_scaling_factor(m.fs.treatment.disposal_SW_mixer.cart_filt_disposal_state[0.0].flow_mass_phase_comp["Liq","H2O"], 1e3)


def init_system(m, permian_fo_config):
    treat = m.fs.treatment

    treat.feed.initialize()
    propagate_state(arc=treat.feed_to_chem_addition)

    init_chem_addition(m, treat.chem_addition)
    propagate_state(arc=treat.chem_addition_to_ec)

    init_ec(m, treat.ec)
    propagate_state(arc=treat.ec_to_cart_filt)
    propagate_state(arc=treat.ec_disposal_to_translator)
    treat.zo_to_sw_ec_disposal.initialize()
    treat.zo_to_sw_ec_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_ec_disposal.outlet.pressure[0].fix(101325)
    propagate_state(arc=treat.ec_disposal_translator_to_SW_mixer)

    init_cart_filt(m, treat.cart_filt)
    propagate_state(arc=treat.cart_filt_to_translator)
    treat.zo_to_sw_feed.initialize()
    propagate_state(arc=treat.cart_filt_disposal_to_translator)
    treat.zo_to_sw_cart_filt_disposal.initialize()
    treat.zo_to_sw_cart_filt_disposal.outlet.temperature[0].fix(25 + 273.15)
    treat.zo_to_sw_cart_filt_disposal.outlet.pressure[0].fix(101325)
    propagate_state(arc=treat.cart_filt_disposal_translator_to_SW_mixer)

    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
    propagate_state(arc=treat.cart_filt_translated_to_fo)
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].fix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].fix()

    fix_dof_and_initialize(
        treat.FO,
        strong_draw_mass_frac=permian_fo_config["strong_draw_mass_frac"],
        product_draw_mass_frac=permian_fo_config["product_draw_mass_frac"],
        RO_recovery_ratio=permian_fo_config["RO_recovery_ratio"],
        NF_recovery_ratio=permian_fo_config["NF_recovery_ratio"],
    )
    # unfix FO fs and set operating point
    treat.FO.fs.HX1.area.unfix()
    treat.FO.fs.HX2.area.unfix()
    treat.FO.fs.HX1.weak_draw_outlet.temperature.fix(
        permian_fo_config["HX1_cold_out_temp"]
    )
    treat.FO.fs.HX1.product_water_outlet.temperature.fix(
        permian_fo_config["HX1_hot_out_temp"]
    )

    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    treat.FO.fs.fo.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    propagate_state(arc=treat.fo_to_translator)
    treat.draw_to_sw.initialize()

    propagate_state(arc=treat.fo_translator_to_product)
    treat.product.pressure[0].fix(101325)
    treat.product.initialize()

    treat.FO.fs.fo.brine.pressure[0].fix(101325)
    propagate_state(arc=treat.fo_disposal_translator_to_SW_mixer)

    treat.disposal_SW_mixer.initialize()
    propagate_state(arc=treat.SW_mixer_to_DWI)

    treat.DWI.unit.properties[0].temperature.fix()
    treat.DWI.unit.properties[0].pressure.fix()
    init_dwi(m, treat.DWI)


def add_treatment_costing(m):
    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )

    add_ec_costing(m, m.fs.treatment.ec, flowsheet_costing_block=m.fs.treatment.costing)

    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.FO.fs.fo.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing
    )

    add_dwi_costing(
        m, m.fs.treatment.DWI, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.costing.cost_process()


def run_permian_FO_DWI(
    operating_condition,
    permian_fo_config,
):
    m = build_permian_FO_DWI(permian_fo_config)
    treat = m.fs.treatment

    set_operating_conditions(m, operating_condition)
    set_permian_scaling(m)

    treat.feed.properties[0].flow_vol

    init_system(m, permian_fo_config)

    print("DOF after init: ", degrees_of_freedom(m))
    results = solver.solve(m)
    assert_optimal_termination(results)

    add_treatment_costing(m)

    flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]

    treat.costing.base_currency = pyunits.USD_2023

    electricity_price = 0.0434618999  # 2018 price $/kWh
    heat_price = 0.0166  # 2023 price $/kWh

    treat.costing.electricity_cost.fix(electricity_price)
    treat.costing.heat_cost.fix(heat_price)

    treat.costing.add_LCOW(flow_vol)
    flow_in = treat.feed.properties[0].flow_vol

    treat.costing.add_specific_energy_consumption(flow_in, name="SEC_in")
    m.fs.treatment.costing.SEC_th_in = Expression(
        expr=pyunits.convert(
            m.fs.treatment.costing.aggregate_flow_heat / flow_in,
            to_units=pyunits.kilowatt * pyunits.hr * pyunits.m**-3,
        )
    )
    m.fs.treatment.costing._add_flow_component_breakdowns(
        "heat", "SEC_th_in", flow_in, period=pyunits.hr
    )
    treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
    treat.costing.add_specific_thermal_energy_consumption(flow_vol, name="SEC_th")
    treat.costing._add_flow_component_breakdowns(
        "heat", "SEC_th", flow_vol, period=pyunits.hr
    )
    treat.costing.initialize()

    print("dof before solving", degrees_of_freedom(m))
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


##################################################################
##################################################################
###################### SWEEP FUNCTIONS BELOW #####################
##################################################################
##################################################################


def run_recovery_ratio_sweep():

    permian_fo_config = {
        "feed_vol_flow": 0.22,  # initial value for fo model setup
        "feed_TDS_mass": 0.119,  # mass fraction, 0.119 is about 130 g/L, 0.092 for 100 g/L, 0.19 for 200 g/L
        "recovery_ratio": 0.485,  # To get 250 g/L brine, select 0.485 for 130g/L, 0.612 for 100g/L, 0.165 for 200g/L
        "RO_recovery_ratio": 1,  # RO recovery ratio
        "NF_recovery_ratio": 0.8,  # Nanofiltration recovery ratio
        "feed_temperature": 25,
        "strong_draw_temp": 25,  # Strong draw solution inlet temperature (C)
        "strong_draw_mass_frac": 0.9,  # Strong draw solution mass fraction
        "product_draw_mass_frac": 0.01,  # FO product draw solution mass fraction
        "HX1_cold_out_temp": 78 + 273.15,  # HX1 coldside outlet temperature
        "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
    }

    operating_condition = {"feed_vol_flow": 5, "feed_tds": 130}  # MGD  # g/L

    m = run_permian_FO_DWI(
        operating_condition,
        permian_fo_config,
    )

    rd = build_results_dict(m, skips=skips)
    rd["fo_recovery_ratio"] = []
    rd["fo_thermal_energy_flow"] = []
    rd["brine_conc"] = []

    recovery_ratios = [
        0.35,
        0.36,
        0.37,
        0.38,
        0.39,
        0.40,
        0.41,
        0.42,
        0.43,
        0.44,
        0.45,
        0.46,
        0.47,
        0.48,
        0.485,
        0.49,
        0.5,
        0.51,
        0.52,
        0.53,
        0.54,
        0.55,
        0.56,
        0.57,
        0.58,
    ]

    for rr in recovery_ratios:
        permian_fo_config["recovery_ratio"] = rr

        try:
            m = run_permian_FO_DWI(
                operating_condition,
                permian_fo_config,
            )
            rd = results_dict_append(m, rd)
            rd["fo_recovery_ratio"].append(rr * 100)
            rd["fo_thermal_energy_flow"].append(
                value(m.fs.treatment.FO.fs.fo.costing.thermal_energy_flow)
            )
            rd["brine_conc"].append(
                value(
                    m.fs.treatment.FO.fs.fo.brine_props[0].conc_mass_phase_comp[
                        "Liq", "TDS"
                    ]
                )
            )

        except:
            pass

    df = pd.DataFrame.from_dict(rd)
    df.to_csv(f"{save_dir}/permian_RPT2_FO_DWI_no_CST_sweep_recovery_ratio.csv")


def run_flow_tds_sweep():

    permian_fo_config = {
        "feed_vol_flow": 0.22,  # initial value for fo model setup
        "feed_TDS_mass": 0.119,  # mass fraction, 0.119 is about 130 g/L, 0.092 for 100 g/L, 0.19 for 200 g/L
        # "recovery_ratio": 0.485,  # To get 250 g/L brine, select 0.485 for 130g/L, 0.612 for 100g/L, 0.165 for 200g/L
        "recovery_ratio": 0.5,  # To get 250 g/L brine, select 0.485 for 130g/L, 0.612 for 100g/L, 0.165 for 200g/L
        "RO_recovery_ratio": 1,  # RO recovery ratio
        "NF_recovery_ratio": 0.8,  # Nanofiltration recovery ratio
        "feed_temperature": 25,
        "strong_draw_temp": 25,  # Strong draw solution inlet temperature (C)
        "strong_draw_mass_frac": 0.9,  # Strong draw solution mass fraction
        "product_draw_mass_frac": 0.01,  # FO product draw solution mass fraction
        "HX1_cold_out_temp": 78 + 273.15,  # HX1 coldside outlet temperature
        "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
    }

    operating_condition = {
        "feed_vol_flow": 5,  # MGD
        "feed_tds": 130,  # g/L
    }
    m = run_permian_FO_DWI(
        operating_condition,
        permian_fo_config,
    )

    rd = build_results_dict(m, skips=skips)
    rd["tds"] = []
    rd["flow_mgd"] = []
    rd["recovery_ratio"] = []
    rd["feed_TDS_mass"] = []
    rd["fo_thermal_energy_flow"] = []
    rd["brine_conc"] = []

    tds = [100, 130, 200]
    mfs = [0.092, 0.119, 0.19]
    # rrs = [0.612, 0.485, 0.165]
    # NOTE: July 16, 2025
    # rrs that are commented out above were found to result 
    # in brine concentrations that were < 250g/L
    # when using m.fs.treatment.FO.fs.fo.brine_props[0].conc_mass_phase_comp for brine concentration.
    # They were re-evaluated to find appropriate values.
    rrs = [0.6395, 0.5305, 0.275]
    # tds = [100]
    # mfs = [0.092]
    # rrs = [0.612]
    qs = [1, 5, 9]

    for mf, salt, rr in zip(mfs, tds, rrs):
        for q in qs:

            permian_fo_config = {
                "feed_vol_flow": 0.22,  # initial value for fo model setup
                "feed_TDS_mass": mf,  # mass fraction, 0.119 is about 130 g/L, 0.092 for 100 g/L, 0.19 for 200 g/L
                "recovery_ratio": rr,  # To get 250 g/L brine, select 0.485 for 130g/L, 0.612 for 100g/L, 0.165 for 200g/L
                "RO_recovery_ratio": 1,  # RO recovery ratio
                "NF_recovery_ratio": 0.8,  # Nanofiltration recovery ratio
                "feed_temperature": 25,
                "strong_draw_temp": 25,  # Strong draw solution inlet temperature (C)
                "strong_draw_mass_frac": 0.9,  # Strong draw solution mass fraction
                "product_draw_mass_frac": 0.01,  # FO product draw solution mass fraction
                "HX1_cold_out_temp": 78 + 273.15,  # HX1 coldside outlet temperature
                "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
            }

            operating_condition = {
                "feed_vol_flow": q,  # MGD
                "feed_tds": salt,  # g/L
            }

            m = run_permian_FO_DWI(
                operating_condition,
                permian_fo_config,
            )

            rd["tds"].append(salt)
            rd["flow_mgd"].append(q)
            rd["recovery_ratio"].append(rr)
            rd["feed_TDS_mass"].append(mf)
            rd["fo_thermal_energy_flow"].append(
                value(m.fs.treatment.FO.fs.fo.costing.thermal_energy_flow)
            )
            rd["brine_conc"].append(
                value(
                    m.fs.treatment.FO.fs.fo.brine_props[0].conc_mass_phase_comp[
                        "Liq", "TDS"
                    ]
                )
            )
            rd = results_dict_append(m, rd)

    df = pd.DataFrame.from_dict(rd)
    df.to_csv(f"{save_dir}/permian_RPT2_FO_DWI_flow_TDS_sweep.csv")


if __name__ == "__main__":

    run_recovery_ratio_sweep()
    run_flow_tds_sweep()

#     permian_fo_config = {
#         "feed_vol_flow": 0.22,  # initial value for fo model setup
#         "feed_TDS_mass": 0.119,  # mass fraction, 0.119 is about 130 g/L, 0.092 for 100 g/L, 0.19 for 200 g/L
#         "recovery_ratio": 0.485,  # To get 250 g/L brine, select 0.485 for 130g/L, 0.612 for 100g/L, 0.165 for 200g/L
#         "RO_recovery_ratio": 1,  # RO recovery ratio
#         "NF_recovery_ratio": 0.8,  # Nanofiltration recovery ratio
#         "feed_temperature": 25,
#         "strong_draw_temp": 25,  # Strong draw solution inlet temperature (C)
#         "strong_draw_mass_frac": 0.9,  # Strong draw solution mass fraction
#         "product_draw_mass_frac": 0.01,  # FO product draw solution mass fraction
#         "HX1_cold_out_temp": 78 + 273.15,  # HX1 coldside outlet temperature
#         "HX1_hot_out_temp": 32 + 273.15,  # HX1 hotside outlet temperature
#     }

#     operating_condition = {"feed_vol_flow": 5, "feed_tds": 130}  # MGD  # g/L

#     m = run_permian_FO_DWI(
#         operating_condition,
#         permian_fo_config,
#     )
#     print(m.fs.treatment.costing.LCOW())
