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
    "build_permian_st1_md",
    "set_operating_conditions_st1_md",
    "add_treatment_costing_st1_md",
    "set_permian_pretreatment_scaling_st1_md",
    "init_system_st1_md",
    "run_permian_st1_md",
]

# TODO:
# Add back CST
# Update membrane type and MD recovery


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


def build_permian_st1_md(Qin=5, Q_md=0.22478, Cin=118, water_recovery=0.2, rho=None):
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

    m.water_recovery = water_recovery

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

    # Build energy block
    m.fs.energy = energy = Block()
    m.fs.energy.cst = FlowsheetBlock()
    build_cst(m.fs.energy.cst)

    # Add treatment costing 
    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    m.fs.energy.costing = EnergyCosting()

    return m


def set_operating_conditions_st1_md(m, rho, Qin=5, tds=130, **kwargs):

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

    set_cst_op_conditions(m.fs.energy.cst,hours_storage=24)


def set_permian_pretreatment_scaling_st1_md(
    m, calclate_m_scaling_factors=False, **kwargs
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

    # ZO to SW disposal translator
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_disposal.properties_in[0].flow_mass_comp["H2O"],
        1,
    )
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_disposal.properties_in[0].flow_mass_comp["tds"],
        1,
    )
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_disposal.properties_out[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1,
    )
    set_scaling_factor(
        m.fs.treatment.zo_to_sw_disposal.properties_out[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        1,
    )

    # ZO DISPOSAL MIXER
    # CF inlet
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.cart_filt_disposal_state[0].flow_mass_comp[
            "H2O"
        ],
        100,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.cart_filt_disposal_state[0].flow_mass_comp[
            "tds"
        ],
        1e8,
    )

    # EC inlet
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.ec_disposal_state[0].flow_mass_comp["H2O"],
        1,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.ec_disposal_state[0].flow_mass_comp["tds"],
        1,
    )

    # mixed state
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.mixed_state[0].flow_mass_comp["H2O"],
        1,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_ZO_mixer.mixed_state[0].flow_mass_comp["tds"],
        1,
    )
    # SW DISPOSAL MIXER
    # ZO mixer inlet
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.zo_mixer_state[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        100,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.zo_mixer_state[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        10,
    )

    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.md_disposal_state[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-3,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.md_disposal_state[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        1e-2,
    )

    # mixed state outlet
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.mixed_state[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-1,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.mixed_state[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        1e-2,
    )

    # DWI
    set_scaling_factor(
        m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-3,
    )
    set_scaling_factor(
        m.fs.treatment.DWI.unit.properties[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        1e-1,
    )

    if calclate_m_scaling_factors:
        print("calclate_m_scaling_factors\n\n\n")
        calculate_scaling_factors(m)


def init_system_st1_md(m, **kwargs):

    treat = m.fs.treatment

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

    init_cst(m.fs.energy.cst)



def add_costing_st1_md(m, heat_price=0.018, electricity_price=0.0626):
   
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
    m.fs.treatment.costing.add_annual_water_production(
        m.fs.treatment.product.properties[0].flow_vol
    )
    m.fs.treatment.costing.add_LCOW(
        m.fs.treatment.product.properties[0].flow_vol
    )

    # Add energy costing

    add_cst_costing(m.fs.energy.cst, m.fs.energy.costing)

    m.fs.energy.costing.cost_process()
    m.fs.energy.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.energy.costing.add_LCOH()

    # Add system costing
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()

    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOH()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")
    
    m.fs.energy.costing.initialize()
    m.fs.treatment.costing.initialize()
    m.fs.costing.initialize()

    print("\n--------- INITIALIZING SYSTEM COSTING COMPLETE---------\n")


def run_permian_st1_md(Qin=5, tds=130, water_recovery = 0.3, **kwargs):
    """
    Run Permian pretreatment flowsheet
    """
    rho = get_stream_density(Qin, tds)

    m_pretreatment = build_and_run_permian_pretreatment(Qin=5)

    print(
        f"Pretreatment Product Flow: {pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)():.4f} m3/s"
    )

    print(
        f"Pretreatment Product Flow: {pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].conc_mass_phase_comp['Liq', 'TDS'],to_units=pyunits.g / pyunits.L,)():.4f} g/L"
    )

    md_flow = pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].flow_vol_phase['Liq'],to_units=pyunits.m**3 / pyunits.s,)
    md_conc = pyunits.convert(m_pretreatment.fs.treatment.product.properties[0].conc_mass_phase_comp['Liq', 'TDS'],to_units=pyunits.g / pyunits.L,)

    # build_permian_st1_md(Qin=5, Q_md=0.22478, Cin=118, water_recovery=0.2, rho=None)

    m = build_permian_st1_md(Q_md=md_flow(), Cin=md_conc(), water_recovery=water_recovery, rho=rho)
    treat = m.fs.treatment

    # set_operating_conditions_st1_md(m, rho, Qin=5, tds=130, **kwargs)
    set_operating_conditions_st1_md(m, rho, Qin, tds)
    set_permian_pretreatment_scaling_st1_md(
        m, calclate_m_scaling_factors=True
    )  # Doesn't solve without this even before costing

    treat.feed.properties[0].flow_vol

    init_system_st1_md(m)
    print(f"DOF = {degrees_of_freedom(m)}")

    # Unfix CST heat
    m.fs.energy.cst.unit.heat_load.unfix() 

    results = solver.solve(m)
    print_infeasible_constraints(m)
    assert_optimal_termination(results)

    print("\n--------- Before costing solve Completed ---------\n")
    report_MD(m, treat.md)

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
            1e6,
        )

    print("\n--------- CST Inputs Completed ---------\n")
    print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
    print('CST Heat:', value(m.fs.energy.cst.unit.heat))
    print("\n")

    # Add costing  
    add_costing_st1_md(m)

    results = solver.solve(m)

    m.fs.energy.cst.unit.heat_load.fix()
    # m.fs.lcot_objective = Objective(expr=m.fs.costing.LCOT)  

    try:
        results = solver.solve(m)
        print_infeasible_constraints(m)
    except ValueError:
        print_infeasible_constraints(m)
    assert_optimal_termination(results)
    print("\n--------- After costing solve Completed ---------\n")

    print('CST Heat load:', value(m.fs.energy.cst.unit.heat_load))
    print('CST Heat:', value(m.fs.energy.cst.unit.heat))

    return m


if __name__ == "__main__":

    tds = 200
    Qin = 9
    water_recovery = 0.5

    m = run_permian_st1_md(Qin=Qin, tds=tds, water_recovery = water_recovery)
    treat = m.fs.treatment
    report_MD(m, treat.md)
    print(f"DOF = {degrees_of_freedom(m)}")

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