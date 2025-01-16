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
from watertap_contrib.reflo.kurby import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1125 * pyunits.kg / pyunits.m**3
rho_water = 995 * pyunits.kg / pyunits.m**3

solver = get_solver()

__all__ = [
    "build_permian_SOA",
    "set_operating_conditions_SOA",
    "add_treatment_costing",
    "set_SOA_scaling",
    "init_SOA_system",
    "run_permian_SOA",
]


def build_permian_SOA(mvc_recovery=0.5):
    """
    Build Permian state-of-the-art flowsheet
    """

    m = ConcreteModel()
    m.mvc_recovery = mvc_recovery
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.properties = ZO(solute_list=["tds"])

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    # Begin building Treatment Block
    m.fs.treatment = treat = Block()

    treat.feed = Feed(property_package=m.fs.properties)
    treat.product = Product(property_package=m.fs.properties_feed)

    # Add translator blocks
    treat.zo_to_sw_feed = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )

    treat.zo_to_sw_disposal = Translator_ZO_to_SW(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.properties_feed,
    )

    # treat.sw_to_zo = Translator_SW_to_ZO(
    #     inlet_property_package=m.fs.properties_feed,
    #     outlet_property_package=m.fs.properties,
    # )

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
        inlet_list=["zo_mixer", "mvc_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.equality,
    )

    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)

    treat.EC = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.EC)

    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.MVC = FlowsheetBlock(dynamic=False)
    build_mvc(m, treat.MVC)

    treat.DWI = FlowsheetBlock(dynamic=False)
    build_dwi(m, treat.DWI, m.fs.properties_feed)

    # BUILD PRODUCT STREAM
    # feed > chem_addition > EC > cart_filt > ZO_to_SW_translator > MVC > product
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

    treat.cart_filt_translated_to_mvc = Arc(
        source=treat.zo_to_sw_feed.outlet, destination=treat.MVC.feed.inlet
    )

    treat.mvc_to_product = Arc(
        source=treat.MVC.product.outlet, destination=treat.product.inlet
    )

    # BUILD DISPOSAL STREAM
    #        EC > ZO_mixer > ZO_to_SW_translator > disposal_mixer > disposal_mixer > DWI
    # cart_filt > ZO_mixer
    #                                        MVC > disposal_mixer

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

    treat.mvc_disposal_to_translator = Arc(
        source=treat.MVC.disposal.outlet,
        destination=treat.disposal_SW_mixer.mvc_disposal,
    )

    treat.disposal_SW_mixer_to_dwi = Arc(
        source=treat.disposal_SW_mixer.outlet, destination=treat.DWI.feed.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    treat.recovery = Var(
        initialize=0.5,
        bounds=(0, 1.0001),
        units=pyunits.dimensionless,
        doc="Overall system recovery",
    )

    # m.fs.treatment.eq_recovery = Constraint(
    #     expr=m.fs.treatment.recovery * m.fs.treatment.feed.properties[0].flow_vol
    #     == m.fs.treatment.product.properties[0].flow_vol_phase["Liq"]
    # )

    add_treatment_costing(m)

    return m


def set_operating_conditions_SOA(m, Qin=5, tds=130, **kwargs):

    global flow_mass_water, flow_mass_tds, flow_in
    m.fs.properties.dens_mass_default = rho

    Qin = Qin * pyunits.Mgallons / pyunits.day
    m.flow_in = flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho_water, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )

    m.fs.treatment.feed.properties[0].flow_mass_comp["H2O"].fix(flow_mass_water)
    m.fs.treatment.feed.properties[0].flow_mass_comp["tds"].fix(flow_mass_tds)
    m.fs.treatment.feed.properties[0].conc_mass_comp[...]

    set_chem_addition_op_conditions(m, m.fs.treatment.chem_addition, **kwargs)
    set_ec_operating_conditions(m, m.fs.treatment.EC, **kwargs)
    set_cart_filt_op_conditions(m, m.fs.treatment.cart_filt, **kwargs)
    set_mvc_operating_conditions(m, m.fs.treatment.MVC, **kwargs)


def add_treatment_costing(m):

    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    add_chem_addition_costing(
        m, m.fs.treatment.chem_addition, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_ec_costing(m, m.fs.treatment.EC, flowsheet_costing_block=m.fs.treatment.costing)
    add_cartridge_filtration_costing(
        m, m.fs.treatment.cart_filt, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_mvc_costing(
        m, m.fs.treatment.MVC, flowsheet_costing_block=m.fs.treatment.costing
    )
    add_dwi_costing(
        m, m.fs.treatment.DWI, flowsheet_costing_block=m.fs.treatment.costing
    )

    m.fs.treatment.costing.cost_process()


def set_SOA_scaling(m, **kwargs):

    set_permian_pretreatment_scaling(m, calculate_m_scaling_factors=False)

    set_mvc_scaling(m, m.fs.treatment.MVC, calc_blk_scaling_factors=True)

    # SW DISPOSAL MIXER
    # MVC disposal inlet
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.mvc_disposal_state[0].flow_mass_phase_comp[
            "Liq", "H2O"
        ],
        1e-2,
    )
    set_scaling_factor(
        m.fs.treatment.disposal_SW_mixer.mvc_disposal_state[0].flow_mass_phase_comp[
            "Liq", "TDS"
        ],
        0.1,
    )

    calculate_scaling_factors(m)


def init_SOA_system(m, mvc_inlet_temp=50, **kwargs):

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

    treat.zo_to_sw_disposal.outlet.temperature[0].fix(293.15)
    treat.zo_to_sw_disposal.outlet.pressure[0].fix(101325)
    treat.zo_to_sw_disposal.initialize()

    treat.zo_to_sw_feed.properties_out[0].temperature.fix(273.15 + mvc_inlet_temp)
    # treat.zo_to_sw_feed.properties_out[0].temperature.set_value(273.15 + mvc_inlet_temp)  # K
    treat.zo_to_sw_feed.properties_out[0].pressure.fix(101325)
    treat.zo_to_sw_feed.initialize()

    propagate_state(treat.cart_filt_translated_to_mvc)
    init_mvc(m, treat.MVC, feed_props=treat.zo_to_sw_feed.properties_out[0])
    propagate_state(treat.mvc_to_product)
    propagate_state(treat.mvc_disposal_to_translator)

    propagate_state(treat.disposal_ZO_mix_translated_to_disposal_SW_mixer)
    treat.disposal_SW_mixer.initialize()
    treat.disposal_SW_mixer.mixed_state[0].pressure.fix(101325)
    # treat.disposal_SW_mixer.mixed_state[0].temperature.fix()

    propagate_state(treat.disposal_SW_mixer_to_dwi)

    init_dwi(m, treat.DWI)

    treat.product.initialize()
    treat.costing.initialize()
    assert False


def run_permian_SOA(recovery=0.5, **kwargs):
    """
    Run Permian state-of-the-art case study
    """

    m = build_permian_SOA()
    treat = m.fs.treatment
    mvc = treat.MVC
    mvc.feed.properties[0].conc_mass_phase_comp

    set_operating_conditions_SOA(m, **kwargs)
    set_SOA_scaling(m)

    # treat.feed.properties[0].flow_vol
    # treat.feed.properties[0].conc_mass_comp
    treat.chem_addition.unit.chemical_flow_vol[0] = 1e-5
    # treat.feed.properties[0].conc_mass_comp.display()
    # assert False

    # m.fs.treatment.costing.initialize()
    # check_jac(m)
    # assert False

    init_SOA_system(m)

    # check_jac(m)
    # mvc.evaporator.display()
    # assert False
    # check_vars(m.fs.treatment.MVC.hx_brine, skip_list=["properties_feed", "properties_vapor", "_ref"])
    # assert False

    flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]
    treat.product.properties[0].flow_vol_phase["Liq"] = value(recovery * m.flow_in)

    treat.costing.add_LCOW(flow_vol)
    treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
    treat.costing.initialize()
    # results = solver.solve(m)
    # assert_optimal_termination(results)

    # First, minimize external heating
    # mvc.external_heating_obj = Objective(expr=mvc.external_heating)
    try:
        print(f"DOF = {degrees_of_freedom(m)}")
        results = solver.solve(mvc, tee=True)
        # results = solver.solve(m, tee=True)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
        print_variables_close_to_bounds(m)
        print("SOLVE FAILED")
        return m
    print(f"DOF = {degrees_of_freedom(m)}")

    # Next, minimize LCOW for MVC design
    # mvc.del_component(mvc.external_heating_obj)
    # treat.LCOW_obj = Objective(expr=treat.costing.LCOW)

    # treat.costing.electricity_cost.fix(0.07)

    try:
        mvc.evaporator.area.unfix()
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
        print_variables_close_to_bounds(m)
        print("SOLVE FAILED")
        return m
    return m
    # Optimize design
    mvc.external_heating.fix(0)  # Remove external heating
    mvc.evaporator.area.unfix()
    mvc.evaporator.outlet_brine.temperature[0].unfix()
    mvc.compressor.pressure_ratio.unfix()
    mvc.hx_distillate.area.unfix()
    mvc.hx_brine.area.unfix()
    # treat.zo_to_sw_feed.properties_out[0].temperature.unfix()

    # Add flowsheet level recovery
    treat.eq_recovery = Constraint(
        expr=treat.recovery * treat.feed.properties[0].flow_vol
        == treat.product.properties[0].flow_vol_phase["Liq"]
    )
    treat.recovery.fix(recovery)

    # Unfix MVC recovery
    mvc.recovery.unfix()

    try:
        results = solver.solve(m)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
        print_variables_close_to_bounds(m)
        print("SOLVE FAILED")
        return m
    print(f"DOF = {degrees_of_freedom(m)}")

    mvc.evaporator.area.fix()
    # mvc.evaporator.outlet_brine.temperature[0].fix()
    mvc.compressor.pressure_ratio.fix()
    mvc.hx_distillate.area.fix()
    # mvc.hx_brine.area.fix()
    mvc.recovery.fix()  # effectively fixes split_fraction on separator

    mvc.compressor.control_volume.properties_out[0.0].temperature.setub(500)

    # try:
    print(f"DOF = {degrees_of_freedom(m)}")
    # mvc.feed.properties[0].conc_mass_phase_comp
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    mvc.hx_brine.area.fix()
    print(f"DOF = {degrees_of_freedom(m)}")
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    return m


# def run_permian_SOA(recovery=0.5):
#     """
#     Run Permian state-of-the-art case study
#     """

#     m = build_permian_SOA()
#     treat = m.fs.treatment
#     mvc = treat.MVC

#     set_operating_conditions_SOA(m)
#     set_SOA_scaling(m)

#     treat.feed.properties[0].flow_vol

#     init_SOA_system(m)

#     flow_vol = treat.product.properties[0].flow_vol_phase["Liq"]
#     treat.product.properties[0].flow_vol_phase["Liq"] = value(recovery * flow_in)

#     treat.costing.add_LCOW(flow_vol)
#     treat.costing.add_specific_energy_consumption(flow_vol, name="SEC")
#     treat.costing.initialize()
#     results = solver.solve(m)
#     assert_optimal_termination(results)
#     # mvc.external_heating.display()
#     # assert False

#     # First, minimize external heating
#     # mvc.external_heating_obj = Objective(expr=mvc.external_heating)

#     # try:
#     # results = solver.solve(m)
#     # assert_optimal_termination(results)
#     # except:
#     #     print_infeasible_constraints(m)
#     #     print_variables_close_to_bounds(m)
#     #     print("SOLVE FAILED")
#     print(f"DOF = {degrees_of_freedom(m)}")
#         # return m
#     # assert False

#     # Next, minimize LCOW for MVC design
#     # mvc.del_component(mvc.external_heating_obj)
#     treat.LCOW_obj = Objective(expr=treat.costing.LCOW)

#     treat.costing.electricity_cost.fix(0.07)

#     try:
#         results = solver.solve(m)
#         assert_optimal_termination(results)
#     except:
#         print_infeasible_constraints(m)
#         print_variables_close_to_bounds(m)
#         print("SOLVE FAILED")
#         return m

#     # Optimize design
#     mvc.external_heating.fix(0)  # Remove external heating
#     mvc.evaporator.area.unfix()
#     mvc.evaporator.outlet_brine.temperature[0].unfix()
#     mvc.compressor.pressure_ratio.unfix()
#     mvc.hx_distillate.area.unfix()
#     mvc.hx_brine.area.unfix()
#     # treat.zo_to_sw_feed.properties_out[0].temperature.unfix()

#     # Add flowsheet level recovery
#     treat.eq_recovery = Constraint(
#         expr=treat.recovery * treat.feed.properties[0].flow_vol
#         == treat.product.properties[0].flow_vol_phase["Liq"]
#     )
#     treat.recovery.fix(recovery)

#     # Unfix MVC recovery
#     mvc.recovery.unfix()

#     results = solver.solve(m)
#     assert_optimal_termination(results)
#     print(f"DOF = {degrees_of_freedom(m)}")

#     mvc.evaporator.area.fix()
#     # mvc.evaporator.outlet_brine.temperature[0].fix()
#     mvc.compressor.pressure_ratio.fix()
#     mvc.hx_distillate.area.fix()
#     mvc.hx_brine.area.fix()
#     mvc.recovery.fix()  # effectively fixes split_fraction on separator

#     # TODO: This results in 1 DOF
#     # Possible additional DOF are
#     # EC conductivity
#     # temperature/pressure on dispoal mixer
#     # mvc.pump_brine.control_volume.deltaP[0].fix()
#     # treat.zo_to_sw_feed.properties_out[0].temperature.fix()
#     # treat.EC.unit.conductivity.fix()

#     mvc.compressor.control_volume.properties_out[0.0].temperature.setub(500)

#     # try:
#     print(f"DOF = {degrees_of_freedom(m)}")
#     mvc.feed.properties[0].conc_mass_phase_comp
#     results = solver.solve(m)
#     assert_optimal_termination(results)
#     #     print(f"termination {results.solver.termination_condition}")
#     #     print(f"DOF = {degrees_of_freedom(m)}")
#     # except:
#     #     print_infeasible_constraints(m)
#     #     print_variables_close_to_bounds(m)
#     #     print("SOLVE FAILED")

#     return m


if __name__ == "__main__":

    m = run_permian_SOA()
    # treat = m.fs.treatment
    # mvc = treat.MVC
    # mvc.evaporator.outlet_brine.display()
    # # treat.disposal_SW_mixer.display()
    # # print(m.fs.treatment.costing.LCOW(), m.fs.treatment.EC.unit.conductivity())
    # # treat.feed.conc_mass_comp.display()
    # mvc.feed.display()
