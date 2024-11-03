import os
import pathlib
import math
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
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.models.unit_models import Product, Feed, Mixer, StateJunction, Separator
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.util.initialization import *

from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.mvc.components import Evaporator
from watertap.unit_models.mvc.components import Compressor
from watertap.unit_models.mvc.components import Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.analysis.case_studies.permian.components.translator_sw_to_water import (
    Translator_SW_to_Water,
)

__all__ = [
    "build_mvc",
    "set_mvc_operating_conditions",
    "set_mvc_scaling",
    "init_mvc",
    "add_mvc_costing",
]


reflo_dir = pathlib.Path(__file__).resolve().parents[4]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3

_log = idaeslog.getLogger("permian_MVC")


def build_system(single_stage_build=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.product = Product(property_package=m.fs.properties_feed)
    m.fs.disposal = Product(property_package=m.fs.properties_feed)

    m.fs.MVC = mvc = FlowsheetBlock(dynamic=False)

    build_mvc(m, mvc, single_stage_build=single_stage_build)

    m.fs.feed_to_mvc = Arc(source=m.fs.feed.outlet, destination=mvc.feed.inlet)

    m.fs.mvc_to_product = Arc(source=mvc.product.outlet, destination=m.fs.product.inlet)

    m.fs.mvc_to_disposal = Arc(
        source=mvc.disposal.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mvc(m, blk, single_stage_build=True):

    blk.feed = StateJunction(property_package=m.fs.properties_feed)
    blk.product = StateJunction(property_package=m.fs.properties_feed)
    blk.disposal = StateJunction(property_package=m.fs.properties_feed)

    blk.recovery = Var(
        initialize=0.5, bounds=(0, 1), units=pyunits.dimensionless, doc="MVC recovery"
    )

    # Evaporator
    blk.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )
    # Compressor
    blk.compressor = Compressor(property_package=m.fs.properties_vapor)

    # Condenser
    blk.condenser = Condenser(property_package=m.fs.properties_vapor)

    # Translator SW to Water
    blk.tb_sw_to_water = Translator_SW_to_Water(
        inlet_property_package=m.fs.properties_vapor,
        outlet_property_package=m.fs.properties_feed,
    )

    add_external_heating(m, blk)

    if single_stage_build:
        blk.pump_feed = Pump(property_package=m.fs.properties_feed)
        blk.pump_brine = Pump(property_package=m.fs.properties_feed)
        blk.pump_distillate = Pump(property_package=m.fs.properties_feed)

        blk.separator = Separator(
            property_package=m.fs.properties_feed,
            outlet_list=["hx_distillate_cold", "hx_brine_cold"],
            split_basis=SplittingType.totalFlow,
        )

        blk.hx_distillate = HeatExchanger(
            hot_side_name="hot",
            cold_side_name="cold",
            hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
            cold={
                "property_package": m.fs.properties_feed,
                "has_pressure_change": True,
            },
            delta_temperature_callback=delta_temperature_chen_callback,
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )
        # Set lower bound of approach temperatures
        blk.hx_distillate.delta_temperature_in.setlb(0)
        blk.hx_distillate.delta_temperature_out.setlb(0)
        blk.hx_distillate.area.setlb(10)

        blk.hx_brine = HeatExchanger(
            hot_side_name="hot",
            cold_side_name="cold",
            hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
            cold={
                "property_package": m.fs.properties_feed,
                "has_pressure_change": True,
            },
            delta_temperature_callback=delta_temperature_chen_callback,
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )
        # Set lower bound of approach temperatures
        blk.hx_brine.delta_temperature_in.setlb(0)
        blk.hx_brine.delta_temperature_out.setlb(0)
        blk.hx_brine.area.setlb(10)

        blk.mixer_feed = Mixer(
            property_package=m.fs.properties_feed,
            momentum_mixing_type=MomentumMixingType.equality,
            inlet_list=["hx_distillate_cold", "hx_brine_cold"],
        )

        blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump_feed.inlet)
        blk.pump_to_separator = Arc(
            source=blk.pump_feed.outlet, destination=blk.separator.inlet
        )
        blk.sep_dist_cold_to_hx_dist_cold = Arc(
            source=blk.separator.hx_distillate_cold,
            destination=blk.hx_distillate.cold_inlet,
        )
        blk.sep_brine_cold_to_hx_brine_cold = Arc(
            source=blk.separator.hx_brine_cold, destination=blk.hx_brine.cold_inlet
        )
        blk.hx_dist_cold_to_mixer = Arc(
            source=blk.hx_distillate.cold_outlet,
            destination=blk.mixer_feed.hx_distillate_cold,
        )
        blk.hx_brine_cold_to_mixer = Arc(
            source=blk.hx_brine.cold_outlet, destination=blk.mixer_feed.hx_brine_cold
        )
        blk.mixer_feed_to_evaporator = Arc(
            source=blk.mixer_feed.outlet, destination=blk.evaporator.inlet_feed
        )
        blk.evaporator_to_compressor = Arc(
            source=blk.evaporator.outlet_vapor, destination=blk.compressor.inlet
        )
        blk.compressor_to_condenser = Arc(
            source=blk.compressor.outlet, destination=blk.condenser.inlet
        )
        blk.evaporator_to_brine_pump = Arc(
            source=blk.evaporator.outlet_brine, destination=blk.pump_brine.inlet
        )
        ##
        blk.brine_pump_to_hx_brine_hot = Arc(
            source=blk.pump_brine.outlet, destination=blk.hx_brine.hot_inlet
        )
        blk.hx_brine_hot_to_disposal = Arc(
            source=blk.hx_brine.hot_outlet, destination=blk.disposal.inlet
        )
        blk.condenser_to_translator = Arc(
            source=blk.condenser.outlet, destination=blk.tb_sw_to_water.inlet
        )
        blk.translated_to_dist_pump = Arc(
            source=blk.tb_sw_to_water.outlet, destination=blk.pump_distillate.inlet
        )
        ##
        blk.dist_pump_to_hx_dist_hot = Arc(
            source=blk.pump_distillate.outlet, destination=blk.hx_distillate.hot_inlet
        )
        blk.hx_dist_hot_to_product = Arc(
            source=blk.hx_distillate.hot_outlet, destination=blk.product.inlet
        )

        blk.eq_recovery = Constraint(
            expr=blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
            == blk.recovery
            * (
                blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
                + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
            )
        )

        blk.eq_separator_split_frac = Constraint(
            expr=blk.separator.split_fraction[0, "hx_distillate_cold"] == blk.recovery
        )

    else:
        blk.feed_to_evaporator = Arc(
            source=blk.feed.outlet, destination=blk.evaporator.inlet_feed
        )

        blk.evaporator_to_compressor = Arc(
            source=blk.evaporator.outlet_vapor, destination=blk.compressor.inlet
        )

        blk.compressor_to_condenser = Arc(
            source=blk.compressor.outlet, destination=blk.condenser.inlet
        )
        blk.condenser_to_translator = Arc(
            source=blk.condenser.outlet, destination=blk.tb_sw_to_water.inlet
        )
        blk.translated_to_product = Arc(
            source=blk.tb_sw_to_water.outlet, destination=blk.product.inlet
        )

        # blk.condenser_to_product = Arc(
        #     source=blk.condenser.outlet, destination=blk.product.inlet
        # )

        blk.evaporator_to_disposal = Arc(
            source=blk.evaporator.outlet_brine, destination=blk.disposal.inlet
        )

    TransformationFactory("network.expand_arcs").apply_to(m)

    blk.evaporator.connect_to_condenser(blk.condenser)
    # _log.info("MVC flowsheet built")


def set_mvc_operating_conditions(
    m,
    blk,
    recovery=0.5,
    inlet_brine_temp_guess=50,  # degC
    outlet_brine_temp=70,  # degC
    steam_temp_ub=75,
    single_stage_build=True,
):

    # blk.feed.properties[0].temperature.fix(273.15 + 50.52)  # K
    # blk.feed.properties[0].pressure.fix(1e5)  # Pa
    # blk.evaporator.inlet_feed.temperature[0] = 273.15 + inlet_brine_temp_guess  # provide guess
    # blk.evaporator.outlet_brine.temperature[0].fix(273.15 + outlet_brine_temp)
    # blk.evaporator.U.fix(1e3)
    # blk.evaporator.area.fix(400)
    # blk.evaporator.area.fix(9125)

    # blk.compressor.pressure_ratio.fix(2)
    # blk.compressor.efficiency.fix(0.8)

    blk.recovery.fix(recovery)

    if single_stage_build:

        # Feed pump
        blk.pump_feed.efficiency_pump[0].fix(0.8)
        blk.pump_feed.control_volume.deltaP[0].fix(7e3)

        # Separator
        blk.separator.split_fraction[0, "hx_distillate_cold"] = blk.recovery.value

        # Distillate HX
        blk.hx_distillate.overall_heat_transfer_coefficient[0].fix(2e3)
        blk.hx_distillate.area.fix(125)
        blk.hx_distillate.cold.deltaP[0].fix(7e3)
        blk.hx_distillate.hot.deltaP[0].fix(7e3)

        # Brine HX
        blk.hx_brine.overall_heat_transfer_coefficient[0].fix(2e3)
        blk.hx_brine.area.fix(115)
        blk.hx_brine.cold.deltaP[0].fix(7e3)
        blk.hx_brine.hot.deltaP[0].fix(7e3)

        # Evaporator
        blk.evaporator.inlet_feed.temperature[0] = (
            inlet_brine_temp_guess + 273.15
        )  # provide guess
        blk.evaporator.outlet_brine.temperature[0].fix(outlet_brine_temp + 273.15)
        blk.evaporator.U.fix(3e3)  # W/K-m^2
        blk.evaporator.area.setub(1e4)  # m^2

        # Compressor
        blk.compressor.pressure_ratio.fix(1.6)
        blk.compressor.efficiency.fix(0.8)

        # Brine pump
        blk.pump_brine.efficiency_pump[0].fix(0.8)
        blk.pump_brine.control_volume.deltaP[0].fix(4e4)

        # Distillate pump
        blk.pump_distillate.efficiency_pump[0].fix(0.8)
        blk.pump_distillate.control_volume.deltaP[0].fix(4e4)

        # Fix 0 TDS
        blk.tb_sw_to_water.properties_out[0].flow_mass_phase_comp["Liq", "TDS"].fix(
            1e-5
        )

        # Temperature bounds
        blk.evaporator.properties_vapor[0].temperature.setub(steam_temp_ub + 273.15)
        blk.compressor.control_volume.properties_out[0].temperature.setub(450)
    else:
        pass
    # _log.info("MVC operating conditions set")
    # print(f"MVC")
    # print(f"\tblock DOF = {degrees_of_freedom(blk)}\n")
    # print(f"\tevaporator DOF = {degrees_of_freedom(blk.evaporator)}\n")
    # print(f"\tcompressor DOF = {degrees_of_freedom(blk.compressor)}\n")
    # print(f"\tcondenser DOF = {degrees_of_freedom(blk.condenser)}\n")
    print("DOF after setting operating conditions: ", degrees_of_freedom(blk))


def set_mvc_scaling(
    m,
    blk,
    properties_feed=None,
    properties_vapor=None,
    calc_blk_scaling_factors=False,
    single_stage_build=True,
):

    if properties_feed is None:
        properties_feed = m.fs.properties_feed

    if properties_vapor is None:
        properties_vapor = m.fs.properties_vapor

    if single_stage_build:

        set_scaling_factor(blk.pump_feed.control_volume.work, 1e-3)
        set_scaling_factor(blk.pump_brine.control_volume.work, 1e-3)
        set_scaling_factor(blk.pump_distillate.control_volume.work, 1e-3)

        # distillate HX
        set_scaling_factor(blk.hx_distillate.hot.heat, 1e-3)
        set_scaling_factor(blk.hx_distillate.cold.heat, 1e-3)
        set_scaling_factor(blk.hx_distillate.overall_heat_transfer_coefficient, 1e-3)

        set_scaling_factor(blk.hx_distillate.area, 1e-1)
        constraint_scaling_transform(
            blk.hx_distillate.cold_side.pressure_balance[0], 1e-5
        )
        constraint_scaling_transform(
            blk.hx_distillate.hot_side.pressure_balance[0], 1e-5
        )

        # brine HX
        set_scaling_factor(blk.hx_brine.hot.heat, 1e-3)
        set_scaling_factor(blk.hx_brine.cold.heat, 1e-3)
        set_scaling_factor(blk.hx_brine.overall_heat_transfer_coefficient, 1e-3)
        set_scaling_factor(blk.hx_brine.area, 1e-1)
        constraint_scaling_transform(blk.hx_brine.cold_side.pressure_balance[0], 1e-5)
        constraint_scaling_transform(blk.hx_brine.hot_side.pressure_balance[0], 1e-5)

        # # evaporator
        # set_scaling_factor(blk.evaporator.area, 1e-3)
        # set_scaling_factor(blk.evaporator.U, 1e-3)
        # set_scaling_factor(blk.evaporator.delta_temperature_in, 1e-1)
        # set_scaling_factor(blk.evaporator.delta_temperature_out, 1e-1)
        # set_scaling_factor(blk.evaporator.lmtd, 1e-1)

        # # compressor
        # set_scaling_factor(blk.compressor.control_volume.work, 1e-6)

        # # condenser
        # set_scaling_factor(blk.condenser.control_volume.heat, 1e-6)

    else:
        pass

    properties_feed.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    properties_feed.set_default_scaling(
        "conc_mass_phase_comp", 1e-2, index=("Liq", "TDS")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    # Evaporator
    set_scaling_factor(blk.evaporator.area, 1e-3)
    set_scaling_factor(blk.evaporator.U, 1e-3)
    set_scaling_factor(blk.evaporator.delta_temperature_in, 1e-1)
    set_scaling_factor(blk.evaporator.delta_temperature_out, 1e-1)
    set_scaling_factor(blk.evaporator.lmtd, 1e-1)

    # Compressor
    set_scaling_factor(blk.compressor.control_volume.work, 1e-6)

    # Condenser
    set_scaling_factor(blk.condenser.control_volume.heat, 1e-6)

    if calc_blk_scaling_factors:
        calculate_scaling_factors(blk)

    else:
        calculate_scaling_factors(m)


def set_system_operating_conditions(m, Qin=5, tds=130, feed_temp=25):

    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.g / pyunits.liter
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(Qin * tds, to_units=pyunits.kg / pyunits.s)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_in,
            ("conc_mass_phase_comp", ("Liq", "TDS")): tds,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )

    set_mvc_operating_conditions(m, m.fs.MVC)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_water)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(flow_mass_tds)
    # m.fs.feed.properties[0].temperature.fix(273.15 + 50.52)  # K
    # m.fs.feed.properties[0].pressure.fix(1e5)  # Pa
    # m.fs.feed.properties[0].conc_mass_phase_comp[...]


def init_system(m, blk, delta_temperature_in=30, delta_temperature_out=None):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_mvc)

    init_mvc(
        m,
        blk,
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
    )

    m.fs.product.initialize()
    propagate_state(m.fs.mvc_to_product)

    m.fs.disposal.initialize()
    propagate_state(m.fs.mvc_to_disposal)


def init_mvc(m, blk, delta_temperature_in=30, delta_temperature_out=None, solver=None):

    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # Touch feed mass fraction property
    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.feed.initialize()

    # results = solver.solve(blk.feed)
    # _log.info(f"{blk.name} feed {results.solver.termination_condition}")

    # Propagate vapor flow rate based on given recovery
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] = (
        blk.recovery
        * (
            blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        )
    )
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"] = 0

    # Propagate brine salinity and flow rate
    blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"] = (
        blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"] / (1 - blk.recovery)
    )
    blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "H2O"] = (
        1 - blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
    )
    blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"] = (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"] = (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        - blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )

    # initialize feed pump
    propagate_state(blk.feed_to_pump)
    blk.pump_feed.initialize(optarg=optarg, solver="ipopt-watertap")

    # initialize separator
    propagate_state(blk.pump_to_separator)
    # Touch property for initialization
    blk.separator.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.split_fraction[0, "hx_distillate_cold"].fix(blk.recovery.value)
    blk.separator.mixed_state.initialize(optarg=optarg, solver="ipopt-watertap")
    # Touch properties for initialization
    blk.separator.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.initialize(optarg=optarg, solver="ipopt-watertap")
    blk.separator.split_fraction[0, "hx_distillate_cold"].unfix()

    # initialize distillate heat exchanger
    propagate_state(blk.sep_dist_cold_to_hx_dist_cold)
    blk.hx_distillate.cold_outlet.temperature[0] = (
        blk.evaporator.inlet_feed.temperature[0].value
    )
    blk.hx_distillate.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[
        0
    ].value
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
    blk.hx_distillate.hot_inlet.temperature[0] = (
        blk.evaporator.outlet_brine.temperature[0].value
    )
    blk.hx_distillate.hot_inlet.pressure[0] = 101325
    blk.hx_distillate.initialize(solver="ipopt-watertap")

    # initialize brine heat exchanger
    propagate_state(blk.sep_brine_cold_to_hx_brine_cold)
    blk.hx_brine.cold_outlet.temperature[0] = blk.evaporator.inlet_feed.temperature[
        0
    ].value
    blk.hx_brine.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[0].value
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = (
        blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    blk.hx_brine.hot_inlet.temperature[0] = blk.evaporator.outlet_brine.temperature[
        0
    ].value
    blk.hx_brine.hot_inlet.pressure[0] = 101325
    blk.hx_brine.initialize(solver="ipopt-watertap")

    # initialize mixer
    propagate_state(blk.hx_dist_cold_to_mixer)
    propagate_state(blk.hx_brine_cold_to_mixer)
    blk.mixer_feed.initialize(solver="ipopt-watertap")
    blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()
    print(f"\ndof @ 1 = {degrees_of_freedom(m)}\n")

    # initialize evaporator
    propagate_state(blk.mixer_feed_to_evaporator)
    blk.external_heating.fix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    # fixes and unfixes those values
    blk.evaporator.initialize(
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
        solver="ipopt-watertap",
    )
    blk.external_heating.unfix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # initialize compressor
    propagate_state(blk.evaporator_to_compressor)
    blk.compressor.initialize(solver="ipopt-watertap")

    # initialize condenser
    propagate_state(blk.compressor_to_condenser)
    blk.condenser.initialize(
        heat=-blk.evaporator.heat_transfer.value, solver="ipopt-watertap"
    )

    # initialize brine pump
    propagate_state(blk.evaporator_to_brine_pump)
    blk.pump_brine.initialize(optarg=optarg, solver="ipopt-watertap")
    # propagate_state(blk.brine_pump_to_hx_brine_hot)

    # initialize distillate pump
    propagate_state(blk.condenser_to_translator)  # to translator block
    propagate_state(blk.translated_to_dist_pump)  # from translator block to pump
    blk.pump_distillate.control_volume.properties_in[0].temperature = (
        blk.condenser.control_volume.properties_out[0].temperature.value
    )
    blk.pump_distillate.control_volume.properties_in[0].pressure = (
        blk.condenser.control_volume.properties_out[0].pressure.value
    )
    blk.pump_distillate.initialize(optarg=optarg, solver="ipopt-watertap")
    # propagate_state(blk.dist_pump_to_hx_dist_hot)

    # propagate brine state
    propagate_state(blk.hx_brine_hot_to_disposal)
    propagate_state(blk.hx_dist_hot_to_product)

    print(f"\ndof @ 2 = {degrees_of_freedom(m)}\n")

    blk.product.initialize()
    blk.disposal.initialize()

    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    blk.pump_brine.control_volume.deltaP[0].unfix()
    blk.disposal.properties[0].pressure.fix(101325)
    print(f"\ndof @ 3 = {degrees_of_freedom(m)}\n")

    print(f"termination {results.solver.termination_condition}")
    # _log.info_high(f"terminatinon {results.solver.termination_condition}")

    print(f"Initialization done, dof = {degrees_of_freedom(m)}")

def run_sequential_decomposition(m, blk, delta_temperature_in=30, delta_temperature_out=None, tear_solver="cbc"):

    seq = SequentialDecomposition(tear_solver=tear_solver)
    seq.options.log_info = True
    seq.options.iterLim = 1000

    def func_initialize(unit):
        print(unit.local_name)
        print(f"dof = {degrees_of_freedom(unit)}\n")
        if unit.local_name == "feed":
            pass
        elif unit.local_name == "condenser":
            unit.initialize(
                heat=-unit.flowsheet().evaporator.heat_transfer.value,
                optarg=solver.options,
                solver="ipopt-watertap",
            )
        elif unit.local_name == "evaporator":
            unit.flowsheet().external_heating.fix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
            unit.initialize(
                delta_temperature_in=delta_temperature_in, solver="ipopt-watertap"
            )
            unit.flowsheet().external_heating.unfix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        elif unit.local_name == "separator":
            unit.split_fraction[0, "hx_distillate_cold"].fix(
                unit.flowsheet().recovery.value
            )
            unit.initialize(solver="ipopt-watertap")
            unit.split_fraction[0, "hx_distillate_cold"].unfix()
        elif unit.local_name == "mixer_feed":
            unit.initialize(solver="ipopt-watertap")
            unit.pressure_equality_constraints[0, 2].deactivate()
        elif unit.local_name == "hx_distillate":
            unit.cold_outlet.temperature[0] = (
                blk.evaporator.inlet_feed.temperature[0].value
            )
            unit.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[
                0
            ].value
            unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
                blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
            )
            unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
            unit.hot_inlet.temperature[0] = (
                blk.evaporator.outlet_brine.temperature[0].value
            )
            unit.hot_inlet.pressure[0] = 101325
            unit.initialize(solver="ipopt-watertap")
        else:
            unit.initialize(solver="ipopt-watertap")


    seq.run(blk, func_initialize)


def add_external_heating(m, blk):
    # Allows additional heat to be added to evaporator so that an initial feasible solution can be found as a starting
    # guess for optimization in case physically infeasible simulation is proposed

    blk.external_heating = Var(
        initialize=0,
        units=pyunits.watt,
        bounds=(0, None),
        doc="External heating for evaporator",
    )

    blk.evaporator.eq_energy_balance.deactivate()
    blk.evaporator.eq_energy_balance_with_external_heat = Constraint(
        expr=blk.evaporator.heat_transfer
        + blk.external_heating
        + blk.evaporator.properties_feed[0].enth_flow
        == blk.evaporator.properties_brine[0].enth_flow
        + blk.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
    )
    iscale.set_scaling_factor(blk.external_heating, 1e-6)


def add_mvc_costing(m, blk, flowsheet_costing_block=None):

    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    blk.evaporator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )

    blk.pump_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.pump_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.pump_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.hx_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.hx_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.mixer_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    # m.fs.costing.TIC.fix(2)
    # m.fs.costing.electricity_cost = 0.1  # 0.15
    flowsheet_costing_block.heat_exchanger.material_factor_cost.fix(5)
    flowsheet_costing_block.evaporator.material_factor_cost.fix(5)
    flowsheet_costing_block.compressor.unit_cost.fix(1 * 7364)


def display_metrics(m, blk):
    print("\nSystem metrics")
    print(
        "Feed flow rate:                           %.2f kg/s"
        % (
            blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
            + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].value
        )
    )
    print(
        "Feed salinity:                            %.2f g/kg"
        % (blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].value * 1e3)
    )
    print(
        "Brine salinity:                           %.2f g/kg"
        % (
            blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
            * 1e3
        )
    )
    print(
        "Product flow rate:                        %.2f kg/s"
        % blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    print(
        "Recovery:                                 %.2f %%"
        % (blk.recovery.value * 100)
    )
    print(
        "Specific energy consumption:              %.2f kWh/m3"
        % value(m.fs.costing.SEC)
    )
    print(
        "Levelized cost of water:                  %.2f $/m3" % value(m.fs.costing.LCOW)
    )
    print(
        "External Q:                               %.2f W" % blk.external_heating.value
    )  # should be 0 for optimization


def display_design(m, blk):
    print("\nState variables")
    print(
        "Preheated feed temperature:               %.2f K"
        % blk.evaporator.properties_feed[0].temperature.value
    )
    print(
        "Evaporator (brine, vapor) temperature:    %.2f K"
        % blk.evaporator.properties_brine[0].temperature.value
    )
    print(
        "Evaporator (brine, vapor) pressure:       %.2f kPa"
        % (blk.evaporator.properties_vapor[0].pressure.value * 1e-3)
    )
    print(
        "Compressed vapor temperature:             %.2f K"
        % blk.compressor.control_volume.properties_out[0].temperature.value
    )
    print(
        "Compressed vapor pressure:                %.2f kPa"
        % (blk.compressor.control_volume.properties_out[0].pressure.value * 1e-3)
    )
    print(
        "Condensed vapor temperature:              %.2f K"
        % blk.condenser.control_volume.properties_out[0].temperature.value
    )

    print("\nDesign variables")
    print(
        "Brine heat exchanger area:                %.2f m2" % blk.hx_brine.area.value
    )
    print(
        "Distillate heat exchanger area:           %.2f m2"
        % blk.hx_distillate.area.value
    )
    print(
        "Compressor pressure ratio:                %.2f"
        % blk.compressor.pressure_ratio.value
    )
    print(
        "Evaporator area:                          %.2f m2" % blk.evaporator.area.value
    )
    print(
        "Evaporator LMTD:                          %.2f K" % blk.evaporator.lmtd.value
    )



def set_up_optimization(m, blk):

    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    blk.external_heating.fix(0)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()


if __name__ == "__main__":
    m = build_system(single_stage_build=True)
    mvc = m.fs.MVC

    set_mvc_scaling(m, mvc, single_stage_build=True)

    set_system_operating_conditions(m)
    mvc.recovery.fix(0.6)
    # set_mvc_operating_conditions(m, mvc)

    print(f"\ndof = {degrees_of_freedom(m)}\n")
    init_system(m, mvc)

    add_mvc_costing(m, mvc)

    flow_vol = mvc.product.properties[0].flow_vol_phase["Liq"]
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(flow_vol)
    m.fs.costing.add_specific_energy_consumption(flow_vol, name="SEC")

    m.fs.objective = Objective(expr=mvc.external_heating)
    print(f"dof = {degrees_of_freedom(m)}")

    # print(f"termination = {results.solver.termination_condition}")
    # print(f"dof = {degrees_of_freedom(m)}")
    # print(f"LCOW = {m.fs.costing.LCOW()}")
    # print(f"SEC = {m.fs.costing.SEC()}")
    print("\n***---First solve - simulation results---***")
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    display_metrics(m, mvc)
    display_design(m, mvc)

    m.fs.del_component(m.fs.objective)
    print("\n***---Second solve - optimization results---***")
    set_up_optimization(m, mvc)
    print(f"dof = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)
    display_metrics(m, mvc)
    display_design(m, mvc)
    print(f"dof = {degrees_of_freedom(m)}")
    mvc.evaporator.area.fix()
    mvc.evaporator.outlet_brine.temperature[0].fix()
    # mvc.compressor.pressure_ratio.fix()
    mvc.hx_distillate.area.fix()
    mvc.hx_brine.area.fix()
    print(f"dof = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)


    # mvc.evaporator.area.display()
