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
from idaes.models.unit_models import (
    Product,
    Feed,
    Mixer,
    StateJunction,
    Separator,
    Heater,
)
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
    "build_and_run_mvc",
    "set_mvc_operating_conditions",
    "set_mvc_scaling",
    "init_mvc",
    "add_mvc_costing",
    "scale_mvc_costs",
    "display_metrics",
    "display_costing",
    "display_design",
    "run_sequential_decomposition",
]

electricity_cost_base = 0.0434618999  # USD_2018/kWh equivalent to 0.0575 USD_2023/kWh
heat_cost_base = 0.00894

solver = get_solver()

reflo_dir = pathlib.Path(__file__).resolve().parents[4]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3

_log = idaeslog.getLogger("permian_MVC")


def build_and_run_mvc(
    recovery=0.5,
    Qin=4.9,
    tds=118,
    electricity_cost=electricity_cost_base,
    heat_cost=heat_cost_base,
    **kwargs,
):
    m = build_mvc_system(recovery=recovery)
    m.fs.costing.electricity_cost.fix(electricity_cost)
    m.fs.costing.heat_cost.fix(heat_cost)
    mvc = m.fs.MVC

    m.fs.optimal_solve_mvc = Var(initialize=0)
    m.fs.optimal_solve_mvc.fix()

    set_system_operating_conditions(m, Qin=Qin, tds=tds, feed_temp=25)
    set_mvc_operating_conditions(m, mvc)
    set_mvc_scaling(m, mvc)
    init_system(m, mvc)

    m.fs.costing.LCOW_obj = Objective(expr=m.fs.costing.LCOW)

    mvc.evaporator.area.unfix()
    mvc.evaporator.outlet_brine.temperature[0].unfix()
    mvc.compressor.pressure_ratio.unfix()
    # mvc.compressor.pressure_ratio.fix(1.6)
    mvc.hx_distillate.area.unfix()
    mvc.hx_brine.area.unfix()
    mvc.recovery_mass.unfix()
    mvc.recovery_vol.fix(m.recovery_vol)
    print(f"\n~~~~~FOURTH MVC SOLVE~~~~")
    try:
        results = solve_mvc(m)
        print(f"termination MVC FOURTH {results.solver.termination_condition}")
    except:
        pass

    display_metrics(m, mvc)
    display_design(m, mvc)
    # del m.fs.costing.LCOW_obj
    # mvc.evaporator.area.fix()
    # mvc.hx_distillate.area.fix()
    # mvc.hx_brine.area.fix()
    mvc.compressor.pressure_ratio.fix(1.6)
    print(f"dof MVC FINAL = {degrees_of_freedom(m)}")
    # assert degrees_of_freedom(m) == 0
    m.fs.product.properties[0].flow_vol_phase
    # m.fs.product.properties[0].conc_mass_phase_comp
    m.fs.product.initialize()
    m.fs.disposal.properties[0].flow_vol_phase
    # m.fs.disposal.properties[0].conc_mass_phase_comp

    print(f"\n~~~~~FINAL MVC SOLVE~~~~")
    try:
        results = solve_mvc(m)
        print(f"termination MVC FINAL {results.solver.termination_condition}")
        if not check_optimal_termination(results):
            m.fs.optimal_solve_mvc.fix(0)
    except:
        pass
    # m.fs.costing.LCOW.display()

    display_metrics(m, mvc)
    display_design(m, mvc)

    return m


def solve_mvc(m):
    solver = get_solver()
    try:
        results = solver.solve(m, tee=False)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)
        m.fs.optimal_solve_mvc.fix(1)
    except:
        results = solver.solve(m, tee=False)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)
        m.fs.optimal_solve_mvc.fix(1)

    return results


def build_mvc_system(recovery=0.5, **kwargs):
    m = ConcreteModel()
    m.recovery_mass = recovery
    m.recovery_vol = recovery

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    m.fs.costing.electricity_cost.fix(electricity_cost_base)
    m.fs.costing.heat_cost.fix(heat_cost_base)

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.product = Product(property_package=m.fs.properties_feed)
    m.fs.disposal = Product(property_package=m.fs.properties_feed)

    m.fs.MVC = mvc = FlowsheetBlock(dynamic=False)

    build_mvc(m, mvc, **kwargs)

    m.fs.feed_to_mvc = Arc(source=m.fs.feed.outlet, destination=mvc.feed.inlet)

    m.fs.mvc_to_product = Arc(source=mvc.product.outlet, destination=m.fs.product.inlet)

    m.fs.mvc_to_disposal = Arc(
        source=mvc.disposal.outlet, destination=m.fs.disposal.inlet
    )

    add_mvc_costing(m, mvc)
    # print(pyunits.get_units(mvc.evaporator.costing.capital_cost))
    # assert False
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(mvc.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(mvc.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        mvc.product.properties[0].flow_vol, name="SEC"
    )

    m.fs.costing.heat_exchanger.material_factor_cost.fix(5)
    m.fs.costing.evaporator.material_factor_cost.fix(5)
    m.fs.costing.compressor.unit_cost.fix(1 * 7364)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mvc(m, blk, external_heating=True):
    print(f'\n{"=======> BUILDING MVC SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties_feed)
    blk.product = StateJunction(property_package=m.fs.properties_feed)
    blk.disposal = StateJunction(property_package=m.fs.properties_feed)

    blk.recovery_mass = Var(
        initialize=0.5, bounds=(0, 1), units=pyunits.dimensionless, doc="MVC recovery"
    )

    blk.recovery_vol = Var(
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="MVC volumetric recovery",
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
    blk.hx_distillate.area.setlb(5)

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
    blk.hx_brine.area.setlb(5)

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

    blk.dist_pump_to_hx_dist_hot = Arc(
        source=blk.pump_distillate.outlet, destination=blk.hx_distillate.hot_inlet
    )
    blk.hx_dist_hot_to_product = Arc(
        source=blk.hx_distillate.hot_outlet, destination=blk.product.inlet
    )

    blk.eq_recovery_mass = Constraint(
        expr=blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
        == blk.recovery_mass
        * (
            blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        )
    )

    blk.eq_recovery_vol = Constraint(
        expr=blk.product.properties[0].flow_vol_phase["Liq"]
        == blk.recovery_vol * (blk.feed.properties[0].flow_vol_phase["Liq"])
    )

    blk.eq_separator_split_frac = Constraint(
        expr=blk.separator.split_fraction[0, "hx_distillate_cold"] == blk.recovery_mass
    )

    if external_heating:
        add_external_heating(m, blk)

    TransformationFactory("network.expand_arcs").apply_to(m)

    blk.evaporator.connect_to_condenser(blk.condenser)
    _log.info("MVC flowsheet built")


def set_mvc_operating_conditions(
    m,
    blk,
    recovery=0.5,
    inlet_brine_temp_guess=50,  # degC
    outlet_brine_temp=70,  # degC
    steam_temp_ub=75,
    **kwargs,
):
    """
    Generic initial point for MVC system.
    """

    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].fix(0.1)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40)
    # Feed pump
    blk.pump_feed.efficiency_pump[0].fix(0.8)
    blk.pump_feed.control_volume.deltaP[0].fix(7e3)

    # Separator
    blk.separator.split_fraction[0, "hx_distillate_cold"] = 0.5

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
    blk.evaporator.properties_vapor[0].temperature.setub(steam_temp_ub + 273.15)
    blk.evaporator.inlet_feed.temperature[0] = (
        inlet_brine_temp_guess + 273.15
    )  # provide guess
    blk.evaporator.outlet_brine.temperature[0].fix(outlet_brine_temp + 273.15)
    blk.evaporator.U.fix(3e3)  # W/K-m^2
    blk.evaporator.area.setub(1e4)  # m^2
    # blk.evaporator.area.set_value(4275)  # m^2

    # Compressor
    blk.compressor.control_volume.properties_out[0].temperature.setub(450)
    blk.compressor.pressure_ratio.fix(1.6)
    blk.compressor.efficiency.fix(0.8)

    # Brine pump
    blk.pump_brine.efficiency_pump[0].fix(0.8)
    blk.pump_brine.control_volume.deltaP[0].fix(4e4)

    # Distillate pump
    blk.pump_distillate.efficiency_pump[0].fix(0.8)
    blk.pump_distillate.control_volume.deltaP[0].fix(4e4)

    # Fix 0 TDS
    blk.tb_sw_to_water.properties_out[0].flow_mass_phase_comp["Liq", "TDS"].fix(1e-5)

    print("DOF after setting operating conditions: ", degrees_of_freedom(blk))


def set_system_operating_conditions(m, Qin=1, tds=130, feed_temp=25):
    global flow_in

    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.g / pyunits.liter
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_in,
            ("conc_mass_phase_comp", ("Liq", "TDS")): tds,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )


def set_mvc_scaling(
    m,
    blk,
    properties_feed=None,
    properties_vapor=None,
    calc_blk_scaling_factors=True,
):
    if properties_feed is None:
        properties_feed = m.fs.properties_feed

    if properties_vapor is None:
        properties_vapor = m.fs.properties_vapor

    properties_feed.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    set_scaling_factor(blk.external_heating, 1e-6)

    # MVC FEED
    # set_scaling_factor(
    #     blk.feed.properties[0.0].conc_mass_phase_comp["Liq", "TDS"], 1e-2
    # )
    # set_scaling_factor(blk.hx_distillate.hot_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"], 1e2)

    # MVC PRODUCT

    # MVC DISPOSAL

    # PUMPS
    set_scaling_factor(blk.pump_feed.control_volume.work, 1e-3)
    set_scaling_factor(blk.pump_brine.control_volume.work, 1e-3)
    set_scaling_factor(blk.pump_distillate.control_volume.work, 1e-3)

    # DISTILLATE HX
    set_scaling_factor(blk.hx_distillate.hot.heat, 1e-3)
    set_scaling_factor(blk.hx_distillate.cold.heat, 1e-3)
    set_scaling_factor(blk.hx_distillate.overall_heat_transfer_coefficient, 1e-3)

    set_scaling_factor(blk.hx_distillate.area, 1e-1)
    constraint_scaling_transform(blk.hx_distillate.cold_side.pressure_balance[0], 1e-5)
    constraint_scaling_transform(blk.hx_distillate.hot_side.pressure_balance[0], 1e-5)

    # BRINE HX
    set_scaling_factor(blk.hx_brine.hot.heat, 1e-3)
    set_scaling_factor(blk.hx_brine.cold.heat, 1e-3)
    set_scaling_factor(blk.hx_brine.overall_heat_transfer_coefficient, 1e-3)
    set_scaling_factor(blk.hx_brine.area, 1e-1)
    constraint_scaling_transform(blk.hx_brine.cold_side.pressure_balance[0], 1e-5)
    constraint_scaling_transform(blk.hx_brine.hot_side.pressure_balance[0], 1e-5)

    # EVAPORATOR
    # set_scaling_factor(blk.evaporator.properties_brine[0].pressure, 1e-4)
    set_scaling_factor(blk.evaporator.area, 1e-3)
    set_scaling_factor(blk.evaporator.U, 1e-3)
    set_scaling_factor(blk.evaporator.delta_temperature_in, 1e-1)
    set_scaling_factor(blk.evaporator.delta_temperature_out, 1e-1)
    set_scaling_factor(blk.evaporator.lmtd, 1e-1)
    # set_scaling_factor(blk.evaporator.heat_transfer, 1e-7)

    # COMPRESSOR
    set_scaling_factor(blk.compressor.control_volume.work, 1e-6)
    # set_scaling_factor(
    #     blk.compressor.control_volume.properties_in[0].enth_flow_phase["Liq"], 1
    # )

    # CONDENSER
    set_scaling_factor(blk.condenser.control_volume.heat, 1e-5)

    # if hasattr(blk.evaporator, "costing"):
    #     scale_mvc_costs(m, blk)

    if calc_blk_scaling_factors:
        calculate_scaling_factors(blk)

    else:
        calculate_scaling_factors(m)


def init_system(m, blk, **kwargs):
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_mvc)

    init_mvc(m, blk, **kwargs)

    m.fs.product.initialize()
    propagate_state(m.fs.mvc_to_product)

    m.fs.disposal.initialize()
    propagate_state(m.fs.mvc_to_disposal)


def init_mvc(
    m,
    blk,
    feed_props=None,
    delta_temperature_in=10,
    delta_temperature_out=None,
    solver=None,
):
    """
    Initialization routine for generic MVC setup.
    To be used with external heating.
    """

    if solver is None:
        solver = get_solver()
    solver.options["halt_on_ampl_error"] = "yes"
    # solver.options["tee"] = True
    optarg = solver.options

    if feed_props is None:
        feed_props = m.fs.feed.properties[0]

    blk.recovery_mass.fix(0.5)

    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]

    blk.feed.properties[0].temperature.fix(value(feed_props.temperature))
    blk.feed.properties[0].pressure.fix(value(feed_props.pressure))

    solver.solve(blk.feed)

    _log.info(f"{blk.name} feed initialization complete.")

    # Propagate vapor flow rate based on given recovery
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp[
        "Vap", "H2O"
    ] = blk.recovery_mass * (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"] = 0

    # Propagate brine salinity and flow rate
    blk.evaporator.properties_brine[0].mass_frac_phase_comp[
        "Liq", "TDS"
    ] = blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"] / (
        1 - blk.recovery_mass
    )
    blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "H2O"] = (
        1 - blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
    )
    blk.evaporator.properties_brine[0].flow_mass_phase_comp[
        "Liq", "TDS"
    ] = blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"] = (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        - blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )

    # Initialize feed pump
    propagate_state(blk.feed_to_pump)
    blk.pump_feed.initialize(optarg=optarg, solver="ipopt-watertap")
    _log.info(f"{blk.name} feed pump initialization complete.")

    # Initialize separator
    propagate_state(blk.pump_to_separator)
    # Touch property for initialization
    blk.separator.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.split_fraction[0, "hx_distillate_cold"].fix(blk.recovery_mass.value)
    blk.separator.mixed_state.initialize(optarg=optarg, solver="ipopt-watertap")

    # Touch properties for initialization
    blk.separator.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.initialize(optarg=optarg, solver="ipopt-watertap")
    blk.separator.split_fraction[0, "hx_distillate_cold"].unfix()
    _log.info(f"{blk.name} separator initialization complete.")

    # Initialize distillate heat exchanger
    propagate_state(blk.sep_dist_cold_to_hx_dist_cold)
    blk.hx_distillate.cold_outlet.temperature[
        0
    ] = blk.evaporator.inlet_feed.temperature[0].value
    blk.hx_distillate.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[
        0
    ].value
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
    blk.hx_distillate.hot_inlet.temperature[
        0
    ] = blk.evaporator.outlet_brine.temperature[0].value
    blk.hx_distillate.hot_inlet.pressure[0] = 101325
    blk.hx_distillate.initialize(solver="ipopt-watertap")
    _log.info(f"{blk.name} Distillate HX initialization complete.")

    # Initialize brine heat exchanger
    propagate_state(blk.sep_brine_cold_to_hx_brine_cold)
    blk.hx_brine.cold_outlet.temperature[0] = blk.evaporator.inlet_feed.temperature[
        0
    ].value
    blk.hx_brine.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[0].value
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ] = blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[
        0, "Liq", "TDS"
    ] = blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
    blk.hx_brine.hot_inlet.temperature[0] = blk.evaporator.outlet_brine.temperature[
        0
    ].value
    blk.hx_brine.hot_inlet.pressure[0] = 101325
    blk.hx_brine.initialize(solver="ipopt-watertap")
    _log.info(f"{blk.name} Brine HX initialization complete.")

    # Initialize mixer
    propagate_state(blk.hx_dist_cold_to_mixer)
    propagate_state(blk.hx_brine_cold_to_mixer)
    blk.mixer_feed.initialize(solver="ipopt-watertap")
    blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()
    _log.info(f"{blk.name} Mixer initialization complete.")
    # print(f"\ndof @ 1 = {degrees_of_freedom(blk)}\n")

    # Initialize evaporator
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

    # Initialize compressor
    propagate_state(blk.evaporator_to_compressor)
    blk.compressor.initialize(solver="ipopt-watertap")
    _log.info(f"{blk.name} Compressor initialization complete.")

    # Initialize condenser
    try:
        propagate_state(blk.compressor_to_condenser)
        blk.condenser.initialize(
            heat=-blk.evaporator.heat_transfer.value, solver="ipopt-watertap"
        )
        _log.info(f"{blk.name} Condenser initialization complete.")
    except:
        print_infeasible_constraints(blk.condenser)
        assert False

    # Initialize brine pump
    propagate_state(blk.evaporator_to_brine_pump)
    blk.pump_brine.initialize(optarg=optarg, solver="ipopt-watertap")
    _log.info(f"{blk.name} Brine Pump initialization complete.")
    # propagate_state(blk.brine_pump_to_hx_brine_hot)

    # Initialize distillate pump
    propagate_state(blk.condenser_to_translator)  # to translator block
    propagate_state(blk.translated_to_dist_pump)  # from translator block to pump
    blk.pump_distillate.control_volume.properties_in[
        0
    ].temperature = blk.condenser.control_volume.properties_out[0].temperature.value
    blk.pump_distillate.control_volume.properties_in[
        0
    ].pressure = blk.condenser.control_volume.properties_out[0].pressure.value
    blk.pump_distillate.initialize(optarg=optarg, solver="ipopt-watertap")
    _log.info(f"{blk.name} Distillate Pump initialization complete.")
    # propagate_state(blk.dist_pump_to_hx_dist_hot)

    # Propagate brine state
    propagate_state(blk.hx_brine_hot_to_disposal)
    propagate_state(blk.hx_dist_hot_to_product)

    # print(f"\ndof @ 2 = {degrees_of_freedom(blk)}\n")

    run_sequential_decomposition(
        m,
        blk,
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
    )

    blk.product.initialize()
    _log.info(f"{blk.name} Product initialization complete.")
    blk.disposal.initialize()
    _log.info(f"{blk.name} Disposal initialization complete.")

    m.fs.costing.initialize()
    results = solver.solve(blk, tee=False)
    print(f"MVC solve termination {results.solver.termination_condition}")
    _log.info(f"MVC solve termination {results.solver.termination_condition}")
    assert_optimal_termination(results)

    print(f"blk dof at end of init = {degrees_of_freedom(blk)}")

    blk.pump_brine.control_volume.deltaP[0].unfix()
    blk.disposal.properties[0].pressure.fix(101325)
    # print(f"\ndof @ 3 blk = {degrees_of_freedom(blk)}\n")
    # print(f"\ndof @ 3 model = {degrees_of_freedom(m)}\n")
    # _log.info_high(f"terminatinon {results.solver.termination_condition}")

    print(f"\n~~~~~FIRST SOLVE~~~~")
    blk.obj = Objective(expr=blk.external_heating)
    print(f"blk dof = {degrees_of_freedom(blk)}")
    print(f"model dof = {degrees_of_freedom(m)}")
    results = solver.solve(blk, tee=False)
    assert_optimal_termination(results)
    print(f"MVC FIRST solve termination {results.solver.termination_condition}")
    _log.info(f"MVC FIRST solve termination {results.solver.termination_condition}")
    display_metrics(m, blk)
    display_design(m, blk)

    blk.external_heating.fix(0)
    del blk.obj
    blk.external_heating.fix(0)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()

    print(f"\n~~~~~SECOND SOLVE~~~~")

    print(f"blk dof = {degrees_of_freedom(blk)}")
    print(f"model dof = {degrees_of_freedom(m)}")
    results = solver.solve(blk, tee=False)
    assert_optimal_termination(results)
    print(f"MVC SECOND solve termination {results.solver.termination_condition}")
    _log.info(f"MVC SECOND solve termination {results.solver.termination_condition}")

    display_metrics(m, blk)
    display_design(m, blk)

    mvc_feed_state_vars = blk.feed.properties[0].define_port_members()
    feed_state_vars = feed_props.define_port_members()
    blk.feed.properties[0].mass_frac_phase_comp.unfix()
    for k, v in mvc_feed_state_vars.items():
        if v.is_indexed():
            for i, vv in v.items():
                vv.fix(value(feed_state_vars[k][i]))
        else:
            v.fix(value(feed_state_vars[k]))
    try:
        print(f"\n~~~~~THIRD SOLVE~~~~")
        print(f"blk dof = {degrees_of_freedom(blk)}")
        print(f"model dof = {degrees_of_freedom(m)}")
        blk.feed.initialize()
        m.fs.costing.initialize()
        results = solver.solve(blk, tee=False)
        # results = solve_mvc(blk)
        assert_optimal_termination(results)
        print(f"MVC THIRD solve termination {results.solver.termination_condition}")
        _log.info(f"MVC THIRD solve termination {results.solver.termination_condition}")

        display_metrics(m, blk)
        display_design(m, blk)
    except:
        print(f"MVC THIRD SOLVE FAILED!!\n")
        pass

    # print(f"\n~~~~~FOURTH SOLVE~~~~")
    # blk.evaporator.area.unfix()
    # blk.evaporator.outlet_brine.temperature[0].unfix()
    # blk.compressor.pressure_ratio.unfix()
    # blk.hx_distillate.area.unfix()
    # blk.hx_brine.area.unfix()
    # blk.recovery_mass.fix(m.recovery_mass)
    # print(f"blk dof = {degrees_of_freedom(blk)}")
    # print(f"model dof = {degrees_of_freedom(m)}")
    # results = solver.solve(blk, tee=False)
    # print(f"termination {results.solver.termination_condition}")
    # assert_optimal_termination(results)

    # display_metrics(m, blk)
    # display_design(m, blk)

    mvc_feed_state_vars = blk.feed.properties[0].define_port_members()
    feed_state_vars = feed_props.define_port_members()
    for k, v in mvc_feed_state_vars.items():
        if v.is_indexed():
            for i, vv in v.items():
                vv.unfix()
        else:
            v.unfix()

    blk.feed.initialize()

    print(f"Initialization done, blk dof = {degrees_of_freedom(blk)}")
    print(f"Initialization done, model dof = {degrees_of_freedom(m)}")


def run_sequential_decomposition(
    m,
    blk,
    delta_temperature_in=60,
    delta_temperature_out=None,
    tear_solver="cbc",
    iterlim=5,
):
    def func_initialize(unit):
        print(unit.local_name)
        print(f"dof = {degrees_of_freedom(unit)}\n")
        # if unit.local_name == "feed":
        if unit.local_name in ["feed", "product", "disposal"]:
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
                delta_temperature_in=delta_temperature_in,
                delta_temperature_out=delta_temperature_out,
                solver="ipopt-watertap",
            )
            unit.flowsheet().external_heating.unfix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        elif unit.local_name == "separator":
            unit.split_fraction[0, "hx_distillate_cold"].fix(
                unit.flowsheet().recovery_mass.value
            )
            unit.initialize(solver="ipopt-watertap")
            unit.split_fraction[0, "hx_distillate_cold"].unfix()
        elif unit.local_name == "mixer_feed":
            unit.initialize(solver="ipopt-watertap")
            unit.pressure_equality_constraints[0, 2].deactivate()
        # elif unit.local_name == "hx_distillate":
        #     unit.cold_outlet.temperature[0] = blk.evaporator.inlet_feed.temperature[
        #         0
        #     ].value
        #     unit.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[0].value
        #     unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        #         blk.evaporator.properties_vapor[0]
        #         .flow_mass_phase_comp["Vap", "H2O"]
        #         .value
        #     )
        #     unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
        #     unit.hot_inlet.temperature[0] = blk.evaporator.outlet_brine.temperature[
        #         0
        #     ].value
        #     unit.hot_inlet.pressure[0] = 101325
        #     unit.initialize(solver="ipopt-watertap")
        else:
            unit.initialize(solver="ipopt-watertap")

    seq = SequentialDecomposition(tear_solver=tear_solver)
    seq.options.log_info = True
    seq.options.iterLim = iterlim
    seq.run(blk, func_initialize)


def add_external_heating(m, blk):
    """
    Add additional heat source so infeasible designs can be used
    as initial guess prior to optimization
    """

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
    set_scaling_factor(blk.external_heating, 1e-6)


def add_mvc_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    flowsheet_costing_block.base_currency = pyunits.USD_2023

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

    blk.costing = Block()
    blk.costing.capital_cost = Expression(
        expr=blk.evaporator.costing.capital_cost
        + blk.evaporator.costing.capital_cost
        + blk.compressor.costing.capital_cost
        + blk.pump_feed.costing.capital_cost
        + blk.pump_distillate.costing.capital_cost
        + blk.pump_brine.costing.capital_cost
        + blk.hx_brine.costing.capital_cost
        + blk.hx_distillate.costing.capital_cost
    )

    # m.fs.costing.TIC.fix(2)
    # m.fs.costing.electricity_cost = 0.1  # 0.15
    flowsheet_costing_block.heat_exchanger.material_factor_cost.fix(5)
    flowsheet_costing_block.evaporator.material_factor_cost.fix(5)
    flowsheet_costing_block.compressor.unit_cost.fix(1 * 7364)


def set_up_optimization(m, blk):
    """
    Add objective and unfix design variables
    """
    if hasattr(blk, "LCOW_obj"):
        blk.del_component(blk.LCOW_obj)
    blk.LCOW_obj = Objective(expr=m.fs.costing.LCOW)
    blk.external_heating.fix(0)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()


def design_mvc(m, blk):
    """
    After system is initialized, full design routine.
    """
    # First minimize external heating required
    # blk.external_heating.unfix()
    blk.external_heating_obj = Objective(expr=blk.external_heating)
    print(f"MVC DOF @ A = {degrees_of_freedom(blk)}")

    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"MVC first solve {results.solver.termination_condition}")

    print(f"\n~~~~~FIRST SOLVE~~~~")
    display_metrics(m, mvc)
    display_design(m, mvc)

    blk.del_component(blk.external_heating_obj)

    # Add LCOW objective; unfix design variables
    set_up_optimization(m, blk)
    # blk.evaporator.eq_energy_balance.activate()
    # blk.evaporator.eq_energy_balance_with_external_heat.deactivate()

    print(f"MVC DOF @ B = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"MVC second solve {results.solver.termination_condition}")

    print(f"\n~~~~~SECOND SOLVE~~~~")
    display_metrics(m, mvc)
    display_design(m, mvc)

    print(f"MVC DOF @ C = {degrees_of_freedom(m)}")

    # blk.evaporator.eq_energy_balance.activate()
    blk.evaporator.area.fix()
    blk.evaporator.outlet_brine.temperature[0].fix()
    blk.hx_distillate.area.fix()
    blk.hx_brine.area.fix()

    print(f"MVC DOF final = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"MVC final solve {results.solver.termination_condition}")
    print(f"\n~~~~~FINAL SOLVE~~~~")
    display_metrics(m, mvc)
    display_design(m, mvc)


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
        % (blk.recovery_mass.value * 100)
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
        "Evaporator (brine, vapor) pressure:       %.5f kPa"
        % (blk.evaporator.properties_vapor[0].pressure.value * 1e-3)
    )
    print(
        "Compressed vapor temperature:             %.5f K"
        % blk.compressor.control_volume.properties_out[0].temperature.value
    )
    print(
        "Compressed vapor pressure:                %.5f kPa"
        % (blk.compressor.control_volume.properties_out[0].pressure.value * 1e-3)
    )
    print(
        "Condensed vapor temperature:              %.5f K"
        % blk.condenser.control_volume.properties_out[0].temperature.value
    )

    print("\nDesign variables")
    print("Brine heat exchanger area:                %.2f m2" % blk.hx_brine.area.value)
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


def display_costing(m, mvc):
    print("\n\n---------MVC COSTING---------")
    print(f"\n{'LCOW':<30} {m.fs.costing.LCOW():<20.2f} $/m3")
    print(f"{'SEC':<30} {m.fs.costing.SEC():<20.2f} kWh/m3")
    print(
        f"{'Aggregate Electricity Flow':<30} {m.fs.costing.aggregate_flow_electricity():<20.2f} kW"
    )
    print(f"{'Total Capital Cost':<30} ${m.fs.costing.total_capital_cost():<20.2f}")

    print(
        f"\n{'Evaporator Capital Cost':<30} ${mvc.evaporator.costing.capital_cost():<20.2f}"
    )
    print(f"{'Evaporator Area':<30} {mvc.evaporator.area():<20.2f} m2")

    print(
        f"\n{'Compressor Capital Cost':<30} ${mvc.compressor.costing.capital_cost():<20.2f}"
    )
    print(f"{'Compressor Efficiency':<30} {mvc.compressor.efficiency():<20.2f} ---")
    print(f"{'Compressor Work':<30} {mvc.compressor.control_volume.work[0]():<20.2f} W")
    print(
        f"{'Compressor Pressure Ratio':<30} {mvc.compressor.pressure_ratio():<20.2f} ---"
    )
    print(
        f"{'Compressor Steam Flow In':<30} {mvc.compressor.control_volume.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']():<20.2f} kg/s"
    )

    print(
        f"\n{'Feed Pump Capital Cost':<30} ${mvc.pump_feed.costing.capital_cost():<20.2f}"
    )
    print(f"{'Feed Pump Work':<30} {mvc.pump_feed.work_mechanical[0]():<20.2f} W")

    print(
        f"\n{'Distillate Pump Capital Cost':<30} ${mvc.pump_distillate.costing.capital_cost():<20.2f}"
    )
    print(
        f"{'Distillate Pump Work':<30} {pyunits.convert(mvc.pump_distillate.work_mechanical[0], to_units=pyunits.watt)():<20.2f} W"
    )

    print(
        f"\n{'Brine Pump Capital Cost':<30} ${mvc.pump_brine.costing.capital_cost():<20.2f}"
    )
    print(f"{'Brine Pump Work':<30} {mvc.pump_brine.work_mechanical[0]():<20.2f} W")

    print(
        f"\n{'HX Distillate Capital Cost':<30} ${mvc.hx_distillate.costing.capital_cost():<20.2f}"
    )
    print(f"{'HX Distillate Area':<30} {mvc.hx_distillate.area():<20.2f} m2")

    print(
        f"\n{'HX Brine Capital Cost':<30} ${mvc.hx_brine.costing.capital_cost():<20.2f}"
    )
    print(f"{'HX Brine Area':<30} {mvc.hx_brine.area():<20.2f} m2")

    print(
        f"\n{'Mixer Capital Cost':<30} ${mvc.mixer_feed.costing.capital_cost():<20.2f}"
    )
    print(
        f"{'Mixer Flow Vol':<30} {mvc.mixer_feed.mixed_state[0].flow_vol():<20.5f} m3/s"
    )


def calculate_cost_sf(cost):
    print(cost.name, cost.value)
    if cost.value in [0, None]:
        sf = 1e-2
    else:
        sf = 10 ** -(math.log10(abs(cost.value)))
    iscale.set_scaling_factor(cost, sf)


def scale_mvc_costs(m, blk):
    calculate_cost_sf(blk.hx_distillate.costing.capital_cost)
    calculate_cost_sf(blk.hx_brine.costing.capital_cost)
    calculate_cost_sf(blk.mixer_feed.costing.capital_cost)
    calculate_cost_sf(blk.evaporator.costing.capital_cost)
    calculate_cost_sf(blk.compressor.costing.capital_cost)
    # calculate_cost_sf(m.fs.costing.aggregate_capital_cost)
    # calculate_cost_sf(m.fs.costing.aggregate_flow_costs["electricity"])
    # calculate_cost_sf(m.fs.costing.total_capital_cost)
    # calculate_cost_sf(m.fs.costing.total_operating_cost)

    calculate_scaling_factors(m)

    print("Scaled costs")


if __name__ == "__main__":
    m = build_and_run_mvc(recovery=0.5, Qin=5.01, tds=130)
    blk = m.fs.MVC
    # blk.evaporator.costing.display()
    # blk.compressor.costing.display()

    # blk.pump_feed.costing.display()
    # blk.pump_distillate.costing.display()
    # blk.pump_brine.costing.display()
    # blk.hx_distillate.costing.display()
    # blk.hx_brine.costing.display()
    # blk.mixer_feed.costing.display()
    blk.costing.capital_cost.display()
    # m.fs.MVC.recovery_vol.display()
    # m.fs.costing.total_capital_cost.display()
    # m.fs.costing.total_operating_cost.display()
    # m.fs.product.properties[0].flow_vol_phase.display()
    # m.fs.MVC.product.properties[0].flow_vol_phase.display()
    # m.fs.costing.LCOW.display()


## USED TO RE-INIT MVC @ THIRD SOLVE
# blk.condenser.initialize(heat=-blk.evaporator.heat_transfer.value)

# blk.external_heating.fix()
# blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
# blk.evaporator.initialize(
#     delta_temperature_in=delta_temperature_in,
#     delta_temperature_out=delta_temperature_out,
#     solver="ipopt-watertap",
# )
# blk.external_heating.unfix()
# blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

# blk.separator.split_fraction[0, "hx_distillate_cold"].fix(
#     blk.recovery_mass.value
# )
# blk.separator.initialize(solver="ipopt-watertap")
# blk.separator.split_fraction[0, "hx_distillate_cold"].unfix()

# blk.mixer_feed.initialize(solver="ipopt-watertap")
# blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()
