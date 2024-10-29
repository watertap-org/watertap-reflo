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
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.util.initialization import *

from watertap.unit_models.mvc.components import Evaporator
from watertap.unit_models.mvc.components import Compressor
from watertap.unit_models.mvc.components import Condenser

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import TreatmentCosting

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


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.product = Product(property_package=m.fs.properties_vapor)
    m.fs.disposal = Product(property_package=m.fs.properties_feed)

    m.fs.MVC = mvc = FlowsheetBlock(dynamic=False)

    build_mvc(m, mvc)

    m.fs.feed_to_mvc = Arc(source=m.fs.feed.outlet, destination=mvc.feed.inlet)

    m.fs.mvc_to_product = Arc(source=mvc.product.outlet, destination=m.fs.product.inlet)

    m.fs.mvc_to_disposal = Arc(
        source=mvc.disposal.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mvc(m, blk):

    blk.feed = StateJunction(property_package=m.fs.properties_feed)
    blk.product = StateJunction(property_package=m.fs.properties_vapor)
    blk.disposal = StateJunction(property_package=m.fs.properties_feed)

    # Evaporator
    blk.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )
    # Compressor
    blk.compressor = Compressor(property_package=m.fs.properties_vapor)

    # Condenser
    blk.condenser = Condenser(property_package=m.fs.properties_vapor)

    blk.feed_to_evaporator = Arc(
        source=blk.feed.outlet, destination=blk.evaporator.inlet_feed
    )

    blk.evaporator_to_compressor = Arc(
        source=blk.evaporator.outlet_vapor, destination=blk.compressor.inlet
    )

    blk.compressor_to_condenser = Arc(
        source=blk.compressor.outlet, destination=blk.condenser.inlet
    )

    blk.condenser_to_product = Arc(
        source=blk.condenser.outlet, destination=blk.product.inlet
    )

    blk.evaporator_to_disposal = Arc(
        source=blk.evaporator.outlet_brine, destination=blk.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    blk.evaporator.connect_to_condenser(blk.condenser)


def set_mvc_operating_conditions(
    m,
    blk,
    outlet_brine_temp=60,  # degC
):

    blk.evaporator.outlet_brine.temperature[0].fix(273.15 + outlet_brine_temp)
    blk.evaporator.U.fix(1e3)
    # blk.evaporator.area.fix(400)
    blk.evaporator.area.fix(9125)

    blk.compressor.pressure_ratio.fix(2)
    blk.compressor.efficiency.fix(0.8)


def set_mvc_scaling(m, blk):

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
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

    calculate_scaling_factors(m)


def set_system_operating_conditions(m, Qin=5, tds=130):

    mvc = m.fs.MVC

    Qin = Qin * pyunits.Mgallons / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)
    flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
    flow_mass_tds = pyunits.convert(
        Qin * tds * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.s
    )

    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(10)
    # m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.05)

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_water)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(flow_mass_tds)
    m.fs.feed.properties[0].temperature.fix(273.15 + 50.52)  # K
    m.fs.feed.properties[0].pressure.fix(1e5)  # Pa
    m.fs.feed.properties[0].conc_mass_phase_comp[...]

    print(f"\nflow_mass_water = {flow_mass_water()}")
    print(f"\nflow_mass_tds = {flow_mass_tds()}\n")

    # inlet_dict = {"tds": tds * pyunits.kg / pyunits.m**3}

    # for "Liq", solute, solute_conc in inlet_dict.items():
    #     flow_mass_solute = pyunits.convert(
    #         Qin * solute_conc, to_units=pyunits.kg / pyunits.s
    #     )
    #     sf = 1 / value(flow_mass_solute)
    #     m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(flow_mass_solute)
    #     m.fs.properties_feed.set_default_scaling(
    #         "flow_mass_phase_comp",
    #         sf,
    #         index=("Liq", solute),
    #     )
    #     m.fs.properties_feed.set_default_scaling(
    #         "conc_mass_comp",
    #         1 / solute_conc(),
    #         index=("Liq", solute),
    #     )

    # m.fs.properties_feed.set_default_scaling(
    #     "flow_mass_phase_comp",
    #     1 / value(flow_mass_water),
    #     index=("H2O"),
    # )


def init_system(m, blk, delta_temperature_in=30, delta_temperature_out=5):

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


def init_mvc(m, blk, solver=None, delta_temperature_in=30, delta_temperature_out=5):

    if solver is None:
        solver = get_solver()

    blk.feed.initialize()
    propagate_state(blk.feed_to_evaporator)

    blk.evaporator.initialize(
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
    )
    propagate_state(blk.evaporator_to_compressor)
    propagate_state(blk.evaporator_to_disposal)

    blk.compressor.initialize()
    propagate_state(blk.compressor_to_condenser)

    blk.condenser.initialize(heat=-blk.evaporator.heat_transfer.value)
    propagate_state(blk.condenser_to_product)

    blk.product.initialize()

    blk.disposal.initialize()


def add_mvc_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    blk.evaporator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )


if __name__ == "__main__":
    m = build_system()
    mvc = m.fs.MVC

    set_system_operating_conditions(m)

    set_mvc_operating_conditions(m, mvc)

    set_mvc_scaling(m, mvc)

    add_mvc_costing(m, mvc)

    flow_vol = mvc.product.properties[0].flow_vol_phase["Liq"]
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(flow_vol)
    m.fs.costing.add_specific_energy_consumption(flow_vol, name="SEC")

    init_system(m, mvc)

    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(f"dof = {degrees_of_freedom(m)}")
    print(f"LCOW = {m.fs.costing.LCOW()}")
    print(f"SEC = {m.fs.costing.SEC()}")

    # mvc.evaporator.area.display()
