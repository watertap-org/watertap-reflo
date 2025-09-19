from pyomo.environ import (
    ConcreteModel,
    Expression,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import Product, Feed, StateJunction
from idaes.core.util.initialization import propagate_state

from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.zero_order import ElectrocoagulationZO

from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve, calc_scale

__all__ = [
    "build_EC",
    "set_EC_operating_conditions",
    "init_EC",
    "add_EC_costing",
    "add_ec_scaling",
    "report_EC",
    "print_EC_costing_breakdown",
]


def build_system():
    m = ConcreteModel()
    m.db = REFLODatabase()

    m.fs = FlowsheetBlock()
    m.fs.costing = TreatmentCosting()
    m.fs.properties = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.EC = FlowsheetBlock()

    build_EC(m.fs.EC)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.EC.feed.inlet,
    )
    m.fs.unit_to_disposal = Arc(
        source=m.fs.EC.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )
    m.fs.unit_to_product = Arc(
        source=m.fs.EC.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_EC(blk, prop_package=None):

    print(f'\n{"=======> BUILDING EC SYSTEM <=======":^60}\n')
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = ElectrocoagulationZO(
        property_package=prop_package,
        database=m.db,
        electrode_material="aluminum",
        reactor_material="pvc",
        overpotential_calculation="calculated",
        process_subtype="kbhdp",
    )

    blk.unit.SEC = Expression(
        expr=pyunits.convert(
            blk.unit.power_required / blk.unit.properties_in[0].flow_vol,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.treated,
        destination=blk.product.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )


def set_system_operating_conditions(m):

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(175.25054)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(2.143156)
    m.fs.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)


def set_EC_operating_conditions(blk):

    blk.unit.load_parameters_from_database(use_default_removal=True)


def add_ec_scaling(blk):

    iscale.set_scaling_factor(blk.unit.charge_loading_rate, 1e1)
    iscale.set_scaling_factor(blk.unit.electrode_volume, 1)
    iscale.set_scaling_factor(blk.unit.properties_in[0].flow_mass_comp["H2O"], 1e1)
    iscale.set_scaling_factor(blk.unit.properties_in[0].flow_mass_comp["tds"], 1)
    iscale.set_scaling_factor(blk.unit.properties_treated[0].flow_mass_comp["H2O"], 1e1)
    iscale.set_scaling_factor(blk.unit.properties_byproduct[0].flow_mass_comp["H2O"], 1)
    iscale.constraint_scaling_transform(blk.unit.eq_power_required, 1e-4)


def init_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)

    init_EC(m.fs.EC)
    propagate_state(m.fs.unit_to_product)
    propagate_state(m.fs.unit_to_disposal)


def init_EC(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    calculate_variable_from_constraint(
        blk.unit.overpotential, blk.unit.eq_overpotential
    )
    calculate_variable_from_constraint(
        blk.unit.applied_current, blk.unit.eq_applied_current
    )
    calculate_variable_from_constraint(
        blk.unit.anode_area, blk.unit.eq_electrode_area_total
    )
    calculate_variable_from_constraint(
        blk.unit.ohmic_resistance, blk.unit.eq_ohmic_resistance
    )

    blk.unit.initialize()
    propagate_state(blk.unit_to_product)
    propagate_state(blk.unit_to_disposal)


def add_EC_costing(blk, costing_block=None):

    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def report_EC(blk, w=25):
    title = "EC Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{f"Stream":<{w}}{f"Mass Flow Water (kg/s)":<{w}}{f"Mass Flow TDS (kg/s)":<{w}}'
    )
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Feed":<{w}}{value(blk.feed.properties[0].flow_mass_comp["H2O"]):<{w}.2f}{value(blk.feed.properties[0].flow_mass_comp["tds"]):<{w}.2f}'
    )
    print(
        f'{"Product":<{w}}{value(blk.product.properties[0].flow_mass_comp["H2O"]):<{w}.2f}{value(blk.product.properties[0].flow_mass_comp["tds"]):<{w}.2f}'
    )
    print(
        f'{"Disposal":<{w}}{value(blk.disposal.properties[0].flow_mass_comp["H2O"]):<{w}.2f}{value(blk.disposal.properties[0].flow_mass_comp["tds"]):<{w}.2f}'
    )
    print(f'\n\n{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Inlet Salinity":<{w}}{value(blk.feed.properties[0].conc_mass_comp["tds"]):<{w}.2f}{pyunits.get_units(blk.feed.properties[0].conc_mass_comp["tds"])}'
    )
    print(
        f'{"Solution Conductivity":<{w}}{value(blk.unit.conductivity):<{w}.2f}{pyunits.get_units(blk.unit.conductivity)}'
    )
    print(
        f'{"Ohmic Resistance":<{w}s}{f"{pyunits.convert(blk.unit.ohmic_resistance, to_units=pyunits.mohm)():<{w},.3e}"}{"mΩ"}'
    )
    print(
        f'{"Total Electrode Area":<{w}}{value(blk.unit.cathode_area)*2:<{w}.2f}{pyunits.get_units(blk.unit.cathode_area)}'
    )
    print(
        f'{"Total Electrode Area":<{w}}{value(blk.unit.cathode_area)*2:<{w}.2f}{pyunits.get_units(blk.unit.cathode_area)}'
    )
    print(
        f'{"Power Required":<{w}s}{f"{pyunits.convert(blk.unit.power_required, to_units=pyunits.kW)():<{w},.1f}"}{"kW"}'
    )


def print_EC_costing_breakdown(blk, w=25):
    title = "EC Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"EC CAPEX":<{w}s}{f"{blk.unit.costing.capital_cost():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost)}'
    )
    print(
        f'{"    Reactor":<{w}s}{f"{blk.unit.costing.capital_cost_reactor():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost_reactor)}'
    )
    print(
        f'{"    Electrodes":<{w}s}{f"{blk.unit.costing.capital_cost_electrodes():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost_electrodes)}'
    )
    print(
        f'{"    Power Supply":<{w}s}{f"{blk.unit.costing.capital_cost_power_supply():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost_power_supply)}'
    )
    print(
        f'{"    Floc Reactor":<{w}s}{f"{blk.unit.costing.capital_cost_floc_reactor():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost_floc_reactor)}'
    )
    print(
        f'{"EC Fixed OPEX":<{w}s}{f"{blk.unit.costing.fixed_operating_cost():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.fixed_operating_cost)}'
    )
    print(
        f'{"EC Power Required":<{w}s}{f"{blk.unit.costing.electricity_flow():<{w},.0f}"}{"kW"}'
    )
    print(f'{"SEC":<{w}s}{f"{blk.unit.SEC():<{w},.1f}"}{"kWh/m3"}')


def main():

    m = build_system()
    set_system_operating_conditions(m)
    set_EC_operating_conditions(m.fs.EC)
    add_ec_scaling(m.fs.EC)

    scale_flow = calc_scale(m.fs.feed.flow_mass_comp[0, "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_comp[0, "tds"].value)
    scale_tss = calc_scale(m.fs.feed.flow_mass_comp[0, "tss"].value)

    m.fs.properties.set_default_scaling(
        "flow_mass_comp", 10**-scale_flow, index=("H2O")
    )
    m.fs.properties.set_default_scaling("flow_mass_comp", 10**-scale_tds, index=("tds"))
    m.fs.properties.set_default_scaling("flow_mass_comp", 10**-scale_tss, index=("tss"))
    iscale.calculate_scaling_factors(m)

    init_system(m)
    add_EC_costing(m.fs.EC)
    m.fs.costing.cost_process()

    m.fs.costing.add_LCOW(m.fs.EC.unit.properties_treated[0].flow_vol)

    m.fs.costing.add_electricity_intensity(
        m.fs.EC.unit.properties_treated[0].flow_vol, name="SEC"
    )
    results = solve(m)
    assert_optimal_termination(results)

    report_EC(m.fs.EC)
    print_EC_costing_breakdown(m.fs.EC)

    return m


if __name__ == "__main__":
    m = main()
