from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import StateJunction

from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.zero_order.ultra_filtration_zo import UltraFiltrationZO

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_UF",
    "init_UF",
    "set_UF_op_conditions",
    "add_UF_costing",
    "add_UF_scaling",
    "report_UF",
    "print_UF_costing_breakdown",
]


def build_UF(blk, prop_package=None):

    print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.UF_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)
    blk.unit = UltraFiltrationZO(property_package=prop_package, database=m.db)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.treated,
        destination=blk.product.inlet,
    )


def init_UF(blk):

    print(
        "\n\n-------------------- INITIALIZING ULTRAFILTRATION --------------------\n\n"
    )
    print(f"UF Degrees of Freedom: {degrees_of_freedom(blk)}")

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()
    propagate_state(blk.unit_to_disposal)
    propagate_state(blk.unit_to_product)

    blk.product.initialize()
    blk.disposal.initialize()


def set_UF_op_conditions(blk):
    print(f"UF Degrees of Freedom: {degrees_of_freedom(blk)}")
    blk.unit.recovery_frac_mass_H2O.fix(0.99)
    blk.unit.removal_frac_mass_comp[0, "tds"].fix(1e-3)
    blk.unit.removal_frac_mass_comp[0, "tss"].fix(0.9)
    blk.unit.energy_electric_flow_vol_inlet.fix(0.05)


def set_system_conditions(blk):
    blk.feed.properties[0.0].flow_mass_comp["H2O"].fix(171.37)
    blk.feed.properties[0.0].flow_mass_comp["tds"].fix(1.96)
    blk.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)


def add_UF_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
    )


def add_UF_scaling(blk):
    iscale.set_scaling_factor(blk.disposal.properties[0.0].flow_mass_comp["tds"], 1e3)
    iscale.set_scaling_factor(
        blk.unit.properties_byproduct[0.0].flow_mass_comp["tds"], 1e3
    )


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    m.db = Database()
    m.fs.params = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.UF = FlowsheetBlock(dynamic=False)
    build_UF(m.fs.UF, prop_package=m.fs.params)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def report_UF(blk, w=25):
    title = "UF Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Inlet Flow Rate":<{w}s}{value(pyunits.convert(blk.feed.properties[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w}.2f}MGD'
    )
    print(
        f'{"Recovery":<{w}s}{100*blk.unit.recovery_frac_mass_H2O[0].value:<{w}.1f}{"%"}'
    )
    print(
        f'{"TDS Removal":<{w}s}{100*blk.unit.removal_frac_mass_comp[0, "tds"].value:<{w}.1f}{"%"}'
    )
    print(
        f'{"TSS Removal":<{w}s}{100*blk.unit.removal_frac_mass_comp[0, "tss"].value:<{w}.1f}{"%"}'
    )
    print(
        f'{"Energy Consumption":<{w}s}{blk.unit.electricity[0].value:<{w}.3f}{pyunits.get_units(blk.unit.electricity[0])}'
    )
    print(
        f'{"Specific Energy Cons.":<{w}s}{value(blk.unit.energy_electric_flow_vol_inlet):<{w}.3f}{pyunits.get_units(blk.unit.energy_electric_flow_vol_inlet)}'
    )


def print_UF_costing_breakdown(blk, w=25):
    title = "UF Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"UF Capital Cost":<{w}s}{f"{blk.unit.costing.capital_cost():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost)}'
    )
    print("\n\n")


def main():

    m = build_system()
    set_UF_op_conditions(m.fs.UF)
    set_system_conditions(m.fs.UF)
    add_UF_costing(m.fs.UF)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.UF.product.properties[0.0].flow_vol)
    init_UF(m.fs.UF)

    results = solve(m)
    assert_optimal_termination(results)

    report_UF(m.fs.UF)
    print_UF_costing_breakdown(m.fs.UF)

    return m


if __name__ == "__main__":
    m = main()
