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
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import solve

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


def add_UF_costing(blk, costing_blk=None):
    if costing_blk is None:
        m = blk.model()
        costing_blk = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_blk,
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


def report_UF(blk):
    print(f"\n\n-------------------- UF Report --------------------\n")
    print("\n")
    print(
        f'{"Inlet Flow Volume":<30s}{value(blk.feed.properties[0.0].flow_vol):<10.3f}{pyunits.get_units(blk.feed.properties[0.0].flow_vol)}'
    )
    print(f'{"UF Performance:":<30s}')
    print(
        f'{"    Recovery":<30s}{100*blk.unit.recovery_frac_mass_H2O[0.0].value:<10.1f}{"%"}'
    )
    print(
        f'{"    TDS Removal":<30s}{100*blk.unit.removal_frac_mass_comp[0.0,"tds"].value:<10.1f}{"%"}'
    )
    print(
        f'{"    TSS Removal":<30s}{100*blk.unit.removal_frac_mass_comp[0.0,"tss"].value:<10.1f}{"%"}'
    )
    print(
        f'{"    Energy Consumption":<30s}{blk.unit.electricity[0.0].value:<10.3f}{pyunits.get_units(blk.unit.electricity[0.0])}'
    )
    print(
        f'{"    Specific Energy Cons.":<30s}{value(blk.unit.energy_electric_flow_vol_inlet):<10.3f}{pyunits.get_units(blk.unit.energy_electric_flow_vol_inlet)}'
    )


def print_UF_costing_breakdown(
    blk,
):
    print(f"\n\n-------------------- UF Costing Breakdown --------------------\n")
    print(f'{"UF Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')


def main():

    m = build_system()
    set_UF_op_conditions(m.fs.UF)
    set_system_conditions(m.fs.UF)
    add_UF_costing(m.fs.UF)
    m.fs.costing.cost_process()
    init_UF(m.fs.UF)

    results = solve(m)
    assert_optimal_termination(results)

    report_UF(m.fs.UF)
    print_UF_costing_breakdown(m.fs.UF)


if __name__ == "__main__":
    main()
