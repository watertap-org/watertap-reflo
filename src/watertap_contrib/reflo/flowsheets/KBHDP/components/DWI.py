from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_DWI",
    "init_DWI",
    "add_DWI_costing",
    "add_DWI_scaling",
    "print_DWI_costing_breakdown",
]


def build_DWI(blk, prop_package=None):
    print(f'\n{"=======> BUILDING DEEP WELL INJECTION SYSTEM <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.unit = DeepWellInjection(property_package=prop_package)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )


def set_system_op_conditions(blk):
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_mgd = 2.08 * pyunits.Mgallons / pyunits.day  # estimated from upstream flows

    blk.unit.properties[0].temperature.fix()
    blk.unit.properties[0].pressure.fix()
    blk.unit.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mgd * rho)

    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        blk.unit.properties[0].flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
        blk.unit.properties[0].set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / mass_flow_solute),
            index=("Liq", solute),
        )
    blk.unit.properties[0].set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / blk.unit.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )


def init_DWI(blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize()


def add_DWI_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
        costing_method_arguments={
            "cost_method": "as_opex"
        },  # could be "as_capex" or "blm"
    )
    # fix value for KBHDP case study
    costing_block.deep_well_injection.dwi_lcow.fix(0.58)


def add_DWI_scaling(blk):
    iscale.set_scaling_factor(
        blk.feed.properties[0.0].mass_frac_phase_comp["Liq", "H2O"], 0.1
    )
    iscale.set_scaling_factor(
        blk.unit.properties[0.0].mass_frac_phase_comp["Liq", "H2O"], 0.1
    )
    iscale.set_scaling_factor(blk.unit.properties[0.0].flow_vol_phase["Liq"], 1e-2)


def print_DWI_costing_breakdown(blk, w=25):
    title = "DWI Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"DWI Capital Cost":<{w}s}{f"${blk.unit.costing.capital_cost():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost)}'
    )
    print(
        f'{"DWI Operating Cost":<{w}s}{f"${blk.unit.costing.variable_operating_cost():<{w},.0f}"}{pyunits.get_units(blk.unit.costing.variable_operating_cost)}'
    )
    print("\n")


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }
    # BUG MCAS for RO flowsheets WaterParameterBlock for LTMED causing issues
    m.fs.properties = MCASParameterBlock(
        solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
    )

    m.fs.DWI = FlowsheetBlock(dynamic=False)
    build_DWI(m.fs.DWI, prop_package=m.fs.properties)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def main():
    m = build_system()
    set_system_op_conditions(m.fs.DWI)

    init_DWI(m.fs.DWI)
    add_DWI_costing(m.fs.DWI)
    m.fs.costing.cost_process()

    m.fs.costing.add_LCOW(m.fs.DWI.unit.properties[0].flow_vol)

    results = solve(m)
    assert_optimal_termination(results)
    print_DWI_costing_breakdown(m.fs.DWI)

    return m


if __name__ == "__main__":
    m = main()
