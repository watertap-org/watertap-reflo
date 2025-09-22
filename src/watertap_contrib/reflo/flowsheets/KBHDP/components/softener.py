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
from idaes.core import MaterialFlowBasis
from idaes.models.unit_models import Product, Feed, StateJunction

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.unit_models.chemical_softening import (
    ChemicalSoftening,
)
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_softener",
    "set_softener_op_conditions",
    "add_softener_costing",
    "report_softener",
    "print_softening_costing_breakdown",
    "init_softener",
]


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()
    m.fs.properties = MCASParameterBlock(
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Mg_2+",
            "SiO2",
            "Na_+",
            "Cl_-",
            "K_+",
            "SO4_2-",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.softener = FlowsheetBlock(dynamic=False)

    build_softener(m.fs.softener)

    m.fs.feed_to_softener = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.softener.feed.inlet,
    )

    m.fs.softener_to_product = Arc(
        source=m.fs.softener.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.softener_to_disposal = Arc(
        source=m.fs.softener.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def build_softener(blk, prop_package=None):
    print(f'\n{"=======> BUILDING CHEMICAL SOFTENING SYSTEM <=======":^60}\n')
    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = ChemicalSoftening(
        property_package=prop_package,
        silica_removal=False,
        softening_procedure_type="excess_lime_soda",
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )
    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )
    blk.unit_to_disposal = Arc(
        source=blk.unit.waste,
        destination=blk.disposal.inlet,
    )


def set_system_operating_conditions(m, Qin=4):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    rho = 1000 * pyunits.kg / pyunits.m**3
    Qin = Qin * pyunits.Mgal / pyunits.day
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)

    flow_mass_phase_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)

    inlet_dict = {
        "Ca_2+": 0.61 * pyunits.kg / pyunits.m**3,
        "Mg_2+": 0.16 * pyunits.kg / pyunits.m**3,
        "Alkalinity_2-": 0.0821 * pyunits.kg / pyunits.m**3,
        "SiO2": 0.13 * pyunits.kg / pyunits.m**3,
        "Cl_-": 5.5 * pyunits.kg / pyunits.m**3,
        "Na_+": 5.5 * pyunits.kg / pyunits.m**3,
        "K_+": 0.016 * pyunits.kg / pyunits.m**3,
        "SO4_2-": 0.23 * pyunits.kg / pyunits.m**3,
    }
    calc_state_dict = {
        ("flow_vol_phase", "Liq"): value(flow_in),
        ("pressure", None): 101325,
        ("temperature", None): 298,
    }

    for solute, solute_conc in inlet_dict.items():
        calc_state_dict[("conc_mass_phase_comp", ("Liq", solute))] = solute_conc
        flow_mass_solute = pyunits.convert(
            flow_in * solute_conc, to_units=pyunits.kg / pyunits.s
        )
        sf = 1 / value(flow_mass_solute)
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].set_value(
            flow_mass_solute
        )
        m.fs.softener.unit.properties_in[0].flow_mass_phase_comp[
            "Liq", solute
        ].set_value(flow_mass_solute)
        m.fs.softener.unit.properties_in[0].conc_mass_phase_comp[
            "Liq", solute
        ].set_value(solute_conc)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            sf,
            index=("Liq", solute),
        )
        m.fs.properties.set_default_scaling(
            "conc_mass_phase_comp",
            1 / solute_conc(),
            index=("Liq", solute),
        )
        m.fs.properties.set_default_scaling(
            "mass_frac_phase_comp",
            1 / value(flow_mass_solute / flow_mass_phase_water),
            index=("Liq", solute),
        )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water),
        index=("Liq", "H2O"),
    )
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_mol_phase_comp
    m.fs.softener.unit.properties_in[0].flow_mol_phase_comp

    iscale.calculate_scaling_factors(m)

    m.fs.feed.properties.calculate_state(var_args=calc_state_dict, hold_state=True)
    m.fs.softener.unit.properties_in.calculate_state(
        var_args=calc_state_dict, hold_state=False
    )


def set_softener_op_conditions(
    blk, ca_effluent=0.01, mg_effluent=0.01, non_important_removals=0.01
):
    print(
        "\n\n-------------------- SETTING SOFTENER REMOVAL EFFICIENCY --------------------\n\n"
    )

    blk.unit.ca_eff_target.fix(ca_effluent)
    blk.unit.mg_eff_target.fix(mg_effluent)
    blk.unit.removal_efficiency.fix(non_important_removals)
    blk.unit.frac_mass_water_recovery.fix(0.99)

    blk.unit.retention_time_mixer.fix(0.4)
    blk.unit.retention_time_floc.fix(25)
    blk.unit.retention_time_sed.fix(120)
    blk.unit.retention_time_recarb.fix(20)
    blk.unit.number_mixers.value = 2
    blk.unit.number_floc.value = 4

    blk.unit.CO2_CaCO3.fix(0.063)
    blk.unit.MgCl2_dosing.fix(0)
    blk.unit.vel_gradient_mix.fix(300)
    blk.unit.vel_gradient_floc.fix(50)


def add_softener_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
    )


def init_system(m):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_softener)

    init_softener(m.fs.softener)

    propagate_state(m.fs.softener_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.softener_to_disposal)
    m.fs.disposal.initialize()


def init_softener(blk):

    print(
        "\n\n-------------------- INITIALIZING WATER SOFTENER --------------------\n\n"
    )
    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize()
    propagate_state(blk.unit_to_product)
    propagate_state(blk.unit_to_disposal)
    blk.product.initialize()
    blk.disposal.initialize()


def report_softener(blk, w=25):

    unit = blk.unit
    comps = unit.config.property_package.solute_set
    title = "Softener Report"
    side = int(((4 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print("Stream Table:\n")

    print(
        f'{"Flow In":<{w}s}{value(pyunits.convert(unit.properties_in[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day)):<{w}.3f}{pyunits.get_units(pyunits.convert(unit.properties_in[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day))}'
    )
    print(
        f'{"Flow Out":<{w}s}{value(pyunits.convert(unit.properties_out[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day)):<{w}.3f}{pyunits.get_units(pyunits.convert(unit.properties_out[0.0].flow_vol_phase["Liq"], to_units=pyunits.Mgal / pyunits.day))}'
    )
    print(
        f'{"Inlet TDS":<{w}s}{sum(value(unit.properties_in[0].conc_mass_phase_comp["Liq", i]) for i in comps):<{w}.3f}{"g/L"}'
    )
    print(
        f'{"Outlet TDS":<{w}s}{sum(value(unit.properties_out[0].conc_mass_phase_comp["Liq", i]) for i in comps):<{w}.3f}{"g/L"}\n'
    )

    print(
        f'{"Component":<{w}s}{"Conc. In (g/L)":<{w}s}{"Conc. Out (g/L)":<{w}s}{"Removal (%)":<{w}s}'
    )
    print(f"{'-' * (4 * w)}")

    for solute in comps:
        conc_in = unit.properties_in[0].conc_mass_phase_comp["Liq", solute].value
        conc_out = unit.properties_out[0].conc_mass_phase_comp["Liq", solute].value
        removal = (conc_in - conc_out) / conc_in * 100
        print(f"{solute:{w}s}{conc_in:<{w}.3f}{conc_out:<{w}.3f}{removal:<{w}.1f}")

    print("\nDosing Details:\n")
    print(f'{"Chemical":<{w}s}{"Dose":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"CaO Dose":<{w}s}{unit.CaO_dosing.value:<{w}.3f}{pyunits.get_units(unit.CaO_dosing)}'
    )
    print(
        f'{"MgCl Dose":<{w}s}{value(unit.MgCl2_dosing):<{w}.3f}{pyunits.get_units(unit.MgCl2_dosing)}'
    )
    print(
        f'{"Na2CO3 Dose":<{w}s}{value(unit.Na2CO3_dosing):<{w}.3f}{pyunits.get_units(unit.Na2CO3_dosing)}'
    )
    print(
        f'{"Excess CaO":<{w}s}{unit.excess_CaO.value:<{w}.3f}{pyunits.get_units(unit.excess_CaO)}'
    )
    print(
        f'{"CO2 CaCO3":<{w}s}{unit.CO2_CaCO3.value:<{w}.3f}{pyunits.get_units(unit.CO2_CaCO3)}'
    )
    print(
        f'{"Sludge Produced":<{w}s}{unit.sludge_prod.value:<{w}.3f}{pyunits.get_units(unit.sludge_prod)}'
    )


def print_softening_costing_breakdown(blk, w=25):
    unit = blk.unit

    title = "Softener Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Capital Cost":<{w}s}{f"${unit.costing.capital_cost.value:<{w - 1},.0f}"}{pyunits.get_units(unit.costing.capital_cost)}'
    )
    print(
        f'{"Operating Cost":<{w}s}{f"${unit.costing.fixed_operating_cost.value:<{w - 1},.0f}"}{pyunits.get_units(unit.costing.fixed_operating_cost)}'
    )
    print(
        f'{"Electric Power Reqd":<{w}s}{f"{unit.costing.electricity_flow.value:<{w},.4f}"}{pyunits.get_units(unit.costing.electricity_flow)}'
    )

    print("\n\n")


def main():
    m = build_system()
    set_system_operating_conditions(m)
    set_softener_op_conditions(m.fs.softener)
    add_softener_costing(m.fs.softener)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.softener.unit.properties_out[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.softener.unit.properties_out[0].flow_vol_phase["Liq"], name="SEC"
    )
    assert degrees_of_freedom(m) == 0

    init_system(m)

    results = solve(m, tee=False)
    assert_optimal_termination(results)
    report_softener(m.fs.softener)
    print_softening_costing_breakdown(m.fs.softener)

    return m


if __name__ == "__main__":
    m = main()
