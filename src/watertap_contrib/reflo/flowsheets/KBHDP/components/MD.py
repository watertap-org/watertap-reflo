from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed, StateJunction

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.VAGMD_batch import (
    VAGMDBatchSurrogate,
    get_n_time_points,
)
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_md",
    "set_md_model_options",
    "init_md",
    "report_MD",
    "print_MD_costing_breakdown",
]


def build_system(Qin=4, Cin=12, water_recovery=0.5):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.inlet_flow_rate = pyunits.convert(
        Qin * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = water_recovery

    # Property package
    m.fs.properties = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.disposal.properties[0].flow_vol_phase

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    build_md(m.fs.md, prop_package=m.fs.properties)
    add_connections(m)

    return m


def add_connections(m):

    m.fs.feed_to_md = Arc(source=m.fs.feed.outlet, destination=m.fs.md.feed.inlet)

    m.fs.md_to_product = Arc(
        source=m.fs.md.permeate.outlet, destination=m.fs.product.inlet
    )

    m.fs.md_to_disposal = Arc(
        source=m.fs.md.concentrate.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_md_model_options(m, blk, n_time_points=None):

    m.system_capacity = m.water_recovery * pyunits.convert(
        m.inlet_flow_rate, to_units=pyunits.m**3 / pyunits.day
    )
    m.feed_salinity = pyunits.convert(m.inlet_salinity, to_units=pyunits.g / pyunits.L)

    model_options = {
        "dt": None,
        "system_capacity": value(m.system_capacity),  # m3/day
        "feed_flow_rate": 750,  # L/h
        "evap_inlet_temp": 80,
        "cond_inlet_temp": 20,
        "feed_temp": 25,
        "feed_salinity": value(m.feed_salinity),  # g/L
        "initial_batch_volume": 50,  # L
        "module_type": "AS26C7.2L",
        "cooling_system_type": "closed",
        "cooling_inlet_temp": 25,
        "recovery_ratio": m.water_recovery,
    }

    if n_time_points == None:
        # Calculate the number of periods to reach target recovery rate by solving the system first
        n_time_points = get_n_time_points(
            dt=model_options["dt"],
            feed_flow_rate=model_options["feed_flow_rate"],
            evap_inlet_temp=model_options["evap_inlet_temp"],
            cond_inlet_temp=model_options["cond_inlet_temp"],
            feed_temp=model_options["feed_temp"],
            feed_salinity=model_options["feed_salinity"],
            recovery_ratio=model_options["recovery_ratio"],
            initial_batch_volume=model_options["initial_batch_volume"],
            module_type=model_options["module_type"],
            cooling_system_type=model_options["cooling_system_type"],
            cooling_inlet_temp=model_options["cooling_inlet_temp"],
        )

    blk.model_input = model_options
    blk.n_time_points = n_time_points


def build_md(blk, prop_package=None):

    print(f'\n{"=======> BUILDING MEMBRANE DISTILLATION SYSTEM <=======":^60}\n')
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    # Build a feed, permeate and brine state function for MD
    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = StateJunction(property_package=prop_package)
    blk.concentrate = StateJunction(property_package=prop_package)

    set_md_model_options(m, blk, n_time_points=None)

    blk.unit = VAGMDBatchSurrogate(model_input=blk.model_input)


def init_md(blk):

    m = blk.model()
    blk.feed.initialize()

    # Build connection to permeate state junction
    blk.permeate.properties[0]._flow_vol_phase

    @blk.Constraint(
        doc="Assign the permeate flow rate to its respective state junction"
    )
    def get_permeate_flow(b):
        # num_modules = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.num_modules

        vagmd = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd
        num_modules = pyunits.convert(
            vagmd.system_capacity
            / sum(
                b.unit.mp.get_active_process_blocks()[i].fs.vagmd.permeate_flux
                for i in range(blk.n_time_points)
            )
            * blk.n_time_points
            / b.unit.mp.get_active_process_blocks()[0].fs.vagmd.module_area,
            to_units=pyunits.dimensionless,
        )

        return b.permeate.properties[0].flow_vol_phase["Liq"] == pyunits.convert(
            (
                num_modules
                * b.unit.mp.get_active_process_blocks()[-1].fs.acc_distillate_volume
                / (
                    b.unit.mp.get_active_process_blocks()[-1].fs.dt
                    * (blk.n_time_points - 1)
                )
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    blk.permeate.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)
    blk.permeate.properties[0].pressure.fix(101325)
    blk.permeate.properties[0].temperature.fix(298.15)
    # blk.permeate.initialize()

    # Build connection to concentrate state junction
    blk.concentrate.properties[0]._flow_vol_phase
    blk.concentrate.properties[0]._conc_mass_phase_comp

    @blk.Constraint(
        doc="Assign the concentrate flow rate to its respective state junction"
    )
    def get_concentrate_flow(b):
        num_modules = b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.num_modules

        return b.concentrate.properties[0].flow_vol_phase["Liq"] == pyunits.convert(
            (
                b.unit.mp.get_active_process_blocks()[-1].fs.vagmd.system_capacity
                * (1 - m.water_recovery)
                / m.water_recovery
                #     num_modules
                #     * blk.model_input["initial_batch_volume"]
                #     * pyunits.L
                #     * (1 - b.unit.mp.get_active_process_blocks()[-1].fs.acc_recovery_ratio)
                #     / (b.unit.mp.get_active_process_blocks()[-1].fs.dt * (blk.n_time_points - 1))
            ),
            to_units=pyunits.m**3 / pyunits.s,
        )

    @blk.Constraint(
        doc="Assign the concentrate concentration to its respective state junction"
    )
    def get_concentrate_conc(b):
        return b.concentrate.properties[0].conc_mass_phase_comp[
            "Liq", "TDS"
        ] == pyunits.convert(
            b.unit.mp.get_active_process_blocks()[-1]
            .fs.vagmd.feed_props[0]
            .conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.kg / pyunits.m**3,
        )

    blk.concentrate.properties[0].pressure.fix(101325)
    blk.concentrate.properties[0].temperature.fix(298.15)
    # blk.concentrate.initialize()


def init_system(m):

    print(
        "\n\n-------------------- INITIALIZING MEMBRANE DISTILLATION --------------------\n\n"
    )

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_md)

    init_md(m.fs.md)

    propagate_state(m.fs.md_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.md_to_disposal)
    m.fs.disposal.initialize()


def set_system_op_conditions(m):

    feed_flow_rate = m.fs.md.model_input["feed_flow_rate"]
    feed_salinity = m.fs.md.model_input["feed_salinity"]
    feed_temp = m.fs.md.model_input["feed_temp"]

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): pyunits.convert(
                feed_flow_rate * pyunits.L / pyunits.h,
                to_units=pyunits.m**3 / pyunits.s,
            ),
            ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
            ("temperature", None): feed_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def report_MD(blk, w=30):

    active_blks = blk.unit.mp.get_active_process_blocks()

    title = "MD Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(f'{"Membrane type":<{w}s}{active_blks[-1].fs.vagmd.config.module_type}')
    print(
        f'{"Number of Modules":<{w}s}{value(active_blks[-1].fs.vagmd.num_modules):<{w}.2f}'
    )
    print(
        f'{"System Capacity":<{w}s}{value(active_blks[0].fs.vagmd.system_capacity):<{w},.2f}{pyunits.get_units(active_blks[0].fs.vagmd.system_capacity)}'
    )
    print(
        f'{"Feed Salinity":<{w}s}{value(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<{w}.2f}{pyunits.get_units(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )
    print(
        f'{"MD Period 1 Feed Salinity":<{w}s}{value(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"]):<{w}.2f}{pyunits.get_units(blk.feed.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )
    print(
        f'{"Accumulated Recovery":<{w}s}{value(active_blks[-1].fs.acc_recovery_ratio):<{w}.2f}{pyunits.get_units(active_blks[-1].fs.acc_recovery_ratio)}'
    )

    perm_flow = pyunits.convert(
        blk.permeate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Permeate flow rate":<{w}s}{value(perm_flow):<{w},.2f}{pyunits.get_units(perm_flow)}'
    )

    conc_flow = pyunits.convert(
        blk.concentrate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.day,
    )

    print(
        f'{"Brine flow rate":<{w}s}{value(conc_flow):<{w},.2f}{pyunits.get_units(conc_flow)}'
    )

    print(
        f'{"Brine Concentration":<{w}s}{value(blk.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"]):<{w}.2f}{pyunits.get_units(blk.concentrate.properties[0].conc_mass_phase_comp["Liq","TDS"])}'
    )
    print(
        f'{"STEC":<{w}s}{value(active_blks[-1].fs.specific_energy_consumption_thermal):<{w}.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_thermal)}'
    )

    print(
        f'{"SEC":<{w}s}{value(active_blks[-1].fs.specific_energy_consumption_electric):<{w}.2f}{pyunits.get_units(active_blks[-1].fs.specific_energy_consumption_electric)}'
    )

    print(
        f'{"Thermal Power Required":<{w}s}{value(blk.unit.overall_thermal_power_requirement):<{w}.2f}{pyunits.get_units(blk.unit.overall_thermal_power_requirement)}'
    )

    print(
        f'{"Electric Power Required":<{w}s}{value(blk.unit.overall_elec_power_requirement):<{w}.2f}{pyunits.get_units(blk.unit.overall_elec_power_requirement)}'
    )


def print_MD_costing_breakdown(blk, w=30):
    active_blks = blk.unit.mp.get_active_process_blocks()
    title = "MD Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Capital Cost":<{w}s}{value(active_blks[-1].fs.vagmd.costing.capital_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.capital_cost)}'
    )
    print(
        f'{"Fixed Operating Cost":<{w}s}{value(active_blks[-1].fs.vagmd.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.fixed_operating_cost)}'
    )
    print(
        f'{"Membrane Cost":<{w}s}{value(active_blks[-1].fs.vagmd.costing.module_cost):<20,.2f}{pyunits.get_units(active_blks[-1].fs.vagmd.costing.module_cost)}'
    )


def check_md_flows(m, w=25):

    print("\n")
    print(
        f'Sys Feed Flow Rate: {value(pyunits.convert(m.fs.feed.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'MD Feed Flow Rate: {value(pyunits.convert(m.fs.md.feed.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'Sys Perm Flow Rate: {value(pyunits.convert(m.fs.product.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'MD Perm Flow Rate: {value(pyunits.convert(m.fs.md.permeate.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'Sys Conc Flow Rate: {value(pyunits.convert(m.fs.disposal.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'MD Conc Flow Rate: {value(pyunits.convert(m.fs.md.concentrate.properties[0].flow_vol_phase["Liq"], pyunits.m**3/pyunits.day)):<{w}.2f} m3/day'
    )
    print(
        f'Calculated Recovery: {value(m.fs.md.permeate.properties[0].flow_vol_phase["Liq"] / (m.fs.md.permeate.properties[0].flow_vol_phase["Liq"] + m.fs.md.concentrate.properties[0].flow_vol_phase["Liq"])):<{w}.2f}'
    )


def add_md_costing(blk, costing_block=None):

    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.add_costing_module(costing_block)


def main():
    m = build_system()
    set_system_op_conditions(m)
    init_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    add_md_costing(m.fs.md)

    m.fs.costing.capital_recovery_factor.fix(0.08764)
    m.fs.costing.wacc.unfix()

    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])

    assert degrees_of_freedom(m) == 0

    results = solve(m)
    assert_optimal_termination(results)

    report_MD(m.fs.md)
    print_MD_costing_breakdown(m.fs.md)

    return m


if __name__ == "__main__":
    m = main()
