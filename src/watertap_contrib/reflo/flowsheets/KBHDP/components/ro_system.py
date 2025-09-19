from pyomo.environ import (
    ConcreteModel,
    Param,
    Var,
    Constraint,
    NonNegativeReals,
    TransformationFactory,
    RangeSet,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Feed, StateJunction
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve, calc_scale

__all__ = [
    "build_ro",
    "build_ro_stage",
    "init_ro_system",
    "init_ro_stage",
    "set_ro_system_operating_conditions",
    "add_ro_costing",
    "add_ro_scaling",
    "display_RO_flow_table",
    "report_RO",
    "print_RO_costing_breakdown",
]


def build_ro(blk, number_of_stages=1, prop_package=None):
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.RO_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)
    blk.numberOfStages = Param(initialize=number_of_stages)
    blk.Stages = RangeSet(blk.numberOfStages)
    blk.booster_pumps = False

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(number_of_stages - 1)

    blk.primary_mixer = Mixer(
        property_package=prop_package,
        has_holdup=False,
        num_inlets=number_of_stages,
        momentum_mixing_type=MomentumMixingType.minimize_and_equality,
    )

    blk.stage = FlowsheetBlock(RangeSet(number_of_stages), dynamic=False)

    for _, stage in blk.stage.items():
        if stage.index() > 1:
            build_ro_stage(
                stage, booster_pump=blk.booster_pumps, prop_package=prop_package
            )
        else:
            build_ro_stage(stage, prop_package=prop_package)

    blk.ro_feed_to_ro = Arc(
        source=blk.feed.outlet,
        destination=blk.stage[1].feed.inlet,
    )

    blk.stage_retentate_to_next_stage = Arc(
        blk.NonFinalStages,
        rule=lambda blk, n: {
            "source": blk.stage[n].retentate.outlet,
            "destination": blk.stage[n + 1].feed.inlet,
        },
    )

    blk.stage_permeate_to_mixer = Arc(
        blk.Stages,
        rule=lambda blk, n: {
            "source": blk.stage[n].permeate.outlet,
            "destination": getattr(blk.primary_mixer, "inlet_" + str(n)),
        },
    )

    blk.primary_mixer_to_product = Arc(
        source=blk.primary_mixer.outlet,
        destination=blk.product.inlet,
    )

    blk.last_stage_retentate_to_ro_retentate = Arc(
        source=blk.stage[number_of_stages].retentate.outlet,
        destination=blk.disposal.inlet,
    )

    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp


def build_ro_stage(blk, booster_pump=False, prop_package=None):

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.RO_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = StateJunction(property_package=prop_package)
    blk.retentate = StateJunction(property_package=prop_package)
    blk.has_booster_pump = booster_pump

    if booster_pump:
        blk.booster_pump = Pump(property_package=prop_package)
        blk.booster_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

    blk.module = ReverseOsmosis1D(
        property_package=prop_package,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        has_full_reporting=True,
    )

    relax_bounds_for_low_salinity_waters(blk.module)

    if booster_pump:
        blk.stage_feed_to_booster_pump = Arc(
            source=blk.feed.outlet,
            destination=blk.booster_pump.inlet,
        )
        blk.stage_booster_pump_to_module = Arc(
            source=blk.booster_pump.outlet,
            destination=blk.module.inlet,
        )
    else:
        blk.stage_feed_to_module = Arc(
            source=blk.feed.outlet,
            destination=blk.module.inlet,
        )

    blk.stage_module_to_permeate = Arc(
        source=blk.module.permeate,
        destination=blk.permeate.inlet,
    )

    blk.stage_module_to_retentate = Arc(
        source=blk.module.retentate,
        destination=blk.retentate.inlet,
    )

    blk.feed.properties[0].conc_mass_phase_comp
    blk.permeate.properties[0].conc_mass_phase_comp
    blk.retentate.properties[0].conc_mass_phase_comp


def init_system(m):

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_ro)

    init_ro_system(m.fs.ro)
    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_disposal)

    m.fs.product.initialize()
    m.fs.disposal.initialize()


def init_ro_system(blk):

    print("\n\n-------------------- INITIALIZING RO SYSTEM --------------------\n\n")

    blk.feed.initialize()

    propagate_state(blk.ro_feed_to_ro)

    for stage in blk.stage.values():
        init_ro_stage(stage)
        if stage.index() < blk.numberOfStages:
            propagate_state(blk.stage_retentate_to_next_stage[stage.index()])
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])
        else:
            propagate_state(blk.last_stage_retentate_to_ro_retentate)
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])

    blk.disposal.initialize()
    blk.primary_mixer.initialize()
    propagate_state(blk.primary_mixer_to_product)
    blk.product.initialize()
    print(
        "\n\n-------------------- RO INITIALIZATION COMPLETE --------------------\n\n"
    )
    display_RO_flow_table(blk)


def init_ro_stage(stage):

    if stage.has_booster_pump:
        stage.feed.initialize()
        propagate_state(stage.stage_feed_to_booster_pump)
        stage.booster_pump.initialize()
        propagate_state(stage.stage_booster_pump_to_module)
    else:
        stage.feed.initialize()
        propagate_state(stage.stage_feed_to_module)

    stage.module.initialize()

    propagate_state(stage.stage_module_to_retentate)
    propagate_state(stage.stage_module_to_permeate)

    stage.permeate.initialize()
    stage.retentate.initialize()


def set_operating_conditions(
    m,
    feed_flow_mass=171.763,
    feed_salinity=11.33,
    ro_pressure=30e5,
    water_recovery=None,
):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_salinity = Var(
        initialize=35,
        bounds=(0, 2000),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.unfix()

    feed_temperature = 273.15 + 20
    supply_pressure = ro_pressure

    # initialize feed
    m.fs.feed.temperature[0].fix(feed_temperature)
    m.fs.feed_flow_mass.fix(feed_flow_mass)
    m.fs.feed_salinity.fix(feed_salinity)
    m.fs.feed.pressure[0].fix(supply_pressure)

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    m.fs.feed_flow_constraint = Constraint(
        expr=m.fs.feed_flow_mass == m.fs.perm_flow_mass / m.fs.water_recovery
    )

    m.fs.nacl_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"] * 1000
        == m.fs.feed_flow_mass * m.fs.feed_salinity
    )

    m.fs.h2o_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
        == m.fs.feed_flow_mass * (1 - m.fs.feed_salinity / 1000)
    )

    m.fs.feed.properties[0].flow_vol_phase["Liq"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)
    )

    iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)
    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    )


def relax_bounds_for_low_salinity_waters(blk):
    blk.feed_side.cp_modulus.setub(5)
    for e in blk.feed_side.K:
        blk.feed_side.K[e].setub(0.01)
        # blk.feed_side.K[e].setlb(1e-7)

    for e in blk.feed_side.cp_modulus:
        blk.feed_side.cp_modulus[e].setlb(1e-5)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.recovery_mass_phase_comp[e].setlb(1e-9)
            blk.recovery_mass_phase_comp[e].setub(1e-1)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.flux_mass_phase_comp[e].setlb(1e-9)
            blk.flux_mass_phase_comp[e].setub(1e-1)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "H2O":
            blk.recovery_mass_phase_comp[e].setlb(1e-4)
            blk.recovery_mass_phase_comp[e].setub(0.999)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "H2O":
            blk.flux_mass_phase_comp[e].setlb(1e-5)
            blk.flux_mass_phase_comp[e].setub(0.999)


# # TODO: move to utils
# def calc_scale(value):
#     return math.floor(math.log(value, 10))


def add_ro_scaling(blk):
    print("Setting RO scaling")

    for _s, stage in blk.stage.items():
        module = stage.module
        iscale.set_scaling_factor(module.area, 1e5)
        iscale.set_scaling_factor(module.feed_side.area, 1)
        iscale.set_scaling_factor(module.width, 1e4)
        iscale.set_scaling_factor(module.length, 1e1)
        iscale.set_scaling_factor(module.feed_side.N_Sh_comp, 1e-4)

        for e in module.feed_side.properties:
            iscale.set_scaling_factor(
                module.feed_side.properties[e].flow_mass_phase_comp["Liq", "NaCl"], 1
            )
            iscale.set_scaling_factor(
                module.feed_side.properties[e].dens_mass_phase["Liq"], 1
            )
            iscale.set_scaling_factor(
                module.feed_side.properties[e].mass_frac_phase_comp["Liq", "NaCl"], 1e1
            )

        for temp_stream in [
            module.eq_permeate_isothermal,
            module.feed_side.eq_equal_temp_interface,
            module.feed_side.eq_feed_isothermal,
            module.eq_permeate_outlet_isothermal,
        ]:
            for e in temp_stream:
                iscale.constraint_scaling_transform(temp_stream[e], 1e-2)
            for pressure_stream in [
                module.eq_permeate_outlet_isobaric,
                module.feed_side.eq_equal_pressure_interface,
            ]:
                for e in pressure_stream:
                    iscale.constraint_scaling_transform(pressure_stream[e], 1e-5)
            for e in module.eq_pressure_drop:
                iscale.constraint_scaling_transform(module.eq_pressure_drop[e], 1e-7)

        for e in module.feed_side.eq_N_Sh_comp:
            if e[-1] == "NaCl":
                iscale.constraint_scaling_transform(
                    module.feed_side.eq_N_Sh_comp[e], 1e-3
                )

        for e in module.feed_side.eq_friction_factor:
            iscale.constraint_scaling_transform(
                module.feed_side.eq_friction_factor[e], 1e-2
            )
        for e in module.feed_side.eq_dP_dx:
            iscale.constraint_scaling_transform(module.feed_side.eq_dP_dx[e], 1e-2)

        iscale.set_scaling_factor(module.mixed_permeate[0.0].dens_mass_phase["Liq"], 1)
        iscale.set_scaling_factor(module.mixed_permeate[0.0].flow_vol_phase["Liq"], 100)
        iscale.constraint_scaling_transform(
            module.mixed_permeate[0.0].eq_flow_vol_phase["Liq"], 100
        )

        for e in module.recovery_mass_phase_comp:
            if e[-1] == "H2O":
                iscale.set_scaling_factor(module.recovery_mass_phase_comp, 1e1)


def set_ro_system_operating_conditions(blk, mem_area=10000):
    print(
        "\n\n-------------------- SETTING RO OPERATING CONDITIONS --------------------\n\n"
    )
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.95  # spacer porosity in membrane stage [-]
    area = mem_area  # membrane area [m^2]
    length = 7  # effective membrane width [m]
    pressure_atm = 101325  # atmospheric pressure [Pa]
    pump_efi = 0.8  # pump efficiency [-]

    for idx, stage in blk.stage.items():
        stage.module.A_comp.fix(mem_A)
        stage.module.B_comp.fix(mem_B)
        stage.module.area.fix(area / idx)
        stage.module.feed_side.velocity[0, 0].fix(0.35)
        stage.module.width.setub(20000)
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        stage.module.feed_side.channel_height.fix(height)
        stage.module.feed_side.spacer_porosity.fix(spacer_porosity)
        stage.module.feed_side.friction_factor_darcy.setub(50)

        for e in stage.module.flux_mass_phase_comp:
            if e[-1] == "H2O":
                stage.module.flux_mass_phase_comp[e].setlb(1e-5)
                stage.module.flux_mass_phase_comp[e].setub(0.99)

    blk.total_membrane_area = Var(
        initialize=10000,
        domain=NonNegativeReals,
        units=pyunits.m**2,
        doc="Total RO System Membrane Area",
    )

    blk.eq_total_membrane_area = Constraint(
        expr=blk.total_membrane_area
        == sum([stage.module.area for idx, stage in blk.stage.items()])
    )

    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")


def add_ro_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    for stage in blk.stage.values():
        stage.module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing_block
        )


def display_RO_flow_table(blk, w=25):
    title = "RO System Flow Table"
    side = int(((5 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{"Unit":<{w}s}{"Mass Flow Water (kg/s)":<{w}s}{"Pressure (bar)":<{w}s}{"Mass Flow NaCl (kg/s)":<{w}s}{"Conc. (g/L)":<{w}s}'
    )
    print(f"{'-' * (5 * w)}")
    print(
        f'{"Feed":<{w}s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<{w}.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Product":<{w}s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<{w}.1f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Disposal":<{w}s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<{w}.1f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )

    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Feed":<{w}s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<{w}.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Permeate":<{w}s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<{w}.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Retentate":<{w}s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<{w}.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
        )
    print("\n\n")


def report_RO(blk, w=25):
    m = blk.model()
    title = "RO Operating Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    try:
        print(f'{"Recovery":<{w}s}{value(100*m.fs.water_recovery):<{w}.1f}{"%"}')
    except:
        print(
            f'{"Recovery":<{w}s}{value(100*blk.parent_block().water_recovery):<{w}.1f}{"%"}'
        )
    print(
        f'{"RO Operating Pressure":<{w}s}{value(pyunits.convert(blk.feed.properties[0].pressure, to_units=pyunits.bar)):<{w}.1f}{"bar":<{w}}'
    )
    print(
        f'{"RO Membrane Area":<{w}s}{value(blk.stage[1].module.area):<{w}.1f}{"m^2":<{w}}'
    )


def print_RO_costing_breakdown(blk, w=25):
    title = "RO Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"RO Capital Cost":<{w}s}{f"${value(blk.stage[1].module.costing.capital_cost):<{w},.0f}{pyunits.get_units(blk.stage[1].module.costing.capital_cost)}"}'
    )
    print(
        f'{"RO Operating Cost":<{w}s}{f"${value(blk.stage[1].module.costing.fixed_operating_cost):<{w},.0f}{pyunits.get_units(blk.stage[1].module.costing.fixed_operating_cost)}"}'
    )
    print("\n\n")


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.RO_properties = NaClParameterBlock()
    m.fs.costing = TreatmentCosting()

    m.fs.feed = Feed(property_package=m.fs.RO_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.disposal = Product(property_package=m.fs.RO_properties)

    m.fs.ro = FlowsheetBlock()
    build_ro(m.fs.ro)

    m.fs.feed_to_ro = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.ro.feed.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.ro.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.ro_to_disposal = Arc(
        source=m.fs.ro.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def main():

    m = build_system()
    set_operating_conditions(m)
    set_ro_system_operating_conditions(m.fs.ro)
    add_ro_scaling(m.fs.ro)
    iscale.calculate_scaling_factors(m)
    init_system(m)
    results = solve(m)
    assert_optimal_termination(results)
    add_ro_costing(m.fs.ro)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    results = solve(m)
    assert_optimal_termination(results)

    display_RO_flow_table(m.fs.ro)
    report_RO(m.fs.ro)
    print_RO_costing_breakdown(m.fs.ro)
    return m


if __name__ == "__main__":
    m = main()
