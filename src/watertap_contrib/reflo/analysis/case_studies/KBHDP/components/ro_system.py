import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

# from idaes.core.solvers import get_solver
from watertap.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

# from analysisWaterTAP.utils.flowsheet_utils import *
from watertap.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)

# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

__all__ = [
    "build_ro",
    "build_ro_stage",
    "init_ro_system",
    "init_ro_stage",
    "set_ro_system_operating_conditions",
    "add_ro_costing",
    "add_ro_scaling",
    "display_ro_system_build",
    "display_dof_breakdown",
    "display_flow_table",
    "report_RO",
    "print_RO_costing_breakdown",
]


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


def _initialize(blk, verbose=False):
    if verbose:
        print("\n")
        print(
            f"{blk.name:<30s}{f'Degrees of Freedom at Initialization = {degrees_of_freedom(blk):<10.0f}'}"
        )
        print("\n")
    try:
        blk.initialize()
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")

        blk.report()
        print_infeasible_bounds(blk)
        print_close_to_bounds(blk)
        assert False


def print_RO_op_pressure_est(blk):
    solver = get_solver()
    operating_pressure = calculate_operating_pressure(
        feed_state_block=blk.feed.properties[0],
        over_pressure=0.15,
        water_recovery=0.8,
        NaCl_passage=0.01,
        solver=solver,
    )

    operating_pressure_psi = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.psi
    )()
    operating_pressure_bar = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.bar
    )()
    print(
        f"\nOperating Pressure Estimate = {round(operating_pressure_bar, 2)} bar = {round(operating_pressure_psi, 2)} psi\n"
    )


_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")


def build_ro(m, blk, number_of_stages=1, prop_package=None) -> None:
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.RO_properties)
    blk.product = StateJunction(property_package=m.fs.RO_properties)
    blk.disposal = StateJunction(property_package=m.fs.RO_properties)
    blk.numberOfStages = Param(initialize=number_of_stages)
    blk.Stages = RangeSet(blk.numberOfStages)
    blk.booster_pumps = False

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(number_of_stages - 1)

    blk.primary_mixer = Mixer(
        property_package=m.fs.RO_properties,
        has_holdup=False,
        num_inlets=number_of_stages,
        momentum_mixing_type=MomentumMixingType.minimize_and_equality,
    )

    blk.stage = FlowsheetBlock(RangeSet(number_of_stages), dynamic=False)

    for idx, stage in blk.stage.items():
        if stage.index() > 1:
            build_ro_stage(m, stage, booster_pump=blk.booster_pumps)
        else:
            build_ro_stage(m, stage)

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


def build_ro_stage(m, blk, booster_pump=False):
    # Define IO
    blk.feed = StateJunction(property_package=m.fs.RO_properties)
    blk.permeate = StateJunction(property_package=m.fs.RO_properties)
    blk.retentate = StateJunction(property_package=m.fs.RO_properties)
    blk.has_booster_pump = booster_pump

    if booster_pump:
        blk.booster_pump = Pump(property_package=m.fs.RO_properties)
        blk.booster_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
        )

    blk.module = ReverseOsmosis1D(
        property_package=m.fs.RO_properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        has_full_reporting=True,
    )

    relax_bounds_for_low_salinity_waters(m, blk.module)

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


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")

    assert_no_degrees_of_freedom(m)
    _initialize(m.fs.feed)
    print(m.fs.feed.report())
    propagate_state(m.fs.feed_to_ro)

    init_ro_system(m, m.fs.ro)
    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_disposal)

    _initialize(m.fs.product)
    _initialize(m.fs.disposal)


def init_ro_system(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING RO SYSTEM --------------------\n\n")

    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    print("\n\n")

    _initialize(blk.feed)
    propagate_state(blk.ro_feed_to_ro)

    for stage in blk.stage.values():
        init_ro_stage(m, stage, solver=solver)
        if stage.index() < blk.numberOfStages:
            propagate_state(blk.stage_retentate_to_next_stage[stage.index()])
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])
        else:
            propagate_state(blk.last_stage_retentate_to_ro_retentate)
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])

    _initialize(blk.disposal)
    _initialize(blk.primary_mixer)
    propagate_state(blk.primary_mixer_to_product)
    _initialize(blk.product)
    print(
        "\n\n-------------------- RO INITIALIZATION COMPLETE --------------------\n\n"
    )
    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    # for stage in blk.stage.values():
    #     print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
    print("\n\n")
    display_flow_table(blk)
    print(blk.report())
    print(
        f'RO Recovery: {100 * (value(blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]) / value(blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"])):<5.2f}%'
    )
    # print(
    #     f'{"Average Flux":<30s}{value(blk.stage[1].module.flux_vol_phase_avg[0, "Liq"]):<10.2f}{pyunits.get_units(blk.stage[1].module.flux_vol_phase_avg[0, "Liq"])}'
    # )


def init_ro_stage(m, stage, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    if stage.has_booster_pump:
        _initialize(stage.feed)
        propagate_state(stage.stage_feed_to_booster_pump)
        _initialize(stage.booster_pump)
        propagate_state(stage.stage_booster_pump_to_module)
    else:
        _initialize(stage.feed)
        propagate_state(stage.stage_feed_to_module)

    # display_inlet_conditions(stage)

    _initialize(stage.module)

    propagate_state(stage.stage_module_to_retentate)
    propagate_state(stage.stage_module_to_permeate)

    _initialize(stage.permeate)
    _initialize(stage.retentate)


def set_operating_conditions(
    m, Qin=None, Qout=None, Cin=None, water_recovery=None, ro_pressure=25e5
):
    print(
        "\n\n-------------------- SETTING SYSTEM OPERATING CONDITIONS --------------------\n\n"
    )
    if Cin is None:
        Cin = 35

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
    pressure_atm = 101325
    supply_pressure = ro_pressure

    #     # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    if Qin is not None:
        m.fs.feed_flow_mass.fix(Qin)

    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    m.fs.feed_salinity.fix(Cin)
    iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

    m.fs.feed_flow_constraint = Constraint(
        expr=m.fs.feed_flow_mass == m.fs.perm_flow_mass / m.fs.water_recovery
    )
    iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)

    m.fs.nacl_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"] * 1000
        == m.fs.feed_flow_mass * m.fs.feed_salinity
    )

    m.fs.h2o_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
        == m.fs.feed_flow_mass * (1 - m.fs.feed_salinity / 1000)
    )

    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    # m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)
    )

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    )

    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # print(f"RO Degrees of Freedom: {degrees_of_freedom(m.fs.ro)}")


def relax_bounds_for_low_salinity_waters(m, blk):
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


def calc_scale(value):
    return math.floor(math.log(value, 10))


def add_ro_scaling(m, blk):
    print("Setting RO scaling")
    for idx, stage in blk.stage.items():
        module = stage.module
        iscale.set_scaling_factor(module.area, 1e5)
        iscale.set_scaling_factor(module.feed_side.area, 1)
        iscale.set_scaling_factor(module.width, 1e4)
        set_scaling_factor(module.length, 1e1)
        # set_scaling_factor(module.feed_side.velocity, 10)
        set_scaling_factor(module.feed_side.N_Sh_comp, 1e-4)

        for e in module.feed_side.properties:
            set_scaling_factor(
                module.feed_side.properties[e].flow_mass_phase_comp["Liq", "NaCl"], 1
            )
            set_scaling_factor(module.feed_side.properties[e].dens_mass_phase["Liq"], 1)
            # set_scaling_factor(module.feed_side.properties[e].dens_mass_phase["Liq"], 1e3)
            set_scaling_factor(module.feed_side.properties[e].mass_frac_phase_comp["Liq", "NaCl"], 1e1)

        for temp_stream in [
            module.eq_permeate_isothermal,
            module.feed_side.eq_equal_temp_interface,
            module.feed_side.eq_feed_isothermal,
            module.eq_permeate_outlet_isothermal,
        ]:
            for e in temp_stream:
                constraint_scaling_transform(temp_stream[e], 1e-2)
            for pressure_stream in [
                module.eq_permeate_outlet_isobaric,
                module.feed_side.eq_equal_pressure_interface,
            ]:
                for e in pressure_stream:
                    constraint_scaling_transform(pressure_stream[e], 1e-5)
            for e in module.eq_pressure_drop:
                constraint_scaling_transform(module.eq_pressure_drop[e], 1e-7)

        for e in module.feed_side.eq_N_Sh_comp:
            if e[-1] == "NaCl":
                constraint_scaling_transform(module.feed_side.eq_N_Sh_comp[e], 1e-3)

        for e in module.feed_side.eq_friction_factor:
            constraint_scaling_transform(module.feed_side.eq_friction_factor[e], 1e-2)
        for e in module.feed_side.eq_dP_dx:
            constraint_scaling_transform(module.feed_side.eq_dP_dx[e], 1e-2)

        set_scaling_factor(module.mixed_permeate[0.0].dens_mass_phase["Liq"], 1)
        set_scaling_factor(module.mixed_permeate[0.0].flow_vol_phase["Liq"], 100)
        constraint_scaling_transform(
            module.mixed_permeate[0.0].eq_flow_vol_phase["Liq"], 100
        )

        for e in module.recovery_mass_phase_comp:
            if e[-1] == "H2O":
                set_scaling_factor(module.recovery_mass_phase_comp, 1e1)


def set_ro_system_operating_conditions(m, blk, mem_area=100, RO_pressure=15e5):
    print(
        "\n\n-------------------- SETTING RO OPERATING CONDITIONS --------------------\n\n"
    )
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    solver = get_solver()
    # mem_A = 2.75 / 3.6e11  # membrane water permeability coefficient [m/s-Pa]
    # mem_B = 0.23 / 1000.0 / 3600.0  # membrane salt permeability coefficient [m/s]
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
        # stage.module.length.fix(length)
        stage.module.width.setub(20000)
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        stage.module.feed_side.channel_height.fix(height)
        stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

        # stage.module.flux_vol_phase_avg[0, "Liq"].setlb(5)
        # stage.module.flux_vol_phase_avg[0, "Liq"].setub(60)

        stage.module.feed_side.friction_factor_darcy.setub(50)

        for e in stage.module.flux_mass_phase_comp:
            if e[-1] == "H2O":
                stage.module.flux_mass_phase_comp[e].setlb(1e-5)
                stage.module.flux_mass_phase_comp[e].setub(0.99)

    # for idx, stage in blk.stage.items():
    #     # stage.module.width.setub(5000)
    #     # stage.module.feed_side.velocity[0, 0].unfix()
    #     # stage.module.feed_side.velocity[0, 1].setlb(0.0)
    #     stage.module.feed_side.K.setlb(1e-6)
    #     stage.module.feed_side.friction_factor_darcy.setub(50)
    #     stage.module.flux_mass_phase_comp.setub(1)
    #     # stage.module.flux_mass_phase_comp.setlb(1e-5)
    #     stage.module.feed_side.cp_modulus.setub(10)
    #     stage.module.rejection_phase_comp.setlb(1e-4)
    #     stage.module.feed_side.N_Re.setlb(1)
    #     stage.module.recovery_mass_phase_comp.setlb(1e-7)

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

    # stage.eq_min_water_flux = Constraint(
    #     expr=pyunits.convert(
    #         stage.module.flux_mass_phase_comp_avg[
    #             0, "Liq", "H2O"
    #         ]
    #         / stage.module.feed_side.properties[0.0,0.0].dens_mass_phase["Liq"],
    #         to_units=pyunits.liter / pyunits.m**2 / pyunits.hr,
    #     ) >= 10 * pyunits.liter / pyunits.m**2 / pyunits.hr
    # )

    # blk.eq_minimum_water_flux = Constraint(
    #     expr=pyunits.convert(
    #         m.fs.treatment.RO.stage[1].module.flux_mass_phase_comp_avg[
    #             0.0, "Liq", "H2O"
    #         ],
    #         to_units=pyunits.kg / pyunits.hr / pyunits.m**2,
    #     )
    #     <= 40 * pyunits.kg / pyunits.hr / pyunits.m**2
    # )
    #     add_ro_scaling(m, stage)

    # iscale.calculate_scaling_factors(m)

    # ---checking model---
    # assert_units_consistent(m)
    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")


def add_ro_costing(m, blk, costing_blk=None):
    # unit equipment capital and operating costs
    if costing_blk is None:
        costing_blk = m.fs.costing

    for stage in blk.stage.values():
        stage.module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing_blk
        )


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results


def get_sub_blocks(block, decend=False, report=False):
    blocks = []
    for v in block.component_data_objects(
        ctype=Block, active=True, descend_into=decend
    ):
        print(v)
        if report:
            try:
                table = v._get_stream_table_contents()
                for item in table:
                    print(table[item])
            except:
                pass


def display_ro_system_build(m):
    get_sub_blocks(m.fs)
    get_sub_blocks(m.fs.ro)
    for stage in m.fs.ro.stage.values():
        get_sub_blocks(stage)
    print("\n")


def display_dof_breakdown(blk, decend=False, report=False):
    print(
        "\n\n-------------------- DEGREE OF FREEDOM BREAKDOWN --------------------\n\n"
    )
    print(f'{"BLOCK":<40s}{"DEGREES OF FREEDOM":<30s}')
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(f"{v.name:<40s}{degrees_of_freedom(v)}")


def display_inlet_conditions(blk):
    print("\n\n")

    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )


def display_flow_table(blk):
    print("\n\n")
    print("RO System Flow Table")
    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Product":<34s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Disposal":<34s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )


def report_RO(m, blk):
    print(f"\n\n-------------------- RO Report --------------------\n")
    print(f'{"Recovery":<30s}{value(100*m.fs.water_recovery):<10.1f}{"%"}')
    print(
        f'{"RO Operating Pressure":<30s}{value(pyunits.convert(blk.feed.properties[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}'
    )
    print(f'{"RO Membrane Area":<30s}{value(blk.stage[1].module.area):<10.1f}{"m^2"}')
    # print(
    #     f'{"Average Flux":<30s}{value(blk.stage[1].module.flux_vol_phase_avg[0, "Liq"]):<10.2f}{pyunits.get_units(blk.stage[1].module.flux_vol_phase_avg[0, "Liq"])}'
    # )
    print(blk.stage[1].module.report())


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.RO_properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.RO_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.disposal = Product(property_package=m.fs.RO_properties)

    m.fs.ro = FlowsheetBlock(dynamic=False)
    build_ro(m, m.fs.ro, prop_package=m.fs.RO_properties)

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


def print_RO_costing_breakdown(blk):
    print(f"\n\n-------------------- RO Costing Breakdown --------------------\n")
    print(
        f'{"RO Capital Cost":<35s}{f"${value(blk.stage[1].module.costing.capital_cost):<25,.0f}"}'
    )
    print(
        f'{"RO Operating Cost":<35s}{f"${value(blk.stage[1].module.costing.fixed_operating_cost):<25,.0f}"}'
    )


def breakdown_dof(blk):
    equalities = [c for c in activated_equalities_generator(blk)]
    active_vars = variables_in_activated_equalities_set(blk)
    fixed_active_vars = fixed_variables_in_activated_equalities_set(blk)
    unfixed_active_vars = unfixed_variables_in_activated_equalities_set(blk)
    print("\n ===============DOF Breakdown================\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(blk)}")
    print(f"Activated Variables: ({len(active_vars)})")
    for v in active_vars:
        print(f"   {v}")
    print(f"Activated Equalities: ({len(equalities)})")
    for c in equalities:
        print(f"   {c}")

    print(f"Fixed Active Vars: ({len(fixed_active_vars)})")
    for v in fixed_active_vars:
        print(f"   {v}")

    print(f"Unfixed Active Vars: ({len(unfixed_active_vars)})")
    for v in unfixed_active_vars:
        print(f"   {v}")
    print("\n")
    print(f" {f' Active Vars':<30s}{len(active_vars)}")
    print(f"{'-'}{f' Fixed Active Vars':<30s}{len(fixed_active_vars)}")
    print(f"{'-'}{f' Activated Equalities':<30s}{len(equalities)}")
    print(f"{'='}{f' Degrees of Freedom':<30s}{degrees_of_freedom(blk)}")
    print("\nSuggested Variables to Fix:")

    if degrees_of_freedom != 0:
        unfixed_vars_without_constraint = [
            v for v in active_vars if v not in unfixed_active_vars
        ]
        for v in unfixed_vars_without_constraint:
            if v.fixed is False:
                print(f"   {v}")


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    display_ro_system_build(m)
    set_operating_conditions(m, Qin=171.763, Cin=3.717, ro_pressure=30e5)
    set_ro_system_operating_conditions(m, m.fs.ro, mem_area=10000)
    add_ro_scaling(m, m.fs.ro)
    iscale.calculate_scaling_factors(m)
    init_system(m)
    solve(m)

    display_flow_table(m.fs.ro)
    report_RO(m, m.fs.ro)
    # print_RO_costing_breakdown(m.fs.ro)