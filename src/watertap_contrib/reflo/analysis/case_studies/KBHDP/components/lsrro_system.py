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
from idaes.core.solvers import get_solver
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
    MixingType,
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
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)

# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "ro_module %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


def check_jac(m, print_extreme_jacobian_values=True):
    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(m, min_scale=1e-8)
    try:
        cond_number = iscale.jacobian_cond(m, jac=jac_scaled) / 1e10
        print("--------------------------")
        print("COND NUMBER:", cond_number)
    except:
        print("Cond number failed")
        cond_number = None
    if print_extreme_jacobian_values:
        print("--------------------------")
        print("Extreme Jacobian entries:")
        extreme_entries = iscale.extreme_jacobian_entries(
            m, jac=jac_scaled, nlp=nlp, zero=1e-20, large=100
        )
        for val, var, con in extreme_entries:
            if val >= 100:
                print(val, var.name, con.name)
        print("--------------------------")
        print("Extreme Jacobian columns:")
        extreme_cols = iscale.extreme_jacobian_columns(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, var in extreme_cols:
            if val >= 100:
                print(val, var.name)
        print("------------------------")
        print("Extreme Jacobian rows:")
        extreme_rows = iscale.extreme_jacobian_rows(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, con in extreme_rows:
            if val >= 100:
                print(val, con.name)
    return cond_number


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print("\n")


_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")


def build_lsrro(m, blk, number_of_stages=1, prop_package=None):

    print(f'\n{"=======> BUILDING LSRRO SYSTEM <=======":^60}\n')
    print(f"Number of Stages: {number_of_stages}")

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.retentate = StateJunction(property_package=m.fs.properties)
    blk.numberOfStages = Param(initialize=number_of_stages)
    blk.Stages = RangeSet(blk.numberOfStages)

    if number_of_stages > 1:
        blk.IntermediateStages = RangeSet(2, blk.numberOfStages - 1)
        blk.LSRRO_Stages = RangeSet(2, blk.numberOfStages)
    else:
        blk.IntermediateStages = RangeSet(0)

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(blk.numberOfStages - 1)

    blk.stage = FlowsheetBlock(RangeSet(number_of_stages), dynamic=False)

    for stage in blk.stage.values():
        if stage.index() == 1:
            build_lsrro_stage(
                m, stage, stage.index(), intermediate_stage=False, non_final_stage=True
            )
        elif (stage.index() > 1) & (stage.index() < number_of_stages):
            build_lsrro_stage(
                m, stage, stage.index(), intermediate_stage=True, non_final_stage=True
            )
        else:
            build_lsrro_stage(
                m, stage, stage.index(), intermediate_stage=False, non_final_stage=False
            )

    blk.feed_to_first_stage = Arc(
        source=blk.feed.outlet, destination=blk.stage[1].feed.inlet
    )

    blk.pump_to_mixer = Arc(
        blk.LSRRO_Stages,
        rule=lambda blk, n: {
            "source": blk.stage[n].permeate.outlet,
            "destination": blk.stage[n - 1].mixer.downstream,
        },
    )

    blk.stage_retentate_to_next_stage = Arc(
        blk.NonFinalStages,
        rule=lambda blk, n: {
            "source": blk.stage[n].retentate.outlet,
            "destination": blk.stage[n + 1].feed.inlet,
        },
    )

    blk.first_stage_permeate_to_product = Arc(
        source=blk.stage[1].permeate.outlet, destination=blk.product.inlet
    )

    blk.last_retentate_to_disposal = Arc(
        source=blk.stage[blk.numberOfStages].retentate.outlet,
        destination=blk.retentate.inlet,
    )


def build_lsrro_stage(
    m, blk, stage_idx, intermediate_stage=False, non_final_stage=False
):
    print(f"Building LSRRO Stage {stage_idx}")
    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.permeate = StateJunction(property_package=m.fs.properties)
    blk.retentate = StateJunction(property_package=m.fs.properties)

    blk.intermediate_stage = intermediate_stage
    blk.non_final_stage = non_final_stage

    blk.stage_pump = Pump(property_package=m.fs.properties)
    blk.stage_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        # Add costing
    )

    if non_final_stage:
        blk.mixer = Mixer(
            property_package=m.fs.properties,
            has_holdup=False,
            momentum_mixing_type=MomentumMixingType.equality,
            energy_mixing_type=MixingType.none,
            inlet_list=["upstream", "downstream"],
        )

    if stage_idx > 1:
        blk.booster_pump = Pump(property_package=m.fs.properties)

    blk.module = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
    )

    # Define Connections
    blk.stage_feed_to_stage_pump = Arc(
        source=blk.feed.outlet,
        destination=blk.stage_pump.inlet,
    )

    if non_final_stage:
        blk.stage_pump_to_mixer = Arc(
            source=blk.stage_pump.outlet,
            destination=blk.mixer.upstream,
        )

        blk.mixer_to_module = Arc(
            source=blk.mixer.outlet,
            destination=blk.module.inlet,
        )
    else:
        blk.stage_pump_to_module = Arc(
            source=blk.stage_pump.outlet,
            destination=blk.module.inlet,
        )

    if stage_idx > 1:
        blk.module_to_booster_pump = Arc(
            source=blk.module.permeate,
            destination=blk.booster_pump.inlet,
        )
        blk.booster_pump_to_permeate = Arc(
            source=blk.booster_pump.outlet,
            destination=blk.permeate.inlet,
        )
    else:
        blk.module_permeate_to_permeate = Arc(
            source=blk.module.permeate,
            destination=blk.permeate.inlet,
        )

    blk.module_retentate_to_retentate = Arc(
        source=blk.module.retentate,
        destination=blk.retentate.inlet,
    )


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n--------- INITIALIZING SYSTEM ---------\n")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_lsrro_feed)

    init_lsrro_system(m, m.fs.lsrro, verbose=verbose, solver=solver)

    propagate_state(m.fs.lsrro_to_product)
    propagate_state(m.fs.lsrro_to_disposal)

    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)


def init_lsrro_system(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n--------- INITIALIZING LSRRO SYSTEM ---------\n")

    for init_pass in range(3):
        print(f"\n--------- INITIALIZATION PASS: {init_pass} ---------\n")
        forward_init_pass(m, blk, verbose=True, solver=None)

    propagate_state(blk.first_stage_permeate_to_product)
    propagate_state(blk.last_retentate_to_disposal)

    blk.product.initialize(optarg=optarg)
    blk.retentate.initialize(optarg=optarg)


def forward_init_pass(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_first_stage)

    for stage in blk.stage.values():
        # init_lsrro_stage(m, stage, solver=solver)
        init_lsrro_stage(m, stage, solver=solver)

        # Handle propogation of the retentate from the stage to the next stage
        if stage.index() < blk.numberOfStages:
            print("Non-final Stage: ", stage)
            propagate_state(blk.stage_retentate_to_next_stage[stage.index()])
        else:
            print("Final Stage: ", stage)
            propagate_state(blk.last_retentate_to_disposal)
            blk.retentate.initialize(optarg=optarg)

        # Handle propogation of the permeate depending on the stage
        if stage.index() > 1:
            propagate_state(blk.pump_to_mixer[stage.index()])
        else:
            propagate_state(blk.first_stage_permeate_to_product)
            m.fs.product.initialize(optarg=optarg)


def backward_init_pass(m, solver=None):
    pass


def init_lsrro_stage(m, stage, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(f"\n--------- INITIALIZING LSRRO STAGE {stage.index()} ---------\n")
    stage.feed.initialize(optarg=optarg)
    propagate_state(stage.stage_feed_to_stage_pump)

    stage.stage_pump.initialize(optarg=optarg)
    if stage.non_final_stage:
        propagate_state(stage.stage_pump_to_mixer)

        if stage.index() == 1:
            stage.mixer.initialize(optarg=optarg)
        else:
            stage.mixer.initialize(optarg=optarg)

        propagate_state(stage.mixer_to_module)
        stage.mixer.initialize(optarg=optarg)
        propagate_state(stage.mixer_to_module)
    else:
        propagate_state(stage.stage_pump_to_module)

    stage.module.initialize(optarg=optarg)
    print(stage.module.report())

    # NOTE HERE MAKE SURE TO HANDLE THE PROPAGATION OF THE BOOSTER PUMP
    if stage.index() > 1:
        propagate_state(stage.module_to_booster_pump)
        stage.booster_pump.initialize(optarg=optarg)
        propagate_state(stage.booster_pump_to_permeate)
    else:
        propagate_state(stage.module_permeate_to_permeate)

    propagate_state(stage.module_retentate_to_retentate)
    stage.permeate.initialize(optarg=optarg)
    stage.retentate.initialize(optarg=optarg)


def _lsrro_mixer_guess_initializer(
    mixer, solvent_multiplier, solute_multiplier, optarg
):
    print("Mixer Guess Initializer")
    for vname in mixer.upstream.vars:
        if vname == "flow_mass_phase_comp":
            for time, phase, comp in mixer.upstream.vars[vname]:
                if comp in mixer.config.property_package.solute_set:
                    mixer.downstream.vars[vname][time, phase, comp].value = (
                        solute_multiplier
                        * mixer.upstream.vars[vname][time, phase, comp].value
                    )
                elif comp in mixer.config.property_package.solvent_set:
                    mixer.downstream.vars[vname][time, phase, comp].value = (
                        solvent_multiplier
                        * mixer.upstream.vars[vname][time, phase, comp].value
                    )
                else:
                    raise RuntimeError(f"Unknown component {comp}")
        else:  # copy the state
            for idx in mixer.upstream.vars[vname]:
                mixer.downstream.vars[vname][idx].value = mixer.upstream.vars[vname][
                    idx
                ].value

    mixer.initialize(optarg=optarg)


def set_operating_conditions(m, Qin=None, Qout=None, Cin=None, water_recovery=None):
    if Qin is None:
        Qin = 1
    if Cin is None:
        Cin = 35

    feed_temperature = 273.15 + 20
    pressure_atm = 101325
    m.fs.feed.pressure[0].fix(pressure_atm)
    m.fs.feed.temperature[0].fix(feed_temperature)

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
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

    m.fs.product_salinity = Var(
        initialize=200e-6,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Water Recovery",
    )

    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery)
    else:
        m.fs.water_recovery.fix(0.5)

    m.fs.feed_flow_mass.fix(Qin)
    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    m.fs.feed_salinity.fix(Cin)
    iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
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
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)
    )

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Liq", "NaCl")
    )

    _logger.info("ro scaling h2o:{} tds:{}".format(scale_flow, scale_tds))

    assert_units_consistent(m)

    print("\n--------- SETTING OPERATING CONDITIONS ---------\n")
    print(f"Feed Flow Rate: {Qin} kg/s")
    print(f"Feed Salinity: {Cin} g/L")


def scale_system(m):
    constraint_scaling_transform(m.fs.eq_water_recovery, 1e-2)
    set_scaling_factor(m.fs.water_recovery, 1)


def scale_lsrro_stage(m, stage):
    module = stage.module

    set_scaling_factor(module.area, 1e-2)
    set_scaling_factor(module.feed_side.area, 1)
    set_scaling_factor(module.width, 1e-2)

    for e in module.feed_side.K:
        set_scaling_factor(module.feed_side.K[e], 1e-2)

    set_scaling_factor(module.feed_side.dh, 1e2)
    # constraint_scaling_transform(stage.feed_side.eq_dh, 100)


def calc_scale(value):
    return math.floor(math.log(value, 10))


def set_lsrro_system_operating_conditions(
    m,
    blk,
    mem_area=10,
    RO_pump_pressure=65e5,
    B_max=None,
):
    # parameters
    mem_A = 2.0 / 3.6e11  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 10 / 1000.0 / 3600.0  # membrane salt permeability coefficient [m/s]
    mem_B_RO = 0.14 / 1000.0 / 3600.0
    mem_A_RO = 2.0 / 3.6e11  # m
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.90  # spacer porosity in membrane stage [-]
    length = 1  # effective membrane width [m]
    area = mem_area  # membrane area [m^2]
    primary_pump_pressure = RO_pump_pressure  # primary pump pressure [Pa]
    pressure_atm = 101325  # atmospheric pressure [Pa]
    pump_efi = 0.8  # pump efficiency [-]

    print("\n--------- SETTING OPERATING CONDITIONS ---------\n")

    # initialize stages
    for idx, stage in m.fs.lsrro.stage.items():
        if idx == m.fs.lsrro.FirstStage:
            stage.module.A_comp.fix(mem_A_RO)
            stage.module.B_comp.fix(mem_B_RO)
        else:
            stage.module.A_comp.fix(mem_A)
            stage.module.B_comp.fix(mem_B)

        stage.module.area.fix(area)
        stage.module.length.fix(length)
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        if (
            (
                stage.module.config.mass_transfer_coefficient
                == MassTransferCoefficient.calculated
            )
            or stage.module.config.pressure_change_type == PressureChangeType.calculated
        ):
            stage.module.feed_side.channel_height.fix(height)
            stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

        stage.stage_pump.control_volume.properties_out[0].pressure.fix(
            primary_pump_pressure
        )
        stage.stage_pump.efficiency_pump.fix(pump_efi)

        if idx > 1:
            # stage.booster_pump.control_volume.properties_out[0].pressure.unfix()
            # stage.booster_pump.control_volume..deltaP.unfix()
            stage.booster_pump.deltaP.unfix()
            stage.booster_pump.efficiency_pump.fix(pump_efi)

        # scale_lsrro_stage(m, stage)

        print(f"Stage {idx} Pump Pressure: {primary_pump_pressure:<5.2e} Pa")

        print(f"Stage {idx} Membrane Area: {area:<5.2f} m^2")
        print(f"Stage {idx} Membrane Length: {length:<5.2f} m")
        print(
            f"Stage {idx} Membrane A: {mem_A:<5.2e} {pyunits.get_units(stage.module.A_comp)}"
        )
        print(
            f"Stage {idx} Membrane B: {mem_B:<5.2e} {pyunits.get_units(stage.module.B_comp)}"
        )
        print("\n")

    print(f"DEGREES OF FREEDOM: {degrees_of_freedom(m)}")


def optimize(m):

    for stage in m.fs.lsrro.stage.values():
        if stage.index() > 1:
            stage.booster_pump.control_volume.properties_out[0].pressure.unfix()

    print(f"DEGREES OF FREEDOM: {degrees_of_freedom(m)}")


def solve(m, solver=None, tee=True, raise_on_failure=False, debug=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            check_jac(m)

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(m)

        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(m)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(m)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(m)

        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            check_jac(m)

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
                # print(v._get_stream_table_contents())
            except:
                pass


def display_lsrro_system_build(m):
    get_sub_blocks(m.fs)
    get_sub_blocks(m.fs.lsrro)
    for stage in m.fs.lsrro.stage.values():
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
    # print(blk.feed.display())

    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    # assert False


def display_flow_table(m):
    print("\n\n")
    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{m.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{m.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Product":<34s}{m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.product.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Disposal":<34s}{m.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.disposal.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    for idx, stage in m.fs.lsrro.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in m.fs.lsrro.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in m.fs.lsrro.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )


def report_LSRRO(m, blk):
    print(f"\n\n-------------------- RO Report --------------------\n")
    print(f'{"Recovery":<30s}{value(100*m.fs.water_recovery):<10.1f}{"%"}')
    print(
        f'{"RO Operating Pressure":<30s}{value(pyunits.convert(blk.pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}'
    )


def build_system(number_of_stages=2):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    m.fs.DOF = []

    m.fs.lsrro = FlowsheetBlock(dynamic=False)

    build_lsrro(m, m.fs.lsrro, number_of_stages)

    m.fs.feed_to_lsrro_feed = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.lsrro.feed.inlet,
    )

    m.fs.lsrro_to_product = Arc(
        source=m.fs.lsrro.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.lsrro_to_disposal = Arc(
        source=m.fs.lsrro.retentate.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    display_lsrro_system_build(m)
    set_operating_conditions(m)
    set_lsrro_system_operating_conditions(m, m.fs.lsrro, mem_area=20)
    init_system(m)
    display_flow_table(m)
    optimize(m)
    solve(m, raise_on_failure=True, debug=True)
    display_flow_table(m)
