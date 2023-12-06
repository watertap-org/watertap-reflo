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
    Block,
    RangeSet,
    check_optimal_termination,
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
# from analysisWaterTAP.flowsheets.lssro_oaro.lsrro_oaro_utils.utils import *
# from analysisWaterTAP.flowsheets.lssro_oaro.lsrro_oaro_utils.debug import *

# from analysisWaterTAP.utils import flowsheet_utils as fsTools

def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')

def _initialize(m, blk, optarg):
    try:
        blk.initialize()
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")
        
        # blk.display()
        blk.report()
        print_infeasible_bounds(m)
        # print_close_to_bounds(m)
        # print_infeasible_constraints(m)
        assert False
        
        print('\n')

_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

def build_ro(m, blk, number_of_stages=3) -> None:
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)
    blk.numberOfStages = Param(initialize=number_of_stages)
    blk.Stages = RangeSet(blk.numberOfStages)
    blk.booster_pumps = False

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(number_of_stages - 1)

    blk.primary_mixer = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = number_of_stages,
    )

    blk.stage = FlowsheetBlock(
        RangeSet(number_of_stages),
        dynamic=False)
    
    for idx, stage in blk.stage.items():
        if stage.index() > 1:
            build_ro_stage(m, stage, booster_pump=blk.booster_pumps)
        else:
            build_ro_stage(m, stage)

    #FIX This needs to be moved up the chain
    m.fs.ro_feed_to_first_stage = Arc(
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
        source=blk.stage[number_of_stages].retentate.outlet, destination=blk.disposal.inlet
    )

    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp

def build_ro_stage(m, blk, booster_pump=False):
    # Define IO
    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.permeate = StateJunction(property_package=m.fs.properties)
    blk.retentate = StateJunction(property_package=m.fs.properties)
    blk.has_booster_pump = booster_pump

    if booster_pump:
        blk.booster_pump = Pump(property_package=m.fs.properties)

    blk.module = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=2,
        has_full_reporting = True
    )

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

def init_ro_system(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING RO SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    for stage in blk.stage.values():
        print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
    print('\n\n')
    # assert_no_degrees_of_freedom(m)

    blk.feed.initialize(optarg=optarg)
    #FIX This needs to be moved up the chain
    propagate_state(m.fs.ro_feed_to_first_stage)

    for stage in blk.stage.values():
        init_ro_stage(m, stage, solver=solver)
        if stage.index() < blk.numberOfStages:
            propagate_state(blk.stage_retentate_to_next_stage[stage.index()])
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])
        else:
            propagate_state(blk.last_stage_retentate_to_ro_retentate)
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])

    blk.disposal.initialize(optarg=optarg)
    blk.primary_mixer.initialize(optarg=optarg)
    propagate_state(blk.primary_mixer_to_product)
    blk.product.initialize(optarg=optarg)

    print("\n\n-------------------- RO INITIALIZATION COMPLETE --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    for stage in blk.stage.values():
        print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
    print('\n\n')

def init_ro_stage(m, stage, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    if stage.has_booster_pump:
        stage.feed.initialize(optarg=optarg)
        propagate_state(stage.stage_feed_to_booster_pump)
        stage.booster_pump.initialize(optarg=optarg)
        propagate_state(stage.stage_booster_pump_to_module)
    else:
        stage.feed.initialize(optarg=optarg)
        propagate_state(stage.stage_feed_to_module)

    # stage.module.initialize(optarg=optarg)
    _initialize(m, stage.module, optarg)
    propagate_state(stage.stage_module_to_retentate)
    propagate_state(stage.stage_module_to_permeate)
    # print(stage.module.report())

    stage.permeate.initialize(optarg=optarg)
    stage.retentate.initialize(optarg=optarg)

def set_operating_conditions(m, Qin=None, Qout=None, Cin=None, water_recovery=None, primary_pump_pressure=80e5):
    # osParams.add_default_operating_vars(m.fs)
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

    # m.fs.product_salinity = Var(
    #     initialize=200e-6,
    #     domain=NonNegativeReals,
    #     units=pyunits.dimensionless,
    # )

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
    supply_pressure = 2.7e5

#     # initialize feed
    m.fs.feed.pressure[0].fix(pressure_atm)
    m.fs.feed.temperature[0].fix(feed_temperature)

    m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(primary_pump_pressure)
    m.fs.primary_pump.efficiency_pump.fix(0.8)
    iscale.set_scaling_factor(m.fs.primary_pump.control_volume.work, 1e-3)

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    if Qout is not None:
        m.fs.perm_flow_mass.fix(Qout)
    if Qin is not None:
        m.fs.feed_flow_mass.fix(Qin)

#     # iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)
    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    m.fs.feed_salinity.fix(Cin)
    iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

#     # m.fs.product_salinity.fix(500e-6)
#     # m.fs.product_salinity.unfix()

    # m.fs.eq_product_quality = Constraint(
    #     expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    #     <= m.fs.product_salinity
    # )

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
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ].value = m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

#     # REVIEW: Make sure this is applied in the right place
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    )

def calc_scale(value):
    return math.floor(math.log(value, 10))

def set_ro_system_operating_conditions(m, blk, mem_area=100, booster_pump_pressure=80e5):
    # parameters
    mem_A = 2.75 / 3.6e11  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 0.23 / 1000.0 / 3600.0  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.95  # spacer porosity in membrane stage [-]
    width = 500  # effective membrane width [m]
    area = mem_area # membrane area [m^2]
    pressure_atm = 101325  # atmospheric pressure [Pa]
    pump_efi = 0.8  # pump efficiency [-]

    blk.stage[1].module.feed_side.velocity[0, 0].fix(0.35)

    for idx, stage in blk.stage.items():
        stage.module.A_comp.fix(mem_A)
        stage.module.B_comp.fix(mem_B)
        stage.module.area.fix(area/idx)
        # stage.module.width.fix(width/idx)
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        if stage.has_booster_pump:
            stage.booster_pump.control_volume.properties_out[0].pressure.fix(booster_pump_pressure)

        # stage.module.feed_side.velocity[0, 0].fix(0.25)
        if (
            stage.module.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ) or stage.module.config.pressure_change_type == PressureChangeType.calculated:
            stage.module.feed_side.channel_height.fix(height)
            stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

        iscale.set_scaling_factor(stage.module.area, 3)
        iscale.set_scaling_factor(stage.module.feed_side.area, 3)
        iscale.set_scaling_factor(stage.module.width, 3)

    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)

def get_sub_blocks(block, decend = False, report=False):
    blocks = []
    for v in block.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(v)
        if report:
            try:
                table = v._get_stream_table_contents()
                for item in table:
                    print(table[item])
                # print(v._get_stream_table_contents())
            except:
                pass
        
def display_ro_system_build(m):
    get_sub_blocks(m.fs)
    # get_sub_blocks(m.fs.ro)
    # for stage in m.fs.ro.stage.values():
    #     get_sub_blocks(stage)
    print('\n')

def display_flow_table(m):
    print('\n\n')
    print(f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}')
    print(f'{"Feed":<34s}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.feed.pressure[0], to_units=pyunits.bar)):<30.1f}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
    print(f'{"Primary Pump Inlet":<34s}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_in[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.primary_pump.control_volume.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
    print(f'{"Primary Pump Outlet":<34s}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.primary_pump.control_volume.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
    
    for idx, unit in enumerate([m.fs.softener, m.fs.degasifier, m.fs.RO, m.fs.ED]):
        print(f'{str(unit.name).split(".")[1] + " Feed":<34s}{unit.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.feed.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
        print(f'{str(unit.name).split(".")[1] + " Product":<34s}{unit.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.product.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.product.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
        print(f'{str(unit.name).split(".")[1] + " Disposal":<34s}{unit.disposal.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.disposal.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.disposal.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')

    for idx, stage in m.fs.RO.stage.items():
        print(f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
    for idx, stage in m.fs.RO.stage.items():
        print(f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
    for idx, stage in m.fs.RO.stage.items():
        print(f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    
    m.fs.ro = FlowsheetBlock(dynamic=False)
    
    build_ro(m,m.fs.ro, number_of_stages=2)
    display_ro_system_build()
