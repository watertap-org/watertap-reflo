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
    TransformationFactory,
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
        
        blk.display()
        blk.report()
        # print_infeasible_bounds(m)
        # print_close_to_bounds(m)
        # print_infeasible_constraints(m)
        assert False
        
        print('\n')

_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

class ROSystem():
    def __init__(self, m, number_of_stages=3, RO="1D", cost_type="EPA", compaction_type="None", **kwargs) -> None:
        self.stages = number_of_stages
        self.RO_type = RO
        self.compaction_type = compaction_type
        self.cost_type = cost_type
        self.blk = FlowsheetBlock(dynamic=False)
        self.properties = m.fs.properties
        self.m = m

    def build(self, blk) -> None:
        blk.feed = StateJunction(property_package=self.properties)
        blk.product = StateJunction(property_package=self.properties)
        blk.retentate = StateJunction(property_package=self.properties)
        blk.numberOfStages = Param(initialize=self.stages)
        blk.Stages = RangeSet(blk.numberOfStages)
        blk.booster_pumps = False

        blk.FirstStage = blk.Stages.first()
        blk.LastStage = blk.Stages.last()
        blk.NonFinalStages = RangeSet(blk.numberOfStages - 1)

        blk.primary_mixer = Mixer(
            property_package=self.properties,
            has_holdup=False,
            num_inlets = self.stages,
        )

        blk.stage = FlowsheetBlock(
            RangeSet(self.stages),
            dynamic=False)
        
        for idx, stage in blk.stage.items():
            if stage.index() > 1:
                self.build_ro_stage(self.m, stage, booster_pump=blk.booster_pumps)
            else:
                self.build_ro_stage(self.m, stage)
            # else:
            #     build_lsrro_stage(m, stage, intermediate_stage=False)
            
        #FIX This needs to be moved up the chain
        # m.fs.ro_feed_to_first_stage = Arc(
        #     source=blk.feed.outlet,
        #     destination=blk.stage[1].feed.inlet,
        # )

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
            source=blk.stage[blk.numberOfStages].retentate.outlet, destination=blk.retentate.inlet
        )

    def build_stage():
        pass

# def build_system(number_of_stages=1, booster_pumps = False):
#     m = ConcreteModel()
#     m.fs = FlowsheetBlock(dynamic=False)
#     m.fs.properties = NaClParameterBlock()
#     m.fs.costing = WaterTAPCosting()
#     m.fs.feed = Feed(property_package=m.fs.properties)
#     m.fs.product = Product(property_package=m.fs.properties)
#     m.fs.disposal = Product(property_package=m.fs.properties)
#     m.fs.DOF = []

#     m.fs.primary_pump = Pump(property_package=m.fs.properties)
#     # m.fs.primary_pump.costing = UnitModelCostingBlock(
#     #     flowsheet_costing_block=m.fs.costing,
#     #     costing_method=cost_high_pressure_pump_lsrro,
#     # )

#     m.fs.ro = FlowsheetBlock(dynamic=False)

#     build_ro_system(m, m.fs.ro, number_of_stages, booster_pumps = booster_pumps)
    
#     m.fs.feed_to_primary_pump = Arc(
#         source=m.fs.feed.outlet,
#         destination=m.fs.primary_pump.inlet,
#     )

#     m.fs.primary_pump_to_ro_feed = Arc(
#         source=m.fs.primary_pump.outlet,
#         destination=m.fs.ro.feed.inlet,
#     )

#     m.fs.ro_to_product = Arc(
#         source=m.fs.ro.product.outlet,
#         destination=m.fs.product.inlet,
#     )
#     m.fs.ro_to_disposal = Arc(
#         source=m.fs.ro.retentate.outlet,
#         destination=m.fs.disposal.inlet,
#     )

#     TransformationFactory("network.expand_arcs").apply_to(m)

#     return m

# def build_ro_system(m,blk,stages, booster_pumps = False):
#     blk.feed = StateJunction(property_package=m.fs.properties)
#     blk.product = StateJunction(property_package=m.fs.properties)
#     blk.retentate = StateJunction(property_package=m.fs.properties)
#     blk.numberOfStages = Param(initialize=stages)
#     blk.Stages = RangeSet(blk.numberOfStages)
#     blk.booster_pumps = booster_pumps

#     blk.FirstStage = blk.Stages.first()
#     blk.LastStage = blk.Stages.last()
#     blk.NonFinalStages = RangeSet(blk.numberOfStages - 1)

#     blk.primary_mixer = Mixer(
#         property_package=m.fs.properties,
#         has_holdup=False,
#         num_inlets = stages,
#     )

#     blk.stage = FlowsheetBlock(
#         RangeSet(stages),
#         dynamic=False)
    
#     for idx, stage in blk.stage.items():
#         if stage.index() > 1:
#             build_ro_stage(m, stage, booster_pump=blk.booster_pumps)
#         else:
#             build_ro_stage(m, stage)
#         # else:
#         #     build_lsrro_stage(m, stage, intermediate_stage=False)
        

#     m.fs.ro_feed_to_first_stage = Arc(
#         source=blk.feed.outlet,
#         destination=blk.stage[1].feed.inlet,
#     )

#     blk.stage_retentate_to_next_stage = Arc(
#         blk.NonFinalStages,
#         rule=lambda blk, n: {
#             "source": blk.stage[n].retentate.outlet,
#             "destination": blk.stage[n + 1].feed.inlet,
#         },
#     )
    
#     blk.stage_permeate_to_mixer = Arc(
#         blk.Stages,
#         rule=lambda blk, n: {
#             "source": blk.stage[n].permeate.outlet,
#             "destination": getattr(blk.primary_mixer, "inlet_" + str(n)),
#         },
#     )

#     blk.primary_mixer_to_product = Arc(
#         source=blk.primary_mixer.outlet,
#         destination=blk.product.inlet,
#     )

#     blk.last_stage_retentate_to_ro_retentate = Arc(
#         source=blk.stage[blk.numberOfStages].retentate.outlet, destination=blk.retentate.inlet
#     )

    def build_ro_stage(self, m, blk, booster_pump=False):
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

#     if booster_pump:
#         blk.stage_feed_to_booster_pump = Arc(
#             source=blk.feed.outlet,
#             destination=blk.booster_pump.inlet,
#         )
#         blk.stage_booster_pump_to_module = Arc(
#             source=blk.booster_pump.outlet,
#             destination=blk.module.inlet,
#         )
#     else:
#         blk.stage_feed_to_module = Arc(
#             source=blk.feed.outlet,
#             destination=blk.module.inlet,
#         )
    
    
#     blk.stage_module_to_permeate = Arc(
#         source=blk.module.permeate,
#         destination=blk.permeate.inlet,
#     )
    
#     blk.stage_module_to_retentate = Arc(
#         source=blk.module.retentate,
#         destination=blk.retentate.inlet,
#     )

# def init_system(m, verbose=True, solver=None):
#     if solver is None:
#         solver = get_solver()

#     optarg = solver.options

#     print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
#     print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
#     # for stage in m.fs.ro.stage.values():
#     #     print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
#     # print('\n\n')
#     assert_no_degrees_of_freedom(m)
#     m.fs.feed.initialize(optarg=optarg)
#     propagate_state(m.fs.feed_to_primary_pump)

#     m.fs.primary_pump.initialize(optarg=optarg)
#     propagate_state(m.fs.primary_pump_to_ro_feed)

#     m.fs.ro.feed.initialize(optarg=optarg)
#     propagate_state(m.fs.ro_feed_to_first_stage)

#     for stage in m.fs.ro.stage.values():
#         init_ro_stage(m, stage, solver=solver)
#         if stage.index() < m.fs.ro.numberOfStages:
#             propagate_state(m.fs.ro.stage_retentate_to_next_stage[stage.index()])
#             propagate_state(m.fs.ro.stage_permeate_to_mixer[stage.index()])
#         else:
#             propagate_state(m.fs.ro.last_stage_retentate_to_ro_retentate)
#             propagate_state(m.fs.ro.stage_permeate_to_mixer[stage.index()])

#     m.fs.ro.retentate.initialize(optarg=optarg)
#     m.fs.ro.primary_mixer.initialize(optarg=optarg)
#     propagate_state(m.fs.ro.primary_mixer_to_product)
#     m.fs.ro.product.initialize(optarg=optarg)

#     propagate_state(m.fs.ro_to_product)
#     propagate_state(m.fs.ro_to_disposal)

#     m.fs.product.initialize(optarg=optarg)
#     m.fs.disposal.initialize(optarg=optarg)

#     print("\n\n-------------------- INITIALIZATION COMPLETE --------------------\n\n")

# def init_ro_stage(m, stage, solver=None):
#     if solver is None:
#         solver = get_solver()

#     optarg = solver.options

#     if stage.has_booster_pump:
#         stage.feed.initialize(optarg=optarg)
#         propagate_state(stage.stage_feed_to_booster_pump)
#         stage.booster_pump.initialize(optarg=optarg)
#         propagate_state(stage.stage_booster_pump_to_module)
#     else:
#         stage.feed.initialize(optarg=optarg)
#         propagate_state(stage.stage_feed_to_module)

#     stage.module.initialize(optarg=optarg)
#     propagate_state(stage.stage_module_to_retentate)
#     propagate_state(stage.stage_module_to_permeate)
#     # print(stage.module.report())

#     stage.permeate.initialize(optarg=optarg)
#     stage.retentate.initialize(optarg=optarg)

# def set_operating_conditions(m, Qin=None, Qout=None, Cin=None, water_recovery=None, primary_pump_pressure=80e5):
#     # osParams.add_default_operating_vars(m.fs)
#     # if Qin is None:
#     #     Qin = 1
#     if Cin is None:
#         Cin = 35

#     m.fs.water_recovery = Var(
#         initialize=0.5,
#         bounds=(0, 0.99),
#         domain=NonNegativeReals,
#         units=pyunits.dimensionless,
#         doc="System Water Recovery",
#     )

#     m.fs.feed_salinity = Var(
#         initialize=35,
#         bounds=(0, 2000),
#         domain=NonNegativeReals,
#         units=pyunits.dimensionless,
#         doc="System Water Recovery",
#     )

#     # m.fs.product_salinity = Var(
#     #     initialize=200e-6,
#     #     domain=NonNegativeReals,
#     #     units=pyunits.dimensionless,
#     # )

#     m.fs.feed_flow_mass = Var(
#         initialize=1,
#         bounds=(0.00001, 1e6),
#         domain=NonNegativeReals,
#         units=pyunits.kg / pyunits.s,
#         doc="System Feed Flowrate",
#     )
    
#     m.fs.perm_flow_mass = Var(
#         initialize=1,
#         bounds=(0.00001, 1e6),
#         domain=NonNegativeReals,
#         units=pyunits.kg / pyunits.s,
#         doc="System Produce Flowrate",
#     )

#     if water_recovery is not None:
#         m.fs.water_recovery.fix(water_recovery)
#     else:
#         m.fs.water_recovery.unfix()

#     feed_temperature = 273.15 + 20
#     pressure_atm = 101325
#     supply_pressure = 2.7e5

#     # initialize feed
#     m.fs.feed.pressure[0].fix(supply_pressure)
#     m.fs.feed.temperature[0].fix(feed_temperature)

#     m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(primary_pump_pressure)
#     m.fs.primary_pump.efficiency_pump.fix(0.8)
#     iscale.set_scaling_factor(m.fs.primary_pump.control_volume.work, 1e-3)

#     m.fs.eq_water_recovery = Constraint(
#         expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
#         == m.fs.product.properties[0].flow_vol
#     )

#     if Qout is not None:
#         m.fs.perm_flow_mass.fix(Qout)
#     if Qin is not None:
#         m.fs.feed_flow_mass.fix(Qin)

#     # iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)
#     iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
#     m.fs.feed_salinity.fix(Cin)
#     iscale.set_scaling_factor(m.fs.feed_salinity, 0.1)

#     # m.fs.product_salinity.fix(500e-6)
#     # m.fs.product_salinity.unfix()

#     # m.fs.eq_product_quality = Constraint(
#     #     expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
#     #     <= m.fs.product_salinity
#     # )

#     m.fs.feed_flow_constraint = Constraint(
#             expr=m.fs.feed_flow_mass == m.fs.perm_flow_mass / m.fs.water_recovery
#         )
#     iscale.set_scaling_factor(m.fs.perm_flow_mass, 1)

#     m.fs.nacl_mass_constraint = Constraint(
#         expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"] * 1000
#         == m.fs.feed_flow_mass * m.fs.feed_salinity
#     )

#     m.fs.h2o_mass_constraint = Constraint(
#         expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
#         == m.fs.feed_flow_mass * (1 - m.fs.feed_salinity / 1000)
#     )

#     m.fs.feed.properties[0].flow_vol_phase["Liq"]
#     m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

#     m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
#         m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
#     )
#     m.fs.feed.flow_mass_phase_comp[
#         0, "Liq", "H2O"
#     ].value = m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)

#     scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
#     scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

#     # REVIEW: Make sure this is applied in the right place
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
#     )

#     assert_units_consistent(m)
#     # m.fs.DOF.append(["Flowsheet", "After Set Sys. Op. Conds.", degrees_of_freedom(m)])
#     # m.fs.DOF.append(
#     #     ["RO", "After Set Sys. Op. Conds.", degrees_of_freedom(m.fs.ro)]
#     # )

# def calc_scale(value):
#     return math.floor(math.log(value, 10))

# def set_ro_scaling_and_params(
#     m,
# ):
#     # additional parameters, variables or expressions ---------------------------------------------------------------------------
#     m.fs.ro_min_pressure = Param(initialize=1e5, units=pyunits.Pa, mutable=True)
#     m.fs.ro_max_pressure = Param(initialize=85e5, units=pyunits.Pa, mutable=True)

#     # explicitly set the costing parameters used
#     m.fs.costing.utilization_factor.fix(0.9)
#     m.fs.costing.factor_total_investment.fix(2)
#     m.fs.costing.factor_maintenance_labor_chemical.fix(0.03)
#     m.fs.costing.factor_capital_annualization.fix(0.1)
#     m.fs.costing.electricity_cost.set_value(0.07)

#     m.fs.costing.cost_process()

#     product_flow_vol_total = m.fs.product.properties[0].flow_vol
#     m.fs.costing.add_annual_water_production(product_flow_vol_total)
#     m.fs.costing.add_specific_energy_consumption(product_flow_vol_total)
#     m.fs.costing.add_LCOW(product_flow_vol_total)
#     m.fs.cross_flow_velocity = Var(
#         initialize=0.25,
#         bounds=(0, 1),
#         domain=NonNegativeReals,
#         units=pyunits.meter / pyunits.second,
#     )
#     m.fs.cross_flow_velocity.unfix()
    
#     @m.fs.Constraint(list(range(len(m.fs.ro.stage))))
#     def eq_cross_flow_vel(fs, stage):
#         return (
#             m.fs.ro.stage[stage + 1].module.feed_side.velocity[0, 0]
#             == m.fs.cross_flow_velocity
#         )

#     m.fs.eq_cross_flow_vel.deactivate()
#     m.fs.ro.total_area = Expression(
#         expr=sum([stage.module.area for idx, stage in m.fs.ro.stage.items()])
#     )

#     # objective
#     m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
#     # m.fs.lcow_objective.deactivate()
#     # m.fs.total_area_objective = Objective(expr=m.fs.ro.total_area)
#     # m.fs.total_area_objective.deactivate()
#     # Expression    s for parameter sweep -----------------------------------------
#     # Final permeate concentration as mass fraction
#     m.fs.product.properties[0].mass_frac_phase_comp

#     # Touch feed concentration as mass concentration
#     m.fs.feed.properties[0].conc_mass_phase_comp
#     m.fs.product.properties[0].conc_mass_phase_comp
#     m.fs.disposal.properties[0].conc_mass_phase_comp
#     m.fs.primary_pump.control_volume.properties_in[0].conc_mass_phase_comp
#     m.fs.primary_pump.control_volume.properties_out[0].conc_mass_phase_comp
#     m.fs.ro.feed.properties[0].conc_mass_phase_comp
#     m.fs.ro.product.properties[0].conc_mass_phase_comp
#     m.fs.ro.retentate.properties[0].conc_mass_phase_comp

#     # Touch final brine concentration as mass fraction
#     m.fs.disposal.properties[0].mass_frac_phase_comp

#     m.fs.DOF.append(["Flowsheet", "After Set RO Scaling", degrees_of_freedom(m)])
#     m.fs.DOF.append(
#         ["RO", "After Set Scaling", degrees_of_freedom(m.fs.ro)]
#     )

# def set_ro_system_operating_conditions(m, booster_pump_pressure=80e5):
#     # parameters
#     mem_A = 2.75 / 3.6e11  # membrane water permeability coefficient [m/s-Pa]
#     mem_B = 0.23 / 1000.0 / 3600.0  # membrane salt permeability coefficient [m/s]
#     height = 1e-3  # channel height in membrane stage [m]
#     spacer_porosity = 0.95  # spacer porosity in membrane stage [-]
#     width = 500  # effective membrane width [m]
#     area = 37.1612*28*7 # membrane area [m^2]
#     pressure_atm = 101325  # atmospheric pressure [Pa]
#     pump_efi = 0.8  # pump efficiency [-]

#     for idx, stage in m.fs.ro.stage.items():
#         stage.module.A_comp.fix(mem_A)
#         stage.module.B_comp.fix(mem_B)
#         stage.module.area.fix(area/idx)
#         stage.module.width.fix(width/idx)
#         stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

#         if stage.has_booster_pump:
#             stage.booster_pump.control_volume.properties_out[0].pressure.fix(booster_pump_pressure)

#         # stage.module.feed_side.velocity[0, 0].fix(0.25)
#         if (
#             stage.module.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
#         ) or stage.module.config.pressure_change_type == PressureChangeType.calculated:
#             stage.module.feed_side.channel_height.fix(height)
#             stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

#         iscale.set_scaling_factor(stage.module.area, 3)
#         iscale.set_scaling_factor(stage.module.feed_side.area, 3)
#         iscale.set_scaling_factor(stage.module.width, 3)

#     iscale.calculate_scaling_factors(m)

#     # ---checking model---
#     assert_units_consistent(m)

#     m.fs.DOF.append(["Flowsheet", "After Set RO Op. Conds.", degrees_of_freedom(m)])
#     m.fs.DOF.append(
#         ["RO", "After Set RO Op. Conds.", degrees_of_freedom(m.fs.ro)]
#     )

# def optimize(
#     m,
#     water_recovery=None,
#     fixed_pressure=None,
#     Qin = None,
#     Qout = None,
#     Qprod=None,
#     mem_width = None,
#     velocity = None,
#     # A_value=None,
#     # B_value=None,
#     # permeate_quality_limit=None,
#     # set_default_bounds_on_module_dimensions=False,
#     # fixed_pressure=None,

# ):

#     print("\n\nDOF before optimization: ", degrees_of_freedom(m))
#     if water_recovery is not None:
#         print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
#         # product mass flow rate fraction of feed [-]
#         m.fs.water_recovery.fix(water_recovery)
#         m.fs.feed_flow_mass.fix(Qprod/water_recovery)
#         print(f"------- Fixed Feed Flow Rate at {Qprod/water_recovery:<3.1f} kg/s -------")  
#     else:
#         m.fs.water_recovery.unfix()
#         m.fs.water_recovery.setlb(0.01)
#         m.fs.water_recovery.setub(0.99)

#     if fixed_pressure is not None:
#         print(f"\n------- Fixed Pressure at {fixed_pressure} -------\n")
#         m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
#     else:
#         print(f"------- Unfixed Pressure -------")
#         m.fs.primary_pump.control_volume.properties_out[0].pressure.unfix()
#         m.fs.primary_pump.control_volume.properties_out[0].pressure.setlb(1e5)
#         m.fs.primary_pump.control_volume.properties_out[0].pressure.setub(80e5)

#     if velocity is not None:
#         m.fs.ro.stage[1].module.feed_side.velocity[0, 0].fix(0.35)
#         for idx, stage in m.fs.ro.stage.items():
#             stage.module.width.unfix()

#     if mem_width is not None:
#         for idx, stage in m.fs.ro.stage.items():
#             stage.module.width.fix(mem_width/idx)
#     # if Qin is not None:
#     #     m.fs.feed_flow_mass.fix(Qin)
#     #     print(f"\n------- Fixed Feed Flow at {Qin} kg/s -------")
#     # else:
#     #     m.fs.feed_flow_mass.unfix()
#     #     print(f"------- Unfixed Feed Flow -------")


#     # if Qout is not None:
#     #     m.fs.perm_flow_mass.fix(Qout)
#     #     print(f"------- Fixed Product Flow at {Qout} kg/s -------")
#     # else:
#     #     m.fs.perm_flow_mass.unfix()
#     #     print(f"------- Unfixed Product Flow -------")

#     print('\n')
#     # m.fs.DOF.append(["Flowsheet", "Before Optimization", degrees_of_freedom(m)])
#     # m.fs.DOF.append(["RO", "Before Optimization", degrees_of_freedom(m.fs.ro)])

#     # if fixed_pressure is not None:
#     #     print(f"\n------- Fixed Pressure at {fixed_pressure} -------\n")
#     #     m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
#     # else:
#     #     print(f"------- Unfixed Pressure -------")
#     #     m.fs.primary_pump.control_volume.properties_out[0].pressure.unfix()
#     #     m.fs.primary_pump.control_volume.properties_out[0].pressure.setlb(1e5)
#     #     m.fs.primary_pump.control_volume.properties_out[0].pressure.setlb(9e6)

#     # for idx, stage in m.fs.ro.stage.items():
#     #     stage.module.area.unfix()
#     #     stage.module.length.unfix()
#     #     stage.module.width.unfix()
#     #     # stage.recovery_vol_phase[0.0, "Liq"].setlb(0.1)

#     #     if set_default_bounds_on_module_dimensions:
#     #         stage.module.area.setlb(0.1)
#     #         stage.module.area.setub(None)
#     #         # stage.width.setlb(1)
#     #         # stage.width.setub(5000)
#     #     else:
#     #         stage.module.area.setlb(0.1)
#     #         stage.module.area.setub(None)
#     #         stage.module.width.setlb(0.1)
#     #         stage.module.width.setub(None)
#     #         stage.module.length.setlb(0.1)
#     #         # initialize stages
#     #     # stage.length.setub(4)
#     #     stage.module.length.setlb(1)
#     #     stage.module.width.setlb(1)
#     #     stage.module.feed_side.velocity[0, 0].unfix()
#     #     stage.module.feed_side.velocity[0, 0].setlb(0.0)
#     #     stage.module.feed_side.velocity[0, 0].setub(0.9)
#     #     stage.module.feed_side.K.setlb(1e-6)
#     #     stage.module.feed_side.friction_factor_darcy.setub(50)
#     #     stage.module.flux_mass_phase_comp.setub(1)
#     #     stage.module.feed_side.cp_modulus.setub(10)
#     #     stage.module.rejection_phase_comp.setlb(1e-4)
#     #     stage.module.feed_side.velocity.setlb(0.001)
#     #     stage.module.feed_side.N_Re.setlb(1)
#     #     stage.module.recovery_mass_phase_comp.setlb(1e-7)

#         # stage.module.length.fix(40)
#     # add upper bound for permeate concentration
#     # if permeate_quality_limit is not None:
#     #     if isinstance(permeate_quality_limit, (int, float)):
#     #         m.fs.ro.stage.module[1].mixed_permeate[0].mass_frac_phase_comp["Liq", "NaCl"].setub(
#     #             permeate_quality_limit
#     #         )
#     #     else:
#     #         raise TypeError("permeate_quality_limit must be None, integer, or float")

#     # # additional constraints
#     # if water_recovery is not None:
#     #     # product mass flow rate fraction of feed [-]
#     #     m.fs.water_recovery.fix(water_recovery)
#     # else:
#     #     m.fs.water_recovery.unfix()
#     #     m.fs.water_recovery.setlb(0.01)
#     #     m.fs.water_recovery.setub(0.99)
#     # # m.fs.eq_cross_flow_vel.activate()
#     assert_units_consistent(m)

#     print("Flowsheet DOF after optimization: ", degrees_of_freedom(m))
#     # m.fs.DOF.append(["Flowsheet", "After Optimization", degrees_of_freedom(m)])
#     # m.fs.DOF.append(["RO", "After Optimization", degrees_of_freedom(m.fs.ro)])

#     # assert_degrees_of_freedom(
#     #     m,
#     #     4 * m.fs.ro.numberOfStages
#     #     - (1 if (water_recovery is not None) else 0)
#     #     - (1 if value(m.fs.ro.numberOfStages) == 1 else 0),
#     # )

# def solve(model, solver=None, tee=True, raise_on_failure=True):
#     # ---solving---
#     if solver is None:
#         solver = get_solver()

#     print("\n--------- SOLVING ---------\n")

#     results = solver.solve(model, tee=tee)

#     if check_optimal_termination(results):
#         print("\n--------- OPTIMAL SOLVE!!! ---------\n")
#         return results
#     msg = (
#         "The current configuration is infeasible. Please adjust the decision variables."
#     )
#     if raise_on_failure:
#         debug(model, solver=solver, automate_rescale=False, resolve=False)
#         # debug(model, solver=solver, automate_rescale=False, resolve=False)
#         # check_jac(model)
#         raise RuntimeError(msg)
#     # else:
#     #     print(msg)
#     #     # debug(model, solver=solver, automate_rescale=False, resolve=False)
#     #     # check_jac(model)
#     #     return results

# def solve_box_problem(m, solver=None, tee=False, check_termination=False):
#     if solver is None:
#         solver = get_solver()

#     result = standard_solve(m, solver=solver, tee=True, debug=True)
#     return result

# def set_prob_for_box_solve(m):
#     m.fs.water_recovery.unfix()
#     # m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(75e5)
#     # m.fs.primary_pump.control_volume.properties_out[0].pressure.unfix()
#     for idx, stage in m.fs.ro.stage.items():
#         # stage.module.feed_side.velocity[0, 0].fix(0.5)
#         stage.module.recovery_vol_phase[0.0, "Liq"].setlb(0.001)
#         # stage.module.length.unfix()
#     #     stage.area.unfix()
#     #     stage.width.unfix()

#     # for idx, stage in m.fs.ro.Splitters.items():
#     #     stage.split_fraction[0, "upstream"].unfix()
#     # for idx, pump in m.fs.ro.RecircPumps.items():
#     #     pump.control_volume.properties_out[0].pressure.unfix()
#     #     pump.control_volume.properties_out[0].flow_vol_phase.fix(0.001)

#     assert_no_degrees_of_freedom(m)

# def box_solve_problem(m):
#     set_prob_for_box_solve(m)
#     fsTools.standard_solve(
#         m,
#         tee=True,
#         check_close_to_bounds=True,
#         check_var_scailing=True,
#         expected_DOFs=0,
#     )

    def get_sub_blocks(self, block, decend = False, report=False):
        blocks = []
        for v in block.component_data_objects(ctype=Block, active=True, descend_into=decend):
            print(v)
            if report:
                try:
                    print(v.report())
                except:
                    pass
            
    def display_system_build(self):
        self.get_sub_blocks(self.m.fs)
        self.get_sub_blocks(self.m.fs.ro)
        # for stage in self.m.fs.ro.stage.values():
        #     self.get_sub_blocks(stage)
        print('\n\n')

# def display_flow_table(m):
#     print('\n\n')
#     print(f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}')
#     print(f'{"Feed":<34s}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.feed.pressure[0], to_units=pyunits.bar)):<30.1f}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"Primary Pump Inlet":<34s}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_in[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.primary_pump.control_volume.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"Primary Pump Outlet":<34s}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.primary_pump.control_volume.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"RO Feed":<34s}{m.fs.ro.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.ro.feed.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.ro.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.ro.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')

#     for idx, stage in m.fs.ro.stage.items():
#         print(f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     for idx, stage in m.fs.ro.stage.items():
#         print(f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     for idx, stage in m.fs.ro.stage.items():
#         print(f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')

#     print(f'{"RO Product":<34s}{m.fs.ro.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.ro.product.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.ro.product.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.ro.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"RO Retentate":<34s}{m.fs.ro.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.ro.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.ro.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.ro.retentate.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"Product":<34s}{m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.product.pressure[0], to_units=pyunits.bar)():<30.1f}{m.fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')
#     print(f'{"Disposal":<34s}{m.fs.disposal.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.disposal.pressure[0], to_units=pyunits.bar)():<30.1f}{m.fs.disposal.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}')

# def display_system_metrics(m):
#     print('\n')
#     print(f'{"STAGE":<15s}{"RECOVERY %":<15s}{"REJECTION %":<15s}{"AREA (SQ M)":<15s}{"WIDTH (M)":<15s}{"LENGTH (M)":<15s}{"INLET VEL (M/S)":<20s}{"EXIT VEL (M/S)":<20s}{"EXIT DRIVING FORCE (BAR)":<15s}')
#     for idx, stage in m.fs.ro.stage.items():
#         del_pi = value(pyunits.convert(stage.module.feed_side.properties_interface[0.0,1.0].pressure_osm_phase["Liq"] - stage.module.permeate_side[0.0,1.0].pressure_osm_phase["Liq"], to_units=pyunits.bar))
#         del_P = value(pyunits.convert(stage.module.feed_side.properties_interface[0.0,1.0].pressure - stage.module.permeate_side[0.0,1.0].pressure, to_units=pyunits.bar))
#         print(f'{"RO Stage " + str(idx):<15s}{100*stage.module.recovery_vol_phase[0.0, "Liq"].value:<15.1f}{stage.module.rejection_phase_comp[0.0, "Liq", "NaCl"].value:<15.3f}{stage.module.area.value:<15.3f}{stage.module.width.value:<15.3f}{stage.module.length.value:<15.3f}{stage.module.feed_side.velocity[0, 0].value:<20.3f}{stage.module.feed_side.velocity[0, 1].value:<20.3f}{del_P-del_pi:<15.3f}')
#         # stage.display()
#     print('\n')
#     print(f'{"Pump Work":<20s}{value(pyunits.convert(m.fs.primary_pump.control_volume.work[0], to_units=pyunits.kW)):<10.3f}{"kW":<20s}')
#     print(f'{"Pump Delta P":<20s}{value(pyunits.convert(m.fs.primary_pump.control_volume.deltaP[0], to_units=pyunits.bar)):<10.3f}{"bar":<20s}')
#     print(f'{"Pump Head":<20s}{value(pyunits.convert(m.fs.primary_pump.control_volume.deltaP[0], to_units=pyunits.bar)/0.0981):<10.3f}{"ft":<20s}')
#     print('\n')
#     print(f'{"System Recovery":<20s}{100*m.fs.water_recovery.value:<10.1f}{str(pyunits.get_units(m.fs.water_recovery)):<20s}')
#     print(f'{"System SEC":<20s}{value(m.fs.costing.specific_energy_consumption):<10.3f}{str(pyunits.get_units(m.fs.costing.specific_energy_consumption)):<20s}')
#     print(f'{"System LCOW":<20s}{value(m.fs.costing.LCOW):<10.3f}{str(pyunits.get_units(m.fs.costing.LCOW)):<20s}')
    
# def main():
#     m = build_system(number_of_stages=2, booster_pumps = False)
#     display_system_build(m)
#     # set_ro_scaling_and_params(m)
#     # set_operating_conditions(m, Qin=100, Qout = None, Cin=0.95, water_recovery=None, primary_pump_pressure=12e5)
#     # set_ro_system_operating_conditions(m)
#     # init_system(m, verbose=True, solver=None)
    
#     # box_solve_problem(m)
#     # display_flow_table(m)
#     # display_system_metrics(m)

#     # optimize(m, water_recovery=0.8, fixed_pressure=None, Qprod = 77.77, mem_width=500)
#     # # optimize(m, water_recovery=0.9, fixed_pressure=None, Qprod = 77.77, velocity = 0.25)
#     # solve(m, solver=None, tee=False, raise_on_failure=True)
    
#     # display_flow_table(m)
#     # display_system_metrics(m)

    

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
    
    ROSystem(m).build(m.fs.ro)
    # m.fs.ro.display_system_build()
    print(m.fs.ro.report())
