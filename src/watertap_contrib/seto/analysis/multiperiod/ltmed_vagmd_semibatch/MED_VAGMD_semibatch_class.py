import logging
import pandas as pd
import numpy as np

# Pyomo imports
from pyomo.environ import (
    Var,
    Constraint,
    value,
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.scaling import (
    set_scaling_factor,
    get_scaling_factor,
    calculate_scaling_factors,
    constraint_scaling_transform,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    UserModelError,
    InitializationError,
)
import idaes.core.util.scaling as iscale
from idaes.core import UnitModelCostingBlock
import idaes.logger as idaeslog
from idaes.core.solvers.get_solver import get_solver

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

# Flowsheet function imports
from watertap_contrib.seto.analysis.multiperiod.ltmed_vagmd_semibatch.MED_VAGMD_flowsheet import (
    build_med_md_flowsheet,
    fix_dof_and_initialize,
)

__author__ = "Zhuoran Zhang"

_log = idaeslog.getLogger(__name__)
solver = get_solver()

def get_variable_pairs(t1, t2):
    """
    This function returns paris of variables that need to be connected across two time periods
    Args:
        t1: current time block
        t2: next time block
    Returns:
        None
    """
    return [
        # Take MD feed properties from last step mixer
        (t1.fs.S2.MD_feed_state[0].temperature, t2.fs.vagmd.feed_props[0].temperature),

        (t1.fs.S2.MD_feed_state[0].pressure, t2.fs.vagmd.feed_props[0].pressure),

        (t1.fs.S2.MD_feed_state[0].flow_mass_phase_comp["Liq","TDS"],
         t2.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq","TDS"]),

        (t1.fs.S2.MD_feed_state[0].flow_mass_phase_comp["Liq","H2O"],
         t2.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq","H2O"]),

        # Inherit mixer properties from last step
        (t1.fs.S2.remained_liquid_state[0].temperature, t2.fs.M1.remained_liquid_state[0].temperature),

        (t1.fs.S2.remained_liquid_state[0].pressure, t2.fs.M1.remained_liquid_state[0].pressure),

        (t1.fs.S2.remained_liquid_state[0].flow_mass_phase_comp["Liq","TDS"],
         t2.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq","TDS"]),
         
        (t1.fs.S2.remained_liquid_state[0].flow_mass_phase_comp["Liq","H2O"],
         t2.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq","H2O"]),

        # Mixer volume
        (t1.fs.volume_in_tank, t2.fs.volume_in_tank_previous)

        ]

def unfix_dof(m):
    """
    This function unfixes a few degrees of freedom for optimization
    Args:
        m: object containing the integrated nuclear plant flowsheet
    Returns:
        None
    """

    m.fs.vagmd.feed_props[0].temperature.unfix()
    m.fs.vagmd.feed_props[0].pressure.unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    m.fs.M1.remained_liquid_state[0].temperature.unfix()
    m.fs.M1.remained_liquid_state[0].pressure.unfix()
    m.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    m.fs.volume_in_tank_previous.unfix()

    return

@declare_process_block_class("MEDVAGMDsemibatch")
class MEDVAGMDsemibatchData(UnitModelBlockData):

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
        ),
    )
    CONFIG.declare(
        "model_input",
        ConfigValue(
            default=useDefault,
            description="VAGMD model input variables",
            doc="""VAGMD model inputs as a dictionary, which includes:
            dt:              Time interval of simulation, None or Positive values
            system_capacity: System capacity, m3/day
            feed_flow_rate:  Feed flow rate, 400 - 1100 L/h
            evap_inlet_temp: Evaporator inlet temperature, 60 - 80 deg C
            cond_inlet_temp: Condenser inlet temperature, 20 - 30 deg C
            feed_temp:       Feed water temperature, 20 - 30 deg C
            feed_salinily:   Feed water salinity, 35 - 292 g/L
            initial_batch_volume: Batch volume of the feed water, > 50 L
            recovery_ratio:  Target recovery ratio of the batch operation
            module_type:     Aquastill MD module type, "AS7C1.5L" or "AS26C7.2L",
                            "AS7C1.5L" yields maximum permeate produtivity,
                            "AS26C7.2L" yields maximum thermal efficiency.
            cooling_system_type:  Cooling system type, "open" or "closed"
            cooling_inlet_temp:   Cooling water temperature, 20-30 deg C
                                only required when cooling system type is "open""",
        ),
    )

    def build(self):
        super().build()

        self.mp = self.create_multiperiod_ltmed_vagmd_model(
                n_time_points= self.config.model_input["n_time_points"],
                med_feed_salinity = self.config.model_input["med_feed_salinity"],
                med_feed_temp     = self.config.model_input["med_feed_temp"],
                med_steam_temp    = self.config.model_input["med_steam_temp"],
                med_capacity      = self.config.model_input["med_capacity"],
                med_recovry_ratio = self.config.model_input["med_recovry_ratio"],
                md_feed_flow_rate = self.config.model_input["md_feed_flow_rate"],
                dt = self.config.model_input["dt"],
                batch_volume= self.config.model_input["batch_volume"],)

        # Setup the stats of the first time period 
        self.add_initial_constraints(batch_volume = self.config.model_input["batch_volume"],
                                     md_feed_flow_rate=self.config.model_input["md_feed_flow_rate"],
                                     dt = self.config.model_input["dt"],)

        # Setup the refilling phase stats
        self.add_refilling_phase()

        self.vagmd_thermal_power_requirement = Var(
            initialize=2e5,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW-th)",
        )

        self.vagmd_elec_power_requirement = Var(
            initialize=300,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electric power requirement (kW-e)",
        )

        # @self.Constraint(
        #     doc="Calculate the overall thermal power requirement for vagmd through all periods"
        # )
        # def eq_thermal_power_requirement(b):
        #     return b.vagmd_thermal_power_requirement == (
        #         self.mp.get_active_process_blocks()[
        #             -1
        #         ].fs.specific_energy_consumption_thermal
        #         * pyunits.convert(
        #             b.system_capacity, to_units=pyunits.m**3 / pyunits.h
        #         )
        #     )

        # @self.Constraint(
        #     doc="Calculate the overall electric power requirement through all periods"
        # )
        # def eq_elec_power_requirement(b):
        #     return b.overall_elec_power_requirement == (
        #         self.mp.get_active_process_blocks()[
        #             -1
        #         ].fs.specific_energy_consumption_electric
        #         * pyunits.convert(
        #             b.system_capacity, to_units=pyunits.m**3 / pyunits.h
        #         )
        #     )

        # super().calculate_scaling_factors()

        # if iscale.get_scaling_factor(self.overall_thermal_power_requirement) is None:
        #     iscale.set_scaling_factor(self.overall_thermal_power_requirement, 1e-5)

        # if iscale.get_scaling_factor(self.overall_elec_power_requirement) is None:
        #     iscale.set_scaling_factor(self.overall_elec_power_requirement, 1e-3)


    def create_multiperiod_ltmed_vagmd_model(
            self,
            n_time_points=2,
            med_feed_salinity = 30,
            med_feed_temp     = 25,
            med_steam_temp    = 80,
            med_capacity      = 1,
            med_recovry_ratio = 0.5,
            md_feed_flow_rate = 600,
            dt = None,
            batch_volume = 50,
        ):
        """
        This function creates a multi-period pv battery flowsheet object. This object contains 
        a pyomo model with a block for each time instance.
        Args:
            n_time_points: Number of time blocks to create
        Returns:
            Object containing multi-period vagmd batch flowsheet model
        """
        mp = MultiPeriodModel(
            n_time_points=n_time_points,
            process_model_func=build_med_md_flowsheet,
            linking_variable_func=get_variable_pairs,
            initialization_func=fix_dof_and_initialize,
            unfix_dof_func=unfix_dof,
            outlvl=logging.WARNING,
        )

        processing_options = {"med_feed_salinity": med_feed_salinity, 
                        "med_feed_temp": med_feed_temp,
                        "med_steam_temp": med_steam_temp, 
                        "med_capacity": med_capacity,
                        "med_recovry_ratio": med_recovry_ratio,
                        'md_feed_flow_rate': md_feed_flow_rate,
                        'dt': dt,
                        'batch_volume': batch_volume,
                        'phase': 'processing',
                        }
        
        model_options = {t: processing_options for t in range(n_time_points)}

        # create the multiperiod object
        mp.build_multi_period_model(
            model_data_kwargs=model_options,
            flowsheet_options=processing_options,
            initialization_options=None,
            unfix_dof_options=None,
            )

        return mp

    def add_refilling_phase(self):
        """
        This function adds the refilling phase of the batch operation
        and initializes it
        """    
        self.dt_discharging_phase = Var(initialize = 60,
                                      bounds = (0, None),
                                      units = pyunits.s,
                                      doc = 'Time interval of discharging the tank')

        self.dt_refilling_phase = Var(initialize = 120,
                                      bounds = (None, None),
                                      units = pyunits.s,
                                      doc = 'Time interval of the refilling phase')

        # Get the last period for refilling phase calculations
        first_blk = self.mp.get_active_process_blocks()[0]
        last_blk = self.mp.get_active_process_blocks()[-1]
        @self.Constraint(doc = 'Calculate time interval of the refilling phase')
        def eq_dt_refilling(b):
            return ((first_blk.fs.med.brine_props[0].flow_vol_phase["Liq"]
                    - first_blk.fs.vagmd.permeate_props[0].flow_vol_phase["Liq"])
                    * b.dt_refilling_phase
                    == pyunits.convert(self.config.model_input["batch_volume"] * pyunits.L,
                                       to_units = pyunits.m**3)
                       - last_blk.fs.volume_in_tank )
        
        return

    def add_initial_constraints(self, batch_volume, md_feed_flow_rate, dt):
        blk = self.mp.get_active_process_blocks()[0]
        # VAGMD feed temperature and concentration are the same as
        # MED brine at the beginning of the process
        blk.fs.vagmd.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    md_feed_flow_rate * pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): 70,
                ("temperature", None): 35 + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        blk.fs.M1.remained_liquid_state.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 50/1000/dt,
                ("conc_mass_phase_comp", ("Liq", "TDS")): 70 * pyunits.kg / pyunits.m**3,
                ("temperature", None): 35+273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,  
        )

        blk.fs.volume_in_tank_previous.fix(0)

        return 

    def get_model_performance(self):
        """
        This function returns the overall performance of the batch operation in a dictionary,
        and the timewise system performance in a pandas dataframe
        """
        blks = self.mp.get_active_process_blocks()
        n_time_points = len(blks)

        recovery_ratio = (1 - blks[0].fs.med.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].value / 
                            blks[-1].fs.vagmd.evaporator_out_props[0].conc_mass_phase_comp["Liq", "TDS"].value)

        # Water production during one batch (processing phase + refilling phase)
        # batch_water_production = 

        overall_performance = {
            "Overall recovery ratio": recovery_ratio,
            # "Production capacity (L/h)": production_ratio,
            # "Production capacity (m3/day)": production_rate / 1000 * 24,
            # "Gain output ratio": value(blks[-1].fs.gain_output_ratio),
            # "Specific thermal energy consumption (kWh/m3)": value(
            #     blks[-1].fs.specific_energy_consumption_thermal
            # ),
            # "Specific electric energy consumption (kWh/m3)": value(
            #     blks[-1].fs.specific_energy_consumption_electric
            # ),
        }

        time_period = [i for i in range(n_time_points)]
        
        t_minutes = []
        for t in time_period:
            if not t_minutes:
                t_minutes.append(value(blks[t].fs.dt) / 60)
            else:
                t_minutes.append(t_minutes[-1] + value(blks[t].fs.dt) / 60)

        feed_salinity = [
            value(blks[i].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"])
            for i in range(n_time_points)
        ]
        feed_temp = [
            value(blks[i].fs.vagmd.feed_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        tank_volume = [
            value(blks[i].fs.volume_in_tank)
            for i in range(n_time_points)
        ]
        med_brine_temp =[
            value(blks[i].fs.med.brine_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        md_brine_temp =[
            value(blks[i].fs.vagmd.evaporator_out_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]

        remained_flow_rate = [
            value(blks[i].fs.M1.remained_liquid_state[0].flow_vol_phase["Liq"] )
            for i in range(n_time_points)
        ]
        remained_volume = [ remained_flow_rate[i] * 
            value(blks[i].fs.dt * 1000)
            for i in range(n_time_points)
        ]
        mixed_volume = [ value(blks[i].fs.M1.mixed_state[0].flow_vol_phase["Liq"] )
                        * value(blks[i].fs.dt * 1000)
            for i in range(n_time_points)
        ]


        time_series = np.array(
            [
                time_period,
                t_minutes,
                feed_salinity,
                med_brine_temp,
                md_brine_temp,
                feed_temp,
                tank_volume,
                remained_volume,
                mixed_volume,
            ]
        )

        data_table = pd.DataFrame(
            data=time_series,
            index=[
                "Step",
                "Operation time (min)",
                "Feed salinity (g/L)",
                "MED brine temperature",
                "MD brine temperature",
                "Feed temperature (C)",
                "Tank volume",
                "remained volume",
                "mixed_volume"
            ],
        )

        return overall_performance, data_table

if __name__ == "__main__":

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    model_input = {
        "n_time_points": 2,
        "med_feed_salinity": 35,
        "med_feed_temp": 24,
        "med_steam_temp": 80,
        "med_capacity": 1,
        "med_recovry_ratio": 0.5,
        "md_feed_flow_rate": 600,
        "dt": 60,
        "batch_volume": 50,
    }
    m.fs.semibatch = VAGMDbatchSurrogate(model_input=model_input)

    mp = create_multiperiod_ltmed_vagmd_model(
        n_time_points= 2,
        med_feed_salinity = 35,
        med_feed_temp     = 25,
        med_steam_temp    = 80,
        med_capacity      = 1,
        med_recovry_ratio = 0.5,
        md_feed_flow_rate = 600,
        dt = 60,
        batch_volume=50,)
    
    # model_debug(mp)

    # check_jac(mp)    
    results = solver.solve(mp)
    assert_optimal_termination(results)

    data_table = get_model_performance(mp)

    # blk = mp.get_active_process_blocks()[-1]
    # print(blk.fs.M1.mixed_state[0].flow_vol_phase["Liq"].value * 1000 * 3600)
    # print(blk.fs.S2.remained_liquid_state[0].flow_vol_phase["Liq"].value * 1000 * 3600)
    # print(blk.fs.M1.remained_liquid_state[0].flow_vol_phase["Liq"].value * 1000 * 3600)
    # print(blk.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].value* 1000 * 3600)
    