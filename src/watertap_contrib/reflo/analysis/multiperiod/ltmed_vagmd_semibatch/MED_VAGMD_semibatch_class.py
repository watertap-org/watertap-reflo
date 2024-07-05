#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import logging
import pandas as pd
import numpy as np

# Pyomo imports
from pyomo.environ import (
    Var,
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# IDAES imports
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    UnitModelCostingBlock,
    useDefault,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
)
from idaes.core.util.scaling import (
    set_scaling_factor,
    get_scaling_factor,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
import idaes.logger as idaeslog

# WaterTAP imports
from watertap.core.solvers import get_solver

# Flowsheet function imports
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.analysis.multiperiod.ltmed_vagmd_semibatch.MED_VAGMD_flowsheet import (
    build_med_md_flowsheet,
    fix_dof_and_initialize,
)

__author__ = "Zhuoran Zhang"

solver = get_solver()
_logger = idaeslog.getLogger(__name__)


def get_variable_pairs(t1, t2):
    """
    This function returns pairs of variables that need to be connected across two time periods
    Args:
        t1: current time block
        t2: next time block
    """
    return [
        # Take MD feed properties from last step mixer
        (t1.fs.S2.MD_feed_state[0].temperature, t2.fs.vagmd.feed_props[0].temperature),
        (t1.fs.S2.MD_feed_state[0].pressure, t2.fs.vagmd.feed_props[0].pressure),
        (
            t1.fs.S2.MD_feed_state[0].flow_mass_phase_comp["Liq", "TDS"],
            t2.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "TDS"],
        ),
        (
            t1.fs.S2.MD_feed_state[0].flow_mass_phase_comp["Liq", "H2O"],
            t2.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "H2O"],
        ),
        # Inherit mixer properties from last step
        (
            t1.fs.S2.remained_liquid_state[0].temperature,
            t2.fs.M1.remained_liquid_state[0].temperature,
        ),
        (
            t1.fs.S2.remained_liquid_state[0].pressure,
            t2.fs.M1.remained_liquid_state[0].pressure,
        ),
        (
            t1.fs.S2.remained_liquid_state[0].flow_mass_phase_comp["Liq", "TDS"],
            t2.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq", "TDS"],
        ),
        (
            t1.fs.S2.remained_liquid_state[0].flow_mass_phase_comp["Liq", "H2O"],
            t2.fs.M1.remained_liquid_state[0].flow_mass_phase_comp["Liq", "H2O"],
        ),
        # Mixer volume
        (t1.fs.volume_in_tank, t2.fs.volume_in_tank_previous),
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
    **default** = False.""",
        ),
    )
    CONFIG.declare(
        "model_input",
        ConfigValue(
            default=useDefault,
            description="System input variables",
            doc="""Model inputs as a dictionary, which includes:
            See details in the model inputs table in documentation.
            """,
        ),
    )

    def build(self):
        super().build()

        self.mp = self.create_multiperiod_ltmed_vagmd_model(
            n_time_points=self.config.model_input["n_time_points"],
            med_feed_salinity=self.config.model_input["med_feed_salinity"],
            med_feed_temp=self.config.model_input["med_feed_temp"],
            med_steam_temp=self.config.model_input["med_steam_temp"],
            med_capacity=self.config.model_input["med_capacity"],
            med_recovery_ratio=self.config.model_input["med_recovery_ratio"],
            md_feed_flow_rate=self.config.model_input["md_feed_flow_rate"],
            md_cond_inlet_temp=self.config.model_input["md_cond_inlet_temp"],
            md_evap_inlet_temp=self.config.model_input["md_evap_inlet_temp"],
            md_module_type=self.config.model_input["md_module_type"],
            md_cooling_system_type=self.config.model_input["md_cooling_system_type"],
            md_cooling_inlet_temp=self.config.model_input["md_cooling_inlet_temp"],
            md_high_brine_salinity=self.config.model_input["md_high_brine_salinity"],
            dt=self.config.model_input["dt"],
            batch_volume=self.config.model_input["batch_volume"],
        )

        # Initialize and unfix dof for each period
        solver = get_solver()
        active_blks = self.mp.get_active_process_blocks()

        for blk in active_blks:
            fix_dof_and_initialize(
                m=blk,
            )
            result = solver.solve(blk)
            assert_optimal_termination(result)
            unfix_dof(m=blk)

        # Setup the stats of the first time period
        initial_md_feed_salinity = self.config.model_input["med_feed_salinity"] / (
            1 - self.config.model_input["med_recovery_ratio"]
        )
        self.add_initial_constraints(
            batch_volume=self.config.model_input["batch_volume"],
            md_feed_flow_rate=self.config.model_input["md_feed_flow_rate"],
            md_feed_salinity=initial_md_feed_salinity,
            dt=self.config.model_input["dt"],
        )

        # Setup the refilling phase stats
        self.add_refilling_phase()

        self.target_system_capacity = Var(
            initialize=2000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.day,
            doc="Target system capacity",
        )

        self.batch_water_production = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.day,
            doc="Water production during one batch",
        )

        ## Create expressions for water production during one batch (processing phase, refilling phase)
        # MED
        self.water_prod_med_process = sum(
            pyunits.convert(
                active_blks[i].fs.med.distillate_props[0].flow_vol_phase["Liq"]
                * active_blks[i].fs.dt,
                to_units=pyunits.m**3,
            )
            for i in range(self.config.model_input["n_time_points"])
        )

        self.water_prod_med_refill = pyunits.convert(
            active_blks[0].fs.med.distillate_props[0].flow_vol_phase["Liq"]
            * self.dt_refilling_phase,
            to_units=pyunits.m**3,
        )

        self.water_prod_med_batch = (
            self.water_prod_med_process + self.water_prod_med_refill
        )

        # VAGMD
        self.water_prod_vagmd_process = sum(
            pyunits.convert(
                active_blks[i].fs.vagmd.permeate_flux
                * active_blks[i].fs.vagmd.module_area
                * active_blks[i].fs.dt,
                to_units=pyunits.m**3,
            )
            for i in range(self.config.model_input["n_time_points"])
        )

        self.water_prod_vagmd_refill = pyunits.convert(
            active_blks[0].fs.vagmd.permeate_flux
            * active_blks[0].fs.vagmd.module_area
            * self.dt_refilling_phase,
            to_units=pyunits.m**3,
        )

        self.water_prod_vagmd_batch = (
            self.water_prod_vagmd_process + self.water_prod_vagmd_refill
        )

        # Total water production during one batch
        self.total_water_prod_batch = (
            self.water_prod_vagmd_batch + self.water_prod_med_batch
        )

        ## Water production capacity of the specified system
        self.module_capacity = pyunits.convert(
            self.total_water_prod_batch
            / (
                self.dt_refilling_phase
                + sum(
                    active_blks[i].fs.dt
                    for i in range(self.config.model_input["n_time_points"])
                )
            ),
            to_units=pyunits.m**3 / pyunits.day,
        )

    def create_multiperiod_ltmed_vagmd_model(
        self,
        n_time_points=2,
        med_feed_salinity=30,
        med_feed_temp=25,
        med_steam_temp=80,
        med_capacity=1,
        med_recovery_ratio=0.5,
        md_feed_flow_rate=600,
        md_cond_inlet_temp=25,
        md_evap_inlet_temp=80,
        md_module_type="AS26C7.2L",
        md_cooling_system_type="closed",
        md_cooling_inlet_temp=25,
        md_high_brine_salinity=False,
        dt=None,
        batch_volume=50,
    ):
        """
        This function creates a multi-period flowsheet object. This object contains
        a pyomo model with a block for each time instance.
        Returns:
            Object containing multi-period vagmd batch flowsheet model
        """

        # Check if the input configurations are valid
        if type(n_time_points) is not int or n_time_points < 1:
            raise ConfigurationError(
                f"The number of time steps '{n_time_points}' is not available."
                f"Please input a positive integer."
            )

        if md_module_type not in ["AS7C1.5L", "AS26C7.2L"]:
            raise ConfigurationError(
                f"The MD module type '{md_module_type}' is not available."
                f"Available options include 'AS7C1.5L' and 'AS26C7.2L'."
            )

        if md_cooling_system_type not in ["open", "closed"]:
            raise ConfigurationError(
                f"The cooling system type '{md_cooling_system_type}' is not available."
                f"Available options include 'open' and 'closed'."
            )

        input_variables = [
            "med_feed_salinity",
            "med_feed_temp",
            "med_steam_temp",
            "med_recovery_ratio",
            "md_feed_flow_rate",
            "md_cond_inlet_temp",
            "md_evap_inlet_temp",
            "batch_volume",
            "dt",
        ]
        input_values = [
            med_feed_salinity,
            med_feed_temp,
            med_steam_temp,
            med_recovery_ratio,
            md_feed_flow_rate,
            md_cond_inlet_temp,
            md_evap_inlet_temp,
            batch_volume,
            dt,
        ]
        input_ranges = [
            [30, 60],
            [15, 35],
            [60, 85],
            [0.3, 0.5],
            [400, 1100],
            [20, 30],
            [60, 80],
            [50, float("inf")],
            [0, float("inf")],
        ]
        for i in range(len(input_variables)):
            if (
                input_values[i] < input_ranges[i][0]
                or input_values[i] > input_ranges[i][1]
            ):
                raise ConfigurationError(
                    f"The input value for '{input_variables[i]}' is not valid."
                    f"The valid range is {input_ranges[i][0]} - {input_ranges[i][1]}."
                )

        if md_module_type == "AS7C1.5L" and md_high_brine_salinity:
            _logger.info(
                f"For module AS7C1.5L, when the final brine salinity is larger than 175.3 g/L,"
                f"the operational parameters will be fixed at nominal condition"
            )

        mp = MultiPeriodModel(
            n_time_points=n_time_points,
            process_model_func=build_med_md_flowsheet,
            linking_variable_func=get_variable_pairs,
            initialization_func=fix_dof_and_initialize,
            unfix_dof_func=unfix_dof,
            outlvl=logging.WARNING,
        )

        processing_options = {
            "med_feed_salinity": med_feed_salinity,
            "med_feed_temp": med_feed_temp,
            "med_steam_temp": med_steam_temp,
            "med_capacity": med_capacity,
            "med_recovery_ratio": med_recovery_ratio,
            "md_feed_flow_rate": md_feed_flow_rate,
            "md_cond_inlet_temp": md_cond_inlet_temp,
            "md_evap_inlet_temp": md_evap_inlet_temp,
            "md_module_type": md_module_type,
            "md_cooling_system_type": md_cooling_system_type,
            "md_cooling_inlet_temp": md_cooling_inlet_temp,
            "md_high_brine_salinity": md_high_brine_salinity,
            "dt": dt,
            "batch_volume": batch_volume,
            "phase": "processing",
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
        self.dt_refilling_phase = Var(
            initialize=120,
            bounds=(None, None),
            units=pyunits.s,
            doc="Time interval of the refilling phase",
        )

        # Get the last period for refilling phase calculations
        first_blk = self.mp.get_active_process_blocks()[0]
        last_blk = self.mp.get_active_process_blocks()[-1]

        @self.Constraint(doc="Calculate time interval of the refilling phase")
        def eq_dt_refilling(b):
            return (
                first_blk.fs.med.brine_props[0].flow_vol_phase["Liq"]
                - first_blk.fs.vagmd.permeate_props[0].flow_vol_phase["Liq"]
            ) * b.dt_refilling_phase == pyunits.convert(
                self.config.model_input["batch_volume"] * pyunits.L,
                to_units=pyunits.m**3,
            ) - last_blk.fs.volume_in_tank

        return

    def add_initial_constraints(
        self, batch_volume, md_feed_flow_rate, md_feed_salinity, dt
    ):
        blk = self.mp.get_active_process_blocks()[0]
        # VAGMD feed temperature and concentration are the same as
        # MED brine at the beginning of the process
        @self.Constraint(
            doc="Initialize VAGMD feed temperature with MED brine temperature"
        )
        def eq_initial_vagmd_temp(b):
            return (
                blk.fs.vagmd.feed_props[0].temperature
                == blk.fs.med.brine_props[0].temperature
            )

        @self.Constraint(
            doc="Initialize VAGMD feed concentration with MED brine concentration"
        )
        def eq_initial_vagmd_conc(b):
            return (
                blk.fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == md_feed_salinity
            )

        @self.Constraint(doc="Initialize VAGMD feed flow rate with defined value")
        def eq_initial_vagmd_flow_rate(b):
            return blk.fs.vagmd.feed_props[0].flow_vol_phase["Liq"] == pyunits.convert(
                md_feed_flow_rate * pyunits.L / pyunits.h,
                to_units=pyunits.m**3 / pyunits.s,
            )

        blk.fs.vagmd.feed_props[0].pressure.fix(101325)

        blk.fs.M1.remained_liquid_state.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): batch_volume
                / 1000
                / blk.fs.dt.value,  # convert to m3/s
                ("conc_mass_phase_comp", ("Liq", "TDS")): md_feed_salinity
                * pyunits.kg
                / pyunits.m**3,
                ("temperature", None): 35 + 273.15,
                # feed flow is at atmospheric pressure
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        blk.fs.volume_in_tank_previous.fix(0)

        return

    def add_costing_packages(self):
        """
        This function scales the med-md semibatch module to a specific target capacity,
        and adds costing packages accordingly
        """
        self.costing = TreatmentCosting()
        self.costing.base_currency = pyunits.USD_2020

        # The costing model is built upon the last time step
        blk = self.mp.get_active_process_blocks()[-1].fs

        # Create cost block for MED
        blk.med.costing = UnitModelCostingBlock(flowsheet_costing_block=self.costing)

        # Update capacity in MED cost functions
        med_costing = blk.med.costing
        lt_med_params = med_costing.costing_package.lt_med_surrogate

        med_costing.capacity = pyunits.convert(
            self.water_prod_med_batch
            / self.total_water_prod_batch
            * self.target_system_capacity,
            to_units=pyunits.m**3 / pyunits.day,
        )

        med_costing.med_specific_cost_constraint.deactivate()
        med_costing.med_specific_cost_constraint = Constraint(
            expr=med_costing.med_specific_cost
            == (
                lt_med_params.med_sys_A_coeff
                * med_costing.capacity**lt_med_params.med_sys_B_coeff
            )
        )

        med_costing.membrane_system_cost_constraint.deactivate()
        med_costing.membrane_system_cost_constraint = Constraint(
            expr=med_costing.membrane_system_cost
            == med_costing.capacity
            * (
                med_costing.med_specific_cost
                * (1 - lt_med_params.cost_fraction_evaporator)
            )
        )

        # Costing package of VAGMD surrogate unit model is applied here,
        # while the technical model used here is VAGMD base model, so the
        # following missing items need to be created:

        # Update capacity in VAGMD cost functions
        blk.vagmd.system_capacity = pyunits.convert(
            self.water_prod_vagmd_batch
            / self.total_water_prod_batch
            * self.target_system_capacity,
            to_units=pyunits.m**3 / pyunits.day,
        )

        # Number of vagmd module required = VAGMD system capacity / VAGMD module capacity
        blk.vagmd.num_modules = blk.vagmd.system_capacity / (
            self.module_capacity
            * self.water_prod_vagmd_batch
            / self.total_water_prod_batch
        )

        # Create energy flow items
        blk.vagmd.thermal_power_requirement = (
            blk.vagmd.thermal_power
            / pyunits.convert(
                blk.vagmd.permeate_flux * blk.vagmd.module_area,
                to_units=pyunits.m**3 / pyunits.h,
            )
            * pyunits.convert(
                blk.vagmd.system_capacity, to_units=pyunits.m**3 / pyunits.h
            )
        )

        blk.vagmd.elec_power_requirement = (
            (blk.vagmd.feed_pump_power_elec + blk.vagmd.cooling_pump_power_elec)
            / pyunits.convert(
                blk.vagmd.permeate_flux * blk.vagmd.module_area,
                to_units=pyunits.m**3 / pyunits.h,
            )
            * pyunits.convert(
                blk.vagmd.system_capacity, to_units=pyunits.m**3 / pyunits.h
            )
        )

        # Create cost block for VAGMD
        blk.vagmd.costing = UnitModelCostingBlock(flowsheet_costing_block=self.costing)

        # Add LCOW component
        self.costing.cost_process()
        self.costing.add_annual_water_production(self.target_system_capacity)
        self.costing.add_LCOW(self.target_system_capacity)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if get_scaling_factor(self.target_system_capacity) is None:
            set_scaling_factor(self.target_system_capacity, 1e-2)

        if get_scaling_factor(self.batch_water_production) is None:
            set_scaling_factor(self.batch_water_production, 1e0)

        if get_scaling_factor(self.dt_refilling_phase) is None:
            set_scaling_factor(self.dt_refilling_phase, 1e-2)

    def get_model_performance(self):
        """
        This function returns the overall performance of the batch operation in a dictionary,
        and the timewise system performance in a pandas dataframe
        """
        blks = self.mp.get_active_process_blocks()
        n_time_points = len(blks)

        recovery_ratio = (
            1
            - blks[0].fs.med.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].value
            / blks[-1]
            .fs.vagmd.evaporator_out_props[0]
            .conc_mass_phase_comp["Liq", "TDS"]
            .value
        )

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
        mixed_temp_in_tank = [
            value(blks[i].fs.M1.mixed_state[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        feed_temp = [
            value(blks[i].fs.vagmd.feed_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        tank_volume = [value(blks[i].fs.volume_in_tank) for i in range(n_time_points)]
        med_brine_temp = [
            value(blks[i].fs.med.brine_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        md_brine_temp = [
            value(blks[i].fs.vagmd.evaporator_out_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        md_brine_conc = [
            value(
                blks[i]
                .fs.vagmd.evaporator_out_props[0]
                .conc_mass_phase_comp["Liq", "TDS"]
            )
            for i in range(n_time_points)
        ]

        remained_flow_rate = [
            value(blks[i].fs.M1.remained_liquid_state[0].flow_vol_phase["Liq"])
            for i in range(n_time_points)
        ]
        remained_volume = [
            remained_flow_rate[i] * value(blks[i].fs.dt * 1000)
            for i in range(n_time_points)
        ]
        mixed_volume = [
            value(blks[i].fs.M1.mixed_state[0].flow_vol_phase["Liq"])
            * value(blks[i].fs.dt * 1000)
            for i in range(n_time_points)
        ]

        STEC_vagmd = [
            value(blks[i].fs.vagmd_specific_energy_consumption_thermal)
            for i in range(n_time_points)
        ]

        STEC_med = [
            value(blks[i].fs.med.specific_energy_consumption_thermal)
            for i in range(n_time_points)
        ]

        SEEC_vagmd = [
            value(blks[i].fs.vagmd_specific_energy_consumption_electric)
            for i in range(n_time_points)
        ]

        water_prod_med = [
            value(
                blks[i].fs.med.distillate_props[0].flow_vol_phase["Liq"] * blks[i].fs.dt
            )  # m3
            for i in range(n_time_points)
        ]

        water_prod_vagmd = [
            value(
                blks[i].fs.vagmd.permeate_flux
                * blks[i].fs.vagmd.module_area
                * blks[i].fs.dt
                / 1000
                / 3600
            )  # m3
            for i in range(n_time_points)
        ]

        vagmd_thermal_power = [
            value(blks[i].fs.vagmd.thermal_power) for i in range(n_time_points)
        ]

        med_thermal_power = [
            value(blks[i].fs.med.thermal_power_requirement)
            for i in range(n_time_points)
        ]

        time_series = np.array(
            [
                time_period,
                t_minutes,
                feed_salinity,
                md_brine_conc,
                med_brine_temp,
                md_brine_temp,
                mixed_temp_in_tank,
                feed_temp,
                tank_volume,
                STEC_vagmd,
                STEC_med,
                SEEC_vagmd,
                water_prod_med,
                water_prod_vagmd,
                vagmd_thermal_power,
                med_thermal_power,
            ]
        )

        data_table = pd.DataFrame(
            data=time_series,
            index=[
                "Step",
                "Operation time (min)",
                "MD Feed salinity (g/L)",
                "MD brine concentration (g/L)",
                "MED brine temperature (C)",
                "MD brine temperature (C)",
                "Mixed temperature in tank (C)",
                "MD Feed temperature (C)",
                "Tank volume (m3)",
                "STEC in VAGMD (kWh-th/m3)",
                "STEC in MED (kWh-th/m3)",
                "SEEC in VAGMD (kWh-e/m3)",
                "MED water production (m3)",
                "VAGMD water production (m3)",
                "VAGMD thermal power requirement (kW)",
                "MED thermal power requirement (kW)",
            ],
        )

        # Aggregated system performance
        time_refilling_phase = self.dt_refilling_phase.value / 60  # min
        operation_time_batch = self.dt_refilling_phase.value / 60 + t_minutes[-1]  # min
        water_prod_refilling = (
            self.water_prod_med_refill() + self.water_prod_vagmd_refill()
        )
        water_prod_batch = (
            water_prod_refilling
            + self.water_prod_med_process()
            + self.water_prod_vagmd_process()
        )
        thermal_power_batch = (
            (sum(vagmd_thermal_power) + sum(med_thermal_power))
            * blks[0].fs.dt.value
            / 60
            + (
                blks[0].fs.vagmd.thermal_power.value
                + blks[0].fs.med.thermal_power_requirement.value
            )
            * time_refilling_phase
        ) / operation_time_batch

        STEC_batch = (
            thermal_power_batch * (operation_time_batch / 60) / water_prod_batch
        )

        overall_performance = {
            "Overall recovery ratio": recovery_ratio,
            "Time interval of the processing phase (min)": t_minutes[-1],
            "Time interval of the refilling phase (min)": time_refilling_phase,
            "Total operation time of one batch (min)": operation_time_batch,
            "Total water production during refilling phase (m3)": water_prod_refilling,
            "Total water production during one batch (m3)": self.total_water_prod_batch(),
            "Production capacity (m3/day)": self.module_capacity(),
            "Average thermal power requirement during one batch (kW)": thermal_power_batch,
            "Average specifc thermal energy consumption during one batch (kWh/m3)": STEC_batch,
        }

        return overall_performance, data_table

    def get_costing_performance(self):
        """
        This function returns the cost performance of a system with a target capacity
        scaled upon the specified MED-MD module in the unit model
        """
        # The costing model is built upon the last time step
        blk = self.mp.get_active_process_blocks()[-1].fs

        cost_performance = {
            "Target system capacity (m3/day)": self.target_system_capacity.value,
            "Scaled MED system capacity (m3/day)": blk.med.costing.capacity(),
            "Scaled VAGMD system capacity (m3/day)": blk.vagmd.system_capacity(),
            "Number of VAGMD modules required": blk.vagmd.num_modules(),
            "CAPEX of MED ($)": blk.med.costing.capital_cost.value,
            "CAPEX of VAGMD ($)": blk.vagmd.costing.capital_cost.value,
            "Fixed OPEX of MED ($)": blk.med.costing.fixed_operating_cost.value,
            "Fixed OPEX of VAGMD ($)": blk.vagmd.costing.fixed_operating_cost.value,
            "Annual heat cost ($)": self.costing.aggregate_flow_costs["heat"].value,
            "Annual electricity cost ($)": self.costing.aggregate_flow_costs[
                "electricity"
            ].value,
            "Overall LCOW ($/m3)": self.costing.LCOW(),
        }

        return cost_performance
