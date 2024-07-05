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
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

# WaterTAP imports
from watertap.core.solvers import get_solver

# Flowsheet function imports
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_design_model import (
    get_n_time_points,
)


__author__ = "Zhuoran Zhang"

# Set up logger
_logger = idaeslog.getLogger(__name__)


def get_vagmd_batch_variable_pairs(t1, t2):
    """
    This function returns pairs of variables that need to be connected across two time periods

    Args:
        t1: current time block
        t2: next time block

    Returns:
        None
    """
    return [
        (
            t1.fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            t2.fs.pre_feed_salinity,
        ),
        (t1.fs.vagmd.feed_props[0].temperature, t2.fs.pre_feed_temperature),
        (t1.fs.vagmd.evaporator_out_props[0].temperature, t2.fs.pre_evap_out_temp),
        (
            t1.fs.vagmd.permeate_props[0].flow_vol_phase["Liq"],
            t2.fs.pre_permeate_flow_rate,
        ),
        (t1.fs.acc_thermal_energy, t2.fs.pre_acc_thermal_energy),
        (t1.fs.acc_cooling_energy, t2.fs.pre_acc_cooling_energy),
        (t1.fs.acc_distillate_volume, t2.fs.pre_acc_distillate_volume),
        (t1.fs.acc_electric_energy, t2.fs.pre_acc_electric_energy),
        (t1.fs.vagmd.thermal_power, t2.fs.pre_thermal_power),
        (t1.fs.vagmd.cooling_power_thermal, t2.fs.pre_cooling_power),
        (t1.fs.vagmd.feed_pump_power_elec, t2.fs.pre_feed_pump_power_elec),
        (t1.fs.vagmd.cooling_pump_power_elec, t2.fs.pre_cooling_pump_power_elec),
    ]


def unfix_dof(m, feed_flow_rate):
    """
    This function unfixes a few degrees of freedom for optimization

    Args:
        feed_flow_rate: Feed flow rate after the adjustment in the model configuration
                        of AS7C1.5L with high brine salinty

    Returns:
        None
    """
    m.fs.pre_feed_temperature.unfix()
    m.fs.pre_permeate_flow_rate.unfix()
    m.fs.acc_distillate_volume.unfix()
    m.fs.acc_thermal_energy.unfix()
    m.fs.pre_thermal_power.unfix()
    m.fs.acc_cooling_energy.unfix()
    m.fs.pre_cooling_power.unfix()
    m.fs.acc_electric_energy.unfix()
    m.fs.pre_feed_pump_power_elec.unfix()
    m.fs.pre_cooling_pump_power_elec.unfix()

    m.fs.vagmd.feed_props[0].temperature.unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].fix(
        pyunits.convert(
            feed_flow_rate * pyunits.L / pyunits.h, to_units=pyunits.m**3 / pyunits.s
        )
    )

    return


@declare_process_block_class("VAGMDbatchSurrogate")
class VAGMDbatchSurrogateData(UnitModelBlockData):

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

        self.system_capacity = (
            self.config.model_input["system_capacity"] * pyunits.m**3 / pyunits.day
        )

        self.mp = self.create_multiperiod_vagmd_batch_model(
            dt=self.config.model_input["dt"],
            system_capacity=self.config.model_input["system_capacity"],
            feed_flow_rate=self.config.model_input["feed_flow_rate"],
            evap_inlet_temp=self.config.model_input["evap_inlet_temp"],
            cond_inlet_temp=self.config.model_input["cond_inlet_temp"],
            feed_temp=self.config.model_input["feed_temp"],
            feed_salinity=self.config.model_input["feed_salinity"],
            initial_batch_volume=self.config.model_input["initial_batch_volume"],
            recovery_ratio=self.config.model_input["recovery_ratio"],
            module_type=self.config.model_input["module_type"],
            cooling_system_type=self.config.model_input["cooling_system_type"],
            cooling_inlet_temp=self.config.model_input["cooling_inlet_temp"],
        )

        self.overall_thermal_power_requirement = Var(
            initialize=2e5,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW-th)",
        )

        self.overall_elec_power_requirement = Var(
            initialize=300,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electric power requirement (kW-e)",
        )

        @self.Constraint(
            doc="Calculate the overall thermal power requirement through all periods"
        )
        def eq_thermal_power_requirement(b):
            return b.overall_thermal_power_requirement == (
                self.mp.get_active_process_blocks()[
                    -1
                ].fs.specific_energy_consumption_thermal
                * pyunits.convert(
                    b.system_capacity, to_units=pyunits.m**3 / pyunits.h
                )
            )

        @self.Constraint(
            doc="Calculate the overall electric power requirement through all periods"
        )
        def eq_elec_power_requirement(b):
            return b.overall_elec_power_requirement == (
                self.mp.get_active_process_blocks()[
                    -1
                ].fs.specific_energy_consumption_electric
                * pyunits.convert(
                    b.system_capacity, to_units=pyunits.m**3 / pyunits.h
                )
            )

        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.overall_thermal_power_requirement) is None:
            iscale.set_scaling_factor(self.overall_thermal_power_requirement, 1e-5)

        if iscale.get_scaling_factor(self.overall_elec_power_requirement) is None:
            iscale.set_scaling_factor(self.overall_elec_power_requirement, 1e-3)

    def create_multiperiod_vagmd_batch_model(
        self,
        dt=None,
        system_capacity=1000,  # m3/day
        feed_flow_rate=600,  # 400 - 1100 L/h
        evap_inlet_temp=80,  # 60 - 80 deg C
        cond_inlet_temp=25,  # 20 - 30 deg C
        feed_temp=25,  # 20 - 30 deg C
        feed_salinity=35,  # 35 - 292 g/L
        initial_batch_volume=50,  # > 50 L
        recovery_ratio=0.5,  # -
        module_type="AS7C1.5L",
        cooling_system_type="closed",
        cooling_inlet_temp=25,  # 20 - feed_temp deg C, not required when cooling system type is "closed"
    ):
        """
        This function creates a multi-period vagmd batch flowsheet object. This object contains
        a pyomo model with a block for each time instance.

        Args:
            dt:              Time interval of simuulation (s)
            system_capacity: System capacity (m3/day)
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
                                only required when cooling system type is "open"
        Returns:
            Object containing multi-period vagmd batch flowsheet model
        """

        # Check if the input configurations are valid

        if module_type not in ["AS7C1.5L", "AS26C7.2L"]:
            raise ConfigurationError(
                f"The MD module type '{module_type}' is not available."
                f"Available options include 'AS7C1.5L' and 'AS26C7.2L'."
            )

        if cooling_system_type not in ["open", "closed"]:
            raise ConfigurationError(
                f"The cooling system type '{cooling_system_type}' is not available."
                f"Available options include 'open' and 'closed'."
            )

        if cooling_system_type == "open" and (
            cooling_inlet_temp > feed_temp or cooling_inlet_temp < 20
        ):
            raise ConfigurationError(
                f"In open circuit cooling system, the valid cooling water temperature is 20 - {feed_temp} deg C"
            )

        max_allowed_brine_salinity = {"AS7C1.5L": 292.2, "AS26C7.2L": 245.5}
        input_variables = [
            "system_capacity",
            "feed_flow_rate",
            "evap_inlet_temp",
            "cond_inlet_temp",
            "feed_temp",
            "feed_salinity",
            "initial_batch_volume",
            "recovery_ratio",
        ]
        input_values = [
            system_capacity,
            feed_flow_rate,
            evap_inlet_temp,
            cond_inlet_temp,
            feed_temp,
            feed_salinity,
            initial_batch_volume,
            recovery_ratio,
        ]
        input_ranges = [
            [0, float("inf")],
            [400, 1100],
            [60, 80],
            [20, 30],
            [20, 30],
            [35, max_allowed_brine_salinity[module_type]],
            [50, float("inf")],
            [0, 1],
        ]
        for i in range(len(input_variables)):
            if (
                input_values[i] < input_ranges[i][0]
                or input_values[i] > input_ranges[i][1]
            ):
                raise ConfigurationError(
                    f"The input variable '{input_variables[i]}' is not valid."
                    f"The valid range is {input_ranges[i][0]} - {input_ranges[i][1]}."
                )

        final_brine_salinity = feed_salinity / (1 - recovery_ratio)  # g/L
        max_allowed_recovery_ratio = (
            1 - feed_salinity / max_allowed_brine_salinity[module_type]
        )

        if recovery_ratio > max_allowed_recovery_ratio:
            raise ConfigurationError(
                f"The maximum recovery ratio allowed for module {module_type} with a feed"
                f"salinity of {feed_salinity} is {max_allowed_recovery_ratio}."
            )

        if module_type == "AS7C1.5L" and final_brine_salinity > 175.3:
            feed_flow_rate = 1100
            cond_inlet_temp = 25
            evap_inlet_temp = 80

        _logger.info(
            f"For module AS7C1.5L, when the final brine salinity is larger than 175.3 g/L,"
            f"the operational parameters will be fixed at nominal condition"
        )

        # Calculate the number of periods to reach target recovery rate by solving the system first
        n_time_points = get_n_time_points(
            dt=dt,
            feed_flow_rate=feed_flow_rate,
            evap_inlet_temp=evap_inlet_temp,
            cond_inlet_temp=cond_inlet_temp,
            feed_temp=feed_temp,
            feed_salinity=feed_salinity,
            recovery_ratio=recovery_ratio,
            initial_batch_volume=initial_batch_volume,
            module_type=module_type,
            cooling_system_type=cooling_system_type,
            cooling_inlet_temp=cooling_inlet_temp,
        )

        mp = MultiPeriodModel(
            n_time_points=n_time_points,
            process_model_func=build_vagmd_flowsheet,
            linking_variable_func=get_vagmd_batch_variable_pairs,
            initialization_func=fix_dof_and_initialize,
            unfix_dof_func=unfix_dof,
            outlvl=logging.WARNING,
        )

        model_options = {
            "dt": dt,
            "system_capacity": system_capacity,
            "feed_flow_rate": feed_flow_rate,
            "evap_inlet_temp": evap_inlet_temp,
            "cond_inlet_temp": cond_inlet_temp,
            "feed_temp": feed_temp,
            "feed_salinity": feed_salinity,
            "recovery_ratio": recovery_ratio,
            "module_type": module_type,
            "cooling_system_type": cooling_system_type,
            "cooling_inlet_temp": cooling_inlet_temp,
        }

        mp.build_multi_period_model(
            model_data_kwargs={t: model_options for t in range(n_time_points)},
            flowsheet_options=model_options,
            initialization_options={
                "feed_flow_rate": feed_flow_rate,
                "feed_salinity": feed_salinity,
                "feed_temp": feed_temp,
            },
            unfix_dof_options={"feed_flow_rate": feed_flow_rate},
        )

        active_blks = mp.get_active_process_blocks()

        # Initialize and unfix dof for each period
        solver = get_solver()
        for blk in active_blks:
            fix_dof_and_initialize(
                m=blk,
                feed_flow_rate=feed_flow_rate,
                feed_salinity=feed_salinity,
                feed_temp=feed_temp,
            )
            result = solver.solve(blk)
            unfix_dof(m=blk, feed_flow_rate=feed_flow_rate)

        # Set-up for the first time period
        active_blks[0].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(
            feed_salinity
        )
        active_blks[0].fs.vagmd.feed_props[0].temperature.fix(feed_temp + 273.15)
        active_blks[0].fs.acc_distillate_volume.fix(0)
        active_blks[0].fs.pre_feed_temperature.fix(feed_temp + 273.15)
        active_blks[0].fs.pre_permeate_flow_rate.fix(0)
        active_blks[0].fs.acc_thermal_energy.fix(0)
        active_blks[0].fs.pre_thermal_power.fix(0)
        active_blks[0].fs.acc_cooling_energy.fix(0)
        active_blks[0].fs.pre_cooling_power.fix(0)
        active_blks[0].fs.acc_electric_energy.fix(0)
        active_blks[0].fs.pre_cooling_pump_power_elec.fix(0)
        active_blks[0].fs.pre_feed_pump_power_elec.fix(0)

        return mp

    def add_costing_module(self, flowsheet_costing_block):
        """
        This function adds costing module to the batch operation,
        by adding the unit costing package to the last time step,
        and overwritting flow values with accumulated ones
        """

        # Specify the last time step
        vagmd = self.mp.get_active_process_blocks()[-1].fs.vagmd
        vagmd.costing = UnitModelCostingBlock(
            flowsheet_costing_block=flowsheet_costing_block
        )

        # Overwrite the thermal and electric energy flow with the accumulated values
        vagmd.costing.costing_package.cost_flow(
            self.overall_thermal_power_requirement - vagmd.thermal_power_requirement,
            "heat",
        )
        vagmd.costing.costing_package.cost_flow(
            self.overall_elec_power_requirement - vagmd.elec_power_requirement,
            "electricity",
        )

        # Recalculate the number of modules required
        vagmd.eqn_num_modules.deactivate()

        blks = self.mp.get_active_process_blocks()
        n_time_points = len(blks)

        vagmd.num_modules_constraint = Constraint(
            expr=vagmd.num_modules
            == pyunits.convert(
                vagmd.system_capacity
                / sum(
                    self.mp.get_active_process_blocks()[i].fs.vagmd.permeate_flux
                    for i in range(n_time_points)
                )
                * n_time_points
                / self.mp.get_active_process_blocks()[0].fs.vagmd.module_area,
                to_units=pyunits.dimensionless,
            )
        )

    def get_model_performance(self):
        """
        This function returns the overall performance of the batch operation in a dictionary,
        and the timewise system performance in a pandas dataframe
        """
        blks = self.mp.get_active_process_blocks()
        n_time_points = len(blks)

        production_rate = value(blks[-1].fs.acc_distillate_volume) / (
            value(blks[-1].fs.dt) * (n_time_points - 1) / 3600
        )
        overall_performance = {
            "Production capacity (L/h)": production_rate,
            "Production capacity (m3/day)": production_rate / 1000 * 24,
            "Gain output ratio": value(blks[-1].fs.gain_output_ratio),
            "Specific thermal energy consumption (kWh/m3)": value(
                blks[-1].fs.specific_energy_consumption_thermal
            ),
            "Specific electric energy consumption (kWh/m3)": value(
                blks[-1].fs.specific_energy_consumption_electric
            ),
        }

        time_period = [i for i in range(n_time_points)]
        t_minutes = [value(blks[i].fs.dt) * i / 60 for i in range(n_time_points)]
        feed_salinity = [
            value(blks[i].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"])
            for i in range(n_time_points)
        ]
        permeate_flux = [
            value(blks[i].fs.vagmd.permeate_flux) for i in range(n_time_points)
        ]
        acc_distillate_volume = [
            value(blks[i].fs.acc_distillate_volume) for i in range(n_time_points)
        ]
        feed_temp = [
            value(blks[i].fs.vagmd.feed_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        evap_out_temp = [
            value(blks[i].fs.vagmd.evaporator_out_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        cond_out_temp = [
            value(blks[i].fs.vagmd.condenser_out_props[0].temperature - 273.15)
            for i in range(n_time_points)
        ]
        recovery_ratio = [
            value(blks[i].fs.acc_recovery_ratio) for i in range(n_time_points)
        ]
        STEC = [
            value(blks[i].fs.specific_energy_consumption_thermal)
            for i in range(n_time_points)
        ]
        SEEC = [
            value(blks[i].fs.specific_energy_consumption_electric)
            for i in range(n_time_points)
        ]
        GOR = [value(blks[i].fs.gain_output_ratio) for i in range(n_time_points)]

        time_series = np.array(
            [
                time_period,
                t_minutes,
                feed_salinity,
                permeate_flux,
                acc_distillate_volume,
                feed_temp,
                evap_out_temp,
                cond_out_temp,
                recovery_ratio,
                STEC,
                SEEC,
                GOR,
            ]
        )

        data_table = pd.DataFrame(
            data=time_series,
            index=[
                "Step",
                "Operation time (min)",
                "Feed salinity (g/L)",
                "Permeate Flux (kg/hr/m2)",
                "Accumulated distillate volume (L)",
                "Feed temperature (C)",
                "Evaporator outlet temperature (C)",
                "Condenser outlet temperature (C)",
                "Accumulated recovery ratio",
                "Accumulated specific thermal energy consumption (kWh/m3)",
                "Accumulated specific electric energy consumption (kWh/m3)",
                "Gain output ratio",
            ],
        )

        return overall_performance, data_table
