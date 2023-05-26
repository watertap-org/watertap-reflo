import logging
import pandas as pd

# Pyomo imports
from pyomo.environ import (
    Constraint,
    value,
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.solvers.get_solver import get_solver

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.util.model_diagnostics.infeasible import *

# Flowsheet function imports
from watertap_contrib.seto.multiperiod.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)
from watertap_contrib.seto.multiperiod.VAGMD_batch_design_model import get_n_time_points


__author__ = "Zhuoran Zhang"

# Set up logger
_log = idaeslog.getLogger(__name__)


def get_vagmd_batch_variable_pairs(t1, t2):
    """
    This function returns paris of variables that need to be connected across two time periods

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
        (t1.fs.vagmd.thermal_power, t2.fs.pre_thermal_power),
        (t1.fs.vagmd.cooling_power_thermal, t2.fs.pre_cooling_power),
    ]


def unfix_dof(m, feed_flow_rate):
    """
    This function unfixes a few degrees of freedom for optimization

    Args:
        m: object containing the integrated nuclear plant flowsheet

    Returns:
        None
    """

    # Initialize 2

    # m.fs.pre_feed_temperature.unfix()
    # m.fs.pre_feed_salinity.unfix()
    # m.fs.pre_evap_out_temp.unfix()
    # m.fs.pre_permeate_flow_rate.unfix()
    # m.fs.pre_acc_distillate_volume.unfix()
    # m.fs.pre_acc_thermal_energy.unfix()
    # m.fs.pre_thermal_power.unfix()
    # m.fs.pre_acc_cooling_energy.unfix()
    # m.fs.pre_cooling_power.unfix()

    # Initialize 1

    m.fs.pre_feed_temperature.unfix()
    m.fs.pre_permeate_flow_rate.unfix()
    m.fs.acc_distillate_volume.unfix()
    m.fs.acc_thermal_energy.unfix()
    m.fs.pre_thermal_power.unfix()
    m.fs.acc_cooling_energy.unfix()
    m.fs.pre_cooling_power.unfix()

    m.fs.vagmd.feed_props[0].temperature.unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
    m.fs.vagmd.feed_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].fix(
        pyunits.convert(
            feed_flow_rate * pyunits.L / pyunits.h, to_units=pyunits.m**3 / pyunits.s
        )
    )

    return


def create_multiperiod_vagmd_batch_model(
    feed_flow_rate = 600,
    evap_inlet_temp = 80,
    cond_inlet_temp = 25,
    feed_temp = 25,
    feed_salinity = 35,
    recovery_ratio = 0.5,
    initial_batch_volume = 50,
    module_type = "AS7C1.5L",
    cooling_system_type = "closed",
    cooling_inlet_temp = 25, # not required if cooling system type is "open"
):
    """
    This function creates a multi-period vagmd batch flowsheet object. This object contains
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """


    # Calculate the number of periods to reach target recovery rate by solving the system first
    n_time_points = get_n_time_points(
                    feed_flow_rate = feed_flow_rate,
                    evap_inlet_temp = evap_inlet_temp,
                    cond_inlet_temp = cond_inlet_temp,
                    feed_temp = feed_temp,
                    feed_salinity = feed_salinity,
                    recovery_ratio = recovery_ratio,
                    initial_batch_volume = initial_batch_volume,
                    module_type = module_type,
                    cooling_system_type = cooling_system_type,
                    cooling_inlet_temp = cooling_inlet_temp, # not required if cooling system type is "open"
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
        "feed_flow_rate": feed_flow_rate,
        "evap_inlet_temp": evap_inlet_temp,
        "cond_inlet_temp": cond_inlet_temp,
        "feed_temp": feed_temp,
        "feed_salinity": feed_salinity,
        "recovery_ratio": recovery_ratio,
        "module_type": module_type,
        "cooling_system_type": cooling_system_type,
        "cooling_inlet_temp": cooling_inlet_temp
    }
    
    # For module AS7C1.5L, when the final brine salinity is larger than 175.3 g/L,
    # the feed flow rate will be fixed at 1100 L/h
    final_brine_salinity = feed_salinity / (1- recovery_ratio) # g/L
    if module_type == "AS7C1.5L" and final_brine_salinity > 175.3:
        feed_flow_rate = 1100    # L/h

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
    print("dof after creating model and initialized:", degrees_of_freedom(mp))

    active_blks = mp.get_active_process_blocks()
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

    return mp

