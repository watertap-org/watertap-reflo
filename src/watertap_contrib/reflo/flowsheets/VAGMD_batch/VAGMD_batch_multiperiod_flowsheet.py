#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
    ConcreteModel,
    Var,
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)

# IDAES imports
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
import idaes.core.util.scaling as iscale
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

# WaterTAP imports
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing import REFLOCosting

# Flowsheet function imports
from watertap_contrib.reflo.flowsheets.VAGMD_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

__author__ = "Zhuoran Zhang"

__all__ = [
    "build_VAGMD_batch_multiperiod_fs",
    "add_costing_module",
    "get_multiperiod_performance",
]


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


def build_VAGMD_batch_multiperiod_fs(
    m,
    n_time_points,
    system_capacity=1000,  # m3/day
    # Arguments below are the initialized condition for multiperiod model construction,
    dt=None,
    feed_flow_rate=600,  # 400 - 1100 L/h
    evap_inlet_temp=80,  # 60 - 80 deg C
    cond_inlet_temp=25,  # 20 - 30 deg C
    feed_temp=25,  # 20 - 30 deg C
    feed_salinity=35,  # 35 - 292 g/L
    # initial_batch_volume=50,  # > 50 L
    # recovery_ratio=0.5,  # -
    module_type="AS7C1.5L",
    cooling_system_type="closed",
    cooling_inlet_temp=25,  # 20 - feed_temp deg C, not required when cooling system type is "closed"
    **kwargs
):
    """
    This function builds a multiperiod flowsheet of the VAGMD system
    for the simulation of a batch operation

    Returns:
        object: A multiperiod module created for the target Pyomo optimization model m
    """

    if m is None:
        m = ConcreteModel()

    # Create the multiperiod model object
    m.mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_vagmd_flowsheet,
        linking_variable_func=get_vagmd_batch_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    mp = m.mp

    """
    Specify the initialization conditions of each period
    """

    # Model options listed in the VAGMD flowsheet of a single period,
    # see arguments from the imported 'build_vagmd_flowsheet' method.
    model_options = {
        "dt": dt,
        "system_capacity": system_capacity,
        "feed_flow_rate": feed_flow_rate,
        "evap_inlet_temp": evap_inlet_temp,
        "cond_inlet_temp": cond_inlet_temp,
        "feed_temp": feed_temp,
        "feed_salinity": feed_salinity,
        "module_type": module_type,
        "cooling_system_type": cooling_system_type,
        "cooling_inlet_temp": cooling_inlet_temp,
    }

    mp.build_multi_period_model(
        model_data_kwargs={t: model_options for t in range(n_time_points)},
        flowsheet_options=model_options,
        initialization_options={
            "feed_temp": feed_temp,
        },
        unfix_dof_options={"feed_flow_rate": feed_flow_rate},
    )

    # Create accumlative energy terms
    mp.system_capacity = Var(
        initialize=system_capacity,
        bounds=(0, None),
        units=pyunits.m**3 / pyunits.day,
        doc="System capacity (m3/day)",
    )

    mp.overall_thermal_power_requirement = Var(
        initialize=2e5,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Thermal power requirement (kW-th)",
    )

    mp.overall_elec_power_requirement = Var(
        initialize=300,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Electric power requirement (kW-e)",
    )

    @mp.Constraint(
        doc="Calculate the overall thermal power requirement through all periods"
    )
    def eq_thermal_power_requirement(b):
        return b.overall_thermal_power_requirement == (
            mp.get_active_process_blocks()[-1].fs.specific_energy_consumption_thermal
            * pyunits.convert(b.system_capacity, to_units=pyunits.m**3 / pyunits.h)
        )

    @mp.Constraint(
        doc="Calculate the overall electric power requirement through all periods"
    )
    def eq_elec_power_requirement(b):
        return b.overall_elec_power_requirement == (
            mp.get_active_process_blocks()[-1].fs.specific_energy_consumption_electric
            * pyunits.convert(b.system_capacity, to_units=pyunits.m**3 / pyunits.h)
        )

    iscale.calculate_scaling_factors(mp)

    if iscale.get_scaling_factor(mp.overall_thermal_power_requirement) is None:
        iscale.set_scaling_factor(mp.overall_thermal_power_requirement, 1e-5)

    if iscale.get_scaling_factor(mp.overall_elec_power_requirement) is None:
        iscale.set_scaling_factor(mp.overall_elec_power_requirement, 1e-3)

    active_blks = m.mp.get_active_process_blocks()

    # Initialize and unfix dof for each period
    solver = get_solver()
    for blk in active_blks:
        fix_dof_and_initialize(
            m=blk,
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

    return m.mp


def add_costing_module(mp, flowsheet_costing_block):
    """
    This function adds costing module to the target multiperiod module mp,
    by adding the unit costing package to the last time step,
    and overwritting flow values with accumulated ones

    Args:
        mp: A pyomo module that describes the vagmd multiperiod model
        flowsheet_costing_block: Costing block which has been created for the whole flowsheet

    Returns:
        object: A costing module associated to the multiperiod module
    """

    # Specify the last time step
    vagmd = mp.get_active_process_blocks()[-1].fs.vagmd
    vagmd.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )

    # Overwrite the thermal and electric energy flow with the accumulated values
    vagmd.costing.costing_package.cost_flow(
        mp.overall_thermal_power_requirement - vagmd.thermal_power_requirement,
        "heat",
    )
    vagmd.costing.costing_package.cost_flow(
        mp.overall_elec_power_requirement - vagmd.elec_power_requirement,
        "electricity",
    )

    # Recalculate the number of modules required
    vagmd.eqn_num_modules.deactivate()

    blks = mp.get_active_process_blocks()
    n_time_points = len(blks)

    vagmd.num_modules_constraint = Constraint(
        expr=vagmd.num_modules
        == pyunits.convert(
            vagmd.system_capacity
            / sum(
                mp.get_active_process_blocks()[i].fs.vagmd.permeate_flux
                for i in range(n_time_points)
            )
            * n_time_points
            / mp.get_active_process_blocks()[0].fs.vagmd.module_area,
            to_units=pyunits.dimensionless,
        )
    )


def get_multiperiod_performance(mp):
    """
    This function returns the overall performance of the batch operation in a dictionary,
    and the timewise system performance in a pandas dataframe
    """
    blks = mp.get_active_process_blocks()
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


if __name__ == "__main__":
    # Create the multiperiod flowsheet as a module
    m = ConcreteModel()
    mp = build_VAGMD_batch_multiperiod_fs(
        m=m,
        n_time_points=70,
        system_capacity=1000,
    )

    """
    Fix the input variables
    """
    # Fix system capacity
    mp.system_capacity.fix()

    active_blks = mp.get_active_process_blocks()

    # Fix feed properties
    feed_conc = 35  # g/L
    feed_temp = 30  # C

    active_blks[0].fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(35)
    active_blks[0].fs.vagmd.feed_props[0].temperature.fix(30 + 273.15)

    # Fix initial conditions
    active_blks[0].fs.acc_distillate_volume.fix(0)
    # pre_feed_temperature should be same as the feed temperature
    active_blks[0].fs.pre_feed_temperature.fix(30 + 273.15)
    active_blks[0].fs.pre_permeate_flow_rate.fix(0)
    active_blks[0].fs.acc_thermal_energy.fix(0)
    active_blks[0].fs.pre_thermal_power.fix(0)
    active_blks[0].fs.acc_cooling_energy.fix(0)
    active_blks[0].fs.pre_cooling_power.fix(0)
    active_blks[0].fs.acc_electric_energy.fix(0)
    active_blks[0].fs.pre_cooling_pump_power_elec.fix(0)
    active_blks[0].fs.pre_feed_pump_power_elec.fix(0)

    """
    Test solution
    """
    assert degrees_of_freedom(m) == 0

    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)

    overall_performance, data_table = get_multiperiod_performance(mp)

    print(
        "Accumulated recovery ratio: ", data_table.loc["Accumulated recovery ratio"][69]
    )

    """
    Test costing soluiton
    """
    # The costing model is built upon the last time step
    m.costing = REFLOCosting()
    m.costing.base_currency = pyunits.USD_2020

    add_costing_module(mp, m.costing)

    # Fix some global costing params for better comparison to Pyomo model
    m.costing.total_investment_factor.fix(1)
    m.costing.maintenance_labor_chemical_factor.fix(0)
    m.costing.capital_recovery_factor.fix(0.08764)
    m.costing.wacc.unfix()

    m.costing.cost_process()
    m.costing.add_annual_water_production(active_blks[-1].fs.vagmd.system_capacity)
    m.costing.add_LCOW(active_blks[-1].fs.vagmd.system_capacity)

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)

    print("Overall LCOW ($/m3): ", value(m.costing.LCOW))


if __name__ == "__main__":
    m = ConcreteModel()
    build_VAGMD_batch_multiperiod_fs(m)
