
# General python imports
import numpy as np
import pandas as pd
import logging
from collections import deque

# Pyomo imports
from pyomo.environ import Set, Expression, value

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

# Flowsheet function imports

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
        (t1.fs.vagmd.feed_props[0].flow_vol_phase["Liq"], t2.fs.vagmd.pre_feed_props[0].flow_vol_phase["Liq"]),
        (t1.fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"], t2.fs.vagmd.pre_feed_props[0].conc_mass_phase_comp["Liq", "TDS"]),
        (t1.fs.vagmd.acc_distillate_volume[0], t2.fs.vagmd.pre_acc_distillate_volume[0]),
        (t1.fs.vagmd.permeate_props[0].flow_vol_phase["Liq"], t2.fs.vagmd.pre_permeate_props[0].flow_vol_phase["Liq"]),
        (t1.fs.vagmd.permeate_props[0].conc_mass_phase_comp["Liq", "TDS"], t2.fs.vagmd.pre_permeate_props[0].conc_mass_phase_comp["Liq", "TDS"]),
        (t1.fs.vagmd.distillate_props[0].flow_vol_phase["Liq"], t2.fs.vagmd.pre_distillate_props[0].flow_vol_phase["Liq"]),
        (t1.fs.vagmd.distillate_props[0].conc_mass_phase_comp["Liq", "TDS"], t2.fs.vagmd.pre_distillate_props[0].conc_mass_phase_comp["Liq", "TDS"]),
    ]

def create_multiperiod_vagmd_batch_model(
        n_time_points=70,
    ):
    """
    This function creates a multi-period vagmd batch flowsheet object. This object contains 
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """
    mp_nuclear = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_ne_flowsheet,
        linking_variable_func=get_nuclear_link_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    flowsheet_options={
        "np_capacity": 500, 
        "include_turbine": False,
        "pem_capacity": 100,
        "tank_capacity": 5000,
    }

    # create the multiperiod object
    mp_nuclear.build_multi_period_model(
        model_data_kwargs={t: flowsheet_options for t in range(n_time_points)},
        flowsheet_options=flowsheet_options,
        initialization_options={
            "split_frac_grid": 0.95,
            "tank_holdup_previous": 0,
            "flow_mol_to_pipeline": 1,
            "flow_mol_to_turbine": 0,
        },
    )

    mw_h2 = 2.016 * 1e-3  # Molecular mass of hydrogen in kg/mol

    # Define conversion factors
    # TODO: Use pyunits convert instead of defining conversion factors
    MW_to_kW = 1000
    kW_to_MW = 1e-3
    hours_to_s = 3600

    # TODO: Need to find a better way to import costing data.
    # Using a VOM cost of $2.3 per MWh for the nuclear power plant
    # Using a VOM cost of $1.3 per MWh for the PEM electrolyzer
    # Fixed O&M cost for the nuclear power plant is $120,000 per MW-year
    # Fixed O&M cost for the PEM electrolyzer is $47,900 per MW-year
    # Normalized FOM = (120,000 / 8760) = $13.7 per MWh
    # Normalized FOM = (47,900 / 8760) = $5.47 per MWh
    npp_fom = 13.7
    npp_vom = 2.3
    pem_fom = 5.47
    pem_vom = 1.3
    tank_vom = 0.01

    active_process_blks = mp_nuclear.get_active_process_blocks()
    for blk in active_process_blks:
        # Hydrogen demand constraint
        if demand_type == "variable":
            blk.fs.h2_tank.outlet_to_pipeline.flow_mol[0].setub(h2_demand / mw_h2)

        elif demand_type == "fixed":
            blk.fs.h2_tank.outlet_to_pipeline.flow_mol[0].fix(h2_demand / mw_h2)

        # Treating the revenue generated from hydrogen as negative cost
        # To avoid degeneracy, we are adding an operating cost for hydrogen storage
        # FIXME: How do we include the FOM contribution?
        blk.fs.operating_cost = Expression(
            expr=blk.fs.np_power_split.electricity[0] * kW_to_MW * npp_vom +
                blk.fs.pem.electricity[0] * kW_to_MW * pem_vom +
                blk.fs.h2_tank.tank_holdup[0] * mw_h2 * tank_vom -
                blk.fs.h2_tank.outlet_to_pipeline.flow_mol[0] * mw_h2 * hours_to_s * h2_price)

    return mp_nuclear