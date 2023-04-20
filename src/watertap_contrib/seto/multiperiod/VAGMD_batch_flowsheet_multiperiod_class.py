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
from watertap_contrib.seto.unit_models.surrogate import VAGMDsurrogate

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
    ]

def create_multiperiod_vagmd_batch_model(
        n_time_points=71,
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
        linking_variable_func=get_vagmd_batch_variable_pairs,
        # initialization_func=fix_dof_and_initialize,
        # unfix_dof_func=unfix_dof,
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


    return mp_nuclear