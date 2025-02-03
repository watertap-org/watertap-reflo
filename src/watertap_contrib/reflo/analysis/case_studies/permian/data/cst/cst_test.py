#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
import os
import numpy as np
import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    Expression,
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
)

from watertap_contrib.reflo.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap_contrib.reflo.costing import EnergyCosting

from watertap.core.solvers import get_solver

# Get default solver for testing
solver = get_solver()

dataset_filename = os.path.join(os.path.dirname(__file__), "trough_permian_data_heat_load_1_50_hours_storage_0_24.pkl")
print(dataset_filename)
test_surrogate_filename = os.path.join(
    os.path.dirname(__file__), "trough_permian_data_heat_load_1_50_hours_storage_0_24.json"
)
print(test_surrogate_filename)
input_bounds = dict(heat_load=[1, 50], hours_storage=[0, 24])
input_units = dict(heat_load="MW", hours_storage="hour")
input_variables = {
    "labels": ["heat_load", "hours_storage"],
    "bounds": input_bounds,
    "units": input_units,
}
output_bounds = dict(heat_annual_scaled=[1, 50], electricity_annual_scaled=[0, 24])
output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
output_variables = {
    "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
    "units": output_units,
}

trough_dict = dict(
    dataset_filename=dataset_filename,
    input_variables=input_variables,
    output_variables=output_variables,
    scale_training_data=True,
)


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.trough = TroughSurrogate(**trough_dict)

m.fs.trough.heat_load.fix(25)
m.fs.trough.hours_storage.fix(12)

calculate_scaling_factors(m)
m.fs.trough.initialize()

results = solver.solve(m)
assert_optimal_termination(results)

# m.fs.costing = EnergyCosting()
# # set heat and electricity costs to be non-zero
# m.fs.costing.heat_cost.set_value(0.01)
# m.fs.costing.electricity_cost.fix(0.07)
# m.fs.costing.base_currency = pyunits.USD_2021
# m.fs.trough.costing = UnitModelCostingBlock(
#     flowsheet_costing_block=m.fs.costing
# )
# m.fs.costing.maintenance_labor_chemical_factor.fix(0)
# m.fs.costing.total_investment_factor.fix(1)
# m.fs.costing.cost_process()
# m.fs.costing.initialize()

# results = solver.solve(m)
# assert_optimal_termination(results)