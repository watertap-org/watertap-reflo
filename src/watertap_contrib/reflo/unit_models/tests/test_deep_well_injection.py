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

import pytest
from pyomo.environ import (
    ConcreteModel,
    Set,
    Var,
    Param,
    Expression,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port

from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    set_scaling_factor,
)
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.property_models import AirWaterEq
from watertap_contrib.reflo.costing import REFLOCosting
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection


# Get default solver for testing
solver = get_solver()

