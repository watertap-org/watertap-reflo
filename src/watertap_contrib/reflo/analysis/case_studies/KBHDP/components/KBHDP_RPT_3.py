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
    TransformationFactory,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
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
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import *

# Flowsheet function imports
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_design_model import (
    get_n_time_points,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel


# TODO: Add translator block after