import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import (
    Product,
    Feed,
    StateJunction,
    Separator,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core import Database
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock as MCAS,
)

# from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.analysis.case_studies.permian import *

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1000 * pyunits.kg / pyunits.m**3


def build_permian():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = REFLODatabase()

    m.fs.properties = ZO(solute_list=["tds"])

    m.fs.treatment = treat = Block()

    treat.feed = Feed(property_package=m.fs.properties)
    treat.product = Product(property_package=m.fs.properties)
    treat.disposal = Product(property_package=m.fs.properties)

    treat.disposal_mixer = Mixer(
        property_package=m.fs.properties,
        num_inlets=2,
        inlet_list=["ec_disposal", "cart_filt_disposal"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m, treat.chem_addition)
    treat.EC = FlowsheetBlock(dynamic=False)
    build_ec(m, treat.EC)
    treat.cart_filt = FlowsheetBlock(dynamic=False)
    build_cartridge_filtration(m, treat.cart_filt)

    treat.feed_to_chem_addition = Arc(
        source=treat.feed.outlet, destination=treat.chem_addition.feed.inlet
    )

    treat.chem_addition_to_ec = Arc(
        source=treat.chem_addition.product.outlet, destination=treat.EC.feed.inlet
    )

    treat.ec_to_cart_filt = Arc(
        source=treat.EC.product.outlet, destination=treat.cart_filt.feed.inlet
    )

    treat.ec_to_disposal_mix = Arc(
        source=treat.EC.disposal.outlet, destination=treat.disposal_mixer.ec_disposal
    )

    treat.cart_filt_to_disposal_mix = Arc(
        source=treat.cart_filt.disposal.outlet,
        destination=treat.disposal_mixer.cart_filt_disposal,
    )

    treat.disposal_mix_to_disposal = Arc(
        source=treat.disposal_mixer.outlet, destination=treat.disposal.inlet
    )

    treat.cart_filt_to_product = Arc(
        source=treat.cart_filt.product.outlet, destination=treat.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def add_treatment_costing(m):
    m.fs.treatment.costing = TreatmentCosting(case_study_definition=case_study_yaml)
    pass


if __name__ == "__main__":
    m = build_permian()

    add_treatment_costing(m)

    print(f"dof = {degrees_of_freedom(m)}")
