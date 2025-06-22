# Import relevant libraries
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap_contrib.reflo.core import REFLODatabase
import idaes.logger as idaeslog
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import *
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.analysis.case_studies.KBHDP import *
from watertap.core.zero_order_properties import WaterParameterBlock
# Import relevant libraries
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap_contrib.reflo.core import REFLODatabase
import idaes.logger as idaeslog
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import *
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.reflo.analysis.case_studies.KBHDP import *
from watertap.core.zero_order_properties import WaterParameterBlock

"""
This module builds a sample reverse osmosis (RO) system using the REFLO framework.
It includes the following components:
- A flowsheet block to hold the system
- A NaCl parameter block for the RO properties
- A treatment block that includes feed, product, and waste streams
"""

# Build the model
def build_system(RE=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Get unit model input data from a database
    m.db = REFLODatabase()

    # Add uni
    m.fs.RO_properties = NaClParameterBlock()

    return m


# Build the treatment system
def build_treatment(m):
    # Create a treatment block
    m.fs.treatment = Block()
    # Add feed, product and waste streams
    m.fs.treatment.feed = Feed(property_package=m.fs.RO_properties)
    m.fs.treatment.product = Product(property_package=m.fs.RO_properties)

    # Add unit models for treatment and disposal
    m.fs.treatment.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.treatment.RO = FlowsheetBlock(dynamic=False)
    m.fs.treatment.DWI = FlowsheetBlock(dynamic=False)

    # Build unit models
    build_ro(m, m.fs.treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(m, m.fs.treatment.DWI, prop_package=m.fs.RO_properties)


# Add connecions between unit models
def add_connections(m):
    # Connect feed to pump
    m.fs.treatment.feed_to_pump = Arc(
        source=m.fs.treatment.feed.outlet,
        destination=m.fs.treatment.pump.inlet,
    )

    # Connect pump to RO
    m.fs.treatment.pump_to_ro = Arc(
        source=m.fs.treatment.pump.outlet,
        destination=m.fs.treatment.RO.inlet,
    )

    # Connect RO to product
    m.fs.treatment.ro_to_product = Arc(
        source=m.fs.treatment.RO.outlet,
        destination=m.fs.treatment.product.inlet,
    )

    # Connect RO to DWI
    m.fs.treatment.ro_to_dwi = Arc(
        source=m.fs.treatment.RO.waste,
        destination=m.fs.treatment.DWI.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

# Set operating conditions for the RO unit model
def set_operating_conditions(m):
    # Set feed flow rate
    m.fs.treatment.feed.flow_vol[0].fix(1000)  # 1000 m^3/day

    # Set feed concentration
    m.fs.treatment.feed.flow_mass_comp[0, "H2O"].fix(1000)  # 1000 kg/day
    m.fs.treatment.feed.conc_mass_comp[0, "NaCl"].fix(10)  

    # Set pump pressure
    m.fs.treatment.pump.pressure[0].fix(300000)  # 300 kPa

def initialize_system(m):
    # Initialize the feed stream
    m.fs.treatment.feed.initialize()

    # Propagate state from feed to pump
    _prop_state(m.fs.treatment.feed_to_pump)

    # Initialize the pump
    m.fs.treatment.pump.initialize()

    # Initialize the RO unit model
    m.fs.treatment.RO.initialize()

    # Initialize the product stream
    m.fs.treatment.product.initialize()

    # Initialize the DWI unit model
    m.fs.treatment.DWI.initialize()



def build_sample_ro():
    # Build the system
    m = build_system()

    # Build the treatment system
    build_treatment(m)

    # Add connections between unit models
    add_connections(m)

    # Set operating conditions
    set_operating_conditions(m)



    return m