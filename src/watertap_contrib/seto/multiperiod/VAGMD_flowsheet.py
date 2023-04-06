# Pyomo imports
from pyomo.environ import (
    Constraint,
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc

# IDAES imports
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    Translator,
    Mixer,
    MomentumMixingType,
)
from idaes.models.properties.modular_properties.base.generic_property \
    import GenericParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
import idaes.logger as idaeslog

# WaterTAP imports
from watertap.property_models.water_prop_pack import WaterParameterBlock

# SETO imports
from watertap_contrib.seto.unit_models.surrogate import VAGMDsurrogate

__author__ = "Zhuoran Zhang"

def build_vagmd_batch_flowsheet(

):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Load VAGMD unit model
    m.fs.vagmd = VAGMDsurrogate()