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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Suffix,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog

# from watertap_contrib.reflo.costing.units.med_tvc_surrogate import (
#     cost_med_tvc_surrogate,
# )

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


@declare_process_block_class("ForwardOsmosis0D")
class ForwardOsmosis0DData(UnitModelBlockData):
    """
    Forward Osmosis - 0D model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
        ),
    )
    CONFIG.declare(
        "property_package_water",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for feed, distillate, brine and cooling water properties",
            doc="""Property parameter object used to define water property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_drawsolution",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for heating and motive steam properties",
            doc="""Property parameter object used to define steasm property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = (
            self.config.property_package_liquid.get_metadata().get_derived_units
        )

        """
        Specify system configurations
        """       
        self.recovery_ratio = Var(
            initialize=0.3,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        self.regeneration_temp = Var(
            initialize = 90 + 273.15,
            bounds=(0, None),
            units=pyunits.K,
            doc='Regeneration temperature of draw solution in the separator'
        )

        self.separator_temp_loss = Var(
            initialize = 1,
            bounds=(0, None),
            units=pyunits.K,
            doc='Temperature loss of draw solution in the separator'
        )

        """
        Add block for feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_water
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict,
        )


        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

        self.brine_props = self.config.property_package_liquid.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        """
        Add block for strong draw solution entering the membrane module
        """
        tmp_dict["defined_state"] = False

        self.strong_draw_props = self.config.property_package_drawsolution.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        """
        Add block for weak draw solution leaving the membrane module
        """
        tmp_dict["defined_state"] = False

        self.weak_draw_props = self.config.property_package_drawsolution.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        """
        Add block for regenerated draw solution from the separator
        """
        tmp_dict["defined_state"] = False

        self.reg_draw_props = self.config.property_package_drawsolution.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        """
        Add block for the product water from the separator
        """
        tmp_dict["defined_state"] = False

        self.product_props = self.config.property_package_drawsolution.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )