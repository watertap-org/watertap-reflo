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

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# Import base model from WaterTAP REFLO
from watertap_contrib.reflo.unit_models.surrogate.vagmd_surrogate_base import (
    VAGMDBaseData,
)


_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


@declare_process_block_class("VAGMDSurrogate")
class VAGMDData(VAGMDBaseData):
    """
    Vacuum air-gapped membrane distillation surrogate model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be True",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package_seawater",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_water",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
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
    CONFIG.declare(
        "high_brine_salinity",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="""A ConfigBlock with arguments to be passed to indicate if the brine
            has a salinity higher than 175.3 g/L""",
        ),
    )
    CONFIG.declare(
        "module_type",
        ConfigValue(
            default="AS7C1.5L",
            domain=In(["AS7C1.5L", "AS26C7.2L"]),
            doc="""Selection of module type (7 for AS7C1.5L and 26 for AS26C7.2L)""",
        ),
    )
    CONFIG.declare(
        "cooling_system_type",
        ConfigValue(
            default="closed",
            domain=In(["open", "closed"]),
            doc="""Selection of cooling system type (open or closed)""",
        ),
    )

    def build(self):
        super().build()

        self.system_capacity = Var(
            initialize=2000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.day,
            doc="System capacity",
        )

        self.num_modules = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of MD module required",
        )

        self.recovery_ratio = Var(
            initialize=0.1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        @self.Constraint(doc="Calculate recovery ratio")
        def eqn_recovery_ratio(b):
            return (
                b.recovery_ratio
                == 1
                - b.evaporator_out_props[0].flow_vol_phase["Liq"]
                / b.feed_props[0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="Calculate number of module required")
        def eqn_num_modules(b):
            return b.num_modules == pyunits.convert(
                b.system_capacity / b.permeate_flux / b.module_area,
                to_units=pyunits.dimensionless,
            )

        @self.Expression(doc="Calculate specific thermal energy consumption (kWh/m3)")
        def specific_energy_consumption_thermal(b):
            return b.thermal_power / pyunits.convert(
                b.permeate_flux * b.module_area, to_units=pyunits.m**3 / pyunits.h
            )

        @self.Expression(doc="Calculate specific electric energy consumption (kWh/m3)")
        def specific_energy_consumption_electric(b):
            return (
                b.feed_pump_power_elec + b.cooling_pump_power_elec
            ) / pyunits.convert(
                b.permeate_flux * b.module_area, to_units=pyunits.m**3 / pyunits.h
            )

        @self.Expression(doc="Calculate thermal power requirement (kW)")
        def thermal_power_requirement(b):
            return b.specific_energy_consumption_thermal * pyunits.convert(
                b.system_capacity, to_units=pyunits.m**3 / pyunits.h
            )

        @self.Expression(doc="Calculate thermal power requirement (kW)")
        def elec_power_requirement(b):
            return b.specific_energy_consumption_electric * pyunits.convert(
                b.system_capacity, to_units=pyunits.m**3 / pyunits.h
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.system_capacity) is None:
            iscale.set_scaling_factor(self.system_capacity, 1e-3)

        if iscale.get_scaling_factor(self.num_modules) is None:
            iscale.set_scaling_factor(self.num_modules, 1e-3)

        if iscale.get_scaling_factor(self.recovery_ratio) is None:
            iscale.set_scaling_factor(self.recovery_ratio, 1e1)

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        super().initialize_build()
