###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains a base class for all solar energy unit models.
"""

from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.misc import StrEnum
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var, units as pyunits, NonNegativeReals

__author__ = "Kurban Sitterley"


class SolarEnergyType(StrEnum):
    PV = "PV"


@declare_process_block_class("SolarEnergyBase")
class SolarEnergyBaseData(UnitModelBlockData):
    """
    Base model for a solar energy source.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Solar energy models are steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Solar energy models do not include holdup""",
        ),
    )
    CONFIG.declare(
        "solar_energy_type",
        ConfigValue(
            default="PV",
            domain=In(SolarEnergyType),
            description="Indicates type of solar energy source",
        ),
    )

    def build(self):
        super().build()

        self._tech_type = None
        self._scaling = None

        self.electricity = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Electricity production of solar process",
        )

        self.heat = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Heat production of solar process",
        )

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        Placeholder, should be overloaded
        by derived classes as necessary.
        """

        if solver is None:
            solver = get_solver()

        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise InitializationError(
                    f"\nWhile initializing {self.name}, the degrees of freedom "
                    "are {dof}, when zero is required. \n"
                )

    def calculate_scaling_factors(self):
        """
        Placeholder scaling routine, should be overloaded by derived classes
        """
        super().calculate_scaling_factors()

        if callable(self._scaling):
            self._scaling(self)
