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
This module contains a zero-order representation of a solar energy
operation.
"""

from idaes.core import declare_process_block_class
from idaes.core.util.misc import StrEnum

from pyomo.common.config import ConfigValue, In

from watertap.core import build_pt, ZeroOrderBaseData
from watertap_contrib.seto.energy import solar_energy

__author__ = "Kurban Sitterley"


class SolarEnergyZOType(StrEnum):
    PV = "PV"

@declare_process_block_class("SolarEnergyZO")
class SolarEnergyZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a solar energy source.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()
    CONFIG.declare(
        "solar_energy_type",
        ConfigValue(
            default="PV",
            domain=In(SolarEnergyZOType),
            description="Indicates type of solar energy source",
        ),
    )

    def build(self):
        super().build()

        self._tech_type = "solar_energy"

        build_pt(self)
        solar_energy(self)
