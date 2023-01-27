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

from watertap.core import build_pt, ZeroOrderBaseData
from watertap_contrib.seto.energy import solar_energy

__author__ = "Kurban Sitterley"


@declare_process_block_class("PhotovoltaicZO")
class SolarEnergyZOData(ZeroOrderBaseData):
    """
    Zero-Order model for photovoltaic system.
    """

    def build(self):
        super().build()

        self._tech_type = "photovoltaic"

        build_pt(self)
        solar_energy(self)
