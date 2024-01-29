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

from idaes.core import declare_process_block_class

from pyomo.environ import Param, units as pyunits

from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap_contrib.reflo.costing.solar.photovoltaic import cost_pv

__author__ = "Kurban Sitterley"


@declare_process_block_class("Photovoltaic")
class PhotovoltaicData(SolarEnergyBaseData):
    """
    Zero-Order model for photovoltaic system.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "photovoltaic"

        self.oversize_factor = Param(
            initialize=1,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Oversize factor for PV system",
        )

    @property
    def default_costing_method(self):
        return cost_pv
