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

from pyomo.environ import Var, Param, Suffix, cos, sin, units as pyunits

from idaes.core import declare_process_block_class

from watertap_contrib.seto.core import SolarEnergyBaseData

__author__ = "Matthew Boyd"


@declare_process_block_class("FlatPlatePhysical")
class FlatPlatePhysicalData(SolarEnergyBaseData):
    """
    Physical model for flat plate.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        self._tech_type = "flat_plate"

        # ==========PARAMETERS==========
        self.x = Param(
            initialize=1, units=pyunits.degrees, doc="example parameter"
        )

        # ==========VARIABLES==========

        self.y = Var(initialize=1, units=pyunits.hours, doc="example variable")

        # ==========CONSTRAINTS==========

        @self.Constraint(doc="example docstring")
        def eq_x(b):
            return 0
