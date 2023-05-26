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

from pyomo.environ import Var, Param, Suffix, cos, sin, log, exp, units as pyunits

from idaes.core import declare_process_block_class

from watertap_contrib.seto.core import SolarEnergyBaseData

__author__ = "Matthew Boyd"


@declare_process_block_class("FlatPlatePhysical")
class FlatPlatePhysicalData(SolarEnergyBaseData):
    """
    Physical model for flat plate
    based on equations in Solar Engineering of Thermal Processes, Duffie and Beckman, 4th ed.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        self._tech_type = "flat_plate"

        # ==========PARAMETERS==========
        self.area_coll = Param(
            initialize=1, units=pyunits.m**2, doc="area of a single collector"
        )

        self.FRta = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="collector heat removal factor * effective transmittance-absorption product",
        )

        self.FRUL = Param(
            initialize=1,
            units=pyunits.W / (pyunits.m**2 * pyunits.C),
            doc="collector heat removal factor * collector heat loss coefficient",
        )

        self.iam = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="incidence angle modifier coefficient",
        )

        self.mdot_test = Param(
            initialize=1,
            units=pyunits.kg / pyunits.s,
            doc="mass flow rate of fluid during characterization test",
        )

        self.cp_test = Param(
            initialize=1,
            units=pyunits.J / (pyunits.kg * pyunits.C),
            doc="specific heat capacity of fluid during characterization test",
        )

        self.cp_use = Param(
            initialize=1,
            units=pyunits.J / (pyunits.kg * pyunits.C),
            doc="specific heat capacity of fluid during use",
        )

        self.ncoll = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="number of collectors in array",
        )

        self.pump_watts = Param(initialize=1, units=pyunits.W, doc="pump power")

        self.pump_eff = Param(
            initialize=1, units=pyunits.dimensionless, doc="pump efficiency"
        )

        self.T_amb = Param(initialize=1, units=pyunits.C, doc="ambient temperature")

        self.T_in = Param(initialize=1, units=pyunits.C, doc="inlet temperature")

        self.G_trans = Param(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="irradiance transmitted through glazing",
        )

        # ==========VARIABLES==========

        self.FprimeUL = Var(
            initialize=1,
            units=pyunits.W / (pyunits.m**2 * pyunits.C),
            doc="corrected collector heat loss coefficient, D&B Eq. 6.20.4",
        )

        self.r = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="ratio of FRta_use to FRta_test, D&B Eq. 6.20.3",
        )

        self.Q_useful = Var(initialize=1, units=pyunits.W, doc="useful net heat gain")

        self.P_pump = Var(initialize=1, units=pyunits.W, doc="pump power")

        # ==========CONSTRAINTS==========

        @self.Constraint(
            doc="corrected collector heat loss coefficient, D&B Eq. 6.20.4"
        )
        def eq_FprimeUL(b):
            return self.FprimeUL == (
                -self.mdot_test
                * self.cp_test
                / self.area_coll
                * log(1 - self.FRUL * self.area_coll / (self.mdot_test * self.cp_test))
            )

        @self.Constraint(doc="ratio of FRta_use to FRta_test, D&B Eq. 6.20.3")
        def eq_r(b):
            return self.r == (
                (
                    self.mdot_test
                    * self.ncoll
                    * self.cp_use
                    / (self.area_coll * self.ncoll)
                    * (
                        1
                        - exp(
                            -self.area_coll
                            * self.ncoll
                            * self.FprimeUL
                            / (self.mdot_test * self.ncoll * self.cp_use)
                        )
                    )
                )
                / self.FRUL
            )

        @self.Constraint(
            doc="useful net heat gain, not accounting for pipe heat losses or a heat exchanger, D&B Eq. 6.8.1"
        )
        def eq_Q_useful(b):
            return self.Q_useful == (
                self.area_coll
                * self.ncoll
                * self.r
                * (self.FRta * self.G_trans - self.FRUL * (self.T_in - self.T_amb))
            )

        @self.Constraint(doc="pump power")
        def eq_P_pump(b):
            return self.P_pump == (self.pump_watts / self.pump_eff)
