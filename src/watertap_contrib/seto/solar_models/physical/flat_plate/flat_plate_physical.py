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

from pyomo.environ import Var, Param, Constraint, Suffix, cos, sin, units as pyunits

from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale

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
        self.phi = Param(
            initialize=1,
            units=pyunits.degrees,
            doc="Latitude of collector"
        )

        self.lon = Param(
            initialize=1,
            units=pyunits.degrees,
            doc="Longitude of collector"
        )

        self.std_meridian = Param(
            initialize=1,
            units=pyunits.degrees,
            doc="Standard meridian of collector"
        )

        self.standard_time = Param(
            initialize=1,
            units=pyunits.hours,
            doc="Standard time"
        )

        self.G_d = Param(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            doc="Diffuse radiation on horizontal surface"
        )

        self.rho_g = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Ground reflectance"
        )

        self.beta = Param(
            initialize=1,
            units=pyunits.degrees,
            doc="Tilt angle of collector"
        )

        self.G_bn = Param(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            doc="Beam radiation normal to direction of propagation (DNI)"
        )

        self.day_of_year = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Day of year"
        )

        self.solar_time = Var(
            initialize=1,
            units=pyunits.hours,
            doc="Solar time"
        )

        self.omega = Var(
            initialize=1,
            units=pyunits.degrees,
            doc="Hour angle of sun"
        )

        self.theta_z = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, 90),
            doc="Solar zenith angle"
        )

        self.delta = Var(
            initialize=1,
            units=pyunits.degrees,
            doc="Declination of sun"
        )

        self.R_b = Var(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            doc="Ratio of beam radiation on tilted surface to horizontal surface in northern hemisphere for collectors facing south"
        )

        self.theta = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, 90),
            doc="Incidence angle of radiation on tilted surface"
        )

        self.G_b = Var(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            doc="Beam radiation on horizontal surface"
        )

        self.G = Var(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            doc="Total radiation on horizontal surface"
        )

        self.G_T = Var(
            initialize=1,
            units=pyunits.W / pyunits.m ** 2,
            bounds=(0, None),
            doc="Total radiation on tilted surface"
        )

        self.B = Var(
            initialize=1,
            units=pyunits.degrees,
            doc="Solar B"
        )

        self.eqn_of_time = Var(
            initialize=1,
            units=pyunits.minutes,
            doc="Equation of time"
        )

        # Equation for B, D&B 4th ed. eqn. 1.4.2
        self.eqn_B = Constraint(
            expr=self.B
            == (self.day_of_year - 1) * 360 * pyunits.deg / 365
        )

        # Solar time, D&B 4th ed. eqn. 1.5.2
        self.eqn_solar_time = Constraint(
            expr=self.solar_time
            == (
                    self.standard_time + (4 * pyunits.minutes / pyunits.deg * (self.std_meridian - self.lon) + self.eqn_of_time)
                    / (60 * pyunits.minutes / pyunits.hour)
            )
        )

        # Equation of time, D&B 4th ed. eqn. 1.5.3
        self.eq_eqn_of_time = Constraint(
            expr=self.eqn_of_time
            == (
                    229.2 * pyunits.minutes
                    * ( 0.000075 + 0.001868 * cos(pyunits.convert(self.B, to_units=pyunits.rad))
                        - 0.032077 * sin(pyunits.convert(self.B, to_units=pyunits.rad))
                        - 0.014615 * cos(pyunits.convert(2 * self.B, to_units=pyunits.rad))
                        - 0.04089 * sin(pyunits.convert(2 * self.B, to_units=pyunits.rad)))
            )
        )

        # Hour angle of sun, D&B 4th ed. pg. 13
        self.eq_omega = Constraint(
            expr=self.omega
            == 15 * pyunits.deg * (self.solar_time - 12 * pyunits.hours)
        )

        # Solar zenith angle, D&B 4th ed. eqn. 1.6.5
        self.eq_theta_z = Constraint(
            expr=cos(pyunits.convert(self.theta_z, to_units=pyunits.rad))
            == (
                    cos(pyunits.convert(self.phi, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.delta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                    + sin(pyunits.convert(self.phi, to_units=pyunits.rad))
                    * sin(pyunits.convert(self.delta, to_units=pyunits.rad))
            )
        )

        # Declination of sun, D&B 4th ed. eqn. 1.6.1a
        self.eq_delta = Constraint(
            expr=self.delta
            == 23.45 * pyunits.deg * sin(pyunits.convert(360 * pyunits.deg * (284 + self.day_of_year) / 365, to_units=pyunits.rad))
        )

        # Ratio of beam radiation on tilted surface to horizontal surface in northern hemisphere for collectors facing south
        # D&B 4th ed. eqn. 1.8.2
        self.eq_R_b = Constraint(
            expr=self.R_b
            == (
                    (cos(pyunits.convert(self.phi - self.beta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.delta, to_units=pyunits.rad)) * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                    + sin(pyunits.convert(self.phi - self.beta, to_units=pyunits.rad)) * sin(pyunits.convert(self.delta, to_units=pyunits.rad)))
                    /
                    (cos(pyunits.convert(self.phi, to_units=pyunits.rad)) * cos(pyunits.convert(self.delta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                    + sin(pyunits.convert(self.phi, to_units=pyunits.rad)) * sin(pyunits.convert(self.delta, to_units=pyunits.rad)))
            )
        )
        
        # Ratio of beam radiation on tilted surface to horizontal surface, D&B 4th ed. eqn. 1.8.1
        self.eq_R_b_2 = Constraint(
            expr=self.R_b
            == (
                    cos(pyunits.convert(self.theta, to_units=pyunits.rad))
                    / cos(pyunits.convert(self.theta_z, to_units=pyunits.rad))
            )
        )

        # Beam radiation on a horizontal surface
        self.eq_G_b = Constraint(
            expr=self.G_b
            == self.G_bn * cos(pyunits.convert(self.theta, to_units=pyunits.rad))
        )

        # Total radiation on a horizontal surface
        self.eq_G = Constraint(
            expr=self.G
            == self.G_b + self.G_d
        )

        # Total radiation on tilted surface via isotropic diffuse sky model
        self.eq_G_T = Constraint(
            expr=self.G_T
            == (
                    self.G_b * self.R_b
                    + self.G_d * (1 + cos(pyunits.convert(self.beta, to_units=pyunits.rad))) / 2
                    + self.G * self.rho_g * (1 - cos(pyunits.convert(self.beta, to_units=pyunits.rad))) / 2
                )
        )

        # NOTE: this is an alternative formulation for G_T that does not currently work
        # @self.Constraint(doc="Total radiation on tilted surface via isotropic diffuse sky model")
        # def eq_G_T():
        #     return (
        #         self.G_T
        #         == (
        #             self.G_b * self.R_b
        #             + self.G_d * (1 + cos(pyunits.convert(self.beta, to_units=pyunits.rad))) / 2
        #             + self.G * self.rho_g * (1 - cos(pyunits.convert(self.beta, to_units=pyunits.rad))) / 2
        #         )
        #     )

    # def initialize_build(self):
    #     pass
