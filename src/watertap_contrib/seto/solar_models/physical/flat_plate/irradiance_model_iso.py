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


@declare_process_block_class("IrradianceModelIsoSky")
class IrradianceModelIsoSkyData(SolarEnergyBaseData):
    """
    Isotropic sky irradiance model for sloped surfaces
    based on equations in Solar Engineering of Thermal Processes, Duffie and Beckman, 4th ed.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        self._tech_type = "irradiance_model_iso_sky"

        # ==========PARAMETERS==========
        self.phi = Param(initialize=1, units=pyunits.degrees, doc="Latitude of surface")

        self.lon = Param(
            initialize=1, units=pyunits.degrees, doc="Longitude of surface"
        )

        self.std_meridian = Param(
            initialize=1, units=pyunits.degrees, doc="Standard meridian of surface"
        )

        self.standard_time = Param(
            initialize=1, units=pyunits.hours, doc="Standard time"
        )

        self.G_d = Param(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="Diffuse radiation on horizontal surface",
        )

        self.rho_g = Param(
            initialize=1, units=pyunits.dimensionless, doc="Ground reflectance"
        )

        self.beta = Param(
            initialize=1, units=pyunits.degrees, doc="Tilt angle of surface"
        )

        self.G_bn = Param(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="Beam radiation normal to direction of propagation (DNI)",
        )

        self.day_of_year = Param(
            initialize=1, units=pyunits.dimensionless, doc="Day of year"
        )

        # ==========VARIABLES==========

        self.solar_time = Var(initialize=1, units=pyunits.hours, doc="Solar time")

        self.omega = Var(initialize=1, units=pyunits.degrees, doc="Hour angle of sun")

        self.theta_z = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, 90),
            doc="Solar zenith angle",
        )

        self.delta = Var(initialize=1, units=pyunits.degrees, doc="Declination of sun")

        self.R_b = Var(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="Ratio of beam radiation on tilted surface to horizontal surface in northern hemisphere for surfaces facing south",
        )

        self.theta = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, 90),
            doc="Incidence angle of radiation on tilted surface",
        )

        self.G_b = Var(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="Beam radiation on horizontal surface",
        )

        self.G = Var(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            doc="Total radiation on horizontal surface",
        )

        self.G_T = Var(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            bounds=(0, None),
            doc="Total radiation on tilted surface",
        )

        self.B = Var(initialize=1, units=pyunits.degrees, doc="Solar B")

        self.eqn_of_time = Var(
            initialize=1, units=pyunits.minutes, doc="Equation of time"
        )

        self.G_trans = Var(
            initialize=1,
            units=pyunits.W / pyunits.m**2,
            bounds=(0, None),
            doc="Total radiation transmitted through glazing of tilted collector",
        )

        self.theta_ground = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, None),
            doc="Effective incidence angle of ground reflected radiation [deg], D&B Eq. 5.4.1, pg. 213",
        )

        self.theta_diffuse = Var(
            initialize=1,
            units=pyunits.degrees,
            bounds=(0, None),
            doc="Effective incidence angle of sky diffuse radiation [deg], D&B Eq. 5.4.2, pg. 213",
        )

        self.Kta_b = Var(
            initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Incidence angle modifier coefficient for beam radiation",
        )

        self.Kta_d = Var(
            initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Incidence angle modifier coefficient for sky diffuse radiation",
        )

        self.Kta_g = Var(
            initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Incidence angle modifier coefficient for ground reflected radiation",
        )

        # ==========CONSTRAINTS==========

        @self.Constraint(doc="Equation for B, D&B eqn. 1.4.2")
        def eq_B(b):
            return self.B == (self.day_of_year - 1) * 360 * pyunits.deg / 365

        @self.Constraint(doc="Solar time, D&B eqn. 1.5.2")
        def eq_solar_time(b):
            return self.solar_time == (
                self.standard_time
                + (
                    4 * pyunits.minutes / pyunits.deg * (self.std_meridian - self.lon)
                    + self.eqn_of_time
                )
                / (60 * pyunits.minutes / pyunits.hour)
            )

        @self.Constraint(doc="Equation of time, D&B eqn. 1.5.3")
        def eq_eqn_of_time(b):
            return self.eqn_of_time == (
                229.2
                * pyunits.minutes
                * (
                    0.000075
                    + 0.001868 * cos(pyunits.convert(self.B, to_units=pyunits.rad))
                    - 0.032077 * sin(pyunits.convert(self.B, to_units=pyunits.rad))
                    - 0.014615 * cos(pyunits.convert(2 * self.B, to_units=pyunits.rad))
                    - 0.04089 * sin(pyunits.convert(2 * self.B, to_units=pyunits.rad))
                )
            )

        @self.Constraint(doc="Hour angle of sun, D&B pg. 13")
        def eq_omega(b):
            return self.omega == 15 * pyunits.deg * (
                self.solar_time - 12 * pyunits.hours
            )

        @self.Constraint(doc="Solar zenith angle, D&B eqn. 1.6.5")
        def eq_theta_z(b):
            return cos(pyunits.convert(self.theta_z, to_units=pyunits.rad)) == (
                cos(pyunits.convert(self.phi, to_units=pyunits.rad))
                * cos(pyunits.convert(self.delta, to_units=pyunits.rad))
                * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                + sin(pyunits.convert(self.phi, to_units=pyunits.rad))
                * sin(pyunits.convert(self.delta, to_units=pyunits.rad))
            )

        @self.Constraint(doc="Declination of sun, D&B eqn. 1.6.1a")
        def eq_delta(b):
            return self.delta == 23.45 * pyunits.deg * sin(
                pyunits.convert(
                    360 * pyunits.deg * (284 + self.day_of_year) / 365,
                    to_units=pyunits.rad,
                )
            )

        @self.Constraint(
            doc="Ratio of beam radiation on tilted surface to horizontal surface in northern hemisphere for surfaces facing south, D&B eqn. 1.8.2"
        )
        def eq_R_b(b):
            return self.R_b == (
                (
                    cos(pyunits.convert(self.phi - self.beta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.delta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                    + sin(pyunits.convert(self.phi - self.beta, to_units=pyunits.rad))
                    * sin(pyunits.convert(self.delta, to_units=pyunits.rad))
                )
                / (
                    cos(pyunits.convert(self.phi, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.delta, to_units=pyunits.rad))
                    * cos(pyunits.convert(self.omega, to_units=pyunits.rad))
                    + sin(pyunits.convert(self.phi, to_units=pyunits.rad))
                    * sin(pyunits.convert(self.delta, to_units=pyunits.rad))
                )
            )

        @self.Constraint(
            doc="Ratio of beam radiation on tilted surface to horizontal surface, D&B eqn. 1.8.1"
        )
        def eq_R_b_2(b):
            return self.R_b == (
                cos(pyunits.convert(self.theta, to_units=pyunits.rad))
                / cos(pyunits.convert(self.theta_z, to_units=pyunits.rad))
            )

        @self.Constraint(doc="Beam radiation on a horizontal surface")
        def eq_G_b(b):
            return self.G_b == self.G_bn * cos(
                pyunits.convert(self.theta, to_units=pyunits.rad)
            )

        @self.Constraint(doc="Total radiation on a horizontal surface")
        def eq_G(b):
            return self.G == self.G_b + self.G_d

        @self.Constraint(
            doc="Total radiation on tilted surface via isotropic diffuse sky model, D&B eqn. 2.15.1"
        )
        def eq_G_T(b):
            return self.G_T == (
                self.G_b * self.R_b
                + self.G_d
                * (1 + cos(pyunits.convert(self.beta, to_units=pyunits.rad)))
                / 2
                + self.G
                * self.rho_g
                * (1 - cos(pyunits.convert(self.beta, to_units=pyunits.rad)))
                / 2
            )

        @self.Constraint(
            doc="Effective incidence angle of ground reflected radiation [deg], D&B Eq. 5.4.1, pg. 213"
        )
        def eq_theta_ground(b):
            return self.theta_ground == (
                90 * pyunits.deg
                - 0.5788 * self.beta
                + 0.002693 * pyunits.deg**-1 * self.beta**2
            )

        @self.Constraint(
            doc="Effective incidence angle of sky diffuse radiation [deg], D&B Eq. 5.4.2, pg. 213"
        )
        def eq_theta_diffuse(b):
            return self.theta_diffuse == (
                59.7 * pyunits.deg
                - 0.1388 * self.beta
                + 0.001497 * pyunits.deg**-1 * self.beta**2
            )

        def Kta(theta_deg):
            return 1 - 0.136 * (
                1 / cos(pyunits.convert(theta_deg, to_units=pyunits.rad)) - 1
            )

        @self.Constraint(
            doc="Incidence angle modifier coefficient for beam radiation [-]"
        )
        def eq_Kta_b(b):
            return self.Kta_b == Kta(self.theta)

        @self.Constraint(
            doc="Incidence angle modifier coefficient for diffuse radiation [-]"
        )
        def eq_Kta_d(b):
            return self.Kta_d == Kta(self.theta_diffuse)

        @self.Constraint(
            doc="Incidence angle modifier coefficient for ground reflected radiation [-]"
        )
        def eq_Kta_g(b):
            return self.Kta_g == Kta(self.theta_ground)

        @self.Constraint(
            doc="Total radiation transmitted through glazing of tilted collector"
        )
        def eq_G_trans(b):
            return self.G_trans == (
                self.G_b * self.R_b * self.Kta_b
                + self.G_d
                * self.Kta_d
                * (1 + cos(pyunits.convert(self.beta, to_units=pyunits.rad)))
                / 2
                + self.G
                * self.rho_g
                * self.Kta_g
                * (1 - cos(pyunits.convert(self.beta, to_units=pyunits.rad)))
                / 2
            )
