#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pandas as pd
import numpy as np

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    value,
    exp,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.math import smooth_min, smooth_max
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.property_models.air_water_equilibrium_properties import (
    SaturationVaporPressureCalculation,
)
from watertap_contrib.reflo.costing.units.evaporation_pond import cost_evaporation_pond

__author__ = "Kurban Sitterley"


"""
References

T. T. Shi, D. X. Guan, J. B. Wu, A. Z. Wang, C. J. Jin and S. J. Han (2008)
"Comparison of methods for estimating evapotranspiration rate of dry forest canopy: 
    Eddy covariance, Bowen ratio energy balance, and Penman-Monteith equation"
Journal of Geophysical Research: Atmospheres 2008 Vol. 113 Issue D19. doi: 10.1029/2008jd010174

S. Ali, N. C. Ghosh and R. Singh (2007).
"Evaluating best evaporation estimate model for water surface evaporation in semi-arid region, India"
Hydrological Processes 2007 Vol. 22 Issue 8 Pages 1093-1106. doi: 10.1002/hyp.6664

I. M. Vardavas and A. Fountoulakis (1996).
"Estimation of lake evaporation from standard meteorological measurements: 
    application to four Australian lakes in different climatic regions"
Ecological Modelling 1996 Vol. 84 Issue 1-3 Pages 139-150. doi: 10.1016/0304-3800(94)00126-x
    
Food and Agriculture Organization of the United Nations (FAO) (1998)
"FAO Irrigation and Drainage Paper No. 56 - Crop Evapotranspiration"
ISBN 92-5-104219-5
https://www.fao.org/4/x0490e/x0490e00.htm

S. W. Hostetler, P. J. Bartlein (1990)
"Simulation of lake evaporation with application to modeling lake level 
    variations of Harney-Malheur Lake, Oregon"
Water Resources Research 1990 Vol. 26 Issue 10 Pages 2603-2612. doi: 10.1029/WR026i010p02603

S. Ali, N. C. Ghosh, P. K. Mishra and R. K. Singh (2015)
"A holistic water depth simulation model for small ponds"
Journal of Hydrology 2015 Vol. 529 Pages 1464-1477. doi: 10.1016/j.jhydrol.2015.08.035

U.S. Dept. of Interior & Michael C. Mickley (2006)
"Membrane Concentrate Disposal: Practices and Regulation"
Desalination and Water Purification Research and Development Program Report No. 123 (Second Edition)
Chapter 10: Evaporation Pond Disposal
"""

area_correction_factor_param_dict = {
    12: (2.6429, -0.202),
    8: (2.0512, -0.152),
    4: (1.5357, -0.092),
}

default_weather_data_column_dict = {
    "pressure": "Pressure",
    "temperature": "Temperature",
    "shortwave_radiation": "GHI",
    "relative_humidity": "Relative Humidity",
}


@declare_process_block_class("EvaporationPond")
class EvaporationPondData(InitializationMixin, UnitModelBlockData):
    """
    Zero order evaporation pond model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
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
        "property_package",
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
        "weather_data_path",
        ConfigValue(
            default=useDefault,
            domain=str,
            description="Path to weather file with a years worth of hourly data.",
            doc="""""",
        ),
    )

    CONFIG.declare(
        "weather_data_column_dict",
        ConfigValue(
            default=default_weather_data_column_dict,
            domain=dict,
            description="dict of key, value pairs used to access the required weather variables in the weather file",
            doc="""Dictionary of key, value pairs used to access the required weather variables,
    **pressure** - Atmospheric pressure, kPa
    **temperature** - Air temperature, C,
    **shortwave_radiation** - Shortwave radiation, W/m2,
    **relative_humidity** - Relative humidity (%)
    """,
        ),
    )

    CONFIG.declare(
        "dike_height",
        ConfigValue(
            default=8,
            domain=PositiveInt,
            description="Height of the dike used for evaporation pond.",
            doc="""Specifies the height of the dike for the evaporation pond. Units are ft. 
            Must be 4, 8, or 12. Determines coefficient for calculating various costing and design parameters.""",
        ),
    )

    CONFIG.declare(
        "add_enhancement",
        ConfigValue(
            default=False,
            domain=bool,
            description="Flag to indicate if evaporation enhancement is added.",
            doc="""Specifies if the evaporation enhancement is added. 
            If True, user should modify evaporation_rate_enhancement_adjustment_factor to be >1.""",
        ),
    )

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )

    def build(self):
        super().build()

        if self.config.dike_height not in [4, 8, 12]:
            raise ConfigurationError(
                f"Dike height must be either 4, 8, or 12 but {self.config.dike_height} was provided."
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        if not "TDS" in self.config.property_package.component_list:
            raise ConfigurationError(
                "TDS must be present as a component in the influent stream."
            )

        weather_data = pd.read_csv(self.config.weather_data_path, skiprows=2)
        weather_data["day_of_year"] = np.arange(len(weather_data)) // 24
        temp_col = self.config.weather_data_column_dict["temperature"]
        min_temp = 0.1  # degC

        weather_daily_min = weather_data.groupby("day_of_year").min()
        weather_daily_min[temp_col] = np.where(
            weather_daily_min[temp_col] < min_temp,
            min_temp,
            weather_daily_min[temp_col],
        )
        weather_daily_max = weather_data.groupby("day_of_year").max()
        weather_daily_max[temp_col] = np.where(
            weather_daily_max[temp_col] < min_temp,
            min_temp,
            weather_daily_max[temp_col],
        )
        weather_daily_mean = weather_data.groupby("day_of_year").mean()
        weather_daily_mean[temp_col] = np.where(
            weather_daily_mean[temp_col] < min_temp,
            min_temp,
            weather_daily_mean[temp_col],
        )

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        self.days_of_year = Set(initialize=range(0, 365))

        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )
        prop_in = self.properties_in[0]
        self.add_inlet_port(name="inlet", block=self.properties_in)

        tmp_dict["defined_state"] = True
        self.weather = self.config.property_package.state_block_class(
            self.days_of_year,
            doc="Daily atmospheric properties",
            **tmp_dict,
        )
        # avoid warnings for days with air temp < 0 degC
        self.weather[:].temperature["Vap"].setlb(200)
        self.weather[:].temperature["Liq"].setlb(200)
        self.weather[:].flow_mass_phase_comp["Liq", :].fix(0)

        for day, row in weather_daily_mean.iterrows():

            pres_Pa = pyunits.convert(
                row[self.config.weather_data_column_dict["pressure"]] * pyunits.mbar,
                to_units=pyunits.Pa,
            )
            self.weather[day].pressure.fix(pres_Pa)

            temp = row[self.config.weather_data_column_dict["temperature"]]
            temp_K = temp + 273.15
            self.weather[day].temperature["Vap"].fix(temp_K)

            rh = row[self.config.weather_data_column_dict["relative_humidity"]]
            self.weather[day].relative_humidity["H2O"].fix(rh / 100)
            self.weather[day].pressure_vap_sat["H2O"]

            for p, j in self.weather[day].phase_component_set:
                if (p, j) == ("Vap", "H2O"):
                    continue
                self.weather[day].flow_mass_phase_comp[(p, j)].fix(0)

        self.rh_min = Param(
            self.days_of_year,
            initialize=weather_daily_min[
                self.config.weather_data_column_dict["relative_humidity"]
            ]
            * 1e-2,
            units=pyunits.dimensionless,
            doc="Daily minimum relative humidity",
        )

        self.rh_max = Param(
            self.days_of_year,
            initialize=weather_daily_max[
                self.config.weather_data_column_dict["relative_humidity"]
            ]
            * 1e-2,
            units=pyunits.dimensionless,
            doc="Daily maximum relative humidity",
        )

        self.air_temp_min = Param(
            self.days_of_year,
            initialize=weather_daily_min[temp_col] + 273.15,
            units=pyunits.degK,
            doc="Daily minimum air temperature",
        )

        self.air_temp_max = Param(
            self.days_of_year,
            initialize=weather_daily_max[temp_col] + 273.15,
            units=pyunits.degK,
            doc="Daily maximum air temperature",
        )

        self.water_temp_param1 = Param(
            initialize=1.167,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Water temperature correlation slope",
        )

        self.water_temp_param2 = Param(
            initialize=0.175,
            mutable=True,
            units=pyunits.degK,
            doc="Water temperature correlation intercept",
        )
        self.water_activity_param1 = Param(
            initialize=-0.00056678,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Water activity correlation slope",
        )

        self.water_activity_param2 = Param(
            initialize=0.9985307,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Water activity correlation intercept",
        )

        self.emissivity_air_param = Param(
            initialize=1.24,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Air emissivity equation parameter",
        )

        self.emissivity_air_exp = Param(
            initialize=1 / 7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Air emissivity equation exponent",
        )

        self.area_correction_factor_base = Param(
            initialize=area_correction_factor_param_dict[self.config.dike_height][0],
            mutable=True,
            units=pyunits.dimensionless,
            doc="Area correction factor base",
        )

        self.area_correction_factor_exp = Param(
            initialize=area_correction_factor_param_dict[self.config.dike_height][1],
            mutable=True,
            units=pyunits.dimensionless,
            doc="Area correction factor exponent",
        )

        # solids precipitation rate is a function of TDS concentration [g / L]
        # solids_precipitation_rate [ft / yr] = a1 * C**2 + a2 * C + intercept
        # solids_precipitation_rate [ft / yr] = 4.12e-6 * C**2 + 1.92e-4 * C + 1.15e-3

        self.solids_precipitation_rate_a1 = Param(
            initialize=4.12e-6,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate a1 coefficient",
        )

        self.solids_precipitation_rate_a2 = Param(
            initialize=1.92e-4,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate a2 coefficient",
        )

        self.solids_precipitation_rate_intercept = Param(
            initialize=1.15e-3,
            mutable=True,
            units=pyunits.feet * pyunits.year**-1,
            doc="Solids precipitation rate intercept",
        )

        self.dens_solids = Param(
            initialize=2.16,
            mutable=True,
            units=pyunits.g / pyunits.cm**3,
            doc="Density of precipitated solids",
        )

        self.shortwave_radiation = Param(
            self.days_of_year,
            initialize=weather_daily_mean[
                self.config.weather_data_column_dict["shortwave_radiation"]
            ]
            * 0.0864,  # GHI column; W/m2 to MJ/day/m2
            units=pyunits.megajoule * pyunits.day**-1 * pyunits.m**-2,
            doc="Shortwave radiation at location (GHI)",
        )

        self.shortwave_albedo = Param(
            initialize=0.05,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Albedo factor for shortwave radiation",
        )

        self.longwave_albedo = Param(
            initialize=0.03,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Albedo factor for longwave radiation",
        )

        self.emissivity_water = Param(
            initialize=0.97,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Emissivity of water",
        )

        self.evaporation_pond_depth = Param(
            initialize=18,
            mutable=True,
            units=pyunits.inches,
            doc="Depth of evaporation pond",
        )

        self.evaporation_rate_salinity_adjustment_factor = Param(
            initialize=0.7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor to reduce evaporation rate for higher salinity",
        )

        # should be >1 for enhancement; value for organic dye around 8% (1.08)
        self.evaporation_rate_enhancement_adjustment_factor = Var(
            initialize=1,
            bounds=(1, 10),
            units=pyunits.dimensionless,
            doc="Factor to increase evaporation rate due to enhancement",
        )

        if not self.config.add_enhancement:
            self.evaporation_rate_enhancement_adjustment_factor.fix(1)

        self.net_radiation = Var(
            self.days_of_year,
            initialize=10,
            bounds=(0, None),
            units=pyunits.megajoule * pyunits.day**-1 * pyunits.m**-2,
            doc="Net radiation for evaporation",
        )

        self.mass_flux_water_vapor = Var(
            self.days_of_year,
            initialize=1e-5,
            bounds=(1e-12, 1e-3),
            units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            doc="Mass flux of water vapor evaporated according using BREB method",
        )

        self.area_correction_factor = Var(
            initialize=1,
            bounds=(0.99, 10),
            units=pyunits.dimensionless,
            doc="Area correction factor",
        )

        self.total_evaporative_area_required = Var(
            initialize=100000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total evaporative area required",
        )

        self.evaporative_area_per_pond = Var(
            initialize=10000,
            bounds=(0, 405000),
            units=pyunits.m**2,
            doc="Evaporative area required per pond",
        )

        self.evaporation_pond_area = Var(
            initialize=10000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Area of single evaporation pond",
        )

        self.number_evaporation_ponds = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of evaporation ponds",
        )

        self.solids_precipitation_rate = Var(
            initialize=0.01,
            bounds=(0, None),
            units=pyunits.feet / pyunits.year,
            doc="Rate at which solids precipitate on bottom of pond",
        )

        if (
            self.config.property_package.config.saturation_vapor_pressure_calculation
            == SaturationVaporPressureCalculation.Huang
        ):
            a = self.config.property_package.huang_coeff_A
            b_ = self.config.property_package.huang_coeff_B
            c = self.config.property_package.huang_coeff_C
            d1 = self.config.property_package.huang_coeff_D1
            d2 = self.config.property_package.huang_coeff_D2

            @self.Expression(
                self.days_of_year,
                doc="Huang saturation vapor pressure at min air temp",
            )
            def huang_press_sat_vap_min_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_min[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                huang_exp = a - (b_ / (air_temp + d1))
                p_vap_sat = (exp(huang_exp) / (air_temp + d2) ** c) * pyunits.Pa
                return p_vap_sat

            @self.Expression(
                self.days_of_year,
                doc="Huang saturation vapor pressure at max air temp",
            )
            def huang_press_sat_vap_max_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_max[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                huang_exp = a - (b_ / (air_temp + d1))
                p_vap_sat = (exp(huang_exp) / (air_temp + d2) ** c) * pyunits.Pa
                return p_vap_sat

            @self.Expression(
                self.days_of_year, doc="Actual vapor pressure from air temperature"
            )
            def actual_vapor_pressure(b, d):
                pressure_sat_vap_min = pyunits.convert(
                    b.huang_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.Pa,
                )
                pressure_sat_vap_max = pyunits.convert(
                    b.huang_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.Pa,
                )
                return (pressure_sat_vap_min + pressure_sat_vap_max) * 0.5

        if (
            self.config.property_package.config.saturation_vapor_pressure_calculation
            == SaturationVaporPressureCalculation.ArdenBuck
        ):

            a = self.config.property_package.arden_buck_coeff_a  # millibars
            b_ = self.config.property_package.arden_buck_coeff_b
            c = self.config.property_package.arden_buck_coeff_c
            d_ = self.config.property_package.arden_buck_coeff_d

            @self.Expression(
                self.days_of_year,
                doc="Arden-Buck saturation vapor pressure at min air temp",
            )
            def arden_buck_press_sat_vap_min_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_min[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                ardenbuck_exp = (b_ - air_temp / d_) * (air_temp / (c + air_temp))
                return pyunits.convert(a * exp(ardenbuck_exp), to_units=pyunits.Pa)

            @self.Expression(
                self.days_of_year,
                doc="Arden-Buck saturation vapor pressure at max air temp",
            )
            def arden_buck_press_sat_vap_max_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_max[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                ardenbuck_exp = (b_ - air_temp / d_) * (air_temp / (c + air_temp))
                return pyunits.convert(a * exp(ardenbuck_exp), to_units=pyunits.Pa)

            @self.Expression(
                self.days_of_year, doc="Actual vapor pressure from air temperature"
            )
            def actual_vapor_pressure(b, d):
                pressure_sat_vap_min = pyunits.convert(
                    b.arden_buck_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.Pa,
                )
                pressure_sat_vap_max = pyunits.convert(
                    b.arden_buck_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.Pa,
                )
                return (pressure_sat_vap_min + pressure_sat_vap_max) * 0.5

        if (
            self.config.property_package.config.saturation_vapor_pressure_calculation
            == SaturationVaporPressureCalculation.Antoine
        ):
            a = self.config.property_package.antoine_A
            b_ = self.config.property_package.antoine_B
            c = self.config.property_package.antoine_C

            @self.Expression(
                self.days_of_year,
                doc="Antoine saturation vapor pressure at min air temp",
            )
            def antoine_press_sat_vap_min_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_min[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                antoine = a - (b / (c + air_temp))
                p_vap_sat = 10 ** (antoine) * pyunits.mmHg
                return pyunits.convert(p_vap_sat, to_units=pyunits.Pa)

            @self.Expression(
                self.days_of_year,
                doc="Antoine saturation vapor pressure at max air temp",
            )
            def antoine_press_sat_vap_max_temp(b, d):
                air_temp = pyunits.convert(
                    (b.air_temp_max[d] - 273.15 * pyunits.degK) * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )
                antoine = a - (b / (c + air_temp))
                p_vap_sat = 10 ** (antoine) * pyunits.mmHg
                return pyunits.convert(p_vap_sat, to_units=pyunits.Pa)

            @self.Expression(
                self.days_of_year, doc="Actual vapor pressure from air temperature"
            )
            def actual_vapor_pressure(b, d):
                pressure_sat_vap_min = pyunits.convert(
                    b.antoine_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.Pa,
                )
                pressure_sat_vap_max = pyunits.convert(
                    b.antoine_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.Pa,
                )
                return (pressure_sat_vap_min + pressure_sat_vap_max) * 0.5

        @self.Expression(doc="Water activity")
        def water_activity(b):
            salinity = pyunits.convert(
                prop_in.conc_mass_phase_comp["Liq", "TDS"]
                * pyunits.kg**-1
                * pyunits.m**3,
                to_units=pyunits.dimensionless,
            )
            return b.water_activity_param1 * salinity + b.water_activity_param2

        @self.Expression(self.days_of_year, doc="Evaporation rate")
        def evaporation_rate(b, d):
            return pyunits.convert(
                b.mass_flux_water_vapor[d] / prop_in.dens_mass_solvent["H2O"],
                to_units=pyunits.m / pyunits.s,
            )

        @self.Expression(doc="Average mass flux of water vapor leaving pond over year")
        def mass_flux_water_vapor_average(b):
            return sum(b.mass_flux_water_vapor[d] for d in self.days_of_year) / 365

        @self.Expression(doc="Evaporative area per pond in acres")
        def evaporative_area_acre(b):
            return pyunits.convert(b.evaporative_area_per_pond, to_units=pyunits.acre)

        @self.Expression(doc="Total pond area in acres")
        def total_pond_area_acre(b):
            return pyunits.convert(
                b.evaporation_pond_area * b.number_evaporation_ponds,
                to_units=pyunits.acre,
            )

        @self.Expression(doc="Mass flow of precipitated solids")
        def mass_flow_precipitate(b):
            return pyunits.convert(
                b.total_pond_area_acre * b.solids_precipitation_rate * b.dens_solids,
                to_units=pyunits.kg / pyunits.year,
            )

        @self.Expression(self.days_of_year, doc="Emissivity of air")
        def emissivity_air(b, d):
            p_sat_kPa = pyunits.convert(
                pyunits.convert(b.actual_vapor_pressure[d], to_units=pyunits.kilopascal)
                * pyunits.kilopascal**-1,
                to_units=pyunits.dimensionless,
            )
            temp_C = pyunits.convert(
                (b.weather[d].temperature["Vap"] - 273.15 * pyunits.degK)
                * pyunits.degK**-1,
                to_units=pyunits.dimensionless,
            )
            air_temp_C = smooth_max(temp_C, 0.1)

            return smooth_min(
                b.emissivity_air_param
                * ((p_sat_kPa / air_temp_C) ** (b.emissivity_air_exp)),
                0.99,
            )

        @self.Expression(self.days_of_year, doc="Incident longwave radiation")
        def longwave_radiation_in(b, d):
            return pyunits.convert(
                b.emissivity_air[d]
                * Constants.stefan_constant
                * b.weather[d].temperature["Vap"] ** 4,
                to_units=pyunits.megajoule * pyunits.day**-1 * pyunits.m**-2,
            )

        @self.Expression(self.days_of_year, doc="Net incident shortwave radiation")
        def net_shortwave_radiation_in(b, d):
            return (1 - b.shortwave_albedo) * b.shortwave_radiation[d]

        @self.Expression(self.days_of_year, doc="Net incident longwave radiation")
        def net_longwave_radiation_in(b, d):
            return (1 - b.longwave_albedo) * b.longwave_radiation_in[d]

        @self.Expression(self.days_of_year, doc="Net outgoing longwave radiation")
        def net_longwave_radiation_out(b, d):
            return pyunits.convert(
                b.emissivity_water
                * Constants.stefan_constant
                * b.weather[d].temperature["Liq"] ** 4,
                to_units=pyunits.megajoule * pyunits.day**-1 * pyunits.m**-2,
            )

        @self.Expression(self.days_of_year, doc="Net solar radiation for evaporation")
        def net_solar_radiation(b, d):
            return (
                b.net_shortwave_radiation_in[d]
                + b.net_longwave_radiation_in[d]
                - b.net_longwave_radiation_out[d]
            )

        @self.Expression(self.days_of_year, doc="Psychrometric constant equation")
        def psychrometric_constant(b, d):
            mw_ratio = pyunits.convert(
                prop_in.mw_comp["H2O"] / prop_in.mw_comp["Air"],
                to_units=pyunits.dimensionless,
            )
            return pyunits.convert(
                (prop_in.cp_air * b.weather[d].pressure)
                / (mw_ratio * b.weather[d].dh_vap_mass_solvent),
                to_units=pyunits.kPa * pyunits.degK**-1,
            )

        @self.Expression(self.days_of_year, doc="Bowen ratio calculation")
        def bowen_ratio(b, d):
            return pyunits.convert(
                b.psychrometric_constant[d]
                * (
                    (b.weather[d].temperature["Liq"] - b.weather[d].temperature["Vap"])
                    / (
                        b.water_activity * b.weather[d].pressure_vap_sat["H2O"]
                        - b.actual_vapor_pressure[d]
                    )
                ),
                to_units=pyunits.dimensionless,
            )

        @self.Expression(self.days_of_year, doc="Daily water temperature change")
        def daily_temperature_change(b, d):
            if d == b.days_of_year.first():
                # For Jan 1, previous day is Dec 31
                return (
                    b.weather[d].temperature["Liq"]
                    - b.weather[b.days_of_year.last()].temperature["Liq"]
                ) * pyunits.day**-1
            else:
                return (
                    b.weather[d].temperature["Liq"]
                    - b.weather[d - 1].temperature["Liq"]
                ) * pyunits.day**-1

        @self.Expression(
            self.days_of_year,
            doc="Net heat flux in pond (to/from water, soil, ecosystem, etc.)",
        )
        def net_heat_flux_pond(b, d):
            return pyunits.convert(
                prop_in.dens_mass_phase["Liq"]
                * prop_in.cp_mass_solvent["Liq"]
                * b.evaporation_pond_depth
                * b.daily_temperature_change[d],
                to_units=pyunits.megajoule * pyunits.day**-1 * pyunits.m**-2,
            )

        @self.Constraint(self.days_of_year, doc="Daily mass flow of water")
        def eq_flow_mass_evap(b, d):
            return b.weather[d].flow_mass_phase_comp["Vap", "H2O"] == pyunits.convert(
                b.mass_flux_water_vapor[d] * b.total_evaporative_area_required,
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(self.days_of_year)
        def eq_water_temp(b, d):
            air_temp_C = b.weather[d].temperature["Vap"] - 273.15 * pyunits.degK
            # Minimum water temp is 0.1 degC
            return b.weather[d].temperature["Liq"] == smooth_max(
                (b.water_temp_param1 * air_temp_C - b.water_temp_param2)
                + 273.15 * pyunits.degK,
                273.25,
            )

        @self.Constraint(self.days_of_year, doc="Net radiation")
        def eq_net_radiation(b, d):
            return b.net_radiation[d] == smooth_max(
                b.net_solar_radiation[d] - b.net_heat_flux_pond[d], 1e-3
            )

        @self.Constraint(
            self.days_of_year, doc="Mass flux water vapor using BREB method"
        )
        def eq_mass_flux_water_vapor(b, d):
            return b.mass_flux_water_vapor[d] == (
                b.evaporation_rate_salinity_adjustment_factor
                * b.evaporation_rate_enhancement_adjustment_factor
            ) * pyunits.convert(
                b.net_radiation[d]
                / (b.weather[d].dh_vap_mass_solvent * (1 + b.bowen_ratio[d])),
                to_units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            )

        @self.Constraint(doc="Total evaporative area required")
        def eq_total_evaporative_area_required(b):
            return (
                b.total_evaporative_area_required * b.mass_flux_water_vapor_average
                == prop_in.flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(doc="Evaporation pond area")
        def eq_evaporation_pond_area(b):
            return (
                b.evaporation_pond_area
                == b.evaporative_area_per_pond * b.area_correction_factor
            )

        @self.Constraint(doc="Total evaporation pond area")
        def eq_evaporative_area_per_pond(b):
            return (
                b.number_evaporation_ponds * b.evaporative_area_per_pond
                == b.total_evaporative_area_required
            )

        @self.Constraint(doc="Area correction factor calculation")
        def eq_area_correction_factor(b):
            evap_per_pond_acre_dim = pyunits.convert(
                b.evaporative_area_per_pond * pyunits.acre**-1,
                to_units=pyunits.dimensionless,
            )
            return (
                b.area_correction_factor
                == b.area_correction_factor_base
                * evap_per_pond_acre_dim**b.area_correction_factor_exp
            )

        @self.Constraint(doc="Solids precipitation rate")
        def eq_solids_precipitation_rate(b):
            tds_in_dim = pyunits.convert(
                prop_in.conc_mass_phase_comp["Liq", "TDS"] * pyunits.g**-1 * pyunits.L,
                to_units=pyunits.dimensionless,
            )
            return (
                b.solids_precipitation_rate
                == b.solids_precipitation_rate_a1 * tds_in_dim**2
                + b.solids_precipitation_rate_a2 * tds_in_dim
                + b.solids_precipitation_rate_intercept
            )

    def initialize(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        self.weather.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=False,
        )
        init_log.info("Initialization of Weather Complete.")

        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}.".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.area_correction_factor) is None:
            iscale.set_scaling_factor(self.area_correction_factor, 1)

        if iscale.get_scaling_factor(self.mass_flux_water_vapor) is None:
            iscale.set_scaling_factor(self.mass_flux_water_vapor, 1e5)

        if iscale.get_scaling_factor(self.solids_precipitation_rate) is None:
            iscale.set_scaling_factor(self.solids_precipitation_rate, 1e2)

        if iscale.get_scaling_factor(self.total_evaporative_area_required) is None:
            iscale.set_scaling_factor(self.total_evaporative_area_required, 1e-2)

        if iscale.get_scaling_factor(self.number_evaporation_ponds) is None:
            iscale.set_scaling_factor(self.number_evaporation_ponds, 1)

        if iscale.get_scaling_factor(self.evaporative_area_per_pond) is None:
            iscale.set_scaling_factor(self.evaporative_area_per_pond, 1e-3)

        if iscale.get_scaling_factor(self.evaporation_pond_area) is None:
            iscale.set_scaling_factor(self.evaporation_pond_area, 1e-3)

        if iscale.get_scaling_factor(self.net_radiation) is None:
            iscale.set_scaling_factor(self.net_radiation, 1)

    def _get_timeseries_results(self):
        """
        Extract the values of all the model parameters
        indexed by the Set days_of_year
        """
        from collections import defaultdict

        pond_vars = [
            "evaporation_rate",
            "net_radiation",
            "net_solar_radiation",
            "longwave_radiation_in",
            "net_shortwave_radiation_in",
            "net_longwave_radiation_in",
            "net_longwave_radiation_out",
            "net_heat_flux_pond",
            "emissivity_air",
            "psychrometric_constant",
            "bowen_ratio",
            "mass_flux_water_vapor",
            "shortwave_radiation",
            "daily_temperature_change",
            "actual_vapor_pressure",
            "rh_min",
            "rh_max",
            "air_temp_min",
            "air_temp_max",
        ]

        pt = defaultdict(list)

        for day in self.days_of_year:

            pt["day_of_year"].append(day)
            for pv in pond_vars:
                v = self.find_component(pv)[day]
                pt[pv].append(value(v))

            w = self.weather[day]
            pt["atmospheric_pressure"].append(value(w.pressure))
            pt["temperature_air"].append(value(w.temperature["Vap"]))
            pt["temperature_water"].append(value(w.temperature["Liq"]))
            pt["pressure_vap"].append(value(w.pressure_vap["H2O"]))
            pt["pressure_vap_sat"].append(value(w.pressure_vap_sat["H2O"]))
            pt["dh_vap_mass_solvent"].append(value(w.dh_vap_mass_solvent))
            pt["relative_humidity"].append(value(w.relative_humidity["H2O"]))

        pond_timeseries = pd.DataFrame.from_dict(pt)

        return pond_timeseries

    def _get_stream_table_contents(self, time_point=0):

        return create_stream_table_dataframe(
            {"Feed Inlet": self.inlet},
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = dict()
        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_evaporation_pond
