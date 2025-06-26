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
    Constraint,
    check_optimal_termination,
    Param,
    Suffix,
    value,
    exp,
    atan,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

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

from watertap_contrib.reflo.costing.units.waiv import cost_waiv

__author__ = "Kurban Sitterley"


"""
References

"""
# Weather data from https://nsrdb.nrel.gov/data-viewer
default_weather_data_column_dict = {
    "pressure": "Pressure",
    "temperature": "Temperature",
    "relative_humidity": "Relative Humidity",
    "wind_speed": "Wind Speed",
}


@declare_process_block_class("WAIV")
class WAIVData(InitializationMixin, UnitModelBlockData):
    """
    Wind-Aided Intensified eVaporation (WAIV) model
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
        "terminal_process",
        ConfigValue(
            default=True,
            domain=bool,
            description="Flag to indicate if WAIV is terminal process in treatment train",
            doc="""Flag to indicate if WAIV is terminal process in treatment train.
            If True, assumed that all inlet water is evaporated and no outlet properties are not added.
            If False, user can set recovery_mass_water to indicate mass fraction of inlet water recovered from WAIV.
            Outlet properties and Port are added to unit model.
    """,
        ),
    )

    def build(self):
        super().build()

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

        for day, row in weather_daily_mean.iterrows():

            pres_Pa = pyunits.convert(
                row[self.config.weather_data_column_dict["pressure"]] * pyunits.mbar,
                to_units=pyunits.Pa,
            )
            self.weather[day].pressure.fix(pres_Pa)

            temp = row[self.config.weather_data_column_dict["temperature"]]
            temp_K = temp + 273.15
            self.weather[day].temperature["Vap"].fix(temp_K)
            self.weather[day].temperature["Liq"].set_value(temp_K)

            rh = row[self.config.weather_data_column_dict["relative_humidity"]]
            self.weather[day].relative_humidity["H2O"].fix(rh / 100)
            self.weather[day].pressure_vap_sat["H2O"]

            for p, j in self.weather[day].phase_component_set:
                if (p, j) == ("Vap", "H2O"):
                    continue
                self.weather[day].flow_mass_phase_comp[(p, j)].fix(0)

        self.wind_speed = Param(
            self.days_of_year,
            initialize=weather_daily_mean[
                self.config.weather_data_column_dict["wind_speed"]
            ],
            units=pyunits.m / pyunits.s,
            doc="Daily average wind speed",
        )

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

        self.temp_wb_a = Param(
            initialize=0.151977,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter A",
        )

        self.temp_wb_b = Param(
            initialize=8.313659,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter B",
        )

        self.temp_wb_c = Param(
            initialize=0.00391838,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter C",
        )

        self.temp_wb_d = Param(
            initialize=0.023101,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter D",
        )

        self.temp_wb_e = Param(
            initialize=1.676331,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter E",
        )

        self.temp_wb_f = Param(
            initialize=4.686035,
            units=pyunits.dimensionless,
            doc="Wet bulb temperature equation, parameter F",
        )

        self.waiv_module_area = Param(
            initialize=5700,
            mutable=True,
            units=pyunits.m**2,
            doc="Wettable surface area of single WAIV module",
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

        self.evaporation_rate_salinity_adjustment_factor = Param(
            initialize=0.7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor to reduce evaporation rate for higher salinity",
        )

        self.harbeck_N_base = Param(
            initialize=3.719e-9,
            units=pyunits.millibar**-1,
            doc="Harbeck N equation, base parameter",
        )

        self.harbeck_N_exp = Param(
            initialize=-0.0459,
            units=pyunits.dimensionless,
            doc="Harbeck N equation, exponent",
        )

        self.geotextile_area_waiv_module = Param(
            initialize=2850,
            mutable=True,
            units=pyunits.m**2,
            doc="Textile area of single WAIV module",
        )

        self.land_area_per_waiv_module = Param(
            initialize=364,
            mutable=True,
            units=pyunits.m**2,
            doc="Land area required for single WAIV system (tank and other process equipment)",
        )

        self.wetted_surface_to_geotextile_area_factor = Param(
            initialize=0.5,
            units=pyunits.m**2 / pyunits.m**2,
            doc="Textile area required per wetted area",
        )

        self.recovery_mass_water = Var(
            initialize=0.2,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Mass-based recovery of water from WAIV system",
        )

        self.number_recirculation_loops = Var(
            initialize=1.02,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of recirculation loops through WAIV modules",
        )

        self.evaporation_rate = Var(
            self.days_of_year,
            initialize=1,
            bounds=(0, None),
            units=pyunits.mm / pyunits.day,
            doc="Evaporation rate as determined by Harbeck equation",
        )

        self.total_wetted_surface_area_required = Var(
            initialize=100000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total wetted surface area required",
        )

        self.number_waiv_modules = Var(
            initialize=1,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of WAIV modules",
        )

        @self.Constraint(self.days_of_year, doc="Wet bulb temperature")
        def eq_temperature_wet_bulb(b, d):
            rh_pct = b.weather[d].relative_humidity["H2O"] * 100
            temp_degC = pyunits.convert(
                (b.weather[d].temperature["Vap"] - 273.15 * pyunits.degK)
                * pyunits.degK**-1,
                to_units=pyunits.dimensionless,
            )
            return (
                b.weather[d].temperature["Liq"]
                == (
                    (
                        temp_degC
                        * atan(b.temp_wb_a * (rh_pct + b.temp_wb_b) ** (1 / 2))
                        + b.temp_wb_c
                        * ((rh_pct) ** 3) ** (1 / 2)
                        * atan(b.temp_wb_d * rh_pct)
                        - atan(rh_pct - b.temp_wb_e)
                        + atan(temp_degC + rh_pct)
                        - b.temp_wb_f
                    )
                    * pyunits.degK
                    * pyunits.rad**-1
                )
                + 273.15 * pyunits.degK
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
                pressure_vap_sat_min = pyunits.convert(
                    b.huang_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.Pa,
                )
                pressure_vap_sat_max = pyunits.convert(
                    b.huang_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.Pa,
                )
                return (pressure_vap_sat_min + pressure_vap_sat_max) * 0.5

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
                pressure_vap_sat_min = pyunits.convert(
                    b.arden_buck_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.millibar,
                )
                pressure_vap_sat_max = pyunits.convert(
                    b.arden_buck_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.millibar,
                )
                return (pressure_vap_sat_min + pressure_vap_sat_max) * 0.5

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
                pressure_vap_sat_min = pyunits.convert(
                    b.antoine_press_sat_vap_min_temp[d] * b.rh_max[d],
                    to_units=pyunits.Pa,
                )
                pressure_vap_sat_max = pyunits.convert(
                    b.antoine_press_sat_vap_max_temp[d] * b.rh_min[d],
                    to_units=pyunits.Pa,
                )
                return (pressure_vap_sat_min + pressure_vap_sat_max) * 0.5

        @self.Expression(doc="Water activity")
        def water_activity(b):
            salinity = pyunits.convert(
                prop_in.conc_mass_phase_comp["Liq", "TDS"]
                * pyunits.kg**-1
                * pyunits.m**3,
                to_units=pyunits.dimensionless,
            )
            return b.water_activity_param1 * salinity + b.water_activity_param2

        @self.Expression(doc="Harbeck N parameter")
        def harbeck_N(b):
            area_tot_dimensionless = pyunits.convert(
                b.total_wetted_surface_area_required * pyunits.m**-2,
                to_units=pyunits.dimensionless,
            )
            return b.harbeck_N_base * area_tot_dimensionless**b.harbeck_N_exp

        @self.Expression(
            self.days_of_year, doc="Mass transfer driving force for evaporation"
        )
        def mass_transfer_driving_force(b, d):
            return pyunits.convert(
                b.weather[d].pressure_vap_sat["H2O"] * b.water_activity,
                to_units=pyunits.millibar,
            ) - pyunits.convert(b.actual_vapor_pressure[d], to_units=pyunits.millibar)

        @self.Expression(doc="Average annual evaporation rate for WAIV")
        def evaporation_rate_avg(b):
            return sum(b.evaporation_rate[d] for d in b.days_of_year) / 365

        @self.Expression(
            doc="Geotextile membrane area required to achieve wetted surface area required"
        )
        def geotextile_area_required(b):
            return pyunits.convert(
                b.total_wetted_surface_area_required
                * b.wetted_surface_to_geotextile_area_factor,
                to_units=pyunits.m**2,
            )

        @self.Expression(doc="Total land area required")
        def total_land_area_required(b):
            return b.land_area_per_waiv_module * b.number_waiv_modules

        @self.Expression(doc="Steady-state volumetric flow evaporated water")
        def flow_vol_water_evap(b):
            return pyunits.convert(
                b.total_wetted_surface_area_required
                * b.number_recirculation_loops
                * b.evaporation_rate_avg,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Expression(doc="Steady-state mass flow evaporated water")
        def flow_mass_water_evap(b):
            return pyunits.convert(
                b.flow_vol_water_evap * prop_in.dens_mass_solvent["H2O"],
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(self.days_of_year, doc="Evaporation rate from Harbeck method")
        def eq_evaporation_rate(b, d):
            wind_speed_mm_day = pyunits.convert(
                b.wind_speed[d], to_units=pyunits.mm / pyunits.day
            )
            return b.evaporation_rate[d] == smooth_max(
                pyunits.convert(
                    b.harbeck_N
                    * wind_speed_mm_day
                    * b.mass_transfer_driving_force[d]
                    * b.evaporation_rate_salinity_adjustment_factor,
                    to_units=pyunits.mm / pyunits.day,
                ),
                1e-3 * pyunits.mm / pyunits.day,
            )

        @self.Constraint(doc="Total wetted surface area required for WAIV system")
        def eq_total_wetted_surface_area_required(b):
            return (
                b.total_wetted_surface_area_required * b.number_recirculation_loops
                == pyunits.convert(
                    (
                        (
                            prop_in.flow_mass_phase_comp["Liq", "H2O"]
                            * (1 - b.recovery_mass_water)
                        )
                        / prop_in.dens_mass_solvent["H2O"]
                    )
                    / b.evaporation_rate_avg,
                    to_units=pyunits.m**2,
                )
            )

        @self.Constraint(doc="Number of WAIV modules required")
        def eq_number_waiv_modules(b):
            return (
                b.number_waiv_modules
                == b.geotextile_area_required / b.geotextile_area_waiv_module
            )

        @self.Constraint(self.days_of_year, doc="Daily mass flow of water evaporated")
        def eq_flow_mass_evap(b, d):
            return b.weather[d].flow_mass_phase_comp["Vap", "H2O"] == pyunits.convert(
                b.evaporation_rate[d]
                * b.total_wetted_surface_area_required
                * prop_in.dens_mass_solvent["H2O"],
                to_units=pyunits.kg / pyunits.s,
            )

        if self.config.terminal_process:
            self.recovery_mass_water.fix(0)

        if not self.config.terminal_process:

            tmp_dict["defined_state"] = False

            self.properties_out = self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of outlet",
                **tmp_dict,
            )
            prop_out = self.properties_out[0]
            self.add_outlet_port(name="outlet", block=self.properties_out)

            self.max_brine_salinity = Param(
                initialize=300,
                mutable=True,
                units=pyunits.kg / pyunits.m**3,
                doc="Maximum salinity of brine",
            )

            self.brine_concentration_factor = Var(
                initialize=2,
                bounds=(1, None),
                units=pyunits.dimensionless,
                doc="Concentration factor of effluent brine relative to influent TDS",
            )

            @self.Constraint(self.config.property_package.phase_list, doc="Isothermal")
            def eq_isothermal(b, p):
                return prop_out.temperature[p] == prop_in.temperature[p]

            @self.Constraint(doc="Isobaric")
            def eq_isobaric(b):
                return prop_out.pressure == prop_in.pressure

            @self.Constraint(doc="Maximum brine salinity")
            def eq_max_brine_salintiy(b):
                return (
                    prop_out.conc_mass_phase_comp["Liq", "TDS"] <= b.max_brine_salinity
                )

            @self.Constraint(
                self.config.property_package.liq_comps,
                doc="Mass balance for liquid phase",
            )
            def eq_liq_comps_mass_balance(b, j):
                if j == "TDS":
                    return (
                        prop_out.conc_mass_phase_comp["Liq", j]
                        == prop_in.conc_mass_phase_comp["Liq", j]
                        * b.brine_concentration_factor
                    )
                if j == "H2O":
                    return (
                        prop_out.flow_mass_phase_comp["Liq", j]
                        == prop_in.flow_mass_phase_comp["Liq", j]
                        * b.recovery_mass_water
                    )
                if j in b.config.property_package.volatile_comps:
                    prop_out.flow_mass_phase_comp["Liq", j].fix(0)
                    return Constraint.Skip
                if j in b.config.property_package.non_volatile_comps:
                    return (
                        prop_out.flow_mass_phase_comp["Liq", j]
                        == prop_in.flow_mass_phase_comp["Liq", j]
                    )

            @self.Constraint(
                self.config.property_package.vap_comps,
                doc="Mass balance for vapor phase",
            )
            def eq_vap_comps_mass_balance(b, j):
                return (
                    prop_out.flow_mass_phase_comp["Vap", j]
                    == prop_in.flow_mass_phase_comp["Vap", j]
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

        for d in self.days_of_year:
            calculate_variable_from_constraint(
                self.weather[d].temperature["Liq"], self.eq_temperature_wet_bulb[d]
            )

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

        if not self.config.terminal_process:

            self.properties_out.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args,
                hold_state=False,
            )
            init_log.info("Initialization Step 1b Complete.")

        for d in self.days_of_year:
            calculate_variable_from_constraint(
                self.evaporation_rate[d], self.eq_evaporation_rate[d]
            )

        calculate_variable_from_constraint(
            self.total_wetted_surface_area_required,
            self.eq_total_wetted_surface_area_required,
        )

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

        if iscale.get_scaling_factor(self.total_wetted_surface_area_required) is None:
            iscale.set_scaling_factor(self.total_wetted_surface_area_required, 1e-2)

        if iscale.get_scaling_factor(self.evaporation_rate) is None:
            iscale.set_scaling_factor(self.evaporation_rate, 1)

        if iscale.get_scaling_factor(self.total_wetted_surface_area_required) is None:
            iscale.set_scaling_factor(self.total_wetted_surface_area_required, 1e-4)

        if iscale.get_scaling_factor(self.number_waiv_modules) is None:
            iscale.set_scaling_factor(self.number_waiv_modules, 1)

        if iscale.get_scaling_factor(self.recovery_mass_water) is None:
            iscale.set_scaling_factor(self.recovery_mass_water, 1)

        if iscale.get_scaling_factor(self.number_recirculation_loops) is None:
            iscale.set_scaling_factor(self.number_recirculation_loops, 1)

        if not self.config.terminal_process:
            if iscale.get_scaling_factor(self.brine_concentration_factor) is None:
                iscale.set_scaling_factor(self.brine_concentration_factor, 1)

    def _get_timeseries_results(self):
        """
        Extract the values of all the model parameters
        indexed by the Set days_of_year
        """
        from collections import defaultdict

        waiv_vars = [
            "evaporation_rate",
            "wind_speed",
            "actual_vapor_pressure",
            "rh_min",
            "rh_max",
            "air_temp_min",
            "air_temp_max",
            "mass_transfer_driving_force",
        ]

        wt = defaultdict(list)

        for day in self.days_of_year:

            wt["day_of_year"].append(day)
            for wv in waiv_vars:
                v = self.find_component(wv)[day]
                wt[wv].append(value(v))

            w = self.weather[day]
            wt["atmospheric_pressure"].append(value(w.pressure))
            wt["temperature_air"].append(value(w.temperature["Vap"]))
            wt["temperature_water"].append(value(w.temperature["Liq"]))
            wt["pressure_vap"].append(value(w.pressure_vap["H2O"]))
            wt["pressure_vap_sat"].append(value(w.pressure_vap_sat["H2O"]))
            wt["dh_vap_mass_solvent"].append(value(w.dh_vap_mass_solvent))
            wt["relative_humidity"].append(value(w.relative_humidity["H2O"]))
            wt["flow_mass_phase_vap_h2o"].append(
                value(w.flow_mass_phase_comp["Vap", "H2O"])
            )

        waiv_timeseries = pd.DataFrame.from_dict(wt)

        return waiv_timeseries

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
        return cost_waiv
