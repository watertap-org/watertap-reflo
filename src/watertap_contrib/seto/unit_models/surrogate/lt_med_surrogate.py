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

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    Param,
    Suffix,
    Constraint,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


@declare_process_block_class("LTMEDSurrogate")
class LTMEDData(UnitModelBlockData):
    """
    Low-temperature multi-effect distillation surrogate model
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
        "property_package_water",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for feed, distillate, brine, and cooling water properties",
            doc="""Property parameter object used to define water property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_steam",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for steam properties",
            doc="""Property parameter object used to define steam property calculations,
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

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        """
        Add system configurations
        """
        self.number_effects = Param(
            initialize=12,
            mutable=True,
            within=Set(initialize=[3, 6, 9, 12, 14]),
            units=pyunits.dimensionless,
            doc="Number of effects",
        )

        self.delta_T_last_effect = Param(
            initialize=10,
            mutable=True,
            units=pyunits.K,
            doc="Temperature increase in last effect",
        )

        self.delta_T_cooling_reject = Param(
            initialize=-3,
            mutable=True,
            units=pyunits.K,
            doc="Temperature decrease in cooling reject water",
        )

        self.recovery_ratio = Var(
            initialize=0.50,
            bounds=(0.30, 0.50),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        self.thermal_loss = Param(
            initialize=0.054,
            units=pyunits.dimensionless,
            doc="System thermal loss",  ## units??
        )

        """
        Add block for feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_water
        tmp_dict["defined_state"] = False

        self.feed_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict
        )

        """
        Add block for distillate
        """
        self.distillate_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of distillate",
            **tmp_dict
        )

        # Distillate temperature is the same as last effect vapor temperature,
        # which is assumed to be 10 deg higher than condenser inlet seawater temperature
        @self.Constraint(doc="Distillate temperature")
        def eq_distillate_temp(b):
            return (
                b.distillate_props[0].temperature
                == b.feed_props[0].temperature + b.delta_T_last_effect
            )

        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

        self.brine_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        # Relationship between brine salinity and feed salinity
        @self.Constraint(doc="Brine salinity")
        def eq_brine_salinity(b):
            return b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"] == b.feed_props[
                0
            ].conc_mass_phase_comp["Liq", "TDS"] / (1 - b.recovery_ratio)

        @self.Constraint(doc="Brine temperature")
        def eq_brine_temp(b):
            return (
                b.brine_props[0].temperature
                == b.distillate_props[0].temperature
                + b.brine_props[0].boiling_point_elevation_phase["Liq"]
            )

        """
        Add block for reject cooling water
        """
        self.cooling_out_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of cooling reject",
            **tmp_dict
        )

        # Use the source water for cooling (same salinity)
        # TODO: Fix this rather than constraint?
        @self.Constraint(doc="Cooling reject salinity")
        def eq_cooling_salinity(b):
            return (
                b.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
            )

        # Assumption: the temperature of cooling reject is 3 degC lower than in the condenser
        @self.Constraint(doc="Cooling reject temperature")
        def eq_cooling_temp(b):
            return (
                b.cooling_out_props[0].temperature
                == b.distillate_props[0].temperature + b.delta_T_cooling_reject
            )

        """
        Add block for heating steam
        """
        tmp_dict["parameters"] = self.config.property_package_steam
        tmp_dict["defined_state"] = True

        self.steam_props = self.config.property_package_steam.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of heating steam",
            **tmp_dict
        )

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="distillate", block=self.distillate_props)
        self.add_port(name="brine", block=self.brine_props)
        self.add_port(name="steam", block=self.steam_props)

        """
        Mass balances
        """
        # Distillate flow rate calculation
        @self.Constraint(doc="Distillate volumetric flow rate")
        def eq_dist_vol_flow(b):
            return (
                b.distillate_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"] * b.recovery_ratio
            )

        # Brine flow rate calculation
        @self.Constraint(doc="Brine volumetric flow rate")
        def eq_brine_vol_flow(b):
            return (
                b.brine_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                - b.distillate_props[0].flow_vol_phase["Liq"]
            )

        """
        Add Vars for model outputs
        """
        self.thermal_power_requirement = Var(
            initialize=5000,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW)",
        )

        self.specific_area = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**2 / (pyunits.m**3 / pyunits.d),
            doc="Specific area (m2/m3/day))",
        )

        self.specific_thermal_energy_consumption = Var(
            initialize=65,
            bounds=(0, None),
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific thermal power consumption (kWh/m3)",
        )

        self.gain_output_ratio = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg / pyunits.kg,
            doc="Gained output ratio (kg of distillate water per kg of heating steam)",
        )

        """
        Add Vars for intermediate model variables
        """
        # TODO: Make Expressions?
        self.feed_cool_mass_flow = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.kg / pyunits.s,
            doc="Feed + cooling water mass flow rate (kg/s)",
        )

        self.feed_cool_vol_flow = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.h,
            doc="Feed + cooling water volume flow rate (m3/h)",
        )

        @self.Expression(doc="Temperature in last effect")
        def temp_last_effect(b):
            return (
                b.feed_props[0].temperature - 273.15 * pyunits.K + b.delta_T_last_effect
            )

        # TODO: Make Expressions?
        feed_conc_ppm = pyunits.convert(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.mg / pyunits.L,
        )
        temp_steam = self.steam_props[0].temperature - 273.15

        # Surrogate equations for calculating gain output ratio
        gain_output_ratio_coeffs = self._get_gain_output_ratio_coeffs()

        @self.Constraint(doc="Gain output ratio surrogate equation")
        def eq_gain_output_ratio(b):
            return (
                b.gain_output_ratio
                == feed_conc_ppm * gain_output_ratio_coeffs[b.number_effects.value][0]
                + b.recovery_ratio * gain_output_ratio_coeffs[b.number_effects.value][1]
                + feed_conc_ppm
                * b.recovery_ratio
                * gain_output_ratio_coeffs[b.number_effects.value][2]
                + b.temp_last_effect
                * gain_output_ratio_coeffs[b.number_effects.value][3]
                + b.temp_last_effect
                * feed_conc_ppm
                * gain_output_ratio_coeffs[b.number_effects.value][4]
                + b.temp_last_effect
                * b.recovery_ratio
                * gain_output_ratio_coeffs[b.number_effects.value][5]
                + temp_steam * gain_output_ratio_coeffs[b.number_effects.value][6]
                + temp_steam
                * feed_conc_ppm
                * gain_output_ratio_coeffs[b.number_effects.value][7]
                + temp_steam
                * b.recovery_ratio
                * gain_output_ratio_coeffs[b.number_effects.value][8]
                + temp_steam
                * b.temp_last_effect
                * gain_output_ratio_coeffs[b.number_effects.value][9]
                + 1 * gain_output_ratio_coeffs[b.number_effects.value][10]
                + temp_steam**2 * gain_output_ratio_coeffs[b.number_effects.value][11]
                + b.temp_last_effect**2
                * gain_output_ratio_coeffs[b.number_effects.value][12]
                + b.recovery_ratio**2
                * gain_output_ratio_coeffs[b.number_effects.value][13]
                + feed_conc_ppm**2
                * gain_output_ratio_coeffs[b.number_effects.value][14]
            )

        specific_area_coeffs = self._get_specific_area_coeffs()

        @self.Constraint(doc="Specific area surrogate equation")
        def eq_specific_area(b):
            if b.number_effects.value in [3, 6, 9]:
                return (
                    b.specific_area
                    == feed_conc_ppm * specific_area_coeffs[b.number_effects.value][0]
                    + feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][1]
                    + b.recovery_ratio * specific_area_coeffs[b.number_effects.value][2]
                    + b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][3]
                    + b.recovery_ratio
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][4]
                    + b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][5]
                    + b.recovery_ratio**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][6]
                    + b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][7]
                    + b.temp_last_effect
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][8]
                    + b.temp_last_effect
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][9]
                    + b.temp_last_effect
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][10]
                    + b.temp_last_effect
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][11]
                    + b.temp_last_effect
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][12]
                    + b.temp_last_effect**2
                    * specific_area_coeffs[b.number_effects.value][13]
                    + b.temp_last_effect**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][14]
                    + b.temp_last_effect**2
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][15]
                    + temp_steam * specific_area_coeffs[b.number_effects.value][16]
                    + temp_steam
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][17]
                    + temp_steam
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][18]
                    + temp_steam
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][19]
                    + temp_steam
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][20]
                    + temp_steam
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][21]
                    + temp_steam
                    * b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][22]
                    + temp_steam
                    * b.temp_last_effect
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][23]
                    + temp_steam
                    * b.temp_last_effect
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][24]
                    + temp_steam
                    * b.temp_last_effect**2
                    * specific_area_coeffs[b.number_effects.value][25]
                    + temp_steam**2 * specific_area_coeffs[b.number_effects.value][26]
                    + temp_steam**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][27]
                    + temp_steam**2
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][28]
                    + temp_steam**2
                    * b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][29]
                    + 1 * specific_area_coeffs[b.number_effects.value][30]
                    + temp_steam**3 * specific_area_coeffs[b.number_effects.value][31]
                    + b.temp_last_effect**3
                    * specific_area_coeffs[b.number_effects.value][32]
                    + b.recovery_ratio**3
                    * specific_area_coeffs[b.number_effects.value][33]
                    + feed_conc_ppm**3
                    * specific_area_coeffs[b.number_effects.value][34]
                )

            else:
                return (
                    b.specific_area
                    == feed_conc_ppm * specific_area_coeffs[b.number_effects.value][0]
                    + feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][1]
                    + feed_conc_ppm**3
                    * specific_area_coeffs[b.number_effects.value][2]
                    + b.recovery_ratio * specific_area_coeffs[b.number_effects.value][3]
                    + b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][4]
                    + b.recovery_ratio
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][5]
                    + b.recovery_ratio
                    * feed_conc_ppm**3
                    * specific_area_coeffs[b.number_effects.value][6]
                    + b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][7]
                    + b.recovery_ratio**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][8]
                    + b.recovery_ratio**2
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][9]
                    + b.recovery_ratio**3
                    * specific_area_coeffs[b.number_effects.value][10]
                    + b.recovery_ratio**3
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][11]
                    + b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][12]
                    + b.temp_last_effect
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][13]
                    + b.temp_last_effect
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][14]
                    + b.temp_last_effect
                    * feed_conc_ppm**3
                    * specific_area_coeffs[b.number_effects.value][15]
                    + b.temp_last_effect
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][16]
                    + b.temp_last_effect
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][17]
                    + b.temp_last_effect
                    * b.recovery_ratio
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][18]
                    + b.temp_last_effect
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][19]
                    + b.temp_last_effect
                    * b.recovery_ratio**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][20]
                    + b.temp_last_effect
                    * b.recovery_ratio**3
                    * specific_area_coeffs[b.number_effects.value][21]
                    + b.temp_last_effect**2
                    * specific_area_coeffs[b.number_effects.value][22]
                    + b.temp_last_effect**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][23]
                    + b.temp_last_effect**2
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][24]
                    + b.temp_last_effect**2
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][25]
                    + b.temp_last_effect**2
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][26]
                    + b.temp_last_effect**2
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][27]
                    + b.temp_last_effect**3
                    * specific_area_coeffs[b.number_effects.value][28]
                    + b.temp_last_effect**3
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][29]
                    + b.temp_last_effect**3
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][30]
                    + temp_steam * specific_area_coeffs[b.number_effects.value][31]
                    + temp_steam
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][32]
                    + temp_steam
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][33]
                    + temp_steam
                    * feed_conc_ppm**3
                    * specific_area_coeffs[b.number_effects.value][34]
                    + temp_steam
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][35]
                    + temp_steam
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][36]
                    + temp_steam
                    * b.recovery_ratio
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][37]
                    + temp_steam
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][38]
                    + temp_steam
                    * b.recovery_ratio**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][39]
                    + temp_steam
                    * b.recovery_ratio**3
                    * specific_area_coeffs[b.number_effects.value][40]
                    + temp_steam
                    * b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][41]
                    + temp_steam
                    * b.temp_last_effect
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][42]
                    + temp_steam
                    * b.temp_last_effect
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][43]
                    + temp_steam
                    * b.temp_last_effect
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][44]
                    + temp_steam
                    * b.temp_last_effect
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][45]
                    + temp_steam
                    * b.temp_last_effect
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][46]
                    + temp_steam
                    * b.temp_last_effect**2
                    * specific_area_coeffs[b.number_effects.value][47]
                    + temp_steam
                    * b.temp_last_effect**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][48]
                    + temp_steam
                    * b.temp_last_effect**2
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][49]
                    + temp_steam
                    * b.temp_last_effect**3
                    * specific_area_coeffs[b.number_effects.value][50]
                    + temp_steam**2 * specific_area_coeffs[b.number_effects.value][51]
                    + temp_steam**2
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][52]
                    + temp_steam**2
                    * feed_conc_ppm**2
                    * specific_area_coeffs[b.number_effects.value][53]
                    + temp_steam**2
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][54]
                    + temp_steam**2
                    * b.recovery_ratio
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][55]
                    + temp_steam**2
                    * b.recovery_ratio**2
                    * specific_area_coeffs[b.number_effects.value][56]
                    + temp_steam**2
                    * b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][57]
                    + temp_steam**2
                    * b.temp_last_effect
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][58]
                    + temp_steam**2
                    * b.temp_last_effect
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][59]
                    + temp_steam**2
                    * b.temp_last_effect**2
                    * specific_area_coeffs[b.number_effects.value][60]
                    + temp_steam**3 * specific_area_coeffs[b.number_effects.value][61]
                    + temp_steam**3
                    * feed_conc_ppm
                    * specific_area_coeffs[b.number_effects.value][62]
                    + temp_steam**3
                    * b.recovery_ratio
                    * specific_area_coeffs[b.number_effects.value][63]
                    + temp_steam**3
                    * b.temp_last_effect
                    * specific_area_coeffs[b.number_effects.value][64]
                    + 1 * specific_area_coeffs[b.number_effects.value][65]
                    + temp_steam**4 * specific_area_coeffs[b.number_effects.value][66]
                    + b.temp_last_effect**4
                    * specific_area_coeffs[b.number_effects.value][67]
                    + b.recovery_ratio**4
                    * specific_area_coeffs[b.number_effects.value][68]
                    + feed_conc_ppm**4
                    * specific_area_coeffs[b.number_effects.value][69]
                )

        # Steam flow rate calculation
        @self.Constraint(doc="Steam flow rate")
        def eq_steam_mass_flow(b):
            return (
                b.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
                == sum(
                    b.distillate_props[0].flow_mass_phase_comp["Liq", j]
                    for j in b.distillate_props.component_list
                )
                / b.gain_output_ratio
            )

        # Energy consumption
        @self.Constraint(doc="specific_thermal_energy_consumption calculation")
        def eq_specific_thermal_energy_consumption(b):
            return b.specific_thermal_energy_consumption == pyunits.convert(
                1
                / b.gain_output_ratio
                * (b.steam_props[0].dh_vap_mass)
                * b.distillate_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(doc="Thermal power requirement calculation")
        def eq_thermal_power_requirement(b):
            return b.thermal_power_requirement == pyunits.convert(
                b.specific_thermal_energy_consumption
                * b.distillate_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.kW,
            )

        # Mass flow rate
        @self.Constraint(doc="Feed and cooling water mass flow rate (kg/s)")
        def eq_feed_cool_mass_flow(b):
            feed = b.feed_props[0]
            cool = b.cooling_out_props[0]
            brine = b.brine_props[0]
            dist = b.distillate_props[0]
            feed_mass_flow_tot = sum(
                feed.flow_mass_phase_comp["Liq", j] for j in feed.component_list
            )
            brine_mass_flow_tot = sum(
                brine.flow_mass_phase_comp["Liq", j] for j in brine.component_list
            )
            dist_mass_flow_tot = sum(
                dist.flow_mass_phase_comp["Liq", j] for j in dist.component_list
            )
            enth_feed = pyunits.convert(
                feed.enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
            )
            enth_brine = pyunits.convert(
                brine.enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
            )
            enth_dist = pyunits.convert(
                dist.enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
            )
            enth_cool = pyunits.convert(
                cool.enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
            )
            return b.feed_cool_mass_flow * (enth_cool - enth_feed) == (
                (1 - b.thermal_loss) * b.thermal_power_requirement
                - brine_mass_flow_tot * enth_brine
                - dist_mass_flow_tot * enth_dist
                + enth_cool * feed_mass_flow_tot
            )

        # Volume flow rates
        @self.Constraint(doc="Feed and cooling water mass flow rate (m3/h)")
        def eq_feed_cool_vol_flow(b):
            return b.feed_cool_vol_flow == pyunits.convert(
                b.feed_cool_mass_flow / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hour,
            )

        @self.Constraint(doc="Cooling water mass flow rate (m3/h)")
        def eq_cool_vol_flow(b):
            return b.feed_cool_vol_flow == pyunits.convert(
                b.cooling_out_props[0].flow_vol_phase["Liq"]
                + b.feed_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        dist_vol_flow = self.distillate_props[0].flow_vol_phase["Liq"]

        if iscale.get_scaling_factor(self.recovery_ratio) is None:
            iscale.set_scaling_factor(self.recovery_ratio, 1e1)

        if iscale.get_scaling_factor(self.specific_thermal_energy_consumption) is None:
            iscale.set_scaling_factor(self.specific_thermal_energy_consumption, 1e-3)

        if iscale.get_scaling_factor(self.thermal_power_requirement) is None:
            iscale.set_scaling_factor(
                self.thermal_power_requirement,
                iscale.get_scaling_factor(dist_vol_flow)
                * iscale.get_scaling_factor(self.specific_thermal_energy_consumption),
            )

        if iscale.get_scaling_factor(self.specific_area) is None:
            iscale.set_scaling_factor(self.specific_area, 0.1)

        if iscale.get_scaling_factor(self.gain_output_ratio) is None:
            iscale.set_scaling_factor(self.gain_output_ratio, 0.1)

        if iscale.get_scaling_factor(self.feed_cool_mass_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_mass_flow, 1e-3)

        if iscale.get_scaling_factor(self.feed_cool_vol_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_vol_flow, 1e-3)

        # Transforming constraints

        sf = iscale.get_scaling_factor(self.distillate_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_distillate_temp, sf)

        sf = iscale.get_scaling_factor(self.cooling_out_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_cooling_temp, sf)

        sf = iscale.get_scaling_factor(
            self.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.eq_cooling_salinity, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_brine_temp, sf)

        sf = iscale.get_scaling_factor(
            self.brine_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.eq_brine_salinity, sf)

        sf = iscale.get_scaling_factor(self.feed_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_dist_vol_flow, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_brine_vol_flow, sf)

        sf = iscale.get_scaling_factor(self.gain_output_ratio)
        iscale.constraint_scaling_transform(self.eq_gain_output_ratio, sf)

        sf = iscale.get_scaling_factor(self.specific_area)
        iscale.constraint_scaling_transform(self.eq_specific_area, sf)

        sf = iscale.get_scaling_factor(
            self.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_steam_mass_flow, sf)

        sf = iscale.get_scaling_factor(self.specific_thermal_energy_consumption)
        iscale.constraint_scaling_transform(
            self.eq_specific_thermal_energy_consumption, sf
        )

        sf = iscale.get_scaling_factor(self.thermal_power_requirement)
        iscale.constraint_scaling_transform(self.eq_thermal_power_requirement, sf)

        sf = (
            iscale.get_scaling_factor(self.feed_cool_mass_flow)
            * iscale.get_scaling_factor(
                self.cooling_out_props[0].enth_mass_phase["Liq"]
            )
            * 1e-3
        )
        iscale.constraint_scaling_transform(self.eq_feed_cool_mass_flow, sf)

        sf = iscale.get_scaling_factor(self.feed_cool_vol_flow)
        iscale.constraint_scaling_transform(self.eq_feed_cool_vol_flow, sf)

        sf = (
            iscale.get_scaling_factor(self.cooling_out_props[0].flow_vol_phase["Liq"])
            / 3600
        )
        iscale.constraint_scaling_transform(self.eq_cool_vol_flow, sf)

    def _get_gain_output_ratio_coeffs(self):
        return {
            3: [
                1.60e-07,
                0.826895712,
                -2.04e-07,
                0.003340838,
                -5.56e-09,
                0.000666667,
                -0.003295958,
                1.17e-10,
                -0.000549708,
                -2.46e-06,
                2.662545127,
                -1.98e-07,
                7.41e-07,
                -0.675925926,
                -4.12e-13,
            ],
            6: [
                5.86e-07,
                2.940942982,
                -1.44e-06,
                0.007985234,
                -2.26e-08,
                0.001472222,
                -0.007157144,
                8.73e-09,
                -0.001540936,
                4.88e-06,
                4.741753363,
                -1.04e-06,
                -1.67e-05,
                -2.333333333,
                -7.41e-13,
            ],
            9: [
                1.67e-06,
                5.9507846,
                -3.94e-06,
                0.012607651,
                -5.50e-08,
                0.002222222,
                -0.010203548,
                2.47e-08,
                -0.002549708,
                2.94e-05,
                6.350104873,
                -1.05e-05,
                -5.09e-05,
                -4.653703704,
                -2.88e-12,
            ],
            12: [
                3.30e-06,
                9.621851852,
                -7.98e-06,
                0.016637037,
                -1.09e-07,
                0.002,
                -0.012637326,
                5.13e-08,
                -0.003277778,
                8.28e-05,
                7.592772368,
                -3.09e-05,
                -9.98e-05,
                -7.425925926,
                -6.09e-12,
            ],
            14: [
                5.27e-06,
                12.44928443,
                -1.15e-05,
                0.019398098,
                -1.59e-07,
                0.001666667,
                -0.013396636,
                7.88e-08,
                -0.003333333,
                0.000121663,
                8.195669495,
                -5.34e-05,
                -0.00013251,
                -9.627577763,
                -1.80e-11,
            ],
        }

    def _get_specific_area_coeffs(self):
        return {
            3: [
                0.000596217,
                -3.66e-09,
                0,
                -2.44e-05,
                1.93e-09,
                0,
                5.60e-05,
                0,
                -2.95e-07,
                5.30e-11,
                0,
                7.14e-06,
                0,
                0.064807392,
                2.06e-07,
                0.00974051,
                0,
                -1.16e-05,
                -3.96e-11,
                0,
                -5.27e-06,
                0,
                -0.05718687,
                -2.61e-07,
                -0.011936049,
                -0.000702529,
                0.013464849,
                1.65e-07,
                0.003686623,
                0.000759933,
                0,
                -0.00019293,
                -0.000182949,
                0,
                3.20e-14,
            ],
            6: [
                0.00040105,
                -6.57e-09,
                0,
                -1.56e-05,
                3.67e-10,
                0,
                2.62e-05,
                0,
                7.08e-07,
                8.73e-12,
                0,
                1.46e-06,
                0,
                0.032775092,
                5.04e-08,
                0.002499309,
                0,
                -3.30e-06,
                -6.53e-12,
                0,
                -1.02e-06,
                0,
                -0.028641745,
                -6.66e-08,
                -0.002735652,
                -0.000301667,
                0.005600544,
                4.10e-08,
                0.000713386,
                0.000329733,
                0,
                -7.31e-05,
                -0.00010089,
                0,
                4.94e-14,
            ],
            9: [
                0.000596217,
                -3.66e-09,
                0,
                -2.44e-05,
                1.93e-09,
                0,
                5.60e-05,
                0,
                -2.95e-07,
                5.30e-11,
                0,
                7.14e-06,
                0,
                0.064807392,
                2.06e-07,
                0.00974051,
                0,
                -1.16e-05,
                -3.96e-11,
                0,
                -5.27e-06,
                0,
                -0.05718687,
                -2.61e-07,
                -0.011936049,
                -0.000702529,
                0.013464849,
                1.65e-07,
                0.003686623,
                0.000759933,
                0,
                -0.00019293,
                -0.000182949,
                0,
                3.20e-14,
            ],
            12: [
                0.000000e00,
                3.304374e-08,
                -6.761157e-13,
                0.000000e00,
                0.000000e00,
                -5.496094e-09,
                1.958695e-13,
                0.000000e00,
                0.000000e00,
                6.105760e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -8.520444e-10,
                2.227182e-14,
                0.000000e00,
                0.000000e00,
                6.610470e-10,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                3.561099e-06,
                9.693068e-12,
                0.000000e00,
                1.148815e-06,
                0.000000e00,
                0.000000e00,
                -1.049434e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                1.734819e-10,
                -1.367980e-14,
                0.000000e00,
                0.000000e00,
                -5.044097e-10,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -2.717657e-06,
                -3.905605e-11,
                0.000000e00,
                -1.397796e-06,
                0.000000e00,
                0.000000e00,
                -4.132341e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                4.380618e-07,
                2.072263e-11,
                0.000000e00,
                4.398758e-07,
                0.000000e00,
                0.000000e00,
                5.991695e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -1.833849e-08,
                0.000000e00,
                -3.505028e-06,
                0.000000e00,
                1.713226e-06,
                0.000000e00,
                0.000000e00,
                6.273434e-18,
            ],
            14: [
                0.000000e00,
                4.368251e-08,
                2.260942e-13,
                0.000000e00,
                0.000000e00,
                -4.762111e-08,
                1.282504e-12,
                0.000000e00,
                0.000000e00,
                3.611544e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -3.197411e-09,
                8.656950e-14,
                0.000000e00,
                0.000000e00,
                3.617982e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                6.034818e-06,
                4.276140e-11,
                0.000000e00,
                4.908627e-06,
                0.000000e00,
                0.000000e00,
                -1.357859e-08,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                1.131144e-10,
                -5.279738e-14,
                0.000000e00,
                0.000000e00,
                -2.836528e-09,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -3.906973e-06,
                -1.642642e-10,
                0.000000e00,
                -6.562415e-06,
                0.000000e00,
                0.000000e00,
                -1.071227e-07,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                7.876487e-07,
                8.941839e-11,
                0.000000e00,
                2.276216e-06,
                0.000000e00,
                0.000000e00,
                1.772933e-07,
                0.000000e00,
                0.000000e00,
                0.000000e00,
                -6.455677e-08,
                0.000000e00,
                -1.423548e-05,
                0.000000e00,
                7.012498e-06,
                0.000000e00,
                0.000000e00,
                1.716278e-19,
            ],
        }
