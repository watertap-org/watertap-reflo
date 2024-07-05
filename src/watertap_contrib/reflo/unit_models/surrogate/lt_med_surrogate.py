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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Suffix,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.units.lt_med_surrogate import cost_lt_med_surrogate


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
        "property_package_liquid",
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
        "property_package_vapor",
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
    CONFIG.declare(
        "number_effects",
        ConfigValue(
            default=12,
            domain=PositiveInt,
            description="Number of effects of the LT_MED system",
            doc="""A ConfigBlock specifying the number of effects, which should be an integer between 3 to 14.""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Check if the number of effects is valid
        if self.config.number_effects not in [i for i in range(3, 15)]:
            raise ConfigurationError(
                f"The number of effects was specified as {self.config.number_effects}. The number of effects should be specified as an integer between 3 to 14."
            )

        """
        Add system configurations
        """
        self.delta_T_last_effect = Var(
            initialize=10,
            units=pyunits.K,
            doc="Temperature increase in last effect",
        )

        self.delta_T_cooling_reject = Var(
            initialize=-3,
            units=pyunits.K,
            doc="Temperature decrease in cooling reject water",
        )

        # These two variables should be fixed with the default values,
        # with which the surrogate model was developed
        self.delta_T_last_effect.fix()
        self.delta_T_cooling_reject.fix()

        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            ["Liq"],
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
        tmp_dict["parameters"] = self.config.property_package_liquid
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package_liquid.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict,
        )

        """
        Add block for distillate
        """
        tmp_dict["defined_state"] = False
        self.distillate_props = self.config.property_package_liquid.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of distillate",
            **tmp_dict,
        )

        # Distillate temperature is the same as last effect vapor temperature,
        # which is assumed to be 10 deg higher than condenser inlet seawater temperature
        @self.Constraint(doc="Distillate temperature")
        def eq_distillate_temp(b):
            return (
                b.distillate_props[0].temperature
                == b.feed_props[0].temperature + b.delta_T_last_effect
            )

        # Salinity in distillate is zero
        self.distillate_props[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)

        """
        Add block for brine
        """
        self.brine_props = self.config.property_package_liquid.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        # Relationship between brine salinity and feed salinity
        @self.Constraint(doc="Brine salinity")
        def eq_brine_salinity(b):
            return b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"] == b.feed_props[
                0
            ].conc_mass_phase_comp["Liq", "TDS"] / (1 - b.recovery_vol_phase[0, "Liq"])

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
        self.cooling_out_props = self.config.property_package_liquid.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of cooling reject",
            **tmp_dict,
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
        tmp_dict["parameters"] = self.config.property_package_vapor
        tmp_dict["defined_state"] = False

        self.steam_props = self.config.property_package_vapor.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of heating steam",
            **tmp_dict,
        )

        self.steam_props[0].flow_mass_phase_comp["Liq", "H2O"].fix(0)

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
                == b.feed_props[0].flow_vol_phase["Liq"]
                * b.recovery_vol_phase[0, "Liq"]
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
        Isobaric
        """

        @self.Constraint(self.flowsheet().config.time, doc="Isobaric")
        def eq_feed_to_distillate_isobaric(b, t):
            return b.feed_props[t].pressure == b.distillate_props[t].pressure

        @self.Constraint(self.flowsheet().config.time, doc="Isobaric")
        def eq_feed_to_brine_isobaric(b, t):
            return b.feed_props[t].pressure == b.brine_props[t].pressure

        @self.Constraint(self.flowsheet().config.time, doc="Isobaric")
        def eq_feed_to_cooling_isobaric(b, t):
            return b.feed_props[t].pressure == b.cooling_out_props[t].pressure

        """
        Add Vars for model outputs
        """
        self.thermal_power_requirement = Var(
            initialize=5000,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW)",
        )

        self.specific_area_per_m3_day = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**2 / (pyunits.m**3 / pyunits.d),
            doc="Specific area (m2/m3/day))",
        )

        self.specific_area_per_kg_s = Var(
            initialize=400,
            bounds=(0, None),
            units=pyunits.m**2 / (pyunits.k / pyunits.s),
            doc="Specific area (m2/kg/s))",
        )

        self.specific_energy_consumption_thermal = Var(
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
        # Create the total flow rate of brackish/seawater for the calculation of cooling water flow rate
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

        # Set alias for use in the surrogate equations
        self.feed_conc_ppm = pyunits.convert(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.mg / pyunits.L,
        )
        self.temp_steam = self.steam_props[0].temperature - 273.15

        # Get coefficients for surrogate equations
        self.gain_output_ratio_coeffs = self._get_gain_output_ratio_coeffs()
        self.specific_area_coeffs = self._get_specific_area_coeffs()

        # Surrogate equations were built for 3,6,9,12,14 effects
        # For intermediate number of effects (4,5,7,8,10,11,13), linear interpolation is adopted
        @self.Constraint(doc="Gain output ratio surrogate equation")
        def eq_gain_output_ratio(b):
            if b.config.number_effects in [3, 6, 9, 12, 14]:
                return b.gain_output_ratio == self._get_gain_output_ratio(
                    b.config.number_effects
                )
            else:  # b.config.number_effects in [4, 5, 7, 8, 10, 11, 13]:
                # find out the closest numbers of effects that have a surrogate equation
                interp_effects, i = [3, 6, 9, 12, 14], 1
                while interp_effects[i] < b.config.number_effects:
                    i += 1
                interp_point1, interp_point2 = interp_effects[i - 1], interp_effects[i]

                # implement linear interpolation using 2 points
                return b.gain_output_ratio * (
                    interp_point2 - interp_point1
                ) == self._get_gain_output_ratio(interp_point1) * (
                    interp_point2 - b.config.number_effects
                ) + self._get_gain_output_ratio(
                    interp_point2
                ) * (
                    b.config.number_effects - interp_point1
                )

        @self.Constraint(doc="Specific area surrogate equation")
        def eq_specific_area_per_m3_day(b):
            if b.config.number_effects in [3, 6, 9, 12, 14]:
                return b.specific_area_per_m3_day == self._get_specific_area(
                    b.config.number_effects
                )
            else:  # b.config.number_effects in [4, 5, 7, 8, 10, 11, 13]:
                # find out the closest numbers of effects that have a surrogate equation
                interp_effects, i = [3, 6, 9, 12, 14], 1
                while interp_effects[i] < b.config.number_effects:
                    i += 1
                interp_point1, interp_point2 = interp_effects[i - 1], interp_effects[i]

                # implement linear interpolation using 2 points
                return b.specific_area_per_m3_day * (
                    interp_point2 - interp_point1
                ) == self._get_specific_area(interp_point1) * (
                    interp_point2 - b.config.number_effects
                ) + self._get_specific_area(
                    interp_point2
                ) * (
                    b.config.number_effects - interp_point1
                )

        @self.Constraint(doc="Convert specific area to m2/kg/s for CAPEX calculation")
        def eq_specific_area_kg_s(b):
            return b.specific_area_per_kg_s == pyunits.convert(
                b.specific_area_per_m3_day / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**2 / pyunits.kg * pyunits.s,
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
        @self.Constraint(doc="specific_energy_consumption_thermal calculation")
        def eq_specific_thermal_energy_consumption(b):
            return b.specific_energy_consumption_thermal == pyunits.convert(
                1
                / b.gain_output_ratio
                * (b.steam_props[0].dh_vap_mass)
                * b.distillate_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(doc="Thermal power requirement calculation")
        def eq_thermal_power_requirement(b):
            return b.thermal_power_requirement == pyunits.convert(
                b.specific_energy_consumption_thermal
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

    def initialize_build(
        blk,
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
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)
        # ---------------------------------------------------------------------
        flags = blk.feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state

        if state_args is None:
            blk.state_args = state_args = {}
            state_dict = blk.feed_props[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args_dist = deepcopy(state_args)
        for p, j in blk.distillate_props.phase_component_set:
            if j == "TDS":
                state_args_dist["flow_mass_phase_comp"][(p, j)] = 0

        blk.distillate_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_dist,
        )

        state_args_steam = {}
        state_dict_steam = blk.steam_props[
            blk.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_steam.keys():
            if state_dict_steam[k].is_indexed():
                state_args_steam[k] = {}
                for m in state_dict_steam[k].keys():
                    state_args_steam[k][m] = state_dict_steam[k][m].value
            else:
                state_args_steam[k] = state_dict_steam[k].value

        for p, j in blk.steam_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_steam["flow_mass_phase_comp"][(p, j)] = 0

        blk.steam_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_steam,
        )

        state_args_brine = deepcopy(state_args)
        for p, j in blk.brine_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_brine["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, j)]
                    * (1 - blk.recovery_vol_phase[0, "Liq"])
                    * pyunits.kg
                    / pyunits.s
                )

        blk.brine_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_brine,
        )

        blk.cooling_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        if degrees_of_freedom(blk) != 0:
            raise InitializationError(
                f"{blk.name} degrees of freedom were not 0 at the beginning of initialization. DoF = {degrees_of_freedom(blk)}."
            )

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_props.release_state(flags, outlvl=outlvl)

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        dist_vol_flow = self.distillate_props[0].flow_vol_phase["Liq"]

        if iscale.get_scaling_factor(self.recovery_vol_phase[0, "Liq"]) is None:
            iscale.set_scaling_factor(self.recovery_vol_phase[0, "Liq"], 1e1)

        if iscale.get_scaling_factor(self.specific_energy_consumption_thermal) is None:
            iscale.set_scaling_factor(self.specific_energy_consumption_thermal, 1e-3)

        if iscale.get_scaling_factor(self.thermal_power_requirement) is None:
            iscale.set_scaling_factor(
                self.thermal_power_requirement,
                iscale.get_scaling_factor(dist_vol_flow)
                * iscale.get_scaling_factor(self.specific_energy_consumption_thermal),
            )

        if iscale.get_scaling_factor(self.specific_area_per_m3_day) is None:
            iscale.set_scaling_factor(self.specific_area_per_m3_day, 0.1)

        if iscale.get_scaling_factor(self.specific_area_per_kg_s) is None:
            iscale.set_scaling_factor(self.specific_area_per_kg_s, 1e-2)

        if iscale.get_scaling_factor(self.gain_output_ratio) is None:
            iscale.set_scaling_factor(self.gain_output_ratio, 0.1)

        if iscale.get_scaling_factor(self.feed_cool_mass_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_mass_flow, 1e-3)

        if iscale.get_scaling_factor(self.feed_cool_vol_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_vol_flow, 1e-3)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Water Inlet": self.feed,
                "Distillate Outlet": self.distillate,
                "Brine Outlet": self.brine,
                "Heating Steam Inlet": self.steam_props,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Gained output ratio"] = self.gain_output_ratio
        var_dict["Thermal power requirement (kW)"] = self.thermal_power_requirement
        var_dict[
            "Specific thermal energy consumption (kWh/m3)"
        ] = self.specific_energy_consumption_thermal
        var_dict["Feed water volumetric flow rate"] = self.feed_props[0].flow_vol_phase[
            "Liq"
        ]
        var_dict["Cooling water volumetric flow rate"] = self.cooling_out_props[
            0
        ].flow_vol_phase["Liq"]
        var_dict["Heating steam mass flow rate"] = self.steam_props[
            0
        ].flow_mass_phase_comp["Vap", "H2O"]
        var_dict["Specific area (m2/m3/day)"] = self.specific_area_per_m3_day

        return {"vars": var_dict}

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
                0.00025156,
                -5.32e-09,
                0,
                -2.27e-05,
                1.13e-10,
                0,
                2.01e-05,
                0,
                1.71e-07,
                7.04e-13,
                0,
                2.01e-07,
                0,
                0.0131458,
                6.57e-09,
                0.00040718,
                0,
                -4.99e-07,
                -5.30e-13,
                0,
                -7.50e-08,
                0,
                -0.01124119,
                -9.02e-09,
                -0.00037472,
                -0.00010574,
                0.00182065,
                5.45e-09,
                5.19e-05,
                0.00011658,
                0,
                -2.24e-05,
                -4.41e-05,
                0,
                3.93e-14,
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

    # Surrogate equations for calculating gain output ratio
    def _get_gain_output_ratio(self, num_effect):
        return (
            self.feed_conc_ppm * self.gain_output_ratio_coeffs[num_effect][0]
            + self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][1]
            + self.feed_conc_ppm
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][2]
            + self.temp_last_effect * self.gain_output_ratio_coeffs[num_effect][3]
            + self.temp_last_effect
            * self.feed_conc_ppm
            * self.gain_output_ratio_coeffs[num_effect][4]
            + self.temp_last_effect
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][5]
            + self.temp_steam * self.gain_output_ratio_coeffs[num_effect][6]
            + self.temp_steam
            * self.feed_conc_ppm
            * self.gain_output_ratio_coeffs[num_effect][7]
            + self.temp_steam
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][8]
            + self.temp_steam
            * self.temp_last_effect
            * self.gain_output_ratio_coeffs[num_effect][9]
            + 1 * self.gain_output_ratio_coeffs[num_effect][10]
            + self.temp_steam**2 * self.gain_output_ratio_coeffs[num_effect][11]
            + self.temp_last_effect**2 * self.gain_output_ratio_coeffs[num_effect][12]
            + self.recovery_vol_phase[0, "Liq"] ** 2
            * self.gain_output_ratio_coeffs[num_effect][13]
            + self.feed_conc_ppm**2 * self.gain_output_ratio_coeffs[num_effect][14]
        )

    # Specific area surrogate equation, as a function of the number of effect
    def _get_specific_area(self, num_effect):
        if num_effect in [3, 6, 9]:
            return (
                self.feed_conc_ppm * self.specific_area_coeffs[num_effect][0]
                + self.feed_conc_ppm**2 * self.specific_area_coeffs[num_effect][1]
                + self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][2]
                + self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][3]
                + self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][4]
                + self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][5]
                + self.recovery_vol_phase[0, "Liq"] ** 2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][6]
                + self.temp_last_effect * self.specific_area_coeffs[num_effect][7]
                + self.temp_last_effect
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][8]
                + self.temp_last_effect
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][9]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][10]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][11]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][12]
                + self.temp_last_effect**2 * self.specific_area_coeffs[num_effect][13]
                + self.temp_last_effect**2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][14]
                + self.temp_last_effect**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][15]
                + self.temp_steam * self.specific_area_coeffs[num_effect][16]
                + self.temp_steam
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][17]
                + self.temp_steam
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][18]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][19]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][20]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][21]
                + self.temp_steam
                * self.temp_last_effect
                * self.specific_area_coeffs[num_effect][22]
                + self.temp_steam
                * self.temp_last_effect
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][23]
                + self.temp_steam
                * self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][24]
                + self.temp_steam
                * self.temp_last_effect**2
                * self.specific_area_coeffs[num_effect][25]
                + self.temp_steam**2 * self.specific_area_coeffs[num_effect][26]
                + self.temp_steam**2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][27]
                + self.temp_steam**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][28]
                + self.temp_steam**2
                * self.temp_last_effect
                * self.specific_area_coeffs[num_effect][29]
                + 1 * self.specific_area_coeffs[num_effect][30]
                + self.temp_steam**3 * self.specific_area_coeffs[num_effect][31]
                + self.temp_last_effect**3 * self.specific_area_coeffs[num_effect][32]
                + self.recovery_vol_phase[0, "Liq"] ** 3
                * self.specific_area_coeffs[num_effect][33]
                + self.feed_conc_ppm**3 * self.specific_area_coeffs[num_effect][34]
            )

        else:
            return (
                self.feed_conc_ppm * self.specific_area_coeffs[num_effect][0]
                + self.feed_conc_ppm**2 * self.specific_area_coeffs[num_effect][1]
                + self.feed_conc_ppm**3 * self.specific_area_coeffs[num_effect][2]
                + self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][3]
                + self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][4]
                + self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][5]
                + self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm**3
                * self.specific_area_coeffs[num_effect][6]
                + self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][7]
                + self.recovery_vol_phase[0, "Liq"] ** 2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][8]
                + self.recovery_vol_phase[0, "Liq"] ** 2
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][9]
                + self.recovery_vol_phase[0, "Liq"] ** 3
                * self.specific_area_coeffs[num_effect][10]
                + self.recovery_vol_phase[0, "Liq"] ** 3
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][11]
                + self.temp_last_effect * self.specific_area_coeffs[num_effect][12]
                + self.temp_last_effect
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][13]
                + self.temp_last_effect
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][14]
                + self.temp_last_effect
                * self.feed_conc_ppm**3
                * self.specific_area_coeffs[num_effect][15]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][16]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][17]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][18]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][19]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][20]
                + self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"] ** 3
                * self.specific_area_coeffs[num_effect][21]
                + self.temp_last_effect**2 * self.specific_area_coeffs[num_effect][22]
                + self.temp_last_effect**2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][23]
                + self.temp_last_effect**2
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][24]
                + self.temp_last_effect**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][25]
                + self.temp_last_effect**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][26]
                + self.temp_last_effect**2
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][27]
                + self.temp_last_effect**3 * self.specific_area_coeffs[num_effect][28]
                + self.temp_last_effect**3
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][29]
                + self.temp_last_effect**3
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][30]
                + self.temp_steam * self.specific_area_coeffs[num_effect][31]
                + self.temp_steam
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][32]
                + self.temp_steam
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][33]
                + self.temp_steam
                * self.feed_conc_ppm**3
                * self.specific_area_coeffs[num_effect][34]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][35]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][36]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][37]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][38]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][39]
                + self.temp_steam
                * self.recovery_vol_phase[0, "Liq"] ** 3
                * self.specific_area_coeffs[num_effect][40]
                + self.temp_steam
                * self.temp_last_effect
                * self.specific_area_coeffs[num_effect][41]
                + self.temp_steam
                * self.temp_last_effect
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][42]
                + self.temp_steam
                * self.temp_last_effect
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][43]
                + self.temp_steam
                * self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][44]
                + self.temp_steam
                * self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][45]
                + self.temp_steam
                * self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][46]
                + self.temp_steam
                * self.temp_last_effect**2
                * self.specific_area_coeffs[num_effect][47]
                + self.temp_steam
                * self.temp_last_effect**2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][48]
                + self.temp_steam
                * self.temp_last_effect**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][49]
                + self.temp_steam
                * self.temp_last_effect**3
                * self.specific_area_coeffs[num_effect][50]
                + self.temp_steam**2 * self.specific_area_coeffs[num_effect][51]
                + self.temp_steam**2
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][52]
                + self.temp_steam**2
                * self.feed_conc_ppm**2
                * self.specific_area_coeffs[num_effect][53]
                + self.temp_steam**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][54]
                + self.temp_steam**2
                * self.recovery_vol_phase[0, "Liq"]
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][55]
                + self.temp_steam**2
                * self.recovery_vol_phase[0, "Liq"] ** 2
                * self.specific_area_coeffs[num_effect][56]
                + self.temp_steam**2
                * self.temp_last_effect
                * self.specific_area_coeffs[num_effect][57]
                + self.temp_steam**2
                * self.temp_last_effect
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][58]
                + self.temp_steam**2
                * self.temp_last_effect
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][59]
                + self.temp_steam**2
                * self.temp_last_effect**2
                * self.specific_area_coeffs[num_effect][60]
                + self.temp_steam**3 * self.specific_area_coeffs[num_effect][61]
                + self.temp_steam**3
                * self.feed_conc_ppm
                * self.specific_area_coeffs[num_effect][62]
                + self.temp_steam**3
                * self.recovery_vol_phase[0, "Liq"]
                * self.specific_area_coeffs[num_effect][63]
                + self.temp_steam**3
                * self.temp_last_effect
                * self.specific_area_coeffs[num_effect][64]
                + 1 * self.specific_area_coeffs[num_effect][65]
                + self.temp_steam**4 * self.specific_area_coeffs[num_effect][66]
                + self.temp_last_effect**4 * self.specific_area_coeffs[num_effect][67]
                + self.recovery_vol_phase[0, "Liq"] ** 4
                * self.specific_area_coeffs[num_effect][68]
                + self.feed_conc_ppm**4 * self.specific_area_coeffs[num_effect][69]
            )

    @property
    def default_costing_method(self):
        return cost_lt_med_surrogate
