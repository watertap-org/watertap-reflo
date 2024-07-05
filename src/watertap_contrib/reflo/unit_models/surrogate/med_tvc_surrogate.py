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
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.units.med_tvc_surrogate import (
    cost_med_tvc_surrogate,
)

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


@declare_process_block_class("MEDTVCSurrogate")
class MEDTVCData(UnitModelBlockData):
    """
    Multi-effect distillation with thermal vapor compressor model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
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
            description="Property package to use for feed, distillate, brine and cooling water properties",
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
            description="Property package to use for heating and motive steam properties",
            doc="""Property parameter object used to define steasm property calculations,
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
            description="Number of effects of the MED_TVC system",
            doc="""A ConfigBlock specifying the number of effects, which should be an integer between 8 to 16.""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = (
            self.config.property_package_liquid.get_metadata().get_derived_units
        )

        # Check if the number of effects is valid
        if self.config.number_effects not in [i for i in range(8, 17)]:
            raise ConfigurationError(
                f"The number of effects was specified as {self.config.number_effects}. The number of effects should be specified as an integer between 8 to 16."
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
            initialize=0.30,
            bounds=(0.30, 0.40),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        self.thermal_loss = Param(
            initialize=0.054, units=pyunits.dimensionless, doc="System thermal loss"
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

        # distillate temperature is the same as last effect vapor temperature, which is 10 deg higher than condenser inlet seawater temperature
        @self.Constraint(doc="distillate temperature")
        def eq_distillate_temp(b):
            return (
                b.distillate_props[0].temperature
                == b.feed_props[0].temperature + b.delta_T_last_effect
            )

        # salinity in distillate is zero
        self.distillate_props[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)

        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

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

        # Assumption: the temperature of cooling reject is 3 degC lower than in the condenser
        @self.Constraint(doc="Cooling reject temperature")
        def eq_cooling_temp(b):
            return (
                b.cooling_out_props[0].temperature
                == b.distillate_props[0].temperature + b.delta_T_cooling_reject
            )

        @self.Constraint(doc="Cooling reject salinity")
        def eq_cooling_salinity(b):
            return (
                b.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
            )

        """
        Add block for heating steam
        """
        tmp_dict["parameters"] = self.config.property_package_vapor
        tmp_dict["defined_state"] = False

        self.heating_steam_props = self.config.property_package_vapor.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of heating steam",
            **tmp_dict,
        )

        @self.Constraint(doc="Flow rate of liquid heating steam is zero")
        def eq_heating_steam_liquid_mass(b):
            return b.heating_steam_props[0].flow_mass_phase_comp["Liq", "H2O"] == 0

        """
        Add block for motive steam
        """
        tmp_dict["parameters"] = self.config.property_package_vapor
        tmp_dict["defined_state"] = False

        self.motive_steam_props = self.config.property_package_vapor.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of motive steam",
            **tmp_dict,
        )

        @self.Constraint(doc="Flow rate of liquid motive steam is zero")
        def eq_motive_steam_liquid_mass(b):
            return b.motive_steam_props[0].flow_mass_phase_comp["Liq", "H2O"] == 0

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="distillate", block=self.distillate_props)
        self.add_port(name="brine", block=self.brine_props)
        self.add_port(name="steam", block=self.heating_steam_props)
        self.add_port(name="motive", block=self.motive_steam_props)

        """
        Mass balances
        """
        # Distillate flow rate calculation
        @self.Constraint(doc="Distallate volumetric flow rate")
        def eq_dist_vol_flow(b):
            return (
                b.distillate_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                * b.recovery_vol_phase[0, "Liq"]
            )

        # Brine flow rate calculation
        @self.Constraint(doc="Brine volume flow rate")
        def eq_brine_vol_flow(b):
            return (
                b.brine_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                - b.distillate_props[0].flow_vol_phase["Liq"]
            )

        """
        Isobaric
        """

        @self.Constraint(doc="Isobaric")
        def eq_feed_to_distillate_isobaric(b):
            return b.feed_props[0].pressure == b.distillate_props[0].pressure

        @self.Constraint(doc="Isobaric")
        def eq_feed_to_brine_isobaric(b):
            return b.feed_props[0].pressure == b.brine_props[0].pressure

        @self.Constraint(doc="Isobaric")
        def eq_feed_to_cooling_isobaric(b):
            return b.feed_props[0].pressure == b.cooling_out_props[0].pressure

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
            doc="Feed and cooling water mass flow rate (kg/s)",
        )

        self.feed_cool_vol_flow = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.h,
            doc="Feed and cooling water volume flow rate (m3/h)",
        )

        self.feed_conc_ppm = pyunits.convert(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.mg / pyunits.L,
        )

        self.capacity = pyunits.convert(
            self.feed_props[0].flow_vol_phase["Liq"]
            * self.recovery_vol_phase[0, "Liq"],
            to_units=pyunits.m**3 / pyunits.day,
        )

        # Set alias for use in the surrogate equations
        self.feed_temperature = self.feed_props[0].temperature - 273.15
        self.motive_pressure = pyunits.convert(
            self.motive_steam_props[0].pressure, to_units=pyunits.bar
        )

        # Get coefficients for surrogate equations
        self.gain_output_ratio_coeffs = self._get_gain_output_ratio_coeffs()
        self.specific_area_coeffs = self._get_specific_area_coeffs()
        self.heating_steam_mass_flow_rate_coeffs = (
            self._get_heating_steam_mass_flow_rate_coeffs()
        )
        self.motive_steam_mass_flow_rate_coeffs = (
            self._get_motive_steam_mass_flow_rate_coeffs()
        )

        # Surrogate equations were built for 8,10,12,14,16 effects
        # For intermediate number of effects (9,11,13,15), linear interpolation is adopted
        @self.Constraint(doc="Gain output ratio surrogate equation")
        def eq_gain_output_ratio(b):
            if b.config.number_effects in [8, 10, 12, 14, 16]:
                return b.gain_output_ratio == self._get_gain_output_ratio(
                    b.config.number_effects
                )
            else:  # b.config.number_effects in [9, 11, 13, 15]:
                return (
                    b.gain_output_ratio
                    == (
                        self._get_gain_output_ratio(b.config.number_effects - 1)
                        + self._get_gain_output_ratio(b.config.number_effects + 1)
                    )
                    / 2
                )

        @self.Constraint(doc="specific area surrogate equation")
        def eq_specific_area(b):
            if b.config.number_effects in [8, 10, 12, 14, 16]:
                return b.specific_area_per_m3_day == self._get_specific_area(
                    b.config.number_effects
                )
            else:  # b.config.number_effects in [9, 11, 13, 15]:
                return (
                    b.specific_area_per_m3_day
                    == (
                        self._get_specific_area(b.config.number_effects - 1)
                        + self._get_specific_area(b.config.number_effects + 1)
                    )
                    / 2
                )

        @self.Constraint(doc="Convert specific area to m2/kg/s for CAPEX calculation")
        def eq_specific_area_kg_s(b):
            return b.specific_area_per_kg_s == pyunits.convert(
                b.specific_area_per_m3_day / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**2 / pyunits.kg * pyunits.s,
            )

        @self.Constraint(doc="heating steam mass flow rate surrogate equation")
        def eq_heating_steam_mass_flow_rate(b):
            if b.config.number_effects in [8, 10, 12, 14, 16]:
                return b.heating_steam_props[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ] == self._get_heating_steam_mass_flow_rate(b.config.number_effects)
            else:  # b.config.number_effects in [9, 11, 13, 15]:
                return (
                    b.heating_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
                    == (
                        self._get_heating_steam_mass_flow_rate(
                            b.config.number_effects - 1
                        )
                        + self._get_heating_steam_mass_flow_rate(
                            b.config.number_effects + 1
                        )
                    )
                    / 2
                )

        @self.Constraint(doc="motive steam mass flow rate surrogate equation")
        def eq_motive_steam_mass_flow_rate(b):
            if b.config.number_effects in [8, 10, 12, 14, 16]:
                return b.motive_steam_props[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ] == self._get_motive_steam_mass_flow_rate(b.config.number_effects)
            else:  # b.config.number_effects in [9, 11, 13, 15]:
                return (
                    b.motive_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
                    == (
                        self._get_motive_steam_mass_flow_rate(
                            b.config.number_effects - 1
                        )
                        + self._get_motive_steam_mass_flow_rate(
                            b.config.number_effects + 1
                        )
                    )
                    / 2
                )

        # Energy consumption
        @self.Constraint(doc="Specific thermal energy consumption calculation")
        def eq_specific_energy_consumption_thermal(b):
            return b.specific_energy_consumption_thermal == pyunits.convert(
                1
                / b.gain_output_ratio
                * (
                    b.motive_steam_props[0].enth_mass_phase["Vap"]
                    - b.heating_steam_props[0].enth_mass_phase["Liq"]
                )
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

        # Enthalpy balance (calculate feed and cooling water mass flow rate (kg/s))
        @self.Constraint(doc="System overall enthalpy balance")
        def eq_system_enthalpy_balance(b):
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
        @self.Constraint(doc="Feed and cooling water volumetric flow rate (m3/h)")
        def eq_feed_cool_vol_flow(b):
            return b.feed_cool_vol_flow == pyunits.convert(
                b.feed_cool_mass_flow / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hour,
            )

        @self.Constraint(doc="Cooling water volumetric flow rate (m3/h)")
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
        state_dict_steam = blk.heating_steam_props[
            blk.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_steam.keys():
            if state_dict_steam[k].is_indexed():
                state_args_steam[k] = {}
                for m in state_dict_steam[k].keys():
                    state_args_steam[k][m] = state_dict_steam[k][m].value
            else:
                state_args_steam[k] = state_dict_steam[k].value

        for p, j in blk.heating_steam_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_steam["flow_mass_phase_comp"][(p, j)] = 0

        blk.heating_steam_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_steam,
        )

        blk.motive_steam_props.initialize(
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
        # Check degree of freedom
        assert degrees_of_freedom(blk) == 0

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_props.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

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
            iscale.set_scaling_factor(self.specific_area_per_m3_day, 1e-1)

        if iscale.get_scaling_factor(self.specific_area_per_kg_s) is None:
            iscale.set_scaling_factor(self.specific_area_per_kg_s, 1e-2)

        if iscale.get_scaling_factor(self.gain_output_ratio) is None:
            iscale.set_scaling_factor(self.gain_output_ratio, 1e-1)

        if iscale.get_scaling_factor(self.feed_cool_mass_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_mass_flow, 1e-3)

        if iscale.get_scaling_factor(self.feed_cool_vol_flow) is None:
            iscale.set_scaling_factor(self.feed_cool_vol_flow, 1e-3)

        # Transforming constraints
        sf = iscale.get_scaling_factor(self.distillate_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_distillate_temp, sf)

        sf = iscale.get_scaling_factor(
            self.heating_steam_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_heating_steam_liquid_mass, sf)

        sf = iscale.get_scaling_factor(
            self.motive_steam_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_motive_steam_liquid_mass, sf)

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

        sf = iscale.get_scaling_factor(self.feed_props[0].pressure)
        iscale.constraint_scaling_transform(self.eq_feed_to_distillate_isobaric, sf)
        iscale.constraint_scaling_transform(self.eq_feed_to_brine_isobaric, sf)
        iscale.constraint_scaling_transform(self.eq_feed_to_cooling_isobaric, sf)

        sf = iscale.get_scaling_factor(self.feed_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_dist_vol_flow, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_brine_vol_flow, sf)

        sf = iscale.get_scaling_factor(self.gain_output_ratio)
        iscale.constraint_scaling_transform(self.eq_gain_output_ratio, sf)

        sf = iscale.get_scaling_factor(self.specific_area_per_m3_day)
        iscale.constraint_scaling_transform(self.eq_specific_area, sf)

        sf = iscale.get_scaling_factor(self.specific_area_per_kg_s)
        iscale.constraint_scaling_transform(self.eq_specific_area_kg_s, sf)

        sf = iscale.get_scaling_factor(
            self.heating_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_heating_steam_mass_flow_rate, sf)

        sf = iscale.get_scaling_factor(
            self.motive_steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_motive_steam_mass_flow_rate, sf)

        sf = iscale.get_scaling_factor(self.specific_energy_consumption_thermal)
        iscale.constraint_scaling_transform(
            self.eq_specific_energy_consumption_thermal, sf
        )

        sf = iscale.get_scaling_factor(self.thermal_power_requirement)
        iscale.constraint_scaling_transform(self.eq_thermal_power_requirement, sf)

        sf = (
            iscale.get_scaling_factor(self.feed_cool_mass_flow)
            * iscale.get_scaling_factor(
                self.cooling_out_props[0].enth_mass_phase["Liq"]
            )
            * 1e3
        )
        iscale.constraint_scaling_transform(self.eq_system_enthalpy_balance, sf)

        sf = iscale.get_scaling_factor(self.feed_cool_vol_flow)
        iscale.constraint_scaling_transform(self.eq_feed_cool_vol_flow, sf)

        sf = (
            iscale.get_scaling_factor(self.cooling_out_props[0].flow_vol_phase["Liq"])
            / 3600
        )
        iscale.constraint_scaling_transform(self.eq_cool_vol_flow, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Water Inlet": self.feed,
                "Distillate Outlet": self.distillate,
                "Brine Outlet": self.brine,
                "Heating Steam Inlet": self.heating,
                "Motive Steam Inlet": self.motive,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Gained output ratio"] = self.gain_output_ratio
        var_dict["Thermal power reqruiement (kW)"] = self.thermal_power_requirement
        var_dict[
            "Specific thermal energy consumption (kWh/m3)"
        ] = self.specific_energy_consumption_thermal
        var_dict["Feed water volumetric flow rate"] = self.feed_props[0].flow_vol_phase[
            "Liq"
        ]
        var_dict["Cooling water volumetric flow rate"] = self.cooling_out_props[
            0
        ].flow_vol_phase["Liq"]
        var_dict["Heating steam mass flow rate"] = self.heating_steam_props[
            0
        ].flow_mass_phase_comp["Vap", "H2O"]
        var_dict["Motive steam mass flow rate"] = self.motive_steam_props[
            0
        ].flow_mass_phase_comp["Vap", "H2O"]
        var_dict["Specific area (m2/m3/day)"] = self.specific_area_per_m3_day

        return {"vars": var_dict}

    def _get_gain_output_ratio_coeffs(self):
        return {
            8: [
                -1.42e-06,
                7.764150858,
                -1.15e-05,
                0.582960027,
                1.24e-07,
                0.056555556,
                -2.57e-19,
                3.63e-25,
                3.23e-19,
                3.04e-21,
                0.104803303,
                -2.08e-09,
                0.012249322,
                -0.000299729,
                6.28e-22,
                -5.098563697,
                -0.001263974,
                -3.56e-25,
                -0.005858519,
                -8.896296296,
                -1.98e-11,
            ],
            10: [
                -2.28e-06,
                10.98911472,
                -1.31e-05,
                0.689378726,
                1.47e-07,
                0.080111111,
                -2.98e-19,
                -6.15e-25,
                4.06e-19,
                4.44e-21,
                0.124314705,
                -1.81e-10,
                0.018184282,
                -0.000334959,
                6.98e-22,
                -6.475922544,
                -0.001524489,
                -5.72e-25,
                -0.006804444,
                -12.28888889,
                -1.68e-11,
            ],
            12: [
                -3.05e-06,
                14.55682317,
                -1.30e-05,
                0.425476523,
                1.61e-07,
                0.097638889,
                6.39e-06,
                2.07e-13,
                -2.09e-06,
                -3.91e-08,
                0.143261656,
                6.78e-10,
                0.027642276,
                -0.00025,
                4.65e-08,
                -2.620877082,
                -0.001921861,
                -5.22e-11,
                -0.001659259,
                -15.9537037,
                -1.83e-11,
            ],
            14: [
                -5.77e-06,
                18.58044758,
                -9.55e-06,
                0.758813856,
                1.78e-07,
                0.148784722,
                2.27e-06,
                -4.93e-13,
                -3.10e-06,
                -3.26e-08,
                0.157319103,
                4.40e-09,
                0.033536585,
                -0.000335027,
                1.13e-08,
                -7.631229505,
                -0.002018089,
                -1.01e-11,
                -0.006560185,
                -21.38310185,
                -7.92e-12,
            ],
            16: [
                -8.38e-06,
                20.76389007,
                -3.04e-06,
                0.967432521,
                1.85e-07,
                0.216269841,
                1.58e-05,
                -6.41e-13,
                -2.27e-06,
                -4.83e-07,
                0.185976144,
                2.71e-09,
                0.046360821,
                -0.000913957,
                -4.66e-08,
                -10.82641732,
                -0.002230429,
                -1.06e-11,
                -0.009318519,
                -25.09448224,
                -4.22e-12,
            ],
        }

    def _get_specific_area_coeffs(self):
        return {
            8: [
                -0.000113631,
                -13.11517461,
                0.000106181,
                -0.494480795,
                2.52e-06,
                0.214458333,
                9.34e-06,
                2.00e-11,
                1.72e-06,
                1.04e-07,
                -0.045572653,
                -6.57e-09,
                -0.000291328,
                6.29e-05,
                -7.80e-13,
                12.94569696,
                0.000878721,
                -1.10e-10,
                0.007633796,
                6.371296296,
                3.41e-10,
            ],
            10: [
                -0.000275164,
                -30.34222168,
                0.000245662,
                -0.964861655,
                5.70e-06,
                0.49675,
                -7.92e-07,
                1.92e-12,
                5.90e-07,
                5.40e-09,
                0.001209732,
                -2.56e-08,
                -0.005934959,
                -1.36e-06,
                1.14e-09,
                25.22049669,
                1.96e-05,
                3.20e-12,
                0.013141389,
                14.14444444,
                8.80e-10,
            ],
            12: [
                -0.00059378,
                -70.68831403,
                0.000537912,
                -2.002884776,
                1.10e-05,
                1.071916667,
                2.61e-05,
                -7.91e-11,
                -8.33e-06,
                -4.78e-07,
                -0.020075596,
                9.86e-08,
                0.014342818,
                0.000726897,
                -3.18e-08,
                52.69512155,
                -0.000184651,
                -6.33e-11,
                0.024886019,
                34.90462963,
                2.02e-09,
            ],
            14: [
                -0.00121164,
                -152.5537433,
                0.001101146,
                -3.461298325,
                2.14e-05,
                2.249461806,
                4.67e-05,
                -1.99e-10,
                -2.11e-05,
                -9.28e-07,
                -0.054534982,
                3.70e-07,
                0.068449356,
                0.001176152,
                -3.36e-08,
                97.79031075,
                -0.000342235,
                -6.33e-11,
                0.038331481,
                78.17853009,
                4.22e-09,
            ],
            16: [
                -0.003702813,
                -566.9002008,
                0.003490278,
                -10.24802298,
                5.87e-05,
                8.029940476,
                0.000180849,
                -1.67e-09,
                -0.000223884,
                -5.19e-06,
                -0.699017368,
                3.49e-06,
                0.438254936,
                0.010243428,
                -2.60e-06,
                312.1661769,
                0.00564047,
                1.39e-09,
                0.098896296,
                317.9402872,
                1.27e-08,
            ],
        }

    def _get_heating_steam_mass_flow_rate_coeffs(self):
        return {
            8: [
                -7.59e-06,
                -34.3908768,
                9.69e-06,
                -0.131363139,
                1.48e-07,
                0.200569444,
                0.002087075,
                -1.31e-10,
                -0.000541286,
                -5.70e-06,
                -0.109758232,
                -3.05e-09,
                0.000630081,
                0.000228591,
                -9.38e-08,
                8.029964774,
                0.002098453,
                -2.59e-10,
                0.000841111,
                40.08888889,
                -1.42e-12,
            ],
            10: [
                1.59e-05,
                -29.09136792,
                -2.04e-05,
                -0.066787419,
                -7.21e-08,
                0.145083333,
                0.001764972,
                -1.33e-10,
                -0.000559573,
                -5.46e-06,
                0.005482395,
                -6.00e-08,
                -0.017926829,
                8.64e-06,
                -7.40e-08,
                5.722781223,
                5.07e-05,
                8.14e-12,
                0.00026963,
                36.86296296,
                -6.77e-11,
            ],
            12: [
                -4.48e-06,
                -38.28037765,
                1.02e-07,
                -0.413517098,
                1.86e-07,
                0.304472222,
                0.001604017,
                -1.50e-10,
                -0.00055003,
                -6.73e-06,
                -0.008556022,
                -3.24e-09,
                -0.002256098,
                0.000298408,
                -9.02e-08,
                12.82816936,
                3.81e-05,
                -3.04e-11,
                0.005071204,
                41.34814815,
                -1.16e-11,
            ],
            14: [
                -7.83e-07,
                -40.03225543,
                -7.58e-06,
                -0.359841547,
                1.82e-07,
                0.30859375,
                0.00148361,
                -1.61e-10,
                -0.000593416,
                -6.66e-06,
                -0.001763606,
                -3.31e-09,
                0.001147527,
                0.00014878,
                -6.89e-08,
                11.96567201,
                -5.86e-05,
                -7.44e-11,
                0.004328148,
                45.25318287,
                -2.16e-11,
            ],
            16: [
                3.35e-06,
                -38.68477421,
                -1.90e-05,
                -0.86243177,
                2.00e-07,
                0.233531746,
                0.00126754,
                -1.68e-10,
                -0.000627918,
                -4.45e-06,
                -0.077815222,
                -4.48e-09,
                -0.046564073,
                0.001702812,
                -8.14e-07,
                20.5725045,
                0.002389423,
                5.92e-10,
                0.012013148,
                49.04383976,
                -3.29e-11,
            ],
        }

    def _get_motive_steam_mass_flow_rate_coeffs(self):
        return {
            8: [
                1.49e-05,
                -29.83707466,
                2.86e-05,
                -2.345398834,
                -9.10e-07,
                0.267902778,
                0.002171651,
                4.00e-10,
                -0.000324519,
                -2.78e-05,
                -0.489303064,
                -7.25e-08,
                0.033783875,
                0.006929505,
                -4.33e-06,
                43.26022001,
                0.005568346,
                -1.65e-10,
                0.035492222,
                28.25555556,
                4.56e-11,
            ],
            10: [
                2.43e-05,
                -27.0360671,
                8.10e-06,
                -1.966013646,
                -8.70e-07,
                0.257819444,
                0.001822688,
                3.13e-10,
                -0.000341778,
                -2.35e-05,
                -0.349886659,
                -9.33e-08,
                0.025067751,
                0.00570271,
                -3.59e-06,
                36.5052648,
                0.003546979,
                4.13e-12,
                0.029556389,
                25.92222222,
                7.53e-12,
            ],
            12: [
                1.12e-05,
                -32.60635053,
                1.54e-05,
                -1.366729008,
                -6.01e-07,
                0.362472222,
                0.001598178,
                2.60e-10,
                -0.000344721,
                -2.04e-05,
                -0.312617943,
                -5.56e-08,
                0.035501355,
                0.004905556,
                -3.20e-06,
                28.5739492,
                0.003185042,
                2.08e-11,
                0.019064444,
                28.63333333,
                3.01e-11,
            ],
            14: [
                1.20e-05,
                -33.97061713,
                7.27e-06,
                -1.345914579,
                -4.99e-07,
                0.39046875,
                0.001457515,
                2.05e-10,
                -0.000376928,
                -1.86e-05,
                -0.283422499,
                -4.47e-08,
                0.04076897,
                0.004456843,
                -2.85e-06,
                28.05931057,
                0.002831752,
                2.13e-11,
                0.018738704,
                30.55700231,
                1.65e-11,
            ],
            16: [
                1.33e-05,
                -33.90192114,
                -1.17e-06,
                -1.544393645,
                -4.29e-07,
                0.364702381,
                0.001274023,
                1.65e-10,
                -0.000401789,
                -1.57e-05,
                -0.282901411,
                -2.79e-08,
                0.014909988,
                0.004355014,
                -3.09e-06,
                31.51108079,
                0.004025387,
                4.05e-10,
                0.021729259,
                33.11413454,
                8.59e-12,
            ],
        }

    # Gain output ratio surrogate equation, as a function of the number of effect
    def _get_gain_output_ratio(self, num_effect):
        return (
            self.feed_conc_ppm * self.gain_output_ratio_coeffs[num_effect][0]
            + self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][1]
            + self.feed_conc_ppm
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][2]
            + self.feed_temperature * self.gain_output_ratio_coeffs[num_effect][3]
            + self.feed_temperature
            * self.feed_conc_ppm
            * self.gain_output_ratio_coeffs[num_effect][4]
            + self.feed_temperature
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][5]
            + self.capacity * self.gain_output_ratio_coeffs[num_effect][6]
            + self.capacity
            * self.feed_conc_ppm
            * self.gain_output_ratio_coeffs[num_effect][7]
            + self.capacity
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][8]
            + self.capacity
            * self.feed_temperature
            * self.gain_output_ratio_coeffs[num_effect][9]
            + self.motive_pressure * self.gain_output_ratio_coeffs[num_effect][10]
            + self.motive_pressure
            * self.feed_conc_ppm
            * self.gain_output_ratio_coeffs[num_effect][11]
            + self.motive_pressure
            * self.recovery_vol_phase[0, "Liq"]
            * self.gain_output_ratio_coeffs[num_effect][12]
            + self.motive_pressure
            * self.feed_temperature
            * self.gain_output_ratio_coeffs[num_effect][13]
            + self.motive_pressure
            * self.capacity
            * self.gain_output_ratio_coeffs[num_effect][14]
            + 1 * self.gain_output_ratio_coeffs[num_effect][15]
            + self.motive_pressure**2 * self.gain_output_ratio_coeffs[num_effect][16]
            + self.capacity**2 * self.gain_output_ratio_coeffs[num_effect][17]
            + self.feed_temperature**2 * self.gain_output_ratio_coeffs[num_effect][18]
            + self.recovery_vol_phase[0, "Liq"] ** 2
            * self.gain_output_ratio_coeffs[num_effect][19]
            + self.feed_conc_ppm**2 * self.gain_output_ratio_coeffs[num_effect][20]
        )

    # Specific area surrogate equation, as a function of the number of effect
    def _get_specific_area(self, num_effect):
        return (
            self.feed_conc_ppm * self.specific_area_coeffs[num_effect][0]
            + self.recovery_vol_phase[0, "Liq"]
            * self.specific_area_coeffs[num_effect][1]
            + self.feed_conc_ppm
            * self.recovery_vol_phase[0, "Liq"]
            * self.specific_area_coeffs[num_effect][2]
            + self.feed_temperature * self.specific_area_coeffs[num_effect][3]
            + self.feed_temperature
            * self.feed_conc_ppm
            * self.specific_area_coeffs[num_effect][4]
            + self.feed_temperature
            * self.recovery_vol_phase[0, "Liq"]
            * self.specific_area_coeffs[num_effect][5]
            + self.capacity * self.specific_area_coeffs[num_effect][6]
            + self.capacity
            * self.feed_conc_ppm
            * self.specific_area_coeffs[num_effect][7]
            + self.capacity
            * self.recovery_vol_phase[0, "Liq"]
            * self.specific_area_coeffs[num_effect][8]
            + self.capacity
            * self.feed_temperature
            * self.specific_area_coeffs[num_effect][9]
            + self.motive_pressure * self.specific_area_coeffs[num_effect][10]
            + self.motive_pressure
            * self.feed_conc_ppm
            * self.specific_area_coeffs[num_effect][11]
            + self.motive_pressure
            * self.recovery_vol_phase[0, "Liq"]
            * self.specific_area_coeffs[num_effect][12]
            + self.motive_pressure
            * self.feed_temperature
            * self.specific_area_coeffs[num_effect][13]
            + self.motive_pressure
            * self.capacity
            * self.specific_area_coeffs[num_effect][14]
            + 1 * self.specific_area_coeffs[num_effect][15]
            + self.motive_pressure**2 * self.specific_area_coeffs[num_effect][16]
            + self.capacity**2 * self.specific_area_coeffs[num_effect][17]
            + self.feed_temperature**2 * self.specific_area_coeffs[num_effect][18]
            + self.recovery_vol_phase[0, "Liq"] ** 2
            * self.specific_area_coeffs[num_effect][19]
            + self.feed_conc_ppm**2 * self.specific_area_coeffs[num_effect][20]
        )

    # Heating steam mass flow rate surrogate equation, as a function of the number of effect
    def _get_heating_steam_mass_flow_rate(self, num_effect):
        return (
            self.feed_conc_ppm * self.heating_steam_mass_flow_rate_coeffs[num_effect][0]
            + self.recovery_vol_phase[0, "Liq"]
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][1]
            + self.feed_conc_ppm
            * self.recovery_vol_phase[0, "Liq"]
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][2]
            + self.feed_temperature
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][3]
            + self.feed_temperature
            * self.feed_conc_ppm
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][4]
            + self.feed_temperature
            * self.recovery_vol_phase[0, "Liq"]
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][5]
            + self.capacity * self.heating_steam_mass_flow_rate_coeffs[num_effect][6]
            + self.capacity
            * self.feed_conc_ppm
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][7]
            + self.capacity
            * self.recovery_vol_phase[0, "Liq"]
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][8]
            + self.capacity
            * self.feed_temperature
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][9]
            + self.motive_pressure
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][10]
            + self.motive_pressure
            * self.feed_conc_ppm
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][11]
            + self.motive_pressure
            * self.recovery_vol_phase[0, "Liq"]
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][12]
            + self.motive_pressure
            * self.feed_temperature
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][13]
            + self.motive_pressure
            * self.capacity
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][14]
            + 1 * self.heating_steam_mass_flow_rate_coeffs[num_effect][15]
            + self.motive_pressure**2
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][16]
            + self.capacity**2
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][17]
            + self.feed_temperature**2
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][18]
            + self.recovery_vol_phase[0, "Liq"] ** 2
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][19]
            + self.feed_conc_ppm**2
            * self.heating_steam_mass_flow_rate_coeffs[num_effect][20]
        )

    # Motive steam mass flow rate surrogate equation, as a function of the number of effect
    def _get_motive_steam_mass_flow_rate(self, num_effect):
        return (
            self.feed_conc_ppm * self.motive_steam_mass_flow_rate_coeffs[num_effect][0]
            + self.recovery_vol_phase[0, "Liq"]
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][1]
            + self.feed_conc_ppm
            * self.recovery_vol_phase[0, "Liq"]
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][2]
            + self.feed_temperature
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][3]
            + self.feed_temperature
            * self.feed_conc_ppm
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][4]
            + self.feed_temperature
            * self.recovery_vol_phase[0, "Liq"]
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][5]
            + self.capacity * self.motive_steam_mass_flow_rate_coeffs[num_effect][6]
            + self.capacity
            * self.feed_conc_ppm
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][7]
            + self.capacity
            * self.recovery_vol_phase[0, "Liq"]
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][8]
            + self.capacity
            * self.feed_temperature
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][9]
            + self.motive_pressure
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][10]
            + self.motive_pressure
            * self.feed_conc_ppm
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][11]
            + self.motive_pressure
            * self.recovery_vol_phase[0, "Liq"]
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][12]
            + self.motive_pressure
            * self.feed_temperature
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][13]
            + self.motive_pressure
            * self.capacity
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][14]
            + 1 * self.motive_steam_mass_flow_rate_coeffs[num_effect][15]
            + self.motive_pressure**2
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][16]
            + self.capacity**2
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][17]
            + self.feed_temperature**2
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][18]
            + self.recovery_vol_phase[0, "Liq"] ** 2
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][19]
            + self.feed_conc_ppm**2
            * self.motive_steam_mass_flow_rate_coeffs[num_effect][20]
        )

    @property
    def default_costing_method(self):
        return cost_med_tvc_surrogate
