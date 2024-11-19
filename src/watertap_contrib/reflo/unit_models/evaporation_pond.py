#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import math
from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)

from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.units.evaporation_pond import cost_evaporation_pond

__author__ = "Kurban Sitterley"


"""
REFERENCES



"""


_log = idaeslog.getLogger(__name__)

area_correction_factor_param_dict = {
    12: (2.6429, -0.202),
    8: (2.0512, -0.152),
    4: (1.5357, -0.092),
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
        "dike_height",
        ConfigValue(
            default=8,
            domain=PositiveInt,
            description="Height of the dike used for evaporation pond.",
            doc="""A ConfigBlock specifying the height of the dike of the evaporation pond. Units are ft. Must be 4, 8, or 12. Determines coefficient for calculating various costing and design parameters.""",
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

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of liquid outlet",
            **tmp_dict,
        )

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_out.flow_mass_phase_comp["Liq", "H2O"].set_value(0)

        self.air_temperature_C = prop_in.temperature["Vap"] - 273.15 * pyunits.degK

        self.add_inlet_port(name="inlet", block=self.properties_in)
        self.add_outlet_port(name="outlet", block=self.properties_out)

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

        ## solids precipitation rate is a function of TDS concentration [g / L]
        ## solids_precipitation_rate [ft / yr] = a1 * C**2 + a2 * C + intercept
        ## solids_precipitation_rate [ft / yr] = 4.12e-6 * C**2 + 1.92e-4 * C + 1.15e-3

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

        self.pressure_atm = Param(
            initialize=101325,
            mutable=True,
            units=pyunits.Pa,
            doc="Atmospheric pressure",
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

        self.evaporation_rate_enhancement_adjustment_factor = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor to increase evaporation rate due to enhancement",
        )

        self.water_temperature_calc_slope = Param(
            initialize=1.04,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Slope of equation to calculate water temperature based on air temperature",
        )

        self.water_temperature_calc_intercept = Param(
            initialize=0.22,
            mutable=True,
            units=pyunits.degK,
            doc="Intercept of equation to calculate water temperature based on air temperature",
        )

        self.net_solar_radiation = Param(
            initialize=150,
            mutable=True,
            units=pyunits.watt / pyunits.m**2,
            doc="Net incident solar radiation",  # net shortwave radiation - net longwave radiation
        )

        self.differential_head = Param(
            initialize=40,
            mutable=True,
            units=pyunits.m,
            doc="Differential head for pumping energy",
        )

        self.net_heat_flux_out = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.watt / pyunits.m**2,
            doc="Net heat flux out of system (water, soil, ecosystem, etc.)",
        )

        self.area_correction_factor = Var(
            initialize=1,
            bounds=(0.99, 10),
            units=pyunits.dimensionless,
            doc="Area correction factor",
        )

        self.bowen_ratio = Var(
            initialize=0.3,
            bounds=(0, 2),
            units=pyunits.dimensionless,
            doc="Bowen ratio for BREB calculation of evaporation rate",
        )

        self.psychrometric_constant = Var(
            initialize=0.06,
            bounds=(0, None),
            units=pyunits.kPa * pyunits.degK**-1,
            doc="Psychrometric constant",
        )

        self.mass_flux_water_vapor = Var(
            initialize=1e-5,
            bounds=(0, 1e-3),
            units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            doc="Mass flux of water vapor evaporated according using BREB method",
        )

        self.evaporation_rate = Var(
            initialize=0.03,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Evaporation rate",
        )

        self.solids_precipitation_rate = Var(
            initialize=0.01,
            bounds=(0, None),
            units=pyunits.feet / pyunits.year,
            doc="Rate at which solids precipitate on bottom of pond",
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

        @self.Expression()
        def evap_rate_mm_d(b):
            return pyunits.convert(
                b.evaporation_rate, to_units=pyunits.mm / pyunits.day
            )

        @self.Expression(doc="Total pond area in acres")
        def total_pond_area_acre(b):
            return pyunits.convert(
                b.evaporation_pond_area * b.number_evaporation_ponds,
                to_units=pyunits.acre,
            )

        @self.Expression(doc="Evaporative area per pond in acres")
        def evaporative_area_acre(b):
            return pyunits.convert(b.evaporative_area_per_pond, to_units=pyunits.acre)

        @self.Expression(doc="Individual pond area in acres")
        def pond_area_acre(b):
            return pyunits.convert(b.evaporation_pond_area, to_units=pyunits.acre)

        @self.Expression(doc="Net radiation available for evaporation")
        def net_radiation(b):
            return b.net_solar_radiation - b.net_heat_flux_out

        @self.Expression(doc="Mass flow of precipitated solids")
        def mass_flow_precipitate(b):
            return pyunits.convert(
                b.total_pond_area_acre * b.solids_precipitation_rate * b.dens_solids,
                to_units=pyunits.kg / pyunits.year,
            )

        @self.Expression(doc="Differential pressure required for pumping")
        def differential_pressure(b):
            return pyunits.convert(
                prop_in.dens_mass_phase["Liq"]
                * Constants.acceleration_gravity
                * b.differential_head,
                to_units=pyunits.Pa,
            )

        # @self.Constraint(doc="Mass balance of liquid water")
        # def eq_mass_balance_liquid_water(b):
        #     return (
        #         prop_out.flow_mass_phase_comp["Liq", "H2O"]
        #         == prop_in.flow_mass_phase_comp["Liq", "H2O"]
        #         - prop_out.flow_mass_phase_comp["Vap", "H2O"]
        #     )

        @self.Constraint(
            # non_volatile_comps,
            doc="Mass transfer term for precipitated solids and non-volatile components",
        )
        def eq_mass_transfer_non_volatile_solutes(b):
            return prop_out.flow_mass_phase_comp["Liq", "TDS"] == pyunits.convert(
                b.mass_flow_precipitate,
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(doc="Mass transfer term for evaporated water")
        def eq_mass_transfer_evaporated_water(b):
            return prop_out.flow_mass_phase_comp["Vap", "H2O"] == pyunits.convert(
                b.mass_flux_water_vapor * b.total_evaporative_area_required,
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(doc="Overall mass balance of water")
        def eq_mass_balance_water(b):
            return (
                prop_in.flow_mass_phase_comp["Liq", "H2O"]
                + prop_in.flow_mass_phase_comp["Vap", "H2O"]
                == prop_out.flow_mass_phase_comp["Liq", "H2O"]
                + prop_out.flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(doc="Mass balance for air")
        def eq_air_mass_balance(b):
            return (
                prop_in.flow_mass_phase_comp["Vap", "Air"]
                == prop_out.flow_mass_phase_comp["Vap", "Air"]
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

        @self.Constraint(doc="Water temperature as function of air temperature")
        def eq_water_temperature_constraint(b):
            return (
                prop_out.temperature["Liq"]
                == (
                    b.water_temperature_calc_slope * b.air_temperature_C
                    + b.water_temperature_calc_intercept
                )
                + 273.15 * pyunits.degK
            )

        @self.Constraint(doc="Net heat flux out of surroundings/ecosystem")
        def eq_net_heat_flux_out(b):
            daily_temperature_change = (
                prop_out.temperature["Liq"] - prop_in.temperature["Vap"]
            ) / pyunits.day
            return b.net_heat_flux_out == pyunits.convert(
                prop_in.dens_mass_phase["Liq"]
                * prop_in.cp_mass_solvent["Liq"]
                * b.evaporation_pond_depth
                * daily_temperature_change,
                to_units=pyunits.watt / pyunits.m**2,
            )

        @self.Constraint(doc="Air temperature is isothermal")
        def eq_isothermal(b):
            return prop_in.temperature["Vap"] == prop_out.temperature["Vap"]

        @self.Constraint()
        def eq_isobaric(b):
            return prop_in.pressure == prop_out.pressure

        @self.Constraint(doc="Psychrometric constant equation")
        def eq_psychrometric_constant(b):
            mw_ratio = pyunits.convert(
                prop_in.mw_comp["H2O"] / prop_in.mw_comp["Air"],
                to_units=pyunits.dimensionless,
            )
            return b.psychrometric_constant == pyunits.convert(
                (prop_in.cp_air * b.pressure_atm)
                / (mw_ratio * prop_out.dh_vap_mass_solvent),
                to_units=pyunits.kPa * pyunits.degK**-1,
            )

        @self.Constraint(doc="Bowen ratio calculation")
        def eq_bowen_ratio(b):
            return b.bowen_ratio == pyunits.convert(
                b.psychrometric_constant
                * (
                    (prop_out.temperature["Liq"] - prop_in.temperature["Vap"])
                    / (prop_in.pressure_vap_sat["H2O"] - prop_in.pressure_vap["H2O"])
                ),
                to_units=pyunits.dimensionless,
            )

        @self.Constraint(doc="Mass flux water vapor using BREB method")
        def eq_mass_flux_water_vapor(b):
            return (
                b.mass_flux_water_vapor
                == b.evaporation_rate_salinity_adjustment_factor
                * b.evaporation_rate_enhancement_adjustment_factor
                * pyunits.convert(
                    b.net_radiation
                    / (prop_out.dh_vap_mass_solvent * (1 + b.bowen_ratio)),
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.s),
                )
            )

        @self.Constraint(doc="Evaporation rate")
        def eq_evaporation_rate(b):
            return b.evaporation_rate == pyunits.convert(
                b.mass_flux_water_vapor / prop_in.dens_mass_phase["Liq"],
                to_units=pyunits.m / pyunits.s,
            )

        @self.Constraint(doc="Total evaporative area required")
        def eq_total_evaporative_area_required(b):
            return (
                b.total_evaporative_area_required * b.mass_flux_water_vapor
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

        # ---------------------------------------------------------------------
        flags = self.properties_in.initialize(
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
            self.state_args = state_args = {}
            state_dict = self.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        calculate_variable_from_constraint(
            self.properties_out[0].temperature["Liq"],
            self.eq_water_temperature_constraint,
        )
        calculate_variable_from_constraint(
            self.net_heat_flux_out, self.eq_net_heat_flux_out
        )
        calculate_variable_from_constraint(
            self.psychrometric_constant, self.eq_psychrometric_constant
        )
        calculate_variable_from_constraint(self.bowen_ratio, self.eq_bowen_ratio)
        calculate_variable_from_constraint(
            self.net_heat_flux_out, self.eq_net_heat_flux_out
        )
        calculate_variable_from_constraint(
            self.mass_flux_water_vapor, self.eq_mass_flux_water_vapor
        )
        calculate_variable_from_constraint(
            self.total_evaporative_area_required,
            self.eq_total_evaporative_area_required,
        )

        state_args_out = deepcopy(state_args)
        for k, v in state_args_out.items():
            if k == "flow_mass_phase_comp":
                for p, j in v.keys():
                    if p == "Liq":
                        if j == "H2O":
                            state_args_out[k][(p, j)] = 0
                        elif j == "TDS":
                            state_args_out[k][(p, j)] = state_args[k][(p, j)] * 0.5
                        else:
                            state_args_out[k][(p, j)] = state_args[k][(p, j)]
                    if p == "Vap":
                        if j == "H2O":
                            state_args_out[k][(p, j)] = state_args[k][("Liq", j)]
                        elif j == "Air":
                            state_args_out[k][(p, j)] = state_args[k][(p, j)]
                        else:
                            state_args_out[k][(p, j)] = state_args[k][(p, j)]

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            nponds = math.floor(self.number_evaporation_ponds.value)
            if nponds < 1:
                nponds = 1
            self.number_evaporation_ponds.fix(nponds)
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        self.number_evaporation_ponds.unfix()

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.net_heat_flux_out) is None:
            iscale.set_scaling_factor(self.net_heat_flux_out, 0.01)

        if iscale.get_scaling_factor(self.area_correction_factor) is None:
            iscale.set_scaling_factor(self.area_correction_factor, 1)

        if iscale.get_scaling_factor(self.bowen_ratio) is None:
            iscale.set_scaling_factor(self.bowen_ratio, 10)

        if iscale.get_scaling_factor(self.psychrometric_constant) is None:
            iscale.set_scaling_factor(self.psychrometric_constant, 100)

        if iscale.get_scaling_factor(self.mass_flux_water_vapor) is None:
            iscale.set_scaling_factor(self.mass_flux_water_vapor, 1e5)

        if iscale.get_scaling_factor(self.evaporation_rate) is None:
            iscale.set_scaling_factor(self.evaporation_rate, 1e8)

        if iscale.get_scaling_factor(self.solids_precipitation_rate) is None:
            iscale.set_scaling_factor(self.solids_precipitation_rate, 1e2)

        if iscale.get_scaling_factor(self.total_evaporative_area_required) is None:
            iscale.set_scaling_factor(self.total_evaporative_area_required, 1e-4)

        if iscale.get_scaling_factor(self.number_evaporation_ponds) is None:
            iscale.set_scaling_factor(self.number_evaporation_ponds, 0.1)

        if iscale.get_scaling_factor(self.evaporative_area_per_pond) is None:
            iscale.set_scaling_factor(self.evaporative_area_per_pond, 1e-3)

        if iscale.get_scaling_factor(self.evaporation_pond_area) is None:
            iscale.set_scaling_factor(self.evaporation_pond_area, 1e-3)

    def _get_stream_table_contents(self, time_point=0):

        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Precipitated Solids Outlet": self.solids,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = dict()
        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_evaporation_pond
