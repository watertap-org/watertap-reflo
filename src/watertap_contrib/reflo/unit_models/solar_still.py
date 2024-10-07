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

from copy import deepcopy
from pyomo.environ import (
    Var,
    Constraint,
    check_optimal_termination,
    Param,
    Suffix,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block

from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap_contrib.reflo.costing.units.solar_still import cost_solar_still
from watertap_contrib.reflo.unit_models.util.water_yield_calculation import (
    get_solar_still_daily_water_yield,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("SolarStill")
class SolarStillData(InitializationMixin, UnitModelBlockData):
    """
    Zero dimensional solar still model
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
        "water_yield_calculation_args",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Dict of arguments to pass to water yield calculation function.",
            doc="""Dictionary of arguments passed to helper function get_solar_still_daily_water_yield""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        unit_log = idaeslog.getModelLogger(self.name, tag="unit")

        if len(self.config.water_yield_calculation_args) > 0:
            required_args = ["input_weather_file_path", "interval_day", "salinity", "water_depth_basin", "length_basin"]
            if not all(a in self.config.water_yield_calculation_args.keys() for a in required_args):
                raise ConfigurationError(
                    f'Water yield calculation dict must contain values for "input_weather_file_path", "interval_day", '
                    f'"salinity", "water_depth_basin", and "length_basin", but only '
                    f'{[*self.config.water_yield_calculation_args.keys()]}  were found.'
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
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_waste = self.properties_waste[0]
        prop_waste.flow_mass_phase_comp["Liq", "H2O"].fix(0)
        comps = self.config.property_package.component_list

        self.water_yield = Param(
            initialize=2.5,
            mutable=True,
            units=pyunits.kg / (pyunits.m**2 * pyunits.day),
            doc="Average daily mass water yield per unit area",
        )

        self.length_basin = Param(
            initialize=0.6,
            mutable=True,
            units=pyunits.m,
            doc="Dimension of one side of solar still",
        )

        self.water_depth_basin = Param(
            initialize=0.02,
            mutable=True,
            units=pyunits.m,
            doc="Depth of water in solar solar still",
        )

        if self.config.water_yield_calculation_args != {}:

            unit_log.info(
                f"Found water yield calculation arguments in {self.name} configuration."
                f" Calculating daily water yield assuming {self.config.water_yield_calculation_args['salinity']}"
                f" g/L TDS with {self.config.water_yield_calculation_args['interval_day']} day intervals..."
            )

            daily_water_yield_mass = get_solar_still_daily_water_yield(
                **self.config.water_yield_calculation_args
            )

            unit_log.info(f"Water yield calculation complete.")

            self.water_yield.set_value(daily_water_yield_mass)
            self.length_basin.set_value(
                self.config.water_yield_calculation_args["length_basin"]
            )
            self.water_depth_basin.set_value(
                self.config.water_yield_calculation_args["water_depth_basin"]
            )

        else:

            unit_log.info(
                f"No water yield calculation arguments found in {self.name} configuration. Assuming default values."
            )

        self.dens_mass_salt = Param(
            initialize=2.16,
            units=pyunits.g * pyunits.cm**-3,
            doc="Density of dried salts",
        )

        self.number_stills = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of solar stills",
        )

        self.total_area = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total area of solar still",
        )

        @self.Expression(doc="Annual water yield per unit area")
        def annual_water_yield(b):
            return pyunits.convert(
                b.water_yield
                / b.properties_in[0].dens_mass_phase["Liq"]
                * 365
                * pyunits.day,
                to_units=pyunits.m**3 / pyunits.m**2,
            )

        @self.Expression(doc="Volumetric flow of salts")
        def flow_vol_salt(b):
            return pyunits.convert(
                sum(
                    prop_waste.flow_mass_phase_comp["Liq", jj]
                    for jj in comps
                    if jj != "H2O"
                )
                / b.dens_mass_salt,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Expression(doc="Deposition rate for salts")
        def deposition_rate(b):
            return pyunits.convert(
                b.flow_vol_salt / b.total_area,
                to_units=pyunits.cm / pyunits.d,
            )

        @self.Expression(doc="Area of single still")
        def area_single_still(b):
            return b.length_basin**2

        @self.Expression(doc="Production per still")
        def yield_per_still(b):
            return pyunits.convert(
                b.water_yield * b.area_single_still, to_units=pyunits.kg / pyunits.s
            )

        @self.Expression(doc="Evaporation rate")
        def evaporation_rate(b):
            evap_rate = b.water_yield / prop_out.dens_mass_phase["Liq"]
            return pyunits.convert(evap_rate, to_units=pyunits.mm / pyunits.day)

        @self.Constraint(doc="Number of solar stills")
        def eq_number_stills(b):
            return b.number_stills == pyunits.convert(
                b.total_area / b.area_single_still, to_units=pyunits.dimensionless
            )

        @self.Constraint(doc="Total area of solar stills")
        def eq_total_area(b):
            return b.total_area == pyunits.convert(
                (prop_out.flow_vol_phase["Liq"] * prop_out.dens_mass_phase["Liq"])
                / b.water_yield,
                to_units=pyunits.m**2,
            )

        @self.Constraint(comps, doc="Mass flow out")
        def eq_flow_mass_out(b, j):
            if j == "H2O":
                return (
                    prop_out.flow_mass_phase_comp["Liq", j]
                    == prop_in.flow_mass_phase_comp["Liq", j]
                )
            else:
                prop_out.flow_mass_phase_comp["Liq", j].fix(0)
                return Constraint.Skip

        @self.Constraint(comps, doc="Mass balance")
        def eq_mass_balance(b, j):
            if j == "H2O":
                # return prop_in.flow_mass_phase_comp["Liq", j] == prop_out.flow_mass_phase_comp["Liq", j]
                return Constraint.Skip
            else:
                return (
                    prop_in.flow_mass_phase_comp["Liq", j]
                    == prop_waste.flow_mass_phase_comp["Liq", j]
                )

    def initialize_build(
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

        if state_args is None:
            self.state_args = state_args = {}
            self.state_dict = state_dict = self.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args_out = deepcopy(state_args)

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_waste = deepcopy(state_args)

        self.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        self.state_args_out = state_args_out
        self.state_args_waste = state_args_waste

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.number_stills) is None:
            iscale.set_scaling_factor(self.number_stills, 0.1)

        if iscale.get_scaling_factor(self.total_area) is None:
            iscale.set_scaling_factor(self.total_area, 0.1)

        if iscale.get_scaling_factor(self.water_yield) is None:
            iscale.set_scaling_factor(self.water_yield, 1e-2)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Waste Outlet": self.waste,
            },
            time_point=time_point,
        )

    @property
    def default_costing_method(self):
        return cost_solar_still
