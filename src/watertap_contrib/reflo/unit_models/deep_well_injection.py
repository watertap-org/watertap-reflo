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

from pyomo.environ import (
    check_optimal_termination,
    Param,
    Suffix,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.math import smooth_bound
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing.units.deep_well_injection import (
    cost_deep_well_injection,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("DeepWellInjection")
class DeepWellInjectionData(InitializationMixin, UnitModelBlockData):
    """
    Zero order deep well injection model.
    This is a terminal unit and has no outlet.
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
        "injection_well_depth",
        ConfigValue(
            default=2500,
            domain=PositiveInt,
            description="Depth of injection well. Costing is available for 2500, 5000, 7500, and 10000 ft depths.",
            doc="""Depth of injection well for use with BLM costing approach.""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        if self.config.injection_well_depth not in [2500, 5000, 7500, 10000]:
            raise ConfigurationError(
                f"The injection well depth was specified as {self.config.injection_well_depth}. The injection well depth must be 2500, 5000, 7500, or 10000."
            )

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        self.add_port(name="inlet", block=self.properties)

        self.pipe_diameter_coeff = Param(
            initialize=5.329,
            units=pyunits.dimensionless,
            doc="Base parameter for pipe diameter equation",
        )

        self.pipe_diameter_exponent = Param(
            initialize=0.4998,
            units=pyunits.dimensionless,
            doc="Exponent parameter for pipe diameter equation",
        )

        self.injection_pressure = Param(
            initialize=5,
            units=pyunits.bar,
            doc="Pressure required for injection",
        )

        self.injection_well_depth = Param(
            initialize=self.config.injection_well_depth,
            units=pyunits.ft,
            doc="Depth of injection well",
        )

        self.monitoring_well_depth = Param(
            initialize=1000,
            mutable=True,
            units=pyunits.ft,
            doc="Depth of monitoring well",
        )

        @self.Expression(doc="Injection pipe diameter in inches")
        def pipe_diameter(b):
            flow_mgd_dimensionless = pyunits.convert(
                pyunits.convert(
                    b.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.Mgallons / pyunits.day,
                )
                * pyunits.day
                * pyunits.Mgallons**-1,
                to_units=pyunits.dimensionless,
            )
            pipe_diameter = (
                b.pipe_diameter_coeff * flow_mgd_dimensionless**b.pipe_diameter_exponent
            )
            return smooth_bound(pipe_diameter, 2, 24) * pyunits.inches

    def initialize(
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

        flags = blk.properties.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        blk.properties.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {"Feed Inlet": self.inlet},
            time_point=time_point,
        )

    @property
    def default_costing_method(self):
        return cost_deep_well_injection
