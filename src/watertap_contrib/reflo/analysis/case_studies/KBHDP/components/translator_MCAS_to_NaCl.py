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
"""
Translator block for converting from MCAS to NaCl
"""

# Import Pyomo libraries
from pyomo.environ import check_optimal_termination
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.translator import TranslatorData
from idaes.core.util.config import (
    is_physical_parameter_block,
)

from watertap.core.solvers import get_solver

__author__ = "Zachary Binger"


@declare_process_block_class("Translator_MCAS_to_NaCl")
class Translator_MCAS_to_NaCl_Data(TranslatorData):
    """
    Translator block for converting from MCAS to NaCl
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Translator blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Translator blocks do not contain holdup.""",
        ),
    )
    CONFIG.declare(
        "outlet_state_defined",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Indicated whether outlet state will be fully defined",
            doc="""Indicates whether unit model will fully define outlet state.
If False, the outlet property package will enforce constraints such as sum
of mole fractions and phase equilibrium.
**default** - True.
**Valid values:** {
**True** - outlet state will be fully defined,
**False** - outlet property package should enforce sumation and equilibrium
constraints.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Indicates whether outlet is in phase equilibrium",
            doc="""Indicates whether outlet property package should enforce
phase equilibrium constraints.
**default** - False.
**Valid values:** {
**True** - outlet property package should calculate phase equilibrium,
**False** - outlet property package should notcalculate phase equilibrium.}
""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for incoming stream",
            doc="""Property parameter object used to define property
calculations for the incoming stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "inlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the incoming stream",
            doc="""A ConfigBlock with arguments to be passed to the property
block associated with the incoming stream,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for outgoing stream",
            doc="""Property parameter object used to define property
calculations for the outgoing stream,
**default** - None.
**Valid values:** {
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "outlet_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the outgoing stream",
            doc="""A ConfigBlock with arguments to be passed to the property
block associated with the outgoing stream,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """
        super(Translator_MCAS_to_NaCl_Data, self).build()
        solute_set = self.config.inlet_property_package.solute_set

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def eq_flow_vol_rule(blk, t):
            return (
                blk.properties_out[t].flow_mass_phase_comp["Liq", "H2O"]
                == blk.properties_in[t].flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality temperature equation",
        )
        def eq_temperature_rule(blk, t):
            return blk.properties_out[t].temperature == blk.properties_in[t].temperature

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality pressure equation",
        )
        def eq_pressure_rule(blk, t):
            return blk.properties_out[t].pressure == blk.properties_in[t].pressure

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality solute equation",
        )
        def eq_solute_mass_flow(blk, t):
            return blk.properties_out[t].flow_mass_phase_comp["Liq", "NaCl"] == sum(
                blk.properties_in[t].flow_mass_phase_comp["Liq", i] for i in solute_set
            )

    def initialize_build(
        self,
        state_args_in=None,
        state_args_out=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        This method calls the initialization method of the state blocks.

        Keyword Arguments:
            state_args_in : a dict of arguments to be passed to the inlet
                            property package (to provide an initial state for
                            initialization (see documentation of the specific
                            property package) (default = None).
            state_args_out : a dict of arguments to be passed to the outlet
                             property package (to provide an initial state for
                             initialization (see documentation of the specific
                             property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_in,
            hold_state=True,
        )

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        self.properties_in.release_state(flags=flags, outlvl=outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
