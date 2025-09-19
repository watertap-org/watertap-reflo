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

from pyomo.environ import check_optimal_termination

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.translator import TranslatorData

from watertap.core.solvers import get_solver

__author__ = "Mukta Hardikar"


@declare_process_block_class("TranslatorNaCltoNaCl")
class TranslatorNaCltoNaClData(TranslatorData):
    """
    Translator block for Seawater to NaCl property packages
    """

    CONFIG = TranslatorData.CONFIG()

    def build(self):

        super().build()

        @self.Constraint(doc="Isothermal")
        def eq_temperature(b):
            return b.properties_in[0].temperature == b.properties_out[0].temperature

        @self.Constraint(doc="Isobaric")
        def eq_pressure(b):
            return b.properties_in[0].pressure == b.properties_out[0].pressure

        @self.Constraint(
            doc="Equality mass flow water equation",
        )
        def eq_flow_mass_water(b):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(
            doc="Equality mass flow TDS equation",
        )
        def eq_flow_mass_tds(b):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
                == b.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"]
            )

        self.properties_out[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        self.properties_out[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)

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
