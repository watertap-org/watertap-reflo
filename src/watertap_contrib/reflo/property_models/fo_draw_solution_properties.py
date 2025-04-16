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

"""
This module contains a property package for a proprietary draw solution used in
Trevi's FO system.
"""
from pyomo.environ import (
    Param,
    Expression,
    units as pyunits,
    Reals,
    NonNegativeReals,
    Var,
    Constraint,
    Suffix,
    value,
    check_optimal_termination,
    Expr_if,
    units as pyunits,
)
from pyomo.gdp import Disjunction
from idaes.core import (
    EnergyBalanceType,
    MaterialBalanceType,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.base.components import Solvent, Solute
from idaes.core.base.phases import LiquidPhase
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
)

from watertap.core.util.scaling import transform_property_constraints
from watertap.core.solvers import get_solver

__author__ = "Zhuoran Zhang"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("FODrawSolutionParameterBlock")
class FODrawSolutionParameterBlockData(PhysicalParameterBlock):
    """
    Property parameters derived from experimental data of the draw
    solution that is used in the Trevi's FO system
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = FODrawSolutionStateBlock

        # components
        self.H2O = Solvent()
        self.DrawSolution = Solute()

        # phases
        self.Liq = LiquidPhase()

        # ---------------------------------------------------------------------
        # mass density parameters, equation derived from experimental data
        dens_units = pyunits.kg / pyunits.m**3

        self.dens_mass_param_A0 = Param(
            initialize=1.000446e3,
            units=dens_units,
            doc="Mass density parameter A0",
        )

        self.dens_mass_param_A1 = Param(
            initialize=1.1849e2,
            units=dens_units,
            doc="Mass density parameter A1",
        )

        self.dens_mass_param_A2 = Param(
            initialize=1.1455e2,
            units=dens_units,
            doc="Mass density parameter A2",
        )

        self.dens_mass_param_A3 = Param(
            initialize=-1.6343e2,
            units=dens_units,
            doc="Mass density parameter A3",
        )

        # osmotic coefficient parameters, equation derived from experimental data
        self.osm_coeff_param_0 = Param(
            initialize=-1.24602e5,
            units=pyunits.Pa,
            doc="Osmotic coefficient parameter 0",
        )
        self.osm_coeff_param_1 = Param(
            initialize=6.39608e6,
            units=pyunits.Pa,
            doc="Osmotic coefficient parameter 1",
        )
        self.osm_coeff_param_2 = Param(
            initialize=-2.81942e7,
            units=pyunits.Pa,
            doc="Osmotic coefficient parameter 2",
        )
        self.osm_coeff_param_3 = Param(
            initialize=7.37766e7,
            units=pyunits.Pa,
            doc="Osmotic coefficient parameter 3",
        )
        # specific heat parameters, derived from experimental data
        cp_units = pyunits.J / (pyunits.kg * pyunits.K)
        self.cp_phase_param_A0 = Param(
            initialize=4.20003e3,
            units=cp_units,
            doc="Specific heat of seawater parameter A0",
        )
        self.cp_phase_param_A1 = Param(
            initialize=-3.49888e3,
            units=cp_units,
            doc="Specific heat of seawater parameter A1",
        )
        self.cp_phase_param_A2 = Param(
            initialize=1.33883e3,
            units=cp_units,
            doc="Specific heat of seawater parameter A2",
        )
        self.separation_heat = Param(
            initialize=105,
            units=pyunits.MJ / pyunits.m**3,
            doc="Separation heat per m3 of the separated water",
        )
        # ---------------------------------------------------------------------
        # Set default scaling factors
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("cp_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("enth_mass_phase", 1e-5, index="Liq")
        self.set_default_scaling("heat_separation_phase", 1e-5, index="Liq")
        self.set_default_scaling("liquid_separation", 1, index="Liq")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )  # Not used

        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "pressure_osm_phase": {"method": "_pressure_osm_phase"},
                "cp_mass_phase": {"method": "_cp_mass_phase"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "heat_separation_phase": {"method": "_heat_separation_phase"},
                "enth_flow": {"method": "_enth_flow"},
            }
        )


class _FODrawSolutionStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.

        Keyword Arguments:
        state_args : Dictionary with initial guesses for the state vars
                     chosen. Note that if this method is triggered
                     through the control volume, and if initial guesses
                     were not provied at the unit model level, the
                     control volume passes the inlet values as initial
                     guess.
        outlvl : sets output level of initialization routine
        state_vars_fixed: Flag to denote if state vars have already been
                          fixed.
                          - True - states have already been fixed and
                                   initialization does not need to worry
                                   about fixing and unfixing variables.
                         - False - states have not been fixed. The state
                                   block will deal with fixing/unfixing.
        optarg : solver options dictionary object (default=None, use
                 default solver options)
        solver : str indicating which solver to use during
                 initialization (default = None, use default solver)
        hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization (default=False).
                     - True - states varaibles are not unfixed, and
                             a dict of returned containing flags for
                             which states were fixed during
                             initialization.
                    - False - state variables are unfixed after
                             initialization by calling the
                             relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix state variables
        flags = fix_state_vars(self, state_args)

        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise InitializationError(
                    "\nWhile initializing {sb_name}, the degrees of freedom "
                    "are {dof}, when zero is required. \nInitialization assumes "
                    "that the state variables should be fixed and that no other "
                    "variables are fixed. \nIf other properties have a "
                    "predetermined value, use the calculate_state method "
                    "before using initialize to determine the values for "
                    "the state variables and avoid fixing the property variables."
                    "".format(sb_name=self.name, dof=dof)
                )

        # # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        if not skip_solve:
            # Initialize properties
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solve_indexed_blocks(opt, [self], tee=slc.tee)
                if not check_optimal_termination(results):
                    raise InitializationError(
                        "The property package failed to solve during initialization."
                    )
            init_log.info_high(
                "Property initialization: {}.".format(idaeslog.condition(results))
            )

        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")

        # Unfix state variables
        revert_state_vars(self, flags)
        init_log.info("State Released.")

    def calculate_state(
        self,
        var_args=None,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Solves state blocks given a set of variables and their values. These variables can
        be state variables or properties. This method is typically used before
        initialization to solve for state variables because non-state variables (i.e. properties)
        cannot be fixed in initialization routines.
        Keyword Arguments:
            var_args : dictionary with variables and their values, they can be state variables or properties
                       {(VAR_NAME, INDEX): VALUE}
            hold_state : flag indicating whether all of the state variables should be fixed after calculate state.
                         True - State variables will be fixed.
                         False - State variables will remain unfixed, unless already fixed.
            outlvl : idaes logger object that sets output level of solve call (default=idaeslog.NOTSET)
            solver : solver name string if None is provided the default solver
                     for IDAES will be used (default = None)
            optarg : solver options dictionary object (default={})
        Returns:
            results object from state block solve
        """
        # Get logger
        solve_log = idaeslog.getSolveLogger(self.name, level=outlvl, tag="properties")

        # Initialize at current state values (not user provided)
        self.initialize(solver=solver, optarg=optarg, outlvl=outlvl)

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix variables and check degrees of freedom
        flags = (
            {}
        )  # dictionary noting which variables were fixed and their previous state
        for k in self.keys():
            sb = self[k]
            for (v_name, ind), val in var_args.items():
                var = getattr(sb, v_name)
                if iscale.get_scaling_factor(var[ind]) is None:
                    _log.warning(
                        "While using the calculate_state method on {sb_name}, variable {v_name} "
                        "was provided as an argument in var_args, but it does not have a scaling "
                        "factor. This suggests that the calculate_scaling_factor method has not been "
                        "used or the variable was created on demand after the scaling factors were "
                        "calculated. It is recommended to touch all relevant variables (i.e. call "
                        "them or set an initial value) before using the calculate_scaling_factor "
                        "method.".format(v_name=v_name, sb_name=sb.name)
                    )
                if var[ind].is_fixed():
                    flags[(k, v_name, ind)] = True
                    if value(var[ind]) != val:
                        raise ConfigurationError(
                            "While using the calculate_state method on {sb_name}, {v_name} was "
                            "fixed to a value {val}, but it was already fixed to value {val_2}. "
                            "Unfix the variable before calling the calculate_state "
                            "method or update var_args."
                            "".format(
                                sb_name=sb.name,
                                v_name=var.name,
                                val=val,
                                val_2=value(var[ind]),
                            )
                        )
                else:
                    flags[(k, v_name, ind)] = False
                    var[ind].fix(val)

            if degrees_of_freedom(sb) != 0:
                raise RuntimeError(
                    "While using the calculate_state method on {sb_name}, the degrees "
                    "of freedom were {dof}, but 0 is required. Check var_args and ensure "
                    "the correct fixed variables are provided."
                    "".format(sb_name=sb.name, dof=degrees_of_freedom(sb))
                )

        # Solve
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            solve_log.info_high(
                "Calculate state: {}.".format(idaeslog.condition(results))
            )

        if not check_optimal_termination(results):
            _log.warning(
                "While using the calculate_state method on {sb_name}, the solver failed "
                "to converge to an optimal solution. This suggests that the user provided "
                "infeasible inputs, or that the model is poorly scaled, poorly initialized, "
                "or degenerate."
            )

        # unfix all variables fixed with var_args
        for (k, v_name, ind), previously_fixed in flags.items():
            if not previously_fixed:
                var = getattr(self[k], v_name)
                var[ind].unfix()

        # fix state variables if hold_state
        if hold_state:
            fix_state_vars(self)

        return results


@declare_process_block_class(
    "FODrawSolutionStateBlock", block_class=_FODrawSolutionStateBlock
)
class FODrawSolutionStateBlockData(StateBlockData):
    """A FO draw solution property package."""

    def build(self):
        """Callable method for Block construction."""
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.2, ("Liq", "DrawSolution"): 0.8},
            bounds=(0.0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 373.15),
            domain=NonNegativeReals,
            units=pyunits.K,
            doc="Temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc="Pressure",
        )

    # -----------------------------------------------------------------------------
    # Property Methods
    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=0.1,
            bounds=(0.0, 1),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] / sum(
                b.flow_mass_phase_comp[p, j] for j in b.params.component_list
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_mass_frac_phase_comp,
        )

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(
            self.params.phase_list,
            initialize=1e3,
            bounds=(1, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density of Trevi's FO draw solution",
        )

        def rule_dens_mass_phase(b, p):  # density, eqn derived from experimental data
            s = b.mass_frac_phase_comp[p, "DrawSolution"]
            dens_mass = (
                b.params.dens_mass_param_A0
                + b.params.dens_mass_param_A1 * s
                + b.params.dens_mass_param_A2 * s**2
                + b.params.dens_mass_param_A3 * s**3
            )
            return b.dens_mass_phase[p] == dens_mass

        self.eq_dens_mass_phase = Constraint(
            self.params.phase_list, rule=rule_dens_mass_phase
        )

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=0.001,
            bounds=(0.0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            return (
                b.flow_vol_phase[p]
                == sum(b.flow_mass_phase_comp[p, j] for j in b.params.component_list)
                / b.dens_mass_phase[p]
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=100,
            bounds=(0.0, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, p, j):
            return (
                b.conc_mass_phase_comp[p, j]
                == b.dens_mass_phase[p] * b.mass_frac_phase_comp[p, j]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_conc_mass_phase_comp,
        )

    def _pressure_osm_phase(self):
        self.pressure_osm_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(-1e6, 1e8),
            units=pyunits.Pa,
            doc="Osmotic pressure",
        )

        def rule_pressure_osm_phase(
            b, p
        ):  # osmotic pressure, derived from experimental data
            s = b.mass_frac_phase_comp[p, "DrawSolution"]
            pressure_osm_phase = (
                b.params.osm_coeff_param_0
                + b.params.osm_coeff_param_1 * s
                + b.params.osm_coeff_param_2 * s**2
                + b.params.osm_coeff_param_3 * s**3
            )
            return b.pressure_osm_phase[p] == pressure_osm_phase

        self.eq_pressure_osm_phase = Constraint(
            self.params.phase_list, rule=rule_pressure_osm_phase
        )

    def _cp_mass_phase(self):
        self.cp_mass_phase = Var(
            self.params.phase_list,
            initialize=4e3,
            bounds=(0.0, 1e8),
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity",
        )

        def rule_cp_mass_phase(
            b, p
        ):  # specific heat, equation derived from experimental data
            s = b.mass_frac_phase_comp[p, "DrawSolution"]
            cp_mass_phase = (
                b.params.cp_phase_param_A0
                + b.params.cp_phase_param_A1 * s
                + b.params.cp_phase_param_A2 * s**2
            )

            return b.cp_mass_phase[p] == cp_mass_phase

        self.eq_cp_mass_phase = Constraint(
            self.params.phase_list, rule=rule_cp_mass_phase
        )

    def _heat_separation_phase(self):
        self.heat_separation_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(-1e10, 1e10),
            units=pyunits.J / pyunits.s,
            doc="Heat of liquid separation",
        )
        self.liquid_separation = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="Indication of whether liquid separation will happen",
        )
        # Fix this flag to default value to not consider separation heat in most cases
        self.liquid_separation.fix()

        self.mass_frac_after_separation = Param(
            initialize=0.8,
            mutable=True,
            units=pyunits.dimensionless,
            doc="The strong draw mass fraction after separation",
        )

        def rule_heat_separation_phase(b, p):
            heat = (
                b.params.separation_heat
                * b.flow_mass_phase_comp[p, "H2O"]
                * 1e-3
                * pyunits.m**3
                / pyunits.kg  # Density of pure water
                * b.mass_frac_phase_comp[p, "H2O"]
                / (1 - b.mass_frac_phase_comp[p, "H2O"])
                * (1 - b.mass_frac_after_separation)
                / b.mass_frac_after_separation
            )

            return b.heat_separation_phase[p] == Expr_if(
                IF_=b.liquid_separation >= 0.5,
                THEN_=pyunits.convert(heat, to_units=pyunits.J / pyunits.s),
                ELSE_=0 * pyunits.J / pyunits.s,
            )

        self.eq_heat_separation_phase = Constraint(
            self.params.phase_list, rule=rule_heat_separation_phase
        )

    def _enth_mass_phase(self):
        self.enth_mass_phase = Var(
            self.params.phase_list,
            initialize=1e6,
            bounds=(1, 1e10),
            units=pyunits.J * pyunits.kg**-1,
            doc="Specific enthalpy",
        )

        def rule_enth_mass_phase(b, p):
            return b.enth_mass_phase[p] == b.cp_mass_phase[p] * (
                b.temperature - 273.15 * pyunits.K
            )

        self.eq_enth_mass_phase = Constraint(
            self.params.phase_list, rule=rule_enth_mass_phase
        )

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (
                sum(b.flow_mass_phase_comp["Liq", j] for j in b.params.component_list)
                * b.enth_mass_phase["Liq"]
                + b.heat_separation_phase["Liq"]
            )

        self.enth_flow = Expression(rule=rule_enth_flow)

    # -----------------------------------------------------------------------------
    # General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.enth_flow

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mass_phase_comp": self.flow_mass_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # these variables should have user input
        if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"], default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)

        if (
            iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "DrawSolution"])
            is None
        ):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "DrawSolution"],
                default=1e0,
                warning=True,
            )
            iscale.set_scaling_factor(
                self.flow_mass_phase_comp["Liq", "DrawSolution"], sf
            )

        if self.is_property_constructed("pressure_osm_phase"):
            if iscale.get_scaling_factor(self.pressure_osm_phase["Liq"]) is None:
                iscale.set_scaling_factor(
                    self.pressure_osm_phase["Liq"],
                    iscale.get_scaling_factor(self.pressure),
                )

        if self.is_property_constructed("mass_frac_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.mass_frac_phase_comp["Liq", j])
                    is None
                ):
                    if j == "DrawSolution":
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp["Liq", j]
                        ) / iscale.get_scaling_factor(
                            self.flow_mass_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], sf
                        )
                    elif j == "H2O":
                        iscale.set_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j], 1
                        )

        if self.is_property_constructed("flow_vol_phase"):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"]
            ) / iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol_phase, sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
                if (
                    iscale.get_scaling_factor(self.conc_mass_phase_comp["Liq", j])
                    is None
                ):
                    if j == "H2O":
                        # solvents typically have a mass fraction between 0.5-1
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j], sf_dens
                        )
                    elif j == "DrawSolution":
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j],
                            sf_dens
                            * iscale.get_scaling_factor(
                                self.mass_frac_phase_comp["Liq", j]
                            ),
                        )

        if self.is_property_constructed("enth_mass_phase"):
            iscale.set_scaling_factor(
                self.enth_mass_phase,
                iscale.get_scaling_factor(self.cp_mass_phase["Liq"])
                * iscale.get_scaling_factor(self.temperature),
            )

        if self.is_property_constructed("enth_flow"):
            iscale.set_scaling_factor(
                self.enth_flow,
                iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                * iscale.get_scaling_factor(self.enth_mass_phase["Liq"]),
            )

        transform_property_constraints(self)
