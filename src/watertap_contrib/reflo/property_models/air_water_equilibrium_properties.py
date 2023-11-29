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
Air-water equilibrium property package
"""

# Import Python libraries
import idaes.logger as idaeslog

from enum import Enum, auto
import itertools

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    Set,
    Suffix,
    NonNegativeReals,
    Var,
    Param,
    exp,
    log,
    log10,
    value,
    check_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.base.components import Solute, Solvent
from idaes.core.base.phases import (
    LiquidPhase,
    VaporPhase,
    PhaseType as PT,
)
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import extract_data
from idaes.core.solvers import get_solver
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
    PropertyPackageError,
)
import idaes.core.util.scaling as iscale

# from watertap.property_models.multicomp_aq_sol_prop_pack import DiffusivityCalculation
boltzmann = pyunits.convert(
    Constants.boltzmann_constant,
    to_units=(pyunits.g * pyunits.cm**2) / (pyunits.second**2 * pyunits.degK),
)
# Set up logger
_log = idaeslog.getLogger(__name__)

__author__ = "Kurban Sitterley"


class MolarVolumeCalculation(Enum):
    none = auto()
    TynCalus = auto()


class LiqDiffusivityCalculation(Enum):
    none = auto()
    HaydukLaudie = auto()


class VapDiffusivityCalculation(Enum):
    none = auto()
    WilkeLee = auto()


@declare_process_block_class("AirWaterEq")
class AirWaterEqData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "solute_list",
        ConfigValue(
            domain=list,
            description="Required argument. List of strings that specify names of solute species.",
        ),
    )

    CONFIG.declare(
        "mw_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Required argument. Dict of component names (keys)and molecular weight data (values)",
        ),
    )
    CONFIG.declare(
        "diffusivity_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of phase, solute species names (keys) and bulk ion diffusivity data (values)",
        ),
    )

    CONFIG.declare(
        "molar_volume_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and molar volume of aqueous species",
        ),
    )

    CONFIG.declare(
        "critical_molar_volume_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and critical molar volume of aqueous species. Used for Tyn-Calus method for calculating molar bolume",
        ),
    )

    CONFIG.declare(
        "density_data",
        ConfigValue(
            default={
                "Liq": 998.2,
                "Vap": 1.204,
            },  # default is for pure water and air at 20C
            domain=dict,
            description="Dict of phases and density of vapor and liquid phases",
        ),
    )

    CONFIG.declare(
        "dynamic_viscosity_data",
        ConfigValue(
            default={
                "Liq": 1e-3,
                "Vap": 1.813e-5,
            },  # default is for pure water and air at 20C
            domain=dict,
            description="Dict of phases and dynamic viscosity of vapor and liquid phases",
        ),
    )

    CONFIG.declare(
        "henry_constant_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and Henry's Law Constant",
        ),
    )

    CONFIG.declare(
        "temp_adjust_henry",
        ConfigValue(
            default=False,
            domain=bool,
            description="Flag to indicate if provided Henry's Law Constant should be adjusted for temperature.",
        ),
    )

    CONFIG.declare(
        "standard_enthalpy_change_data",
        ConfigValue(
            default={},  # default is for pure water and air at 20C
            domain=dict,
            description="Dict of solute species names and standard enthalpy change of dissolution in water. Used to temperature adjust Henry's Law Constant.",
        ),
    )

    CONFIG.declare(
        "temperature_boiling_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and boiling temperature data",
        ),
    )

    CONFIG.declare(
        "charge_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and ion charge",
        ),
    )

    CONFIG.declare(
        "liq_diffus_calculation",
        ConfigValue(
            default=LiqDiffusivityCalculation.HaydukLaudie,
            domain=In(LiqDiffusivityCalculation),
            description="Liquid diffusivity calculation flag",
            doc="""
           Options to account for ionic or molecular diffusivity.

           **default** - ``LiqDiffusivityCalculation.HaydukLaudie``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``LiqDiffusivityCalculation.none``", "Users provide data via the diffusivity_data configuration"
           "``LiqDiffusivityCalculation.HaydukLaudie``", "Allow the nonelectrolyte (neutral) species to get diffusivity from the Hayduk-Laudie equation"
       """,
        ),
    )

    CONFIG.declare(
        "vap_diffus_calculation",
        ConfigValue(
            default=VapDiffusivityCalculation.WilkeLee,
            domain=In(VapDiffusivityCalculation),
            description="Vapor diffusivity calculation flag",
            doc="""
           Options to account for ionic or molecular diffusivity.

           **default** - ``VapDiffusivityCalculation.WilkeLee``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``VapDiffusivityCalculation.none``", "Users provide data via the diffusivity_data configuration"
           "``VapDiffusivityCalculation.WilkeLee``", "Allow the nonelectrolyte (neutral) species to get diffusivity from the Wilke-Lee equation"
       """,
        ),
    )

    CONFIG.declare(
        "molar_volume_calculation",
        ConfigValue(
            default=MolarVolumeCalculation.TynCalus,
            domain=In(MolarVolumeCalculation),
            description="Molar volume calculation flag",
            doc="""
           Options to estimate molar volume for a component.

           **default** - ``MolarVolumeCalculation.TynCalus``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``MolarVolumeCalculation.none``", "Users provide data via the molar_volume_data configuration"
           "``MolarVolumeCalculation.TynCalus``", "Allow the component species to calculate molar volume from the Tyn-Calus equation"
       """,
        ),
    )

    def build(self):
        super().build()

        self._state_block_class = AirWaterEqStateBlock

        self.H2O = Solvent(valid_phase_types=[PT.liquidPhase, PT.vaporPhase])

        # Phases
        self.Liq = LiquidPhase()
        self.Vap = VaporPhase()

        # self.component_list = Set(dimen=1)
        self.solute_set = Set()

        for j in self.config.solute_list:
            self.add_component(j, Solute())

        self.phase_solute_set = self.phase_list * self.solute_set

        self.vap_comps = self.solute_set | Set(initialize=["Air"])
        self.vap_idx = list(itertools.product(["Vap"], self.vap_comps))
        self.vap_set = Set(initialize=self.vap_idx, doc="Set for all solutes plus air")

        mw_dict = {"H2O": 18e-3, "Air": 29e-3}
        mw_dict.update(self.config.mw_data)
        self.mw_comp = Param(
            mw_dict.keys() | self.component_list,
            mutable=False,
            initialize=mw_dict,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        self.dens_mass_phase = Var(
            self.phase_list,
            initialize=self.config.density_data,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density for each phase",
        )

        self.visc_d_phase = Var(
            self.phase_list,
            initialize=self.config.dynamic_viscosity_data,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity",
        )

        self.molar_volume_comp = Var(
            self.solute_set,
            initialize=self.config.molar_volume_data,
            units=pyunits.m**3 / pyunits.mol,
            doc="Molar volume of solutes",
        )

        self.critical_molar_volume_comp = Var(
            self.solute_set,
            initialize=self.config.critical_molar_volume_data,
            units=pyunits.m**3 / pyunits.mol,
            doc="Molar volume of solutes",
        )

        self.henry_constant_comp = Var(
            self.solute_set,
            initialize=self.config.henry_constant_data,
            units=pyunits.dimensionless,
            doc="Henry's constant",
        )

        if self.config.temp_adjust_henry:

            self.enth_change_dissolution_comp = Var(
                self.solute_set,
                initialize=self.config.standard_enthalpy_change_data,
                units=pyunits.joule / pyunits.mol,
                doc="Standard enthalpy change of dissolution in water for compound",
            )

        self.temperature_boiling_comp = Var(
            self.solute_set,
            initialize=self.config.temperature_boiling_data,
            units=pyunits.degK,
            doc="Boiling point temperature",
        )

        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)

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
        )

        obj.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                # "solubility_mass_phase_comp": {"method": "_solubility_mass_phase_comp"},
                # "solubility_mass_frac_phase_comp": {
                #     "method": "_solubility_mass_frac_phase_comp"
                # },
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                # "cp_mass_phase": {"method": "_cp_mass_phase"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_vol": {"method": "_flow_vol"},
                # "pressure_sat": {"method": "_pressure_sat"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
            }
        )

        obj.define_custom_properties(
            {
                "mw_comp": {"method": "_mw_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "debye_huckel_constant": {"method": "_debye_huckel_constant"},
                "ionic_strength_molal": {"method": "_ionic_strength_molal"},
                "molar_volume_comp": {"method": "_molar_volume_comp"},
                "critical_molar_volume_comp": {"method": "_critical_molar_volume_comp"},
                "energy_molecular_attraction_phase_comp": {
                    "method": "_energy_molecular_attraction_phase_comp"
                },
                "energy_molecular_attraction": {
                    "method": "_energy_molecular_attraction"
                },
                "collision_molecular_separation_phase_comp": {
                    "method": "_collision_molecular_separation_phase_comp"
                },
                "collision_molecular_separation": {
                    "method": "_collision_molecular_separation"
                },
                "temperature_boiling_comp": {"method": "_temperature_boiling_comp"},
                "henry_constant_comp": {"method": "_henry_constant_comp"},
                "collision_function_comp": {"method": "_collision_function_comp"},
                "collision_function_zeta_comp": {
                    "method": "_collision_function_zeta_comp"
                },
                "collision_function_ee_comp": {"method": "_collision_function_ee_comp"},
                "enth_change_dissolution_comp": {
                    "method": "_enth_change_dissolution_comp"
                },
                # "dh_vap_mass_solvent": {"method": "_dh_vap_mass_solvent"},
                # "cp_mass_solvent": {"method": "_cp_mass_solvent"},
                # "cp_mass_solute": {"method": "_cp_mass_solute"},
                # "temperature_sat_solvent": {"method": "_temperature_sat_solvent"},
                # "enth_mass_solvent": {"method": "_enth_mass_solvent"},
                # "enth_mass_solute": {"method": "_enth_mass_solute"},
                # "enth_flow": {"method": "_enth_flow"},
            }
        )


class _AirWaterEqStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Constraint on water concentration at outlet - unfix in these cases
        for b in self.values():
            if b.config.defined_state is False:
                b.conc_mol_comp["H2O"].unfix()

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
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:
                         flow_mass_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine (default=idaeslog.NOTSET)
            optarg : solver options dictionary object (default=None)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : Solver object to use during initialization if None is provided
                     it will use the default solver for IDAES (default = None)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states variables are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 release_state method
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        # Get loggers
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
                raise PropertyPackageError(
                    "\nWhile initializing {sb_name}, the degrees of freedom "
                    "are {dof}, when zero is required. \nInitialization assumes "
                    "that the state variables should be fixed and that no other "
                    "variables are fixed. \nIf other properties have a "
                    "predetermined value, use the calculate_state method "
                    "before using initialize to determine the values for "
                    "the state variables and avoid fixing the property variables."
                    "".format(sb_name=self.name, dof=dof)
                )

        # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        if not skip_solve:
            # Initialize properties
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            init_log.info_high(
                "Property initialization: {}.".format(idaeslog.condition(results))
            )

        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

        if (not skip_solve) and (not check_optimal_termination(results)):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please "
                f"check the output logs for more information."
            )

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialisation.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        # Unfix state variables
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        revert_state_vars(self, flags)
        init_log.info_high("{} State Released.".format(self.name))

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


@declare_process_block_class("AirWaterEqStateBlock", block_class=_AirWaterEqStateBlock)
class AirWaterEqStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.pressure = Var(
            domain=NonNegativeReals,
            initialize=101325,
            units=pyunits.Pa,
            doc="State pressure [Pa]",
        )

        self.temperature = Var(
            self.phase_list,
            domain=NonNegativeReals,
            initialize=298.15,
            bounds=(273.15, 373.15),
            units=pyunits.degK,
            doc="State temperature [K]",
        )

        self.flow_mass_phase_comp = Var(
            self.phase_component_set,
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

    # Property Methods

    def _mass_frac_phase_comp(self):
        self.mass_frac_phase_comp = Var(
            self.phase_component_set,
            domain=NonNegativeReals,
            initialize=0.5,
            bounds=(0, 1.0001),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        def rule_mass_frac_phase_comp(b, p, j):
            phase_comp_list = [
                (p, j)
                for j in self.params.component_list
                if (p, j) in b.phase_component_set
            ]
            if len(phase_comp_list) == 1:  # one component in this phase
                return b.mass_frac_phase_comp[p, j] == 1
            else:
                return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[
                    p, j
                ] / sum(b.flow_mass_phase_comp[p_j] for p_j in phase_comp_list)

        self.eq_mass_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mass_frac_phase_comp
        )

    def _diffus_phase_comp(self):

        self.diffus_phase_comp = Var(
            self.params.phase_solute_set,
            initialize=1e-9,
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Mass diffusivity of solute components in liquid and vapor",
        )
        self.hl_diffus_cont = Param(
            mutable=True,
            default=13.26e-9,
            initialize=13.26e-9,
            units=pyunits.dimensionless,
            doc="Hayduk-Laudie correlation constant",
        )
        self.hl_visc_coeff = Param(
            mutable=True,
            default=1.14,
            initialize=1.14,
            units=pyunits.dimensionless,
            doc="Hayduk-Laudie viscosity coefficient",
        )
        self.hl_molar_volume_coeff = Param(
            mutable=True,
            default=0.589,
            initialize=0.589,
            units=pyunits.dimensionless,
            doc="Hayduk-Laudie molar volume coefficient",
        )

        self.wilke_lee_param_A = Param(
            within=Reals,
            initialize=1.084,
            units=pyunits.cm**2 * pyunits.degK**-1.5,
            doc="Wilke-Lee parameter A",
        )

        self.wilke_lee_param_B = Param(
            initialize=0.249,
            units=pyunits.cm**2 * pyunits.degK**-1.5,
            doc="Wilke-Lee parameter B",
        )
        #  if (
        #             self.params.config.vap_diffus_calculation
        #             == VapDiffusivityCalculation.WilkeLee
        #         ):
        def rule_diffus_phase_comp(b, p, j):
            diffus_coeff_inv_units = pyunits.s * pyunits.m**-2
            visc_solvent_inv_units = pyunits.cP**-1
            molar_volume_inv_units = pyunits.mol * pyunits.cm**-3
            mw_air = pyunits.convert(b.mw_comp["Air"], to_units=pyunits.g / pyunits.mol)
            mw_j = pyunits.convert(b.mw_comp[j], to_units=pyunits.g / pyunits.mol)
            wilke_lee_denom_units = (
                pyunits.s**3 * pyunits.m * pyunits.kg**-1 * pyunits.nm**-2
            )
            if p == "Liq":
                if (
                    b.params.config.liq_diffus_calculation
                    == LiqDiffusivityCalculation.HaydukLaudie
                ):
                    return (b.diffus_phase_comp[p, j] * diffus_coeff_inv_units) * (
                        (
                            pyunits.convert(b.visc_d_phase[p], to_units=pyunits.cP)
                            * visc_solvent_inv_units
                        )
                        ** b.hl_visc_coeff
                    ) * (
                        (
                            pyunits.convert(
                                b.molar_volume_comp[j],
                                to_units=pyunits.cm**3 * pyunits.mol**-1,
                            )
                            * molar_volume_inv_units
                        )
                        ** b.hl_molar_volume_coeff
                    ) == b.hl_diffus_cont
                elif (
                    b.params.config.liq_diffus_calculation
                    == LiqDiffusivityCalculation.none
                ):
                    if (p, j) not in self.params.config.diffusivity_data.keys():
                        raise ConfigurationError(
                            f"\nThere is no {p} diffusivity provided for {j} in configuration argument 'diffusivity_data' for {b.config.parameters.name}.\n"
                            "Please provide that data or use LiqDiffusivityCalculation.HaydukLaudie configuration."
                        )
                    self.diffus_phase_comp[p, j].fix(
                        self.params.config.diffusivity_data[p, j]
                    )
                    return Constraint.Skip
            if p == "Vap":
                if (
                    self.params.config.vap_diffus_calculation
                    == VapDiffusivityCalculation.WilkeLee
                ):

                    sqrt_term = (1 / mw_j + 1 / mw_air) ** 0.5 * (
                        pyunits.g / pyunits.mol
                    ) ** 0.5  # units = dimensionless
                    numerator = (
                        (b.wilke_lee_param_A - b.wilke_lee_param_B * sqrt_term)
                        * (b.temperature[p] ** 1.5)
                        * sqrt_term
                    )  # units = cm2
                    denominator = (
                        b.pressure
                        * (b.collision_molecular_separation[j]) ** 2
                        * b.collision_function_comp[j]
                    ) * wilke_lee_denom_units  # units = s
                    diffusivity_cm2_s = numerator / denominator
                    # return (
                    #     b.diffus_phase_comp[p, j] * diffus_coeff_inv_units
                    # ) == numerator / denominator
                    # return b.diffus_phase_comp[p, j] == pyunits.convert(diffusivity_cm2_s, to_units=pyunits.m**2/pyunits.s)
                    return b.diffus_phase_comp[p, j] == pyunits.convert(
                        numerator / denominator, to_units=pyunits.m**2 / pyunits.s
                    )
                elif (
                    b.params.config.vap_diffus_calculation
                    == VapDiffusivityCalculation.none
                ):
                    if (p, j) not in self.params.config.diffusivity_data.keys():
                        raise ConfigurationError(
                            f"\nThere is no {p} diffusivity provided for {j} in configuration argument 'diffusivity_data' for {b.config.parameters.name}.\n"
                            "Please provide that data or use VapDiffusivityCalculation.WilkeLee configuration."
                        )
                    self.diffus_phase_comp[p, j].fix(
                        self.params.config.diffusivity_data[p, j]
                    )
                    return Constraint.Skip

        self.eq_diffus_phase_comp = Constraint(
            self.params.phase_solute_set,
            rule=rule_diffus_phase_comp,
        )

        # elif (
        #     self.params.config.liq_diffus_calculation == LiqDiffusivityCalculation.none
        # ):

        #     self.diffus_phase_comp["Liq"].fix(self.params.config.diffusivity_data["Liq"])

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Volumetric flow rate",
        )

        def rule_flow_vol_phase(b, p):
            if p == "Liq":
                return (
                    b.flow_vol_phase[p]
                    == sum(
                        b.flow_mass_phase_comp[p, j]
                        for j in self.params.component_list
                        if (p, j) in self.phase_component_set
                    )
                    / b.dens_mass_phase[p]
                )
            elif p == "Vap":
                return (
                    b.flow_vol_phase[p]
                    == sum(
                        b.flow_mass_phase_comp[p, j]
                        for j in self.params.component_list
                        if (p, j) in self.phase_component_set
                    )
                    / b.dens_mass_phase["Vap"]
                )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    # 12. Total volumetric flow rate
    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    # 13. Vapour pressure of the NaCl solution based on the boiling temperature
    def _pressure_sat(self):
        self.pressure_sat = Var(
            initialize=1e3,
            bounds=(0.001, 1e6),
            units=pyunits.Pa,
            doc="Vapor pressure of NaCl solution",
        )

        def rule_pressure_sat(b):  # vapor pressure, eq6 in Sparrow (2003)
            t = b.temperature - 273.15 * pyunits.K
            x = b.mass_frac_phase_comp["Liq", "NaCl"]

            ps_a = (
                b.params.pressure_sat_param_A1
                + (b.params.pressure_sat_param_A2 * x)
                + (b.params.pressure_sat_param_A3 * x**2)
                + (b.params.pressure_sat_param_A4 * x**3)
                + (b.params.pressure_sat_param_A5 * x**4)
            )

            ps_b = (
                b.params.pressure_sat_param_B1
                + (b.params.pressure_sat_param_B2 * x)
                + (b.params.pressure_sat_param_B3 * x**2)
                + (b.params.pressure_sat_param_B4 * x**3)
                + (b.params.pressure_sat_param_B5 * x**4)
            )

            ps_c = (
                b.params.pressure_sat_param_C1
                + (b.params.pressure_sat_param_C2 * x)
                + (b.params.pressure_sat_param_C3 * x**2)
                + (b.params.pressure_sat_param_C4 * x**3)
                + (b.params.pressure_sat_param_C5 * x**4)
            )

            ps_d = (
                b.params.pressure_sat_param_D1
                + (b.params.pressure_sat_param_D2 * x)
                + (b.params.pressure_sat_param_D3 * x**2)
                + (b.params.pressure_sat_param_D4 * x**3)
                + (b.params.pressure_sat_param_D5 * x**4)
            )

            ps_e = (
                b.params.pressure_sat_param_E1
                + (b.params.pressure_sat_param_E2 * x)
                + (b.params.pressure_sat_param_E3 * x**2)
                + (b.params.pressure_sat_param_E4 * x**3)
                + (b.params.pressure_sat_param_E5 * x**4)
            )

            p_sat = (
                ps_a + (ps_b * t) + (ps_c * t**2) + (ps_d * t**3) + (ps_e * t**4)
            )
            return b.pressure_sat == pyunits.convert(p_sat, to_units=pyunits.Pa)

        self.eq_pressure_sat = Constraint(rule=rule_pressure_sat)

    # 14. Saturation temperature for water vapour at calculated boiling pressure
    def _temperature_sat_solvent(self):
        self.temperature_sat_solvent = Var(
            initialize=298.15,
            bounds=(273.15, 1000.15),
            units=pyunits.K,
            doc="Vapour (saturation) temperature of pure solvent at boiling (i.e. crystallization) pressure",
        )

        def rule_temperature_sat_solvent(b):
            psat = pyunits.convert(b.pressure_sat, to_units=pyunits.kPa)
            return (
                b.temperature_sat_solvent
                == b.params.temp_sat_solvent_A1
                + b.params.temp_sat_solvent_A2
                / (
                    log(psat / b.params.temp_sat_solvent_A3)
                    + b.params.temp_sat_solvent_A4
                )
            )

        self.eq_temperature_sat_solvent = Constraint(rule=rule_temperature_sat_solvent)

    # 15. Mass concentration
    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            ["Liq"],
            self.params.component_list,
            initialize=10,
            bounds=(0, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, j):
            return (
                b.conc_mass_phase_comp["Liq", j]
                == b.dens_mass_phase["Liq"] * b.mass_frac_phase_comp["Liq", j]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.params.component_list, rule=rule_conc_mass_phase_comp
        )

    def _enth_mass_solvent(self):
        self.enth_mass_solvent = Var(
            ["Liq", "Vap"],
            initialize=1e3,
            bounds=(1, 1e4),
            units=pyunits.kJ * pyunits.kg**-1,
            doc="Specific saturated enthalpy of pure solvent",
        )

        def rule_enth_mass_solvent(b, p):
            t = b.temperature[p] - 273.15 * pyunits.K
            h_w = (
                b.params.enth_mass_solvent_param_A1
                + b.params.enth_mass_solvent_param_A2 * t
                + b.params.enth_mass_solvent_param_A3 * t**2
                + b.params.enth_mass_solvent_param_A4 * t**3
            )
            if p == "Liq":  # enthalpy, eq. 55 in Sharqawy
                return b.enth_mass_solvent[p] == pyunits.convert(
                    h_w, to_units=pyunits.kJ * pyunits.kg**-1
                )
            elif p == "Vap":

                return (
                    b.enth_mass_solvent[p]
                    == pyunits.convert(h_w, to_units=pyunits.kJ * pyunits.kg**-1)
                    + +b.dh_vap_mass_solvent
                )

        self.eq_enth_mass_solvent = Constraint(
            ["Liq", "Vap"], rule=rule_enth_mass_solvent
        )

    def _enth_flow(self):
        # enthalpy flow expression for get_enthalpy_flow_terms method

        def rule_enth_flow(b):  # enthalpy flow [J/s]
            return (
                sum(b.flow_mass_phase_comp["Liq", j] for j in b.params.component_list)
                * b.enth_mass_phase["Liq"]
                + b.flow_mass_phase_comp["Vap", "H2O"] * b.enth_mass_solvent["Vap"]
                + b.flow_mass_phase_comp["Sol", "NaCl"] * b.enth_mass_solute["Sol"]
            )

        self.enth_flow = Expression(rule=rule_enth_flow)

    def _flow_mol_phase_comp(self):
        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            initialize=100,
            bounds=(None, None),
            domain=NonNegativeReals,
            units=pyunits.mol / pyunits.s,
            doc="Molar flowrate",
        )

        def rule_flow_mol_phase_comp(b, p, j):
            return (
                b.flow_mol_phase_comp[p, j]
                == b.flow_mass_phase_comp[p, j] / b.params.mw_comp[j]
            )

        self.eq_flow_mol_phase_comp = Constraint(
            self.phase_component_set, rule=rule_flow_mol_phase_comp
        )

    def _mole_frac_phase_comp(self):
        self.mole_frac_phase_comp = Var(
            self.phase_component_set,
            initialize=0.1,
            bounds=(0, 1.0001),
            units=pyunits.dimensionless,
            doc="Mole fraction",
        )

        def rule_mole_frac_phase_comp(b, p, j):
            phase_comp_list = [
                (p, j)
                for j in self.params.component_list
                if (p, j) in b.phase_component_set
            ]
            if len(phase_comp_list) == 1:  # one component in this phase
                return b.mole_frac_phase_comp[p, j] == 1
            else:
                return b.mole_frac_phase_comp[p, j] == b.flow_mol_phase_comp[
                    p, j
                ] / sum(b.flow_mol_phase_comp[p_j] for (p_j) in phase_comp_list)

        self.eq_mole_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mole_frac_phase_comp
        )

    def _energy_molecular_attraction_phase_comp(self):

        self.energy_molecular_attraction_comp_param = Param(
            initialize=1.21,
            units=pyunits.dimensionless,
            doc="Energy of molecular attraction equation parameter for non-air component",
        )

        self.energy_molecular_attraction_air_param = Param(
            initialize=78.6,
            units=pyunits.degK,
            doc="Energy of molecular attraction equation parameter for air",
        )

        self.energy_molecular_attraction_phase_comp = Var(
            self.params.vap_set,
            within=NonNegativeReals,
            initialize=1,
            units=pyunits.erg,
            doc="Energy of molecular attraction for individual components",
        )

        # self.energy_molecular_attraction = Var(
        #     ["Air"],
        #     self.params.solute_set,
        #     within=NonNegativeReals,
        #     initialize=1,
        #     units=pyunits.erg,
        #     doc="Energy of molecular attraction between air and components",
        # )

        def rule_energy_molecular_attraction_phase_comp(b, p, j):
            # boltzmann_units = pyunits.degK / pyunits.joule
            # k = pyunits.convert(
            #     Constants.boltzmann_constant,
            #     to_units=(pyunits.g * pyunits.cm**2)
            #     / (pyunits.second**2 * pyunits.degK),
            # )
            # k = Constants.boltzmann_constant
            if j == "Air":
                return b.energy_molecular_attraction_phase_comp[
                    p, j
                ] == pyunits.convert(
                    boltzmann * b.energy_molecular_attraction_air_param,
                    to_units=pyunits.erg,
                )
            else:
                return b.energy_molecular_attraction_phase_comp[
                    p, j
                ] == pyunits.convert(
                    boltzmann
                    * b.energy_molecular_attraction_comp_param
                    * b.temperature_boiling_comp[j],
                    to_units=pyunits.erg,
                )

        self.eq_energy_molecular_attraction_phase_comp = Constraint(
            self.params.vap_set, rule=rule_energy_molecular_attraction_phase_comp
        )

    def _energy_molecular_attraction(self):
        def rule_energy_molecular_attraction(b, a, j):
            # k = Constants.boltzmann_constant
            return pyunits.convert(
                boltzmann
                * (
                    b.energy_molecular_attraction_phase_comp["Vap", j]
                    / boltzmann
                    * b.energy_molecular_attraction_phase_comp["Vap", "Air"]
                    / boltzmann
                )
                ** 0.5,
                to_units=pyunits.erg,
            )

        self.energy_molecular_attraction = Expression(
            ["Air"], self.params.solute_set, rule=rule_energy_molecular_attraction
        )

    def _collision_molecular_separation_phase_comp(self):

        self.collision_molecular_separation_param = Param(
            initialize=1.18,
            units=(pyunits.nm * pyunits.mol)
            / (pyunits.liter)
            * (pyunits.mol ** (1 / 3)),
            doc="Molecular separation at collision parameter",
        )
        self.collision_molecular_separation_exp = Param(
            initialize=1 / 3,
            units=pyunits.dimensionless,
            doc="Molecular separation at collision exponent",
        )

        self.collision_molecular_separation_phase_comp = Var(
            self.params.vap_set,
            within=NonNegativeReals,
            initialize=0.3711,  # default value is for Air
            units=pyunits.nanometer,
            doc="Molecular separation at collision for components",
        )

        def rule_collision_molecular_separation_phase_comp(b, p, j):
            if j == "Air":
                cms_air = 0.3711 * pyunits.nanometer
                return b.collision_molecular_separation_phase_comp[p, j] == cms_air
            else:
                molar_vol = pyunits.convert(
                    b.molar_volume_comp[j],
                    to_units=pyunits.liter / pyunits.mol,
                )
                return b.collision_molecular_separation_phase_comp[
                    p, j
                ] == b.collision_molecular_separation_param * molar_vol ** (
                    b.collision_molecular_separation_exp
                )

        self.eq_collision_molecular_separation_phase_comp = Constraint(
            self.params.vap_set, rule=rule_collision_molecular_separation_phase_comp
        )

    def _collision_molecular_separation(self):
        def rule_collision_molecular_separation(b, j):
            return 0.5 * (
                b.collision_molecular_separation_phase_comp["Vap", j]
                + b.collision_molecular_separation_phase_comp["Vap", "Air"]
            )

        self.collision_molecular_separation = Expression(
            self.params.solute_set, rule=rule_collision_molecular_separation
        )

    def _collision_function_comp(self):

        self.collision_function_comp = Var(
            self.params.solute_set,
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="Collision function",
        )

        def rule_collision_function_comp(b, j):
            return b.collision_function_comp[j] == 10 ** (
                b.collision_function_zeta_comp[j]
            )

        self.eq_collision_function_comp = Constraint(
            self.params.solute_set, rule=rule_collision_function_comp
        )

    def _collision_function_zeta_comp(self):

        self.collision_function_zeta_param_A = a = Param(
            initialize=-0.14329,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - A parameter",
        )

        self.collision_function_zeta_param_B = b_ = Param(
            initialize=-0.48343,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - B parameter",
        )

        self.collision_function_zeta_param_C = c = Param(
            initialize=0.1939,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - C parameter",
        )

        self.collision_function_zeta_param_D = d = Param(
            initialize=0.13612,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - D parameter",
        )

        self.collision_function_zeta_param_E = e = Param(
            initialize=-0.20578,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - E parameter",
        )

        self.collision_function_zeta_param_F = f = Param(
            initialize=0.083899,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - F parameter",
        )

        self.collision_function_zeta_param_G = g = Param(
            initialize=-0.011491,
            units=pyunits.dimensionless,
            doc="Collision function zeta equation - G parameter",
        )

        self.collision_function_zeta_comp = Var(
            self.params.solute_set,
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="Collision function zeta",
        )

        def rule_collision_function_zeta_comp(b, j):
            ee = b.collision_function_ee_comp[j]
            return b.collision_function_zeta_comp[j] == (
                a
                + b_ * ee
                + c * ee**2
                + d * ee**3
                + e * ee**4
                + f * ee**5
                + g * ee**6
            )

        self.eq_collision_function_zeta_comp = Constraint(
            self.params.solute_set, rule=rule_collision_function_zeta_comp
        )

    def _collision_function_ee_comp(self):

        self.collision_function_ee_comp = Var(
            self.params.solute_set,
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="Collision function ee",
        )

        def rule_collision_function_ee_comp(b, j):
            # k = Constants.boltzmann_constant
            t = b.temperature["Vap"]
            return b.collision_function_ee_comp[j] == log10(
                (boltzmann * t) / b.energy_molecular_attraction["Air", j]
            )

        self.eq_collision_function_ee_comp = Constraint(
            self.params.solute_set, rule=rule_collision_function_ee_comp
        )

    def _temperature_boiling_comp(self):
        add_object_reference(
            self, "temperature_boiling_comp", self.params.temperature_boiling_comp
        )

    def _molar_volume_comp(self):

        if (
            self.params.config.molar_volume_calculation
            == MolarVolumeCalculation.TynCalus
        ):
            self.molar_volume_comp = Var(
                self.params.solute_set,
                initialize=self.params.config.molar_volume_data,
                units=pyunits.m**3 / pyunits.mol,
                doc="Molar volume of solutes",
            )

            self.tyn_calus_param = Param(
                initialize=0.285,
                units=pyunits.dimensionless,
                doc="Tyn-Calus equation parameter",
            )

            self.tyn_calus_exponent = Param(
                initialize=1.048,
                units=pyunits.dimensionless,
                doc="Tyn-Calus equation exponent",
            )

            def rule_molar_volume_comp(b, j):
                return (
                    b.molar_volume_comp[j]
                    == b.tyn_calus_param
                    * b.critical_molar_volume_comp[j] ** b.tyn_calus_exponent
                )

            self.eq_molar_volume_comp = Constraint(
                self.params.solute_set, rule=rule_molar_volume_comp
            )

        else:
            add_object_reference(
                self, "molar_volume_comp", self.params.molar_volume_comp
            )

    def _henry_constant_comp(self):
        if self.params.config.temp_adjust_henry:

            add_object_reference(
                self, "henry_constant_std_comp", self.params.henry_constant_comp
            )  # assume the data provided to config.henry_constant_comp is @20C

            self.henry_constant_comp = Var(
                self.params.solute_set,
                initialize=self.params.config.henry_constant_data,
                units=pyunits.dimensionless,
                doc="Temperature adjusted dimensionless Henry's constant",
            )

            def rule_henry_constant_comp(b, j, doc="van't Hoff equation"):
                t0 = 298 * pyunits.degK
                exponential_term = pyunits.convert(
                    (-b.enth_change_dissolution_comp[j] / Constants.gas_constant)
                    * (1 / b.temperature["Liq"] - 1 / t0),
                    to_units=pyunits.dimensionless,
                )
                return b.henry_constant_comp[j] == b.henry_constant_std_comp[j] * exp(
                    exponential_term
                )

            self.eq_henry_constant_comp = Constraint(
                self.params.solute_set, rule=rule_henry_constant_comp
            )

        else:

            add_object_reference(
                self, "henry_constant_comp", self.params.henry_constant_comp
            )

    def _critical_molar_volume_comp(self):
        add_object_reference(
            self,
            "critical_molar_volume_comp",
            self.params.critical_molar_volume_comp,
        )

    def _mw_comp(self):
        add_object_reference(self, "mw_comp", self.params.mw_comp)

    def _dens_mass_phase(self):
        add_object_reference(self, "dens_mass_phase", self.params.dens_mass_phase)

    def _visc_d_phase(self):
        add_object_reference(self, "visc_d_phase", self.params.visc_d_phase)

    def _enth_change_dissolution_comp(self):

        add_object_reference(
            self,
            "enth_change_dissolution_comp",
            self.params.enth_change_dissolution_comp,
        )

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.flow_mass_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.enth_flow

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

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

        # default scaling factors have already been set with idaes.core.property_base.calculate_scaling_factors()
        # for the following variables: flow_mass_phase_comp, pressure, temperature, dens_mass_phase, enth_mass_phase

        # These variables should have user input
        # if iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"]) is None:
        #     sf = iscale.get_scaling_factor(
        #         self.flow_mass_phase_comp["Liq", "H2O"], default=1e0, warning=True
        #     )
        #     iscale.set_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"], sf)
