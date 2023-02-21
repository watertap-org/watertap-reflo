# Import Python libraries
import idaes.logger as idaeslog

from enum import Enum, auto

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Reals,
    NonNegativeReals,
    log,
    Var,
    Param,
    Set,
    Suffix,
    value,
    check_optimal_termination,
    units as pyunits,
)
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
from idaes.core.base.components import Solute, Solvent, Cation, Anion
from idaes.core.base.phases import AqueousPhase
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
from watertap.core.util.scaling import transform_property_constraints

# Set up logger
_log = idaeslog.getLogger(__name__)

__author__ = "Mukta Hardikar"

@declare_process_block_class("ChemSofteningParameterBlock")
class ChemSofteningParameterData(PhysicalParameterBlock):

    """
    Property Parameter Block Class

    Defines component lists, along with base units and constant
    parameters.
    """
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "solute_list",
        ConfigValue(domain=list, description="List of solute species names"),
    )

    CONFIG.declare(
        "mw_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of component names and molecular weight data",
        ),
    )

    CONFIG.declare(
        "charge", ConfigValue(default={}, domain=dict, description="Ion charge")
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = ChemSofteningStateBlock

        # phases
        self.Liq = AqueousPhase()

        # components
        self.H2O = Solvent()

        # list to hold all species (including water)
        self.component_set = Set()

        self.hardness_set = Set()

        # Add solutes to hardness set and component list
        for j in self.config.solute_list:
            if str(j) in ["Ca_2+", "Mg_2+"]:
                self.hardness_set.add(str(j))

            self.add_component(str(j), Solute())
            self.component_set.add(str(j))

        solute_list = self.config.solute_list

        #self.component_set.add("H2O")
        self.mw_comp = Param(
            # self.ion_set,
            self.component_set,
            mutable=True,
            default=18e-3,
            initialize=self.config.mw_data,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight",
        )

        self.charge_comp = Param(
            # self.ion_set,
            self.component_set,
            mutable=True,
            default=1,
            initialize=self.config.charge,
            units=pyunits.dimensionless,
            doc="Charge",
        )
      
        # for v in self.component_objects(Var):
        #     v.fix()

        # self.pKw = Param(
        #     initialize=14,
        #     units=pyunits.dimensionless,
        #     doc="pKw"
        # )
       
        # Define default value for mass density of solution
        self.dens_mass_default = 1000 * pyunits.kg / pyunits.m**3
        # Define default value for dynamic viscosity of solution
        self.visc_d_default = 0.001 * pyunits.kg / pyunits.m / pyunits.s

        # ---------------------------------------------------------------------
        # Set default scaling factors
        self.default_scaling_factor = {
            ("temperature"): 1e-3,
            ("pressure"): 1e-5,
            ("dens_mass"): 1e-3,
            ("visc_d"): 1e3,
        }

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
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "alkalinity": {"method": None},
                "pH": {"method": None},
                "pOH": {"method": "_pOH"},
                "pKw": {"method": "_pKw"},
                "flow_mass_comp": {"method": "_flow_mass_comp"},
                "temperature": {"method": "_temperature"},
                "pressure": {"method": "_pressure"},
                "dens_mass": {"method": "_dens_mass"},
                "visc_d": {"method": "_visc_d"},
                "mw_comp": {"method": "_mw_comp"},
                "conc_mass_caco3_comp": {"method": "_conc_mass_caco3_comp"},
                "total_hardness": {"method": "_total_hardness"},
            }
        )


class _ChemSofteningStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a whole, rather
    than individual elements of indexed Property Blocks.
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
                         chosen. Note that if this method is triggered through the control
                         volume, and if initial guesses were not provided at the unit model
                         level, the control volume passes the inlet values as initial guess.
                         The keys for the state_args dictionary are:
                         flow_mol_phase_comp : value to initialize phase component flows;
                         pressure : value at which to initialize pressure;
                         temperature : value at which to initialize temperature.
            outlvl : sets output level of initialization routine (default=idaeslog.NOTSET)
            optarg : solver options dictionary object (default=None)
            state_vars_fixed : Flag to denote if state vars have already
                               been fixed.
                               - True - states have already been fixed by the control volume
                               1D. Control volume 0D does not fix the state vars, so will be
                               False if this state block is used with 0D blocks.
                               - False - states have not been fixed. The state block will deal
                               with fixing/unfixing.
            solver : Solver object to use during initialization. If None
                     is provided, it will use the default solver for IDAES (default = None)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - state variables are not unfixed, and a dict of returned
                         containing flags for which states were fixed during initialization.
                         - False - state variables are unfixed after initialization by calling
                         the release_state method.

        Returns:
            If hold_states is True, returns a dict containing flags for which states were fixed
            during initialization.
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix state variables
        flags = fix_state_vars(self, state_args)

        # initialize vars calculated from state vars
        for k in self.keys():
            for j in self[k].params.component_list:
                if self[k].is_property_constructed("flow_mass_comp"):
                    if j == "H2O":
                        self[k].flow_mass_comp[j].set_value(
                            self[k].flow_vol * self[k].dens_mass
                        )
                    else:
                        self[k].flow_mass_comp[j].set_value(
                            self[k].flow_vol * self[k].conc_mass_comp[j]
                        )

        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                # print(f'DOF = {dof}')
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

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)
        init_log.info("State Released.")



@declare_process_block_class("ChemSofteningStateBlock", block_class=_ChemSofteningStateBlock)
class ChemSofteningStateBlockData(StateBlockData):
    def build(self):
        """Callable method for Block construction."""
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_vol = Var(
            initialize=1,
            domain=NonNegativeReals,
            doc="Volumetric flow rate",
            units=pyunits.m**3 / pyunits.s,
        )

        self.conc_mass_comp = Var(
            self.params.component_set,
            initialize=1,
            domain=NonNegativeReals,
            doc="Mass concentration of each solute",
            units=pyunits.kg / pyunits.m**3,
        )

        self.pH = Var(
            initialize=7,
            bounds=(0, 14),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="State pH",
        )

        self.alkalinity = Var(
            initialize=100,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.m**3,
            doc="State alkalinity",
        )

# -----------------------------------------------------------------------------
    # Property Methods
    def _flow_mass_comp(self):

        self.flow_mass_comp = Var(
            self.params.component_list,
            initialize=1,
            domain=NonNegativeReals,
            doc="Mass flowrate of each component",
            units=pyunits.kg / pyunits.s,
        )
        
        def rule_flow_mass_comp(b, j):
            if j == "H2O":
                return b.flow_mass_comp[j] == b.flow_vol * b.dens_mass
            else:
                return b.flow_mass_comp[j] == b.flow_vol * b.conc_mass_comp[j]

        self.eq_flow_mass_comp = Constraint(
            self.params.component_list, rule=rule_flow_mass_comp
        )
        
    def _pOH(self):
        self.pOH = Var(
            initialize=7,
            bounds=(0, 16),
            units=pyunits.dimensionless,
            doc="pOH",
        )
    
        def rule_pOH(b):
            return b.pOH == b.pKw - b.pH
        
        self.eq_pOH = Constraint(rule=rule_pOH)


    def _pKw(self):

        self.pKw = Var(
            initialize=14,
            # domain=(10, 16),
            units=pyunits.dimensionless,
            doc="pKw"
        )

        self.pKw_coeff_A = Param(
            initialize=4470.99
        )

        self.pKw_coeff_B = Param(
            initialize=0.017060
        )

        self.pKw_coeff_C = Param(
            initialize=6.0875
        )

        def rule_pKw(b):
            return b.pKw == b.pKw_coeff_A / b.temperature + b.pKw_coeff_B * b.temperature - b.pKw_coeff_C
        
        self.eq_pKw = Constraint(rule=rule_pKw)


    def _temperature(self):
        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 373.15),
            units=pyunits.K,
            doc="Temperature",
        )

    def _pressure(self):
        self.pressure = Var(
            initialize=101325,
            bounds=(1e5, None),
            units=pyunits.Pa,
            doc="Pressure",
        )

    def _dens_mass(self):
        self.dens_mass = Param(
            initialize=self.params.dens_mass_default,
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Mass density of flow",
        )

    def _visc_d(self):
        self.visc_d = Param(
            initialize=self.params.visc_d_default,
            units=pyunits.kg / pyunits.m / pyunits.s,
            mutable=True,
            doc="Dynamic viscosity of solution",
        )

    def _conc_mass_caco3_comp(self):
        self.conc_mass_caco3_comp = Var(
            self.params.hardness_set,
            initialize=100,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Mass concentration in CaCO3 equivalents",
        )

        self.equivalent_wt_caco3 = Param(
            initialize=0.050, 
            units=pyunits.kg/pyunits.mol, 
            doc="Equivalent weight of CaCO3"
            )

        def rule_conc_mass_caco3_comp(b, j):
            return b.conc_mass_caco3_comp[j] == b.equivalent_wt_caco3 / (b.params.mw_comp[j] / abs(b.params.charge_comp[j])) * b.conc_mass_comp[j]

        self.eq_conc_mass_caco3_comp = Constraint(
            self.params.hardness_set,
            rule=rule_conc_mass_caco3_comp,
        )


    def _total_hardness(self):
        self.total_hardness = Var(
            initialize=1e3,
            bounds=(5e2, 2e3),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass density",
        )
    def _total_hardness(self):
        def rule_total_hardness(b):
            return sum(b.conc_mass_caco3_comp[p] for p in b.params.hardness_set)

        self.total_hardness = Expression(rule=rule_total_hardness)

    def _mw_comp(self):
        add_object_reference(self, "mw_comp", self.params.mw_comp)

    
    
    def get_material_flow_terms(blk, p, j):
        return blk.flow_mass_comp[j]
        
    def get_enthalpy_flow_terms(blk, p):
        raise NotImplementedError

    def get_material_density_terms(blk, p, j):
        return blk.conc_mass_comp[j]

    def get_energy_density_terms(self, p):
        raise NotImplementedError

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.none

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    
    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}

    def define_display_vars(self):
        return {
            "Volumetric Flowrate": self.flow_vol,
            "Mass Concentration": self.conc_mass_comp,
            "Temperature": self.temperature,
        }


    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.flow_vol) is None:
            sf_Q = iscale.get_scaling_factor(self.flow_vol, default=1, warning=True)
            iscale.set_scaling_factor(self.flow_vol, sf_Q)

        for j, v in self.conc_mass_comp.items():
            sf_c = iscale.get_scaling_factor(self.conc_mass_comp[j])
            if sf_c is None:
                try:
                    sf_c = self.params.default_scaling_factor[("conc_mass_comp", j)]
                except KeyError:
                    iscale.get_scaling_factor(
                        self.conc_mass_comp[j], default=1, warning=True
                    )

        if self.is_property_constructed("flow_mass_comp"):
            for j, v in self.flow_mass_comp.items():
                if iscale.get_scaling_factor(v) is None:
                    if j == "H2O":
                        sf = value(self.flow_vol * self.dens_mass) ** -1
                    else:
                        sf = value(self.flow_vol * self.conc_mass_comp[j]) ** -1
                    iscale.set_scaling_factor(v, sf)