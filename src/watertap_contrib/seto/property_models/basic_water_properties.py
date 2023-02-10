###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains the general purpose property package for zero-order
unit models. Zero-order models do not track temperature and pressure, or any
form of energy flow.
"""
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
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError

from pyomo.environ import Expression, Param, PositiveReals, units as pyunits, Var, Constraint
from pyomo.common.config import ConfigValue

# Some more inforation about this module
__author__ = "Kurban Sitterley"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("BasicWaterParameterBlock")
class BasicWaterParameterBlockData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Defines component lists, along with base units and constant
    parameters.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "solute_list",
        ConfigValue(
            domain=list,
            description="List of solute species in the water source",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = BasicWaterStateBlock

        self.Liq = LiquidPhase()

        self.H2O = Solvent()

        # Get component set from database if provided
        self.component_set = None

        # Check definition of solute list
        solute_list = self.config.solute_list

        for j in solute_list:
            self.add_component(str(j), Solute())

        # Define default value for mass density of solution
        self.dens_mass_default = 1000 * pyunits.kg / pyunits.m**3
        # Define default value for dynamic viscosity of solution
        self.visc_d_default = 0.001 * pyunits.kg / pyunits.m / pyunits.s

        # ---------------------------------------------------------------------
        # Set default scaling factors
        self.default_scaling_factor = {
            ("flow_vol"): 1e3,
            ("conc_mass_comp"): 1e2,
            ("temperature"): 1e-3,
            ("mass_frac_comp"): 1e2,
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
                "flow_mass_comp": {"method": None},
                "temperature": {"method": "_temperature"},
                "flow_vol": {"method": "_flow_vol"},
                "conc_mass_comp": {"method": "_conc_mass_comp"},
                "mass_frac_comp": {"method": "_mass_frac_comp"},
                "dens_mass": {"method": "_dens_mass"},
                "visc_d": {"method": "_visc_d"},
            }
        )


class _BasicWaterStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        blk,
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
                     guess.The keys for the state_args dictionary are:

                     flow_mol_comp : value at which to initialize component
                                     flows (default=None)
                     pressure : value at which to initialize pressure
                                (default=None)
                     temperature : value at which to initialize temperature
                                  (default=None)
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
        # For now, there are no constraints in the property package, so only
        # fix state variables if required
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        init_log.info("Initialization Complete.")

        if hold_state is True:
            flags = fix_state_vars(blk, state_args)
            return flags
        else:
            return

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


@declare_process_block_class("BasicWaterStateBlock", block_class=_BasicWaterStateBlock)
class BasicWaterStateBlockData(StateBlockData):
    """
    General purpose StateBlock for Zero-Order unit models.
    """

    def build(self):
        super().build()

        # Create state variables
        self.flow_mass_comp = Var(
            self.params.component_list,
            initialize=1,
            domain=PositiveReals,
            doc="Mass flowrate of each component",
            units=pyunits.kg / pyunits.s,
        )

    # -------------------------------------------------------------------------
    # Other properties
    def _temperature(self):
        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 373.15),
            units=pyunits.K,
            doc="Temperature",
        )

    def _mass_frac_comp(self):
        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.5,
            bounds=(0, 1.001),
            units=pyunits.kg / pyunits.kg,
            doc="Mass fraction",
        )

        def rule_mass_frac_comp(b, j):
            return b.mass_frac_comp[j] == b.flow_mass_comp[j] / sum(
                b.flow_mass_comp[j] for j in self.params.component_list
            )

        self.eq_mass_frac_comp = Constraint(
            self.params.component_list,
            rule=rule_mass_frac_comp,
        )

    def _conc_mass_comp(self):
        self.conc_mass_comp = Var(
            self.params.component_list,
            initialize=1,
            bounds=(0, 2e3),
            units=pyunits.kg / pyunits.m**3,
            doc="Mass concentration",
        )
        def rule_conc_mass_comp(b, j):
            return (b.conc_mass_comp[j] == b.dens_mass * b.mass_frac_comp[j]
            )

        self.eq_conc_mass_comp = Constraint(self.component_list, rule=rule_conc_mass_comp)

    def _dens_mass(self):
        self.dens_mass = Param(
            initialize=self.params.dens_mass_default,
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Mass density of flow",
        )

    def _flow_vol(self):
        self.flow_vol = Expression(
            expr=sum(self.flow_mass_comp[j] for j in self.component_list)
            / self.dens_mass
        )

    def _visc_d(self):
        self.visc_d = Param(
            initialize=self.params.visc_d_default,
            units=pyunits.kg / pyunits.m / pyunits.s,
            mutable=True,
            doc="Dynamic viscosity of solution",
        )

    def get_material_flow_terms(blk, p, j):
        return blk.flow_mass_comp[j]

    def get_enthalpy_flow_terms(blk, p):
        raise NotImplementedError

    def get_material_density_terms(blk, p, j):
        return blk.conc_mass_comp[j]

    def get_energy_density_terms(blk, p):
        raise NotImplementedError

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.none

    def define_state_vars(blk):
        return {"flow_mass_comp": blk.flow_mass_comp}

    def define_display_vars(blk):
        return {
            "Volumetric Flowrate": blk.flow_vol,
            "Mass Concentration": blk.conc_mass_comp,
            "Temperature": blk.temperature
        }

    def get_material_flow_basis(blk):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        d_sf_Q = self.params.default_scaling_factor["flow_vol"]
        d_sf_c = self.params.default_scaling_factor["conc_mass_comp"]

        for j, v in self.flow_mass_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, d_sf_Q * d_sf_c)

        if self.is_property_constructed("flow_vol"):
            if iscale.get_scaling_factor(self.flow_vol) is None:
                iscale.set_scaling_factor(self.flow_vol, d_sf_Q)

        if self.is_property_constructed("conc_mass_comp"):
            for j, v in self.conc_mass_comp.items():
                sf_c = iscale.get_scaling_factor(self.conc_mass_comp[j])
                if sf_c is None:
                    try:
                        sf_c = self.params.default_scaling_factor[("conc_mass_comp", j)]
                    except KeyError:
                        sf_c = d_sf_c
                    iscale.set_scaling_factor(self.conc_mass_comp[j], sf_c)
