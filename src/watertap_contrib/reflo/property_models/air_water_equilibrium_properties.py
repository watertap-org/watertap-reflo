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

import itertools
from enum import Enum, auto

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
    log10,
    value,
    check_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

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
)
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.solvers import get_solver
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
)
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core.util.scaling import transform_property_constraints

boltzmann = pyunits.convert(
    Constants.boltzmann_constant,
    to_units=(pyunits.g * pyunits.cm**2) / (pyunits.second**2 * pyunits.degK),
)

# Set up logger
_log = idaeslog.getLogger(__name__)

__author__ = "Kurban Sitterley"


"""
REFERENCES: 

Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). 
Chapter 7, 14. MWH's Water Treatment: Principles and Design (3rd ed.). doi:10.1002/9781118131473

Aniceto, J. P. S., Zêzere, B., & Silva, C. M. (2021).
Predictive Models for the Binary Diffusion Coefficient at Infinite Dilution in Polar and Nonpolar Fluids. 
Materials (Basel), 14(3). doi.org/10.3390/ma14030542

Wilke, C. R., & Lee, C. Y. (2002).
Estimation of Diffusion Coefficients for Gases and Vapors.
Industrial & Engineering Chemistry, 47(6), 1253-1257. doi:10.1021/ie50546a056

Huang, J. (2018).
A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice.
Journal of Applied Meteorology and Climatology, 57(6), 1265-1272. doi:10.1175/jamc-d-17-0334.1

Sharqawy, M. H., Lienhard, J. H. & Zubair S. M. (2012).
Thermophysical properties of seawater: a review of existing correlations and data
Desalination and Water Treatment, 16(1-3), 354-380. doi: 10.5004/dwt.2010.1079

Buck, A. L. (1981).
New equations for computing vapor pressure and enhancement factor.
J. Appl. Meteorol., 20 (12): 1527-1532. 
https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2

Antoine, C. (1888).
"Tensions des vapeurs; nouvelle relation entre les tensions et les températures"
Comptes Rendus des Séances de l'Académie des Sciences (in French), 107: 681-684, 778-780, 836-837.
https://en.wikipedia.org/wiki/Antoine_equation

"""


class MolarVolumeCalculation(Enum):
    """
    Approach to determine component molar volume. TynCalus is default.
    """

    none = auto()
    TynCalus = auto()


class LiqDiffusivityCalculation(Enum):
    """
    Approach to determine component liquid diffusivity. HaydukLaudie is default.
    """

    none = auto()
    HaydukLaudie = auto()


class VapDiffusivityCalculation(Enum):
    """
    Approach to determine component vapor diffusivity. WilkeLee is default.
    """

    none = auto()
    WilkeLee = auto()


class SaturationVaporPressureCalculation(Enum):
    """
    Approach to determine saturation vapor presssure. ArdenBuck is default.
    """

    none = auto()
    ArdenBuck = auto()
    Huang = auto()


class VaporPressureCalculation(Enum):
    """
    Approach to determine vapor presssure. FromRelativeHumidity is default.
    """

    none = auto()
    FromRelativeHumidity = auto()
    Antoine = auto()


class RelativeHumidityCalculation(Enum):
    """
    Approach to determine relative humidity. none (user-provided) is default.
    """

    none = auto()
    VaporPressureRatio = auto()


class LatentHeatVaporizationCalculation(Enum):
    """
    Approach to determine latent heat of vaporization. Sharqawy is default.
    """

    none = auto()
    Sharqawy = auto()


class SpecificHeatWaterCalculation(Enum):
    """
    Approach to determine specific heat of water in liquid and vapor phase. Sharqawy is default.
    """

    none = auto()
    Sharqawy = auto()


@declare_process_block_class("AirWaterEq")
class AirWaterEqData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "volatile_solute_list",  # volatile components; components that can be in liquid or vapor phase
        ConfigValue(
            default=[],
            domain=list,
            description="List of strings that specify names of volatile solute species.",
        ),
    )

    CONFIG.declare(
        "non_volatile_solute_list",  # non-volatile components; components that can only be in liquid phase
        ConfigValue(
            default=[],
            domain=list,
            description="List of strings that specify names of non-volatile solute species.",
        ),
    )

    CONFIG.declare(
        "mw_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Required argument. Dict of component names (keys) and molecular weight data (values)",
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
        "pressure_vap_sat_data",
        ConfigValue(
            default=0,
            domain=float,
            description="User provided value for saturation vapor pressure of water in air in Pa",
        ),
    )

    CONFIG.declare(
        "pressure_vap_data",
        ConfigValue(
            default=0,
            domain=float,
            description="User provided value for vapor pressure of water in air in Pa",
        ),
    )

    CONFIG.declare(
        "relative_humidity_data",
        ConfigValue(
            default=0,
            domain=float,
            description="User provided value for relative humidity of air",
        ),
    )

    CONFIG.declare(
        "latent_heat_of_vaporization_data",
        ConfigValue(
            default=0,
            domain=float,
            description="User provided value for latent heat of vaporization in kJ / kg",
        ),
    )

    CONFIG.declare(
        "specific_heat_of_water_data",
        ConfigValue(
            default={
                "Liq": 4184,
                "Vap": 1996,
            },
            domain=dict,
            description="Dict for specific heat of water for liquid and vapor phases in J / kg / K",
        ),
    )

    CONFIG.declare(
        "critical_molar_volume_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict of solute species names and critical molar volume of aqueous species."
            "Used for Tyn-Calus method for calculating molar volume",
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
            default=True,
            domain=bool,
            description="Flag to indicate if provided Henry's Law Constant should be adjusted for temperature.",
        ),
    )

    CONFIG.declare(
        "standard_enthalpy_change_data",
        ConfigValue(
            default={},
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

           "``LiqDiffusivityCalculation.none``", "Users provide liquid diffusivity data via the diffusivity_data configuration"
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

           "``VapDiffusivityCalculation.none``", "Users provide vapor diffusivity data via the diffusivity_data configuration"
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

    CONFIG.declare(
        "saturation_vapor_pressure_calculation",
        ConfigValue(
            default=SaturationVaporPressureCalculation.ArdenBuck,
            domain=In(SaturationVaporPressureCalculation),
            description="Saturation vapor pressure calculation flag",
            doc="""
           Options to estimate saturation vapor pressure for water in air.

           **default** - ``SaturationVaporPressureCalculation.ArdenBuck``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``SaturationVaporPressureCalculation.none``", "Users provide data via the pressure_vap_sat_data configuration"
           "``SaturationVaporPressureCalculation.ArdenBuck``", "Saturation vapor pressure is calcualted from Arden Buck equation"
           "``SaturationVaporPressureCalculation.Huang``", "Saturation vapor pressure is calcualted from Huang equation"
       """,
        ),
    )

    CONFIG.declare(
        "vapor_pressure_calculation",
        ConfigValue(
            default=VaporPressureCalculation.FromRelativeHumidity,
            domain=In(VaporPressureCalculation),
            description="Vapor pressure calculation flag",
            doc="""
           Options to estimate vapor pressure for water in air.

           **default** - ``VaporPressureCalculation.FromRelativeHumidity``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``VaporPressureCalculation.none``", "Users provide data via the pressure_vap_data configuration"
           "``VaporPressureCalculation.FromRelativeHumidity``", "Vapor pressure is calcualted from relative humidity and saturation vapor pressure"
           "``VaporPressureCalculation.Antoine``", "Vapor pressure is calcualted from Antoine equation"
       """,
        ),
    )

    CONFIG.declare(
        "relative_humidity_calculation",
        ConfigValue(
            default=RelativeHumidityCalculation.none,
            domain=In(RelativeHumidityCalculation),
            description="Relative humidity calculation flag",
            doc="""
           Options to estimate relative humidity of air-water system.

           **default** - ``RelativeHumidityCalculation.none``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``RelativeHumidityCalculation.none``", "Users provide data via the relative_humidity_data configuration"
           "``RelativeHumidityCalculation.VaporPressureRatio``", "Relative humidity is calcualted as ratio of vapor pressure to saturation vapor pressure"
       """,
        ),
    )

    CONFIG.declare(
        "latent_heat_of_vaporization_calculation",
        ConfigValue(
            default=LatentHeatVaporizationCalculation.Sharqawy,
            domain=In(LatentHeatVaporizationCalculation),
            description="Latent heat of vaporization calculation flag",
            doc="""
           Options to estimate latent heat of vaporization.

           **default** - ``LatentHeatVaporizationCalculation.Sharqawy``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``LatentHeatVaporizationCalculation.none``", "Users provide data via the latent_heat_of_vaporization_data configuration"
           "``LatentHeatVaporizationCalculation.Sharqawy``", "Latent heat of vaporization is calcualted according to Sharqawy (2010)"
       """,
        ),
    )

    CONFIG.declare(
        "specific_heat_of_water_calculation",
        ConfigValue(
            default=SpecificHeatWaterCalculation.Sharqawy,
            domain=In(SpecificHeatWaterCalculation),
            description="Specific heat of water calculation flag",
            doc="""
           Options to estimate specific heat of water.

           **default** - ``SpecificHeatWaterCalculation.Sharqawy``

       .. csv-table::
           :header: "Configuration Options", "Description"

           "``SpecificHeatWaterCalculation.none``", "Users provide data via the specific_heat_of_water_data configuration"
           "``SpecificHeatWaterCalculation.Sharqawy``", "Specific heat of water is calcualted according to Sharqawy (2010)"
       """,
        ),
    )

    def build(self):
        super().build()

        self._state_block_class = AirWaterEqStateBlock

        # Nine Sets are created in total:
        # In liquid phase only
        # 1. All components in Liq phase indexed by component; includes H2O
        # 2. All components in Liq phase indexed by ('Liq', component); includes H2O
        # 3. All non-volatile components in Liq phase indexed by component
        # 4. All non-volatile components in Liq phase indexed by ('Liq', component)

        # In vapor phase only
        # 5. All components in Vap phase indexed by component
        # 6. All components in Vap phase indexed by ('Vap', component)

        # In liquid and vapor phase
        # 7. All volatile components in both phases indexed by component
        # 8. All volatile components in both phases indexed by ('Phase', component)

        # 9. All components in either phase

        # Liq phase indexes and Sets
        self.liq_comp_list = (
            self.config.volatile_solute_list
            + self.config.non_volatile_solute_list
            + ["H2O"]
        )  # assumed all volatile and non-volatile solutes can be in liquid phase except Air
        liq_phase_comps_idx = list(itertools.product(["Liq"], self.liq_comp_list))
        non_volatile_phase_comps_idx = list(
            itertools.product(["Liq"], self.config.non_volatile_solute_list)
        )
        self.liq_comps = Set(
            initialize=self.liq_comp_list, doc="Set for all components in liquid phase"
        )  # component; includes H2O
        self.liq_comp_set = Set(
            initialize=liq_phase_comps_idx,
            doc="Set for liquid phase components with phase",
        )  # (Liq, component); includes H2O
        self.non_volatile_comps = Set(
            initialize=self.config.non_volatile_solute_list,
            doc="Set for all non-volatile components",
        )  # component
        self.non_volatile_comp_set = Set(
            initialize=non_volatile_phase_comps_idx,
            doc="Set for non-volatile components with phase",
        )  # (Liq, non_volatile_component)

        # Vap phase indexes and Sets
        self.vap_comp_list = self.config.volatile_solute_list + ["H2O", "Air"]
        vap_phase_comps_idx = list(itertools.product(["Vap"], self.vap_comp_list))
        volatile_phase_comps_idx = list(
            itertools.product(["Vap", "Liq"], self.config.volatile_solute_list)
        )
        self.vap_comps = Set(
            initialize=self.vap_comp_list, doc="Set for all components in vapor phase"
        )  # component; includes H2O + Air
        self.vap_comp_set = Set(
            initialize=vap_phase_comps_idx,
            doc="Set for vapor phase components with phase",
        )  # (Vap, component); includes H2O + Air
        self.volatile_comps = Set(
            initialize=self.config.volatile_solute_list,
            doc="Set for all volatile components",
        )
        self.volatile_comp_set = Set(
            initialize=volatile_phase_comps_idx,
            doc="Set for volatile components with phase",
        )

        self.component_set = Set(
            initialize=self.vap_comps | self.liq_comps,
            doc="Set for all components in package",
        )

        for j in self.component_set:
            self.add_component(j, Solute())

        self.Liq = LiquidPhase(component_list=self.liq_comps)
        self.Vap = VaporPhase(component_list=self.vap_comps)

        # Dict to store all components for each phase
        self.phase_comp_dict = {"Liq": self.liq_comp_list, "Vap": self.vap_comp_list}

        # Set containing both phases and solutes
        self.phase_solute_set = self.liq_comp_set | self.vap_comp_set

        mw_dict = {"H2O": 18e-3, "Air": 29e-3}
        mw_dict.update(self.config.mw_data)
        self.mw_comp = Param(
            mw_dict.keys() | self.component_set,
            mutable=True,
            initialize=mw_dict,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight kg/mol",
        )

        self.cp_air = Var(
            initialize=1013,
            units=pyunits.joule / (pyunits.kg * pyunits.degK),
            doc="Specific heat capacity of dry air",  # doi:10.1029/2008JD010174
        )

        self.cp_mass_solvent = Var(
            self.phase_list,
            initialize=self.config.specific_heat_of_water_data,
            units=pyunits.joule / (pyunits.kg * pyunits.degK),
            doc="Specific heat capacity of water in liquid and vapor phase",
        )

        self.visc_d_phase = Var(
            self.phase_list,
            initialize=self.config.dynamic_viscosity_data,
            units=pyunits.Pa * pyunits.s,
            doc="Dynamic viscosity",
        )

        self.dens_mass_phase = Var(
            self.phase_list,
            initialize=self.config.density_data,
            units=pyunits.kg / pyunits.m**3,
            doc="Density of liquid and vapor phases",
        )

        self.dh_vap_mass_solvent = Var(
            initialize=self.config.latent_heat_of_vaporization_data,
            units=pyunits.kJ / pyunits.kg,
            doc="Latent heat of vaporization of water",
        )

        self.molar_volume_comp = Var(
            self.volatile_comps,
            initialize=self.config.molar_volume_data,
            units=pyunits.m**3 / pyunits.mol,
            doc="Molar volume of solutes",
        )

        self.molar_volume_comp_crit = Var(
            self.volatile_comps,
            initialize=self.config.critical_molar_volume_data,
            units=pyunits.m**3 / pyunits.mol,
            doc="Critical molar volume of solutes",
        )

        self.henry_comp = Var(
            self.volatile_comps,
            initialize=self.config.henry_constant_data,
            units=pyunits.dimensionless,
            doc="Henry's constant",
        )

        self.pressure_vap_sat = Var(
            ["H2O"],
            initialize=self.config.pressure_vap_sat_data,
            units=pyunits.Pa,
            doc="User provided saturation vapor pressure for water",
        )

        self.pressure_vap = Var(
            ["H2O"],
            initialize=self.config.pressure_vap_data,
            units=pyunits.Pa,
            doc="User provided vapor pressure for water",
        )

        self.relative_humidity = Var(
            ["H2O"],
            initialize=self.config.relative_humidity_data,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="User provided relative humidity",
        )

        if self.config.temp_adjust_henry:

            self.enth_change_dissolution_comp = Var(
                self.volatile_comps,
                initialize=self.config.standard_enthalpy_change_data,
                units=pyunits.joule / pyunits.mol,
                doc="Standard enthalpy change of dissolution in water for compound",
            )

        self.temperature_boiling_comp = Var(
            self.volatile_comps,
            initialize=self.config.temperature_boiling_data,
            units=pyunits.degK,
            doc="Boiling point temperature",
        )

        for v in self.component_objects(Var):
            v.fix()

        # ---default scaling---
        self.set_default_scaling("pressure", 1e-5)
        self.set_default_scaling("temperature", 1e-2, index="Liq")
        self.set_default_scaling("temperature", 1e-2, index="Vap")
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1, index="Vap")
        self.set_default_scaling("visc_d_phase", 1e3, index="Liq")
        self.set_default_scaling("visc_d_phase", 1e5, index="Vap")
        self.set_default_scaling("diffus_phase_comp", 1e10, index="Liq")
        self.set_default_scaling("diffus_phase_comp", 1e6, index="Vap")

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
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
                "conc_mol_phase_comp": {"method": "_conc_mol_phase_comp"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_mass_phase": {"method": "_flow_mass_phase"},
                "flow_vol": {"method": "_flow_vol"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "mw_comp": {"method": "_mw_comp"},
            }
        )

        obj.define_custom_properties(
            {
                "henry_comp": {"method": "_henry_comp"},
                "molar_volume_comp": {"method": "_molar_volume_comp"},
                "pressure_vap": {"method": "_pressure_vap"},
                "pressure_vap_sat": {"method": "_pressure_vap_sat"},
                "relative_humidity": {"method": "_relative_humidity"},
                "molar_volume_comp_crit": {"method": "_molar_volume_comp_crit"},
                "enth_change_dissolution_comp": {
                    "method": "_enth_change_dissolution_comp"
                },
                "energy_molecular_attraction_comp": {
                    "method": "_energy_molecular_attraction_comp"
                },
                "energy_molecular_attraction": {
                    "method": "_energy_molecular_attraction"
                },
                "collision_molecular_separation_comp": {
                    "method": "_collision_molecular_separation_comp"
                },
                "collision_molecular_separation": {
                    "method": "_collision_molecular_separation"
                },
                "collision_function_comp": {"method": "_collision_function_comp"},
                "collision_function_zeta_comp": {
                    "method": "_collision_function_zeta_comp"
                },
                "collision_function_ee_comp": {"method": "_collision_function_ee_comp"},
                "temperature_boiling_comp": {"method": "_temperature_boiling_comp"},
                "dh_vap_mass_solvent": {"method": "_dh_vap_mass_solvent"},
                "cp_mass_solvent": {"method": "_cp_mass_solvent"},
                "cp_air": {"method": "_cp_air"},
            }
        )


class _AirWaterEqStateBlock(StateBlock):
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
            b = self[k]
            for j in b.params.volatile_comps:
                if b.is_property_constructed("molar_volume_comp"):
                    if (
                        b.params.config.molar_volume_calculation
                        == MolarVolumeCalculation.TynCalus
                    ):
                        calculate_variable_from_constraint(
                            b.molar_volume_comp[j], b.eq_molar_volume_comp[j]
                        )
                    else:
                        b.molar_volume_comp[j].set_value(
                            b.params.config.molar_volume_data[j]
                        )
                if b.is_property_constructed("henry_comp"):
                    if b.params.config.temp_adjust_henry:
                        t0 = 298 * pyunits.degK
                        exponential_term = pyunits.convert(
                            (
                                -b.enth_change_dissolution_comp[j]
                                / Constants.gas_constant
                            )
                            * (1 / b.temperature["Liq"] - 1 / t0),
                            to_units=pyunits.dimensionless,
                        )
                        b.henry_comp[j].set_value(
                            b.henry_comp_ref[j] * exp(exponential_term)
                        )
                    else:
                        b.henry_comp[j].set_value(
                            b.params.config.henry_constant_data[j]
                        )

            for p, j in b.params.volatile_comp_set:
                if b.is_property_constructed("diffus_phase_comp"):
                    if p == "Liq":
                        if (
                            b.params.config.liq_diffus_calculation
                            == LiqDiffusivityCalculation.HaydukLaudie
                        ):
                            liq_diff = value(
                                b.hl_diffus_cont
                                / (
                                    (
                                        (
                                            pyunits.convert(
                                                b.visc_d_phase[p], to_units=pyunits.cP
                                            )
                                            * pyunits.cP**-1
                                        )
                                        ** b.hl_visc_coeff
                                    )
                                    * (
                                        (
                                            pyunits.convert(
                                                b.molar_volume_comp[j],
                                                to_units=pyunits.cm**3
                                                * pyunits.mol**-1,
                                            )
                                            * (pyunits.mol * pyunits.cm**-3)
                                        )
                                        ** b.hl_molar_volume_coeff
                                    )
                                )
                            )
                            b.diffus_phase_comp[p, j].set_value(
                                liq_diff * pyunits.m**2 / pyunits.s
                            )
                        else:
                            b.diffus_phase_comp[p, j].set_value(
                                b.params.config.diffusivity_data[(p, j)]
                            )
                    if p == "Vap":
                        if (
                            b.params.config.vap_diffus_calculation
                            == VapDiffusivityCalculation.WilkeLee
                        ):
                            calculate_variable_from_constraint(
                                b.diffus_phase_comp[p, j],
                                b.eq_diffus_phase_comp[p, j],
                            )

            for p in b.phase_list:
                for j in b.params.phase_comp_dict[p]:
                    if b.is_property_constructed("flow_mol_phase_comp"):
                        b.flow_mol_phase_comp[p, j].set_value(
                            b.flow_mass_phase_comp[p, j] / b.params.mw_comp[j]
                        )
                    if b.is_property_constructed("mass_frac_phase_comp"):
                        b.mass_frac_phase_comp[p, j].set_value(
                            b.flow_mass_phase_comp[p, j]
                            / sum(
                                b.flow_mass_phase_comp[p, _j]
                                for _j in b.params.phase_comp_dict[p]
                            )
                        )
                    if b.is_property_constructed("conc_mass_phase_comp"):
                        b.conc_mass_phase_comp[p, j].set_value(
                            b.dens_mass_phase[p] * b.mass_frac_phase_comp[p, j]
                        )

                    if b.is_property_constructed("conc_mol_phase_comp"):
                        b.conc_mol_phase_comp[p, j].set_value(
                            b.conc_mass_phase_comp[p, j] / b.params.mw_comp[j]
                        )

                    if b.is_property_constructed("mol_frac_phase_comp"):
                        b.mol_frac_phase_comp[p, j].set_value(
                            b.flow_mol_phase_comp[p, j]
                            / sum(
                                b.flow_mol_phase_comp[p, _j]
                                for _j in b.params.phase_comp_dict[p]
                            )
                        )

                if b.is_property_constructed("flow_vol_phase"):
                    b.flow_vol_phase[p].set_value(
                        sum(
                            b.flow_mass_phase_comp[p, _j]
                            for _j in b.params.phase_comp_dict[p]
                        )
                        / b.dens_mass_phase[p]
                    )

                if b.is_property_constructed("flow_mass_phase"):
                    b.flow_mass_phase[p].set_value(
                        b.flow_vol_phase[p] * b.dens_mass_phase[p]
                    )

                if b.is_property_constructed("cp_mass_solvent"):
                    if (
                        b.params.config.specific_heat_of_water_calculation
                        == SpecificHeatWaterCalculation.Sharqawy
                    ):
                        calculate_variable_from_constraint(
                            b.cp_mass_solvent[p], b.eq_cp_mass_solvent[p]
                        )
                    else:
                        b.cp_mass_solvent[p].set_value(
                            b.params.config.specific_heat_of_water_data[p]
                        )

            if b.is_property_constructed("dh_vap_mass_solvent"):
                if (
                    b.params.config.latent_heat_of_vaporization_calculation
                    == LatentHeatVaporizationCalculation.Sharqawy
                ):
                    calculate_variable_from_constraint(
                        b.dh_vap_mass_solvent, b.eq_dh_vap_mass_solvent
                    )

            if b.is_property_constructed("pressure_vap_sat"):
                if (
                    b.params.config.saturation_vapor_pressure_calculation
                    != SaturationVaporPressureCalculation.none
                ):
                    calculate_variable_from_constraint(
                        b.pressure_vap_sat["H2O"],
                        b.eq_pressure_vap_sat["H2O"],
                    )
                else:
                    b.pressure_vap_sat["H2O"].set_value(
                        b.params.config.pressure_vap_sat_data
                    )

            if b.is_property_constructed("pressure_vap"):
                if (
                    b.params.config.vapor_pressure_calculation
                    != VaporPressureCalculation.none
                ):
                    if (
                        b.params.config.vapor_pressure_calculation
                        == VaporPressureCalculation.FromRelativeHumidity
                    ):
                        b.pressure_vap["H2O"].set_value(
                            b.relative_humidity["H2O"] * b.pressure_vap_sat["H2O"]
                        )
                    if (
                        b.params.config.vapor_pressure_calculation
                        == VaporPressureCalculation.Antoine
                    ):
                        calculate_variable_from_constraint(
                            b.pressure_vap["H2O"],
                            b.eq_pressure_vap["H2O"],
                        )
                else:
                    b.pressure_vap["H2O"].set_value(b.params.config.pressure_vap_data)

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
                f"Property initialization: {idaeslog.condition(results)}."
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
                        f"While using the calculate_state method on {sb.name}, variable {v_name} "
                        "was provided as an argument in var_args, but it does not have a scaling "
                        "factor. This suggests that the calculate_scaling_factor method has not been "
                        "used or the variable was created on demand after the scaling factors were "
                        "calculated. It is recommended to touch all relevant variables (i.e. call "
                        "them or set an initial value) before using the calculate_scaling_factor "
                        "method."
                    )
                if var[ind].is_fixed():
                    flags[(k, v_name, ind)] = True
                    if value(var[ind]) != val:
                        raise ConfigurationError(
                            f"While using the calculate_state method on {sb.name}, {var.name} was "
                            f"fixed to a value {val}, but it was already fixed to value {value(var[ind])}. "
                            f"Unfix the variable before calling the calculate_state "
                            "method or update var_args."
                            ""
                        )

                else:
                    flags[(k, v_name, ind)] = False
                    var[ind].fix(val)

            if degrees_of_freedom(sb) != 0:
                raise RuntimeError(
                    f"While using the calculate_state method on {sb.name}, the degrees "
                    f"of freedom were {degrees_of_freedom(sb)}, but 0 is required. Check var_args and ensure "
                    "the correct fixed variables are provided."
                )

        # Solve
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            solve_log.info_high(f"Calculate state: {idaeslog.condition(results)}.")

        if not check_optimal_termination(results):
            _log.warning(
                f"While using the calculate_state method on {self.name}, the solver failed "
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
    def _dens_mass_phase(self):

        add_object_reference(self, "dens_mass_phase", self.params.dens_mass_phase)

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
                (p, _j)
                for _j in self.params.phase_comp_dict[p]
                if (p, _j) in self.phase_component_set
            ]
            return b.mass_frac_phase_comp[p, j] == b.flow_mass_phase_comp[p, j] / sum(
                b.flow_mass_phase_comp[p, _j] for p, _j in phase_comp_list
            )

        self.eq_mass_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mass_frac_phase_comp
        )

    def _conc_mass_phase_comp(self):
        self.conc_mass_phase_comp = Var(
            self.phase_component_set,
            initialize=10,
            bounds=(0, 1e6),
            units=pyunits.kg * pyunits.m**-3,
            doc="Mass concentration",
        )

        def rule_conc_mass_phase_comp(b, p, j):
            return (
                b.conc_mass_phase_comp[p, j]
                == b.dens_mass_phase[p] * b.mass_frac_phase_comp[p, j]
            )

        self.eq_conc_mass_phase_comp = Constraint(
            self.phase_component_set, rule=rule_conc_mass_phase_comp
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

            return b.flow_mol_phase_comp[p, j] == b.mole_frac_phase_comp[p, j] * sum(
                b.flow_mol_phase_comp[p, _j] for _j in b.params.phase_comp_dict[p]
            )

        self.eq_mole_frac_phase_comp = Constraint(
            self.phase_component_set, rule=rule_mole_frac_phase_comp
        )

    def _conc_mol_phase_comp(self):

        self.conc_mol_phase_comp = Var(
            self.phase_component_set,
            initialize=500,
            bounds=(0, None),
            units=pyunits.mol * pyunits.m**-3,
            doc="Molar concentration",
        )

        def rule_conc_mol_phase_comp(b, p, j):
            return (
                b.conc_mol_phase_comp[p, j] * b.params.mw_comp[j]
                == b.conc_mass_phase_comp[p, j]
            )

        self.eq_conc_mol_phase_comp = Constraint(
            self.phase_component_set,
            rule=rule_conc_mol_phase_comp,
        )

    def _flow_vol_phase(self):
        self.flow_vol_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.s,
            doc="Total volumetric flow rate for phase",
        )

        def rule_flow_vol_phase(b, p):
            return (
                b.flow_vol_phase[p]
                == sum(
                    b.flow_mass_phase_comp[p, j]
                    for j in self.params.component_set
                    if (p, j) in self.phase_component_set
                )
                / b.dens_mass_phase[p]
            )

        self.eq_flow_vol_phase = Constraint(
            self.params.phase_list, rule=rule_flow_vol_phase
        )

    def _flow_mass_phase(self):
        self.flow_mass_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kg / pyunits.s,
            doc="Total mass flow rate for phase",
        )

        def rule_flow_mass_phase(b, p):
            return b.flow_mass_phase[p] == b.flow_vol_phase[p] * b.dens_mass_phase[p]

        self.eq_flow_mass_phase = Constraint(
            self.params.phase_list, rule=rule_flow_mass_phase
        )

    def _flow_vol(self):
        def rule_flow_vol(b):
            return sum(b.flow_vol_phase[p] for p in self.params.phase_list)

        self.flow_vol = Expression(rule=rule_flow_vol)

    def _diffus_phase_comp(self):

        self.diffus_phase_comp = Var(
            self.params.volatile_comp_set,
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
                    pressure_amb = 101325 * pyunits.Pa
                    sqrt_term = (1 / mw_j + 1 / mw_air) ** 0.5 * (
                        pyunits.g / pyunits.mol
                    ) ** 0.5  # units = dimensionless
                    numerator = (
                        (b.wilke_lee_param_A - b.wilke_lee_param_B * sqrt_term)
                        * (b.temperature[p] ** 1.5)
                        * sqrt_term
                    )  # units = cm2
                    denominator = (
                        pressure_amb
                        * (b.collision_molecular_separation[j]) ** 2
                        * b.collision_function_comp[j]
                    ) * wilke_lee_denom_units  # units = s
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
            self.params.volatile_comp_set,
            rule=rule_diffus_phase_comp,
        )

    def _energy_molecular_attraction_comp(self):

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

        self.energy_molecular_attraction_air = Expression(
            expr=pyunits.convert(
                boltzmann * self.energy_molecular_attraction_air_param,
                to_units=pyunits.erg,
            )
        )

        def expr_energy_molecular_attraction_phase_comp(b, j):
            return pyunits.convert(
                boltzmann
                * b.energy_molecular_attraction_comp_param
                * b.temperature_boiling_comp[j],
                to_units=pyunits.erg,
            )

        self.energy_molecular_attraction_comp = Expression(
            self.params.volatile_comps,
            expr=expr_energy_molecular_attraction_phase_comp,
        )

    def _energy_molecular_attraction(self):
        def expr_energy_molecular_attraction(b, j):
            return pyunits.convert(
                boltzmann
                * (
                    b.energy_molecular_attraction_comp[j]
                    / boltzmann
                    * b.energy_molecular_attraction_air
                    / boltzmann
                )
                ** 0.5,
                to_units=pyunits.erg,
            )

        self.energy_molecular_attraction = Expression(
            self.params.volatile_comps, expr=expr_energy_molecular_attraction
        )

    def _collision_molecular_separation_comp(self):

        self.collision_molecular_separation_param = Param(
            initialize=1.18,
            units=(pyunits.nm * pyunits.mol ** (1 / 3)) / (pyunits.liter ** (1 / 3)),
            doc="Molecular separation at collision equation parameter",
        )
        self.collision_molecular_separation_exp = Param(
            initialize=1 / 3,
            units=pyunits.dimensionless,
            doc="Molecular separation at collision equation exponent",
        )
        self.collision_molecular_separation_air = Param(
            initialize=0.3711,
            units=pyunits.nanometer,
            doc="Molecular separation at collision for air",
        )

        def expr_collision_molecular_separation_comp(b, j):
            molar_vol = pyunits.convert(
                b.molar_volume_comp[j],
                to_units=pyunits.liter / pyunits.mol,
            )
            return b.collision_molecular_separation_param * molar_vol ** (
                b.collision_molecular_separation_exp
            )

        self.collision_molecular_separation_comp = Expression(
            self.params.volatile_comps, expr=expr_collision_molecular_separation_comp
        )

    def _collision_molecular_separation(self):
        def expr_collision_molecular_separation(b, j):
            return 0.5 * (
                b.collision_molecular_separation_comp[j]
                + b.collision_molecular_separation_air
            )

        self.collision_molecular_separation = Expression(
            self.params.volatile_comps, expr=expr_collision_molecular_separation
        )

    def _collision_function_comp(self):
        def expr_collision_function_comp(b, j):
            return 10 ** (b.collision_function_zeta_comp[j])

        self.collision_function_comp = Expression(
            self.params.volatile_comps, expr=expr_collision_function_comp
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

        def expr_collision_function_zeta_comp(b, j):
            ee = b.collision_function_ee_comp[j]
            return (
                a + b_ * ee + c * ee**2 + d * ee**3 + e * ee**4 + f * ee**5 + g * ee**6
            )

        self.collision_function_zeta_comp = Expression(
            self.params.volatile_comps, expr=expr_collision_function_zeta_comp
        )

    def _collision_function_ee_comp(self):
        def expr_collision_function_ee_comp(b, j):
            t = b.temperature["Vap"]
            return log10((boltzmann * t) / b.energy_molecular_attraction[j])

        self.collision_function_ee_comp = Expression(
            self.params.volatile_comps, expr=expr_collision_function_ee_comp
        )

    def _molar_volume_comp(self):

        if (
            self.params.config.molar_volume_calculation
            == MolarVolumeCalculation.TynCalus
        ):
            self.molar_volume_comp = Var(
                self.params.volatile_comps,
                initialize=self.params.config.molar_volume_data,
                units=pyunits.m**3 / pyunits.mol,
                doc="Molar volume of solutes",
            )

            self.tyn_calus_param = Param(
                initialize=0.285,
                units=pyunits.mol**0.048 * pyunits.cm**-0.144,
                doc="Tyn-Calus equation parameter",
            )

            self.tyn_calus_exponent = Param(
                initialize=1.048,
                units=pyunits.dimensionless,
                doc="Tyn-Calus equation exponent",
            )

            def rule_molar_volume_comp(b, j):
                return b.molar_volume_comp[j] == pyunits.convert(
                    b.tyn_calus_param
                    * pyunits.convert(
                        b.molar_volume_comp_crit[j],
                        to_units=pyunits.cm**3 / pyunits.mol,
                    )
                    ** b.tyn_calus_exponent,
                    to_units=pyunits.m**3 / pyunits.mol,
                )

            self.eq_molar_volume_comp = Constraint(
                self.params.volatile_comps, rule=rule_molar_volume_comp
            )

        else:
            add_object_reference(
                self, "molar_volume_comp", self.params.molar_volume_comp
            )

    def _henry_comp(self):
        if self.params.config.temp_adjust_henry:

            add_object_reference(
                self, "henry_comp_ref", self.params.henry_comp
            )  # assume the data provided to config.henry_constant_data is @25C

            self.henry_comp = Var(
                self.params.volatile_comps,
                initialize=self.params.config.henry_constant_data,
                units=pyunits.dimensionless,
                doc="Temperature adjusted dimensionless Henry's constant",
            )

            def rule_henry_comp(b, j, doc="van't Hoff equation"):
                t0 = 298 * pyunits.degK
                exponential_term = pyunits.convert(
                    (-b.enth_change_dissolution_comp[j] / Constants.gas_constant)
                    * (1 / b.temperature["Liq"] - 1 / t0),
                    to_units=pyunits.dimensionless,
                )
                return b.henry_comp[j] == b.henry_comp_ref[j] * exp(exponential_term)

            self.eq_henry_comp = Constraint(
                self.params.volatile_comps, rule=rule_henry_comp
            )

        else:

            add_object_reference(self, "henry_comp", self.params.henry_comp)

    def _pressure_vap_sat(self):

        # Huang, J. (2018). A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice.
        # Journal of Applied Meteorology and Climatology, 57(6), 1265-1272. doi:10.1175/jamc-d-17-0334.1

        if (
            self.params.config.saturation_vapor_pressure_calculation
            != SaturationVaporPressureCalculation.none
        ):
            self.pressure_vap_sat = Var(
                ["H2O"],
                initialize=101325,
                units=pyunits.Pa,
                doc="Saturation vapor pressure of water",
            )

            if (
                self.params.config.saturation_vapor_pressure_calculation
                == SaturationVaporPressureCalculation.Huang
            ):
                self.sat_vap_press_A = a = Param(
                    initialize=34.494,
                    units=pyunits.dimensionless,
                    doc="Huang correlation: A parameter",
                )
                self.sat_vap_press_B = b_ = Param(
                    initialize=4924.99,
                    units=pyunits.dimensionless,
                    doc="Huang correlation: B parameter",
                )
                self.sat_vap_press_C = c = Param(
                    initialize=1.57,
                    units=pyunits.dimensionless,
                    doc="Huang correlation: C parameter",
                )
                self.sat_vap_press_D1 = d1 = Param(
                    initialize=237.1,
                    units=pyunits.dimensionless,
                    doc="Huang correlation: D1 parameter",
                )
                self.sat_vap_press_D2 = d2 = Param(
                    initialize=105,
                    units=pyunits.dimensionless,
                    doc="Huang correlation: D2 parameter",
                )

                def rule_pressure_vap_sat(b, h2o):
                    t = b.temperature["Vap"] - 273.15 * pyunits.degK
                    return (
                        b.pressure_vap_sat[h2o]
                        == exp(a - (b_ / (t + d1))) / (t + d2) ** c
                    )

            if (
                self.params.config.saturation_vapor_pressure_calculation
                == SaturationVaporPressureCalculation.ArdenBuck
            ):
                # https://en.wikipedia.org/wiki/Arden_Buck_equation
                self.arden_buck_coeff_a = a = Param(
                    initialize=6.1121,
                    units=pyunits.millibar,
                    doc="Arden Buck equation for saturation vapor pressure, a coefficient",
                )

                self.arden_buck_coeff_b = b_ = Param(
                    initialize=18.678,
                    units=pyunits.dimensionless,
                    doc="Arden Buck equation for saturation vapor pressure, b coefficient",
                )

                self.arden_buck_coeff_c = c = Param(
                    initialize=257.14,
                    units=pyunits.dimensionless,
                    doc="Arden Buck equation for saturation vapor pressure, c coefficient",
                )

                self.arden_buck_coeff_d = d = Param(
                    initialize=234.5,
                    units=pyunits.dimensionless,
                    doc="Arden Buck equation for saturation vapor pressure, d coefficient",
                )

                air_temp_degC_dimensionless = pyunits.convert(
                    (self.temperature["Vap"] - 273.15 * pyunits.degK)
                    * pyunits.degK**-1,
                    to_units=pyunits.dimensionless,
                )

                self.arden_buck_exponential_term = Expression(
                    expr=(b_ - air_temp_degC_dimensionless / d)
                    * (air_temp_degC_dimensionless / (c + air_temp_degC_dimensionless))
                )

                def rule_pressure_vap_sat(b, h2o):
                    return b.pressure_vap_sat[h2o] == pyunits.convert(
                        a * exp(self.arden_buck_exponential_term), to_units=pyunits.Pa
                    )

            self.eq_pressure_vap_sat = Constraint(["H2O"], rule=rule_pressure_vap_sat)

        else:

            add_object_reference(self, "pressure_vap_sat", self.params.pressure_vap_sat)

    def _pressure_vap(self):
        if (
            self.params.config.vapor_pressure_calculation
            != VaporPressureCalculation.none
        ):

            self.pressure_vap = Var(
                ["H2O"],
                initialize=1000,
                units=pyunits.Pa,
                doc="Vapor pressure of water",
            )

            if (
                self.params.config.vapor_pressure_calculation
                == VaporPressureCalculation.FromRelativeHumidity
            ):

                def rule_pressure_vap(b, h2o):
                    return (
                        b.pressure_vap[h2o]
                        == b.relative_humidity[h2o] * b.pressure_vap_sat[h2o]
                    )

            if (
                self.params.config.vapor_pressure_calculation
                == VaporPressureCalculation.Antoine
            ):
                # Antoine Eq
                # log10(P_vap) = A - B / (C + T)
                # coefficients valid for 1-100 degC

                self.antoine_A = a = Param(
                    initialize=8.07131,
                    units=pyunits.dimensionless,
                    doc="Antoine correlation: A parameter",
                )
                self.antoine_B = b_ = Param(
                    initialize=1730.63,
                    units=pyunits.dimensionless,
                    doc="Antoine correlation: B parameter",
                )
                self.antoine_C = c = Param(
                    initialize=233.426,
                    units=pyunits.dimensionless,
                    doc="Antoine correlation: C parameter",
                )

                def rule_pressure_vap(b, h2o):
                    t = pyunits.convert(
                        (b.temperature["Liq"] - 273.15 * pyunits.degK)
                        * pyunits.degK**-1,
                        to_units=pyunits.dimensionless,
                    )
                    antoine = a - (b_ / (c + t))
                    p_vap = 10 ** (antoine) * pyunits.mmHg
                    return b.pressure_vap[h2o] == pyunits.convert(
                        p_vap, to_units=pyunits.Pa
                    )

            self.eq_pressure_vap = Constraint(["H2O"], rule=rule_pressure_vap)

        else:
            add_object_reference(self, "pressure_vap", self.params.pressure_vap)

    def _relative_humidity(self):

        if (
            self.params.config.relative_humidity_calculation
            == RelativeHumidityCalculation.VaporPressureRatio
        ):

            self.relative_humidity = Var(
                ["H2O"],
                initialize=0.5,
                bounds=(0, 1),
                units=pyunits.dimensionless,
                doc="Relative humidity of air-water system",
            )

            def rule_relative_humidity(b, h2o):
                return (
                    b.relative_humidity[h2o]
                    == b.pressure_vap[h2o] / b.pressure_vap_sat[h2o]
                )

            self.eq_relative_humidity = Constraint(["H2O"], rule=rule_relative_humidity)

        else:
            add_object_reference(
                self, "relative_humidity", self.params.relative_humidity
            )

    def _dh_vap_mass_solvent(self):

        if (
            self.params.config.latent_heat_of_vaporization_calculation
            == LatentHeatVaporizationCalculation.Sharqawy
        ):

            # Latent heat of evaporation of pure water: Parameters from Sharqawy et al. (2010), eq. 54
            enth_mass_units = pyunits.J / pyunits.kg
            t_inv_units = pyunits.K**-1

            self.dh_vap_mass_solvent = Var(
                initialize=2.4e3,
                bounds=(1, 1e5),
                units=pyunits.kJ / pyunits.kg,
                doc="Latent heat of vaporization of pure water",
            )

            self.dh_vap_w_param_0 = a = Param(
                initialize=2.501e6,
                units=enth_mass_units,
                doc="Latent heat of pure water parameter 0",
            )
            self.dh_vap_w_param_1 = b_ = Param(
                initialize=-2.369e3,
                units=enth_mass_units * t_inv_units**1,
                doc="Latent heat of pure water parameter 1",
            )
            self.dh_vap_w_param_2 = c = Param(
                initialize=2.678e-1,
                units=enth_mass_units * t_inv_units**2,
                doc="Latent heat of pure water parameter 2",
            )
            self.dh_vap_w_param_3 = d = Param(
                initialize=-8.103e-3,
                units=enth_mass_units * t_inv_units**3,
                doc="Latent heat of pure water parameter 3",
            )
            self.dh_vap_w_param_4 = e = Param(
                initialize=-2.079e-5,
                units=enth_mass_units * t_inv_units**4,
                doc="Latent heat of pure water parameter 4",
            )

            def rule_dh_vap_mass_solvent(b):
                t = b.temperature["Liq"] - 273.15 * pyunits.K
                dh_vap_liq = a + b_ * t + c * t**2 + d * t**3 + e * t**4
                return b.dh_vap_mass_solvent == pyunits.convert(
                    dh_vap_liq, to_units=pyunits.kJ / pyunits.kg
                )

            self.eq_dh_vap_mass_solvent = Constraint(rule=rule_dh_vap_mass_solvent)

        else:
            add_object_reference(
                self, "dh_vap_mass_solvent", self.params.dh_vap_mass_solvent
            )

    def _cp_mass_solvent(self):
        if (
            self.params.config.specific_heat_of_water_calculation
            == SpecificHeatWaterCalculation.Sharqawy
        ):

            cp_units = pyunits.J / (pyunits.kg * pyunits.K)
            t_inv_units = pyunits.K**-1

            # Specific heat parameters for pure water from eq (9) in Sharqawy et al. (2010)
            self.cp_mass_solvent = Var(
                self.params.phase_list,
                initialize=4e3,
                bounds=(1e-5, 1e5),
                units=pyunits.J / pyunits.kg / pyunits.K,
                doc="Specific heat capacity of pure solvent",
            )

            self.cp_phase_param_A1 = Param(
                initialize=5.328,
                units=cp_units,
                doc="Specific heat of seawater parameter A1",
            )
            self.cp_phase_param_B1 = Param(
                initialize=-6.913e-3,
                units=cp_units * t_inv_units,
                doc="Specific heat of seawater parameter B1",
            )
            self.cp_phase_param_C1 = Param(
                initialize=9.6e-6,
                units=cp_units * t_inv_units**2,
                doc="Specific heat of seawater parameter C1",
            )
            self.cp_phase_param_D1 = Param(
                initialize=2.5e-9,
                units=cp_units * t_inv_units**3,
                doc="Specific heat of seawater parameter D1",
            )

            # Specific heat parameters for Cp vapor from NIST Webbook - Chase, M.W., Jr., NIST-JANAF Themochemical Tables
            self.cp_vap_param_A = Param(
                initialize=30.09200 / 18.01528e-3,
                units=cp_units,
                doc="Specific heat of water vapor parameter A",
            )
            self.cp_vap_param_B = Param(
                initialize=6.832514 / 18.01528e-3,
                units=cp_units * t_inv_units,
                doc="Specific heat of water vapor parameter B",
            )
            self.cp_vap_param_C = Param(
                initialize=6.793435 / 18.01528e-3,
                units=cp_units * t_inv_units**2,
                doc="Specific heat of water vapor parameter C",
            )
            self.cp_vap_param_D = Param(
                initialize=-2.534480 / 18.01528e-3,
                units=cp_units * t_inv_units**3,
                doc="Specific heat of water vapor parameter D",
            )
            self.cp_vap_param_E = Param(
                initialize=0.082139 / 18.01528e-3,
                units=cp_units * t_inv_units**-2,
                doc="Specific heat of water vapor parameter E",
            )

            def rule_cp_mass_solvent(b, p):
                if p == "Liq":
                    # specific heat, eq. 9 in Sharqawy et al. (2010)
                    # Convert T90 to T68, eq. 4 in Sharqawy et al. (2010); primary reference from Rusby (1991)
                    t = (b.temperature["Liq"] - 0.00025 * 273.15 * pyunits.K) / (
                        1 - 0.00025
                    )
                    A = b.cp_phase_param_A1
                    B = b.cp_phase_param_B1
                    C = b.cp_phase_param_C1
                    D = b.cp_phase_param_D1
                    return (
                        b.cp_mass_solvent[p] == (A + B * t + C * t**2 + D * t**3) * 1000
                    )
                if p == "Vap":
                    t = b.temperature["Liq"] / 1000
                    return (
                        b.cp_mass_solvent[p]
                        == b.cp_vap_param_A
                        + b.cp_vap_param_B * t
                        + b.cp_vap_param_C * t**2
                        + b.cp_vap_param_D * t**3
                        + b.cp_vap_param_E / t**2
                    )

            self.eq_cp_mass_solvent = Constraint(
                self.params.phase_list, rule=rule_cp_mass_solvent
            )
        else:
            add_object_reference(self, "cp_mass_solvent", self.params.cp_mass_solvent)

    def _cp_air(self):
        add_object_reference(self, "cp_air", self.params.cp_air)

    def _molar_volume_comp_crit(self):
        add_object_reference(
            self,
            "molar_volume_comp_crit",
            self.params.molar_volume_comp_crit,
        )

    def _temperature_boiling_comp(self):
        add_object_reference(
            self, "temperature_boiling_comp", self.params.temperature_boiling_comp
        )

    def _mw_comp(self):
        add_object_reference(self, "mw_comp", self.params.mw_comp)

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

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def get_material_density_terms(self):
        return self.dens_mass_phase

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
        for j, v in self.mw_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(self.mw_comp[j], value(v) ** -1)

        for p, j in self.phase_component_set:
            if iscale.get_scaling_factor(self.flow_mass_phase_comp[p, j]) is None:
                try:
                    sf = value(self.flow_mass_phase_comp[p, j]) ** -1
                except:
                    sf = 1
                iscale.set_scaling_factor(self.flow_mass_phase_comp[p, j], sf)

        if self.is_property_constructed("flow_mol_phase_comp"):
            for p in self.phase_list:
                for j in self.params.phase_comp_dict[p]:
                    if (
                        iscale.get_scaling_factor(self.flow_mol_phase_comp[p, j])
                        is None
                    ):
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, j]
                        ) / iscale.get_scaling_factor(self.params.mw_comp[j])
                        iscale.set_scaling_factor(self.flow_mol_phase_comp[p, j], sf)

        if self.is_property_constructed("diffus_phase_comp"):
            for p, j in self.params.volatile_comp_set:
                if iscale.get_scaling_factor(self.diffus_phase_comp[p, j]) is None:
                    if (p, j) in self.params.config.diffusivity_data.keys():
                        sf = value(self.params.config.diffusivity_data[p, j]) ** -1
                    else:
                        if p == "Liq":
                            sf = 1e10
                        else:
                            sf = 1e6
                    iscale.set_scaling_factor(self.diffus_phase_comp[p, j], sf)

        if self.is_property_constructed("dens_mass_phase"):
            for p in self.params.phase_list:
                if iscale.get_scaling_factor(self.dens_mass_phase[p]) is None:
                    sf = value(self.params.config.density_data[p]) ** -1
                    iscale.set_scaling_factor(self.dens_mass_phase[p], sf)

        if self.is_property_constructed("visc_d_phase"):
            for p in self.params.phase_list:
                if iscale.get_scaling_factor(self.visc_d_phase[p]) is None:
                    sf = value(self.params.config.dynamic_viscosity_data[p]) ** -1
                    iscale.set_scaling_factor(self.visc_d_phase[p], sf)

        if self.is_property_constructed("flow_vol_phase"):
            for p in self.params.phase_list:
                if iscale.get_scaling_factor(self.flow_vol_phase[p]) is None:
                    flow_vol_phase = sum(
                        value(self.flow_mass_phase_comp[p, j])
                        for j in self.params.component_set
                        if j in self.params.phase_comp_dict[p]
                    ) / value(self.dens_mass_phase[p])
                    sf = flow_vol_phase**-1
                    iscale.set_scaling_factor(self.flow_vol_phase[p], sf)

        if self.is_property_constructed("flow_mass_phase"):
            for p in self.params.phase_list:
                if iscale.get_scaling_factor(self.flow_mass_phase[p]) is None:
                    sf = iscale.get_scaling_factor(
                        self.flow_vol_phase[p]
                    ) * iscale.get_scaling_factor(self.dens_mass_phase[p])
                    iscale.set_scaling_factor(self.flow_mass_phase[p], sf)

        if self.is_property_constructed("cp_mass_solvent"):
            for p in self.params.phase_list:
                if iscale.get_scaling_factor(self.cp_mass_solvent[p]) is None:
                    iscale.set_scaling_factor(self.cp_mass_solvent[p], 1e-3)

        if self.is_property_constructed("mass_frac_phase_comp"):
            for p, j in self.phase_component_set:
                if iscale.get_scaling_factor(self.mass_frac_phase_comp[p, j]) is None:
                    if j in ["Air", "H2O"]:
                        iscale.set_scaling_factor(self.mass_frac_phase_comp[p, j], 1)
                    else:
                        if p == "Vap":
                            j_solv = "Air"
                        if p == "Liq":
                            j_solv = "H2O"
                        sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, j]
                        ) / iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, j_solv]
                        )
                        iscale.set_scaling_factor(self.mass_frac_phase_comp[p, j], sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            for p, j in self.phase_component_set:
                sf_dens = iscale.get_scaling_factor(self.dens_mass_phase[p])
                if iscale.get_scaling_factor(self.conc_mass_phase_comp[p, j]) is None:
                    sf = sf_dens * iscale.get_scaling_factor(
                        self.mass_frac_phase_comp[p, j]
                    )
                    iscale.set_scaling_factor(self.conc_mass_phase_comp[p, j], sf)

        if self.is_property_constructed("conc_mol_phase_comp"):
            for p, j in self.phase_component_set:
                mw = self.params.mw_comp[j]
                if iscale.get_scaling_factor(self.conc_mol_phase_comp[p, j]) is None:
                    sf = iscale.get_scaling_factor(
                        self.conc_mass_phase_comp[p, j]
                    ) / iscale.get_scaling_factor(mw)
                    iscale.set_scaling_factor(self.conc_mol_phase_comp[p, j], sf)

        if self.is_property_constructed("mol_frac_phase_comp"):
            for p, j in self.phase_component_set:
                mw = self.params.mw_comp[j]
                if iscale.get_scaling_factor(self.mol_frac_phase_comp[p, j]) is None:
                    if j in ["Air", "H2O"]:
                        iscale.set_scaling_factor(self.mol_frac_phase_comp[p, j], 1)
                    else:
                        if p == "Vap":
                            j_solv = "Air"
                        if p == "Liq":
                            j_solv = "H2O"
                        flow_mol_j_sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, j]
                        ) / iscale.get_scaling_factor(mw)
                        flow_mol_solv_sf = iscale.get_scaling_factor(
                            self.flow_mass_phase_comp[p, j_solv]
                        ) / iscale.get_scaling_factor(self.params.mw_comp[j_solv])
                        sf = flow_mol_j_sf / flow_mol_solv_sf
                        iscale.set_scaling_factor(self.mol_frac_phase_comp[p, j], sf)

        for j in self.params.volatile_comps:
            if self.is_property_constructed("collision_function_zeta_comp"):
                if (
                    iscale.get_scaling_factor(self.collision_function_zeta_comp[j])
                    is None
                ):
                    iscale.set_scaling_factor(self.collision_function_zeta_comp[j], 10)

            if self.is_property_constructed("collision_function_ee_comp"):
                if (
                    iscale.get_scaling_factor(self.collision_function_ee_comp[j])
                    is None
                ):
                    iscale.set_scaling_factor(self.collision_function_ee_comp[j], 10)

            if self.is_property_constructed("collision_function_comp"):
                if iscale.get_scaling_factor(self.collision_function_comp[j]) is None:
                    iscale.set_scaling_factor(self.collision_function_comp[j], 10)

            if self.is_property_constructed("molar_volume_comp_crit"):
                if iscale.get_scaling_factor(self.molar_volume_comp_crit[j]) is None:
                    sf = value(self.params.config.critical_molar_volume_data[j]) ** -1
                    iscale.set_scaling_factor(self.molar_volume_comp_crit[j], sf)

            if self.is_property_constructed("molar_volume_comp"):
                if iscale.get_scaling_factor(self.molar_volume_comp[j]) is None:
                    if (
                        self.params.config.molar_volume_calculation
                        is MolarVolumeCalculation.TynCalus
                    ):
                        sf = iscale.get_scaling_factor(self.molar_volume_comp_crit[j])
                    else:
                        sf = value(self.params.config.molar_volume_data[j]) ** -1
                    iscale.set_scaling_factor(self.molar_volume_comp[j], sf)

            if self.is_property_constructed("henry_comp"):
                if iscale.get_scaling_factor(self.henry_comp[j]) is None:
                    iscale.set_scaling_factor(self.henry_comp[j], 1)

            if self.is_property_constructed("temperature_boiling_comp"):
                if iscale.get_scaling_factor(self.temperature_boiling_comp[j]) is None:
                    iscale.set_scaling_factor(self.temperature_boiling_comp[j], 0.1)

            if self.is_property_constructed("collision_molecular_separation_comp"):
                if (
                    iscale.get_scaling_factor(
                        self.collision_molecular_separation_comp[j]
                    )
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.collision_molecular_separation_comp[j], 10
                    )
            if self.is_property_constructed("energy_molecular_attraction_comp"):
                if (
                    iscale.get_scaling_factor(self.energy_molecular_attraction_comp[j])
                    is None
                ):
                    iscale.set_scaling_factor(
                        self.energy_molecular_attraction_comp[j], 1e14
                    )

        if self.is_property_constructed("dh_vap_mass_solvent"):
            if iscale.get_scaling_factor(self.dh_vap_mass_solvent) is None:
                iscale.set_scaling_factor(self.dh_vap_mass_solvent, 1e-3)

        if self.is_property_constructed("pressure_vap_sat"):
            if iscale.get_scaling_factor(self.pressure_vap_sat["H2O"]) is None:
                iscale.set_scaling_factor(self.pressure_vap_sat["H2O"], 1e-3)

        if self.is_property_constructed("pressure_vap"):
            if iscale.get_scaling_factor(self.pressure_vap["H2O"]) is None:
                iscale.set_scaling_factor(self.pressure_vap["H2O"], 1e-3)

        if self.is_property_constructed("relative_humidity"):
            if iscale.get_scaling_factor(self.relative_humidity["H2O"]) is None:
                iscale.set_scaling_factor(self.relative_humidity["H2O"], 10)

        transform_property_constraints(self)
