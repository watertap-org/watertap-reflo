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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Suffix,
    value,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import WaterTAP cores
from watertap.core import InitializationMixin

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"

@declare_process_block_class("MgCrystallizerZO")
class MgCrystallizerZOData(InitializationMixin, UnitModelBlockData):
    """
    Zero order Mg Crystallizer model
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

    def build(self):
        super().build()

        #TODO: add references

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add state blocks for inlet, outlet, and waste
        # These include the state variables and any other properties on demand
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add outlet block
        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of liquid outlet",
            **tmp_dict,
        )
        self.properties_out[0].pressure.fix(101325)
        self.properties_out[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        self.properties_out[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        # Initial guess 
        self.properties_out[0].flow_vol_phase["Liq"].value = 1e-5
        self.properties_out[0].flow_vol_phase["Sol"].value = 1e-5
        @self.Constraint(doc="Outlet temperature remains the same")
        def eq_outlet_temperature(b):
            return (b.properties_out[0].temperature == b.properties_in[0].temperature)

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)

        # Model parameters
        self.hydromagnesite_solubility = Param(
            initialize = 0.76071,
            units = pyunits.g / pyunits.L,
            doc = "Solubility of hydromagnesite [Mg5(CO3)4(OH)2 4H2O] at 25 deg C"
        )

        self.calcite_solubility = Param(
            initialize = 0.17631,
            units = pyunits.g / pyunits.L,
            doc = "Solubility of calcite [CaCO3] at 25 deg C"
        )

        self.Mg_molar_weight = Param(
            initialize = 24.305,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of Mg'
        )

        self.Ca_molar_weight = Param(
            initialize = 40.087,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of Mg'
        )

        self.hydromagnesite_molar_weight = Param(
            initialize = 467.64,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of hydromagnesite [Mg5(CO3)4(OH)2 4H2O]'
        )

        self.calcite_molar_weight = Param(
            initialize = 100.0869,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of calcite [CaCO3]'
        )

        self.calcium_hydroxide_molar_weight = Param(
            initialize = 74.095,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of calcium hydroxide [Ca(OH)2]'
        )

        self.co2_molar_weight = Param(
            initialize = 44.01,
            units = pyunits.g / pyunits.mol,
            doc = 'Molar weight of carbon dioxide [CO2]'
        )

        # Model input variables
        self.conversion_factor_mg = Var(
            initialize = 0.9,
            units = pyunits.dimensionless,
            doc = 'Conversion factor of Mg2+ precipitation reaction'
        )

        self.conversion_factor_ca = Var(
            initialize = 1,
            units = pyunits.dimensionless,
            doc = 'Conversion factor of Ca2+ precipitation reaction'
        )

        self.conc_chemical_dose = Var(
            initialize = 0.5,
            units = pyunits.mol / pyunits.L,
            doc = 'Concentration of calcium hydroxide addition'
        )        
        

        # Add Mg and Ca in addition to the WaterTAP crystallizer property package (cryst_prop_pack)
        self.Mg_conc_inlet = Var(
            initialize = 35.49,
            units = pyunits.g / pyunits.L,
            doc = 'Mg concentration in the feed flow'
        )

        self.Ca_conc_inlet = Var(
            initialize = 10.60,
            units = pyunits.g / pyunits.L,
            doc = 'Ca concentration in the feed flow'
        )

        # Add output variables
        self.calcium_hydroxide_vol = Var(
            initialize = 15.53,
            units = pyunits.m**3 / pyunits.day,
            doc = "Volumetric chemical addition of calcium hydroxide"
        )

        self.calcium_hydroxide_mass = Var(
            initialize = 575.26,
            units = pyunits.kg / pyunits.day,
            doc = "Chemical addition mass of calcium hydroxide"
        )

        # Add intermediate variables / paramters
        @self.Expression(doc="Mg concentration in the feed flow after diluted by chemical additions (g/L)")
        def Mg_conc_diluted(b):
            return (b.Mg_conc_inlet * b.properties_in[0].flow_vol_phase["Liq"]
                    / (b.properties_in[0].flow_vol_phase["Liq"] +
                       pyunits.convert(b.calcium_hydroxide_vol,
                                       to_units = pyunits.m**3 / pyunits.s)))

        @self.Expression(doc="Ca concentration in the feed flow after diluted by chemical additions (g/L)")
        def Ca_conc_diluted(b):
            return (b.Ca_conc_inlet * b.properties_in[0].flow_vol_phase["Liq"]
                    / (b.properties_in[0].flow_vol_phase["Liq"] +
                       pyunits.convert(b.calcium_hydroxide_vol,
                                       to_units = pyunits.m**3 / pyunits.s)))

        # self.Mg_solubility = Param(initialize = 0.198, units = pyunits.g/pyunits.L)
        @self.Expression(doc="Mg solubililty at 25 deg C (g/L)")
        def Mg_solubility(b):
            return (b.hydromagnesite_solubility * 5 * b.Mg_molar_weight / b.hydromagnesite_molar_weight)

        @self.Expression(doc="Mg solubililty at 25 deg C (g/L)")
        def Ca_solubility(b):
            return (b.calcite_solubility * b.Ca_molar_weight / b.calcite_molar_weight)

        # Add equations as constraints / expressions
        @self.Expression(doc="Hydromagnesite precipitation rate (kg/day)")
        def hydromagnesite_precipitated(b):
            return ((b.Mg_conc_diluted - b.Mg_solubility)
                    * (pyunits.convert(b.properties_in[0].flow_vol_phase["Liq"],
                                       to_units = pyunits.m**3 / pyunits.day)
                       + b.calcium_hydroxide_vol)
                    * b.conversion_factor_mg
                    * b.hydromagnesite_molar_weight / (5 * b.Mg_molar_weight))

        @self.Expression(doc="Calcite precipitation rate from feed flow (kg/day)")
        def calcite_precipitated_feed(b):
            return ((b.Ca_conc_diluted - b.Ca_solubility)
                    * (pyunits.convert(b.properties_in[0].flow_vol_phase["Liq"],
                                       to_units = pyunits.m**3 / pyunits.day)
                       + b.calcium_hydroxide_vol)
                    * b.conversion_factor_ca
                    * b.calcite_molar_weight / b.Ca_molar_weight)

        @self.Expression(doc="Calcite precipitation rate from chemical dose (kg/day)")
        def calcite_precipitated_dose(b):
            return (b.calcium_hydroxide_mass
                    * b.calcite_molar_weight
                    / b.calcium_hydroxide_molar_weight)

        @self.Expression(doc="Recovery rate of Mg")
        def Mg_recovery_rate(b):
            return (b.hydromagnesite_precipitated
                    * 5 * b.Mg_molar_weight / b.hydromagnesite_molar_weight
                    / (b.Mg_conc_inlet
                       * pyunits.convert(b.properties_in[0].flow_vol_phase["Liq"],
                                         to_units = pyunits.m**3 / pyunits.day)))

        @self.Expression(doc="Recovery rate of Ca")
        def Ca_recovery_rate(b):
            return (b.calcite_precipitated_feed
                    * b.Ca_molar_weight / b.calcite_molar_weight
                    / (b.Ca_conc_inlet
                       * pyunits.convert(b.properties_in[0].flow_vol_phase["Liq"],
                                         to_units = pyunits.m**3 / pyunits.day)))

        @self.Constraint(doc='Calculate NaCl concentration in the waste flow')
        def eqn_conc_NaCl_out(b):
            return (b.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"] ==
                    b.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
                    * b.properties_in[0].flow_vol_phase["Liq"]
                    / (b.properties_in[0].flow_vol_phase["Liq"] +
                       pyunits.convert(b.calcium_hydroxide_vol,
                                       to_units = pyunits.m**3 / pyunits.s)))

        @self.Constraint(doc='Calculate mass flow rate of calcium hydroxide addition')
        def eqn_CaOH2_mass(b):
            return (b.calcium_hydroxide_mass ==
                    b.hydromagnesite_precipitated
                    * 5 * b.calcium_hydroxide_molar_weight
                    / b.hydromagnesite_molar_weight
                    + b.calcite_precipitated_feed
                    * b.calcium_hydroxide_molar_weight
                    / b.calcite_molar_weight
                    )

        @self.Constraint(doc='Calculate volumetric flow rate of calcium hydroxide addition')
        def eqn_CaOH2_vol(b):
            return (b.calcium_hydroxide_vol ==
                    b.calcium_hydroxide_mass
                    / b.calcium_hydroxide_molar_weight
                    / b.conc_chemical_dose)

        @self.Constraint(doc='Calculate the volumetric flow rate of waste flow')
        def eqn_outflow_vol(b):
            return (b.properties_out[0].flow_vol_phase["Liq"]
                    == b.properties_in[0].flow_vol_phase["Liq"]
                    + pyunits.convert(b.calcium_hydroxide_vol,
                                      to_units = pyunits.m**3 / pyunits.s))

        @self.Expression(doc='CO2 mass flow rate required')
        def co2_mass_rate(b):
            return (b.hydromagnesite_precipitated
                    * 4 * b.co2_molar_weight
                    / b.hydromagnesite_molar_weight
                    + b.calcite_precipitated_feed
                    * b.co2_molar_weight
                    / b.calcite_molar_weight)


    def initialize_build(
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
        # ---------------------------------------------------------------------
        flags = blk.properties_in.initialize(
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
            blk.state_args = state_args = {}
            state_dict = blk.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        # Check degree of freedom
        if degrees_of_freedom(blk) != 0:
            raise InitializationError(
                f"{blk.name} degrees of freedom were not 0 at the beginning of initialization. DoF = {degrees_of_freedom(blk)}."
            )

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=False)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.Mg_conc_inlet) is None:
            iscale.set_scaling_factor(self.Mg_conc_inlet, 1e-1)

        if iscale.get_scaling_factor(self.Ca_conc_inlet) is None:
            iscale.set_scaling_factor(self.Ca_conc_inlet, 1e-1)

        if iscale.get_scaling_factor(self.conversion_factor_mg) is None:
            iscale.set_scaling_factor(self.conversion_factor_mg, 1e0)

        if iscale.get_scaling_factor(self.conversion_factor_ca) is None:
            iscale.set_scaling_factor(self.conversion_factor_ca, 1e0)

        if iscale.get_scaling_factor(self.conc_chemical_dose) is None:
            iscale.set_scaling_factor(self.conc_chemical_dose, 1e0)

        if iscale.get_scaling_factor(self.calcium_hydroxide_vol) is None:
            iscale.set_scaling_factor(self.calcium_hydroxide_vol, 1e-2)

        if iscale.get_scaling_factor(self.calcium_hydroxide_mass) is None:
            iscale.set_scaling_factor(self.calcium_hydroxide_mass, 1e-3)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Water Inlet": self.properties_in,
                "Waste Outlet": self.properties_out
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["NaCl concentration in the waste flow (g/L)"] = value(self.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"])
        var_dict["Hydromagnesite precipitated (kg/day)"] = value(self.hydromagnesite_precipitated)
        var_dict["Calcite precipitated (kg/day)"] = value(self.calcite_precipitated_feed + self.calcite_precipitated_dose)
        var_dict["Mg recovery rate"] = value(self.Mg_recovery_rate)
        var_dict["Ca recovery rate"] = value(self.Ca_recovery_rate)
        var_dict["Ca(OH)2 dose (kg/day)"] = value(self.calcium_hydroxide_mass)
        var_dict["Ca(OH)2 dose (m3/day)"] = value(self.calcium_hydroxide_vol)
        var_dict["CO2 dose (kg/day)"] = value(self.co2_mass_rate)

        return {"vars": var_dict}




