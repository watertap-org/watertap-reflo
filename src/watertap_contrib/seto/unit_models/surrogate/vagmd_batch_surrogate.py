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

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    value,
    Suffix,
    ConcreteModel,
    NonNegativeReals,
    Reference,
    Constraint,
    units as pyunits,
    exp,
    log
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

from pyomo.dae import ContinuousSet

# Import IDAES cores
from idaes.core import (
    FlowsheetBlock,
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

# Import Watertap packages
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
import numpy as np


@declare_process_block_class("VAGMDbatch_surrogate")
class VAGMDbatchData(UnitModelBlockData):
    """
    Membrane distillation (air-gapped) - batch operation model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False, True]),
            default=True,
            description="Dynamic model flag - must be True",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "time",
        ConfigValue(
            description="",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "time_units",
        ConfigValue(
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "time_set",
        ConfigValue(
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=True,
            domain=In([False, True]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
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
        "high_brine_salinity",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="""A ConfigBlock with arguments to be passed to indicate if the brine
            has a salinity higher than 175.3 g/L""",
        ),
    )
    CONFIG.declare(
        "module_type",
        ConfigValue(
            default=7,
            doc="""'Selection of module type (7 for AS7C1.5L and 26 for AS26C7.2L)'""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units


        '''
        MD Module type and corresponding area
        '''
        if self.config.module_type == 7:
            module_area = 7.2
        else: # module_type = 26
            module_area = 25.92


        '''
        Model parameters
        '''
        self.module_type = Set(initialize = (7, 26))

        self.Vdisch = Param(initialize = 3.2175,
                        units = pyunits.L,
                        doc = 'Discharge volume (L)')

        self.Vmin = Param(initialize = 50,
                        units = pyunits.L,
                        doc = 'Void volume of an Ags module')


        self.minS = Param(initialize = 35,
                        units = pyunits.g / pyunits.L,
                        doc = 'Minimum feed salinity (g/L)')

        self.maxS = Param(self.module_type,
                    initialize = {7: 292.2, 26: 245.5},
                    units = pyunits.g / pyunits.L,
                    doc = 'Maximum feed salinity for the selected module (g/L)')

        self.U    = Param(initialize = 3168,
                    units = pyunits.W / pyunits.m**2 /pyunits.C,
                    doc = 'Overall heat transfer coefficient (W/m2/C)')

        self.AHX  = Param(initialize = 1.34,
                    units = pyunits.m**2,
                    doc = 'Effective heat transfer surface area (m2)')

        self.CFR  = Param(initialize = 1265,
                    units = pyunits.L / pyunits.h,
                    doc = 'Cooling flow rate (L/h)')
        # Cooling flow rate is fixed to 1265 to maintain vacuum pressure inside the MD module around 200 mbar.
    
        self.APdropCFR = Param(initialize = 170,
                               units = pyunits.mbar,
                               doc = 'Pressure drop of cooling water')

        self.EffPumpCFR  = Param(initialize = 0.6,
                            doc = 'Cooling water pump efficiency')

        '''
        Input variables
        '''
        self.RR = Var(initialize = 0.5,
                bounds = (0, 1), # the maximum value is dynamic to the feed salinity and is constrained later
                doc = 'Water recovery rate')

        self.TEI = Var(
                    initialize = 80,
                    bounds = (60, 80),
                    units = pyunits.K,
                    doc = 'Evaporator inlet temperature (C)')

        self.TCI = Var(
                    initialize = 20,
                    bounds = (20, 30),
                    units = pyunits.K,
                    doc = 'Condenser inlet temperature (C)')

        # self.Ttank = Var(
        #             initialize = 25,
        #             bounds = (self.TCI(), 30),
        #             units = pyunits.K,
        #             doc = 'Iniital temperature of the saline feed (C)')            

        self.j = Var(initialize = 0,
                    doc = 'Type of cooling system (0 for closed system and 1 for open system)')

        self.TCoolIn = Var(
                    initialize = 25,
                    bounds = (20, 30),
                    units = pyunits.K,
                    doc = 'Initial temperature of the cooling water (C)')     
        # TCoolIn is only applicable for open cooling system (j = 1)

        self.V0 =  Var(
                    initialize = 50,
                    bounds = (self.Vmin, None),
                    units = pyunits.L,
                    doc = 'Initial batch volume (L))')

        self.tR = Var(initialize = 30,
                units = pyunits.s,
                doc = 'Time step of the simulation (s)')                     

        '''
        Intermediate variables
        '''

        self.maxR = Var(initialize = 1,
                    bounds = (0,1),
                    doc = 'Maximum recovery rate (%)')

        self.PFlux = Var(self.flowsheet().config.time, 
                         initialize = 10,
                         units = pyunits.L / pyunits.h / pyunits.m**2,
                         doc = 'Permeate flux (L/h/m2)')

        self.permeate_flow_rate = Var(self.flowsheet().config.time, 
                         initialize = 10,
                         units = pyunits.L / pyunits.h ,
                         doc = 'Permeate flow rate (L/h)')

        self.AccVd = Var(self.flowsheet().config.time, 
                         initialize = 0,
                         units = pyunits.L,
                         doc = 'Accumulated volume of distillate(L)')

        self.accumulated_recovery_rate = Var(self.flowsheet().config.time, 
                         initialize = 0,
                         units = pyunits.dimensionless,
                         doc = 'Accumulated recovery rate after each recirculation')

        self.log_mean_temp_dif = Var(self.flowsheet().config.time, 
                         initialize = 0,
                         units = pyunits.C,
                         doc = 'Log mean temperature difference in the MD module')

        '''
        Add variables for output
        '''
        self.thermal_power = Var(self.flowsheet().config.time, 
                         initialize = 0,
                         units = pyunits.kW,
                         doc = 'Thermal power during each timestep (kW)')

        self.thermal_energy = Var(self.flowsheet().config.time, 
                         initialize = 0,
                         units = pyunits.kW,
                         doc = 'Thermal energy consumption during each timestep (kWh)')

        self.specific_thermal_energy_consumption = Var(
                         self.flowsheet().config.time, 
                         initialize = 100,
                         units = pyunits.kWh / pyunits.m**3,
                         doc = 'Specific thermal power consumption (kWh/m3)')


        self.gain_output_ratio = Var(
                         self.flowsheet().config.time, 
                         initialize = 5,
                         units = pyunits.kW / pyunits.kW,
                         doc = 'Gained output ratio')

        """
        Add block for the batched feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict
        )


        """
        Add block for water at the condenser outlet
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.condenser_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of condenser outlet water",
            **tmp_dict
        )

        """
        Add block for water at the evaporator outlet
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.evaporator_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        """
        Add block for inlet cooling water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.cooling_in_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of inlet cooling water",
            **tmp_dict
        )

        """
        Add block for outlet cooling water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.cooling_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet cooling water",
            **tmp_dict
        )

        """
        Add block for the average status of the feed flow
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.avg_feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        """
        Add block for the average status of the cooling flow
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.avg_cool_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        '''
        System constraint
        '''
        # Maximum recovery rate (calculated based on the max allowed salinity)
        @self.Constraint(doc="Maximum recovery rate")
        def maxR_cal(b):
            return b.maxR == (1 - (b.feed_props[0].conc_mass_phase_comp["Liq","TDS"] / b.maxS[b.config.module_type]) ) 


        '''
        Set up for the beginning of the first circulation
        '''
        S0 = self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        FFR0 = pyunits.convert(self.feed_props[0].flow_vol_phase["Liq"], to_units = pyunits.L / pyunits.h)
        Ttank0 = self.feed_props[0].temperature - 273.15 * pyunits.K
        @self.Constraint(doc='initial permeate flux')
        def PFlux_1(b):
            return b.PFlux[0] == b.VAGMD_Models_NaCl(b.TEI, FFR0, b.TCI, S0, Ttank0)[0]

        @self.Constraint(doc='initial permeate flow rate')
        def PFR_1(b):
            return b.permeate_flow_rate[0] == b.PFlux[0] * module_area

        @self.Constraint(doc='initial condenser outlet temperature')
        def TCO_1(b):
            return b.condenser_out_props[0].temperature == b.VAGMD_Models_NaCl(b.TEI, FFR0, b.TCI, S0, Ttank0)[1] + 273.15

        @self.Constraint(doc='initial permeate flux')
        def TEO_1(b):
            return b.evaporator_out_props[0].temperature == b.VAGMD_Models_NaCl(b.TEI, FFR0, b.TCI, S0, Ttank0)[2] + 273.15

        @self.Constraint(doc='initial accumulated volume of distillate')
        def AccVd_1(b):
            return b.AccVd[0] == 0

        @self.Constraint(doc='initial recovery rate')
        def recovery_rate_1(b):
            return b.accumulated_recovery_rate[0] == 0    

        @self.Constraint(doc='initial log mean temperature difference')
        def ATml_1(b):
            TCO = b.condenser_out_props[0].temperature - 273.15 * pyunits.K
            TEO = b.evaporator_out_props[0].temperature - 273.15 * pyunits.K
            return b.log_mean_temp_dif[0] == ((b.TEI - TCO) - (TEO - b.TCI)) / log((b.TEI - TCO)/(TEO - b.TCI + 1e-6))
    
        @self.Constraint(doc='initial average temperature of the feed flow in the heat source')
        def T_F_avg_1(b):
            return b.avg_feed_props[0].temperature == (b.TEI + 273.15*pyunits.K + b.condenser_out_props[0].temperature) / 2
        
        @self.Constraint(doc='initial average salinity of the feed flow')
        def S_F_avg_1(b):
            return (b.avg_feed_props[0].mass_frac_phase_comp["Liq","TDS"] == 
                    b.VAGMD_Models_NaCl(b.TEI, FFR0, b.TCI, S0, Ttank0)[3] /1000
                    )

        @self.Constraint(doc='initial flowrate of the feed flow')
        def V_F_avg_1(b):
            return b.avg_feed_props[0].flow_vol_phase["Liq"] == b.feed_props[0].flow_vol_phase["Liq"]
        
        @self.Constraint(doc='initial thermal power')
        def ThPower_1(b):        
            CpF = b.avg_feed_props[0].cp_mass_phase["Liq"]
            RhoF = b.avg_feed_props[0].dens_mass_phase["Liq"]
            TCO = b.condenser_out_props[0].temperature - 273.15 * pyunits.K
            return b.thermal_power[0] == pyunits.convert((FFR0 * CpF * (b.TEI - TCO)) * (RhoF),
                                                          to_units = pyunits.kW)

        @self.Constraint(doc='initial thermal energy consumption')
        def ThEnergy_1(b):        
            return b.thermal_energy[0] == pyunits.convert(b.thermal_power[0] * b.tR,
                                                          to_units = pyunits.kWh)

        @self.Constraint(doc='initial inlet cooling water temperature')
        def TCoolIn_1(b):        
            return b.cooling_in_props[0].temperature == Ttank0 + 273.15 * pyunits.K

        @self.Constraint(doc='initial outlet cooling water temperature')
        def TCoolOut_1(b):        
            return b.cooling_out_props[0].temperature == Ttank0 + 273.15 * pyunits.K

        @self.Constraint(doc='initial specific thermal energy consumption')
        def SETC_1(b):        
            return b.specific_thermal_energy_consumption[0] == 0

        @self.Constraint(doc='initial gain output ratio')
        def GOR_1(b):        
            return b.gain_output_ratio[0] == 0


        '''
        Calculation for each circulation (time step)
        '''
        @self.Constraint(self.flowsheet().config.time,
                         doc='Accumulated volume of distillate of each time step',
                         )
        def AccVd_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return b.AccVd[t] == b.AccVd[t-1] + pyunits.convert(b.permeate_flow_rate[t-1] * b.tR,
                                                                    to_units = pyunits.L)  
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Salinity of the batched feed water at each time step',
                         )
        def S_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.feed_props[t].conc_mass_phase_comp["Liq","TDS"] * (b.V0 - b.AccVd[t])
                    ==  b.feed_props[t-1].conc_mass_phase_comp["Liq","TDS"] * (b.V0 - b.AccVd[t-1])
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Feed flow rate remains the same during the recirculations',
                         )
        def FFR_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.feed_props[t].flow_vol_phase["Liq"] == b.feed_props[t-1].flow_vol_phase["Liq"]
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Temperature of the batched feed water at each time step',
                         )
        def Ttank_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.feed_props[t].temperature == (
                        pyunits.convert(b.feed_props[t-1].flow_vol_phase["Liq"] * b.tR, 
                                        to_units = pyunits.L)
                        * b.evaporator_out_props[t-1].temperature
                        + (b.V0 - b.AccVd[t-1]) * b.feed_props[t-1].temperature )
                        / (pyunits.convert(b.feed_props[t-1].flow_vol_phase["Liq"] * b.tR, 
                                        to_units = pyunits.L)
                        + b.V0 - b.AccVd[t-1])
                        )
            else:
                return Constraint.Skip
        
        @self.Constraint(self.flowsheet().config.time,
                         doc='Accumulated recovery rate after each circulation',
                         )
        def RR_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.accumulated_recovery_rate[t] == 1  
                        - b.feed_props[0].conc_mass_phase_comp["Liq","TDS"]
                        / b.feed_props[t].conc_mass_phase_comp["Liq","TDS"]
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Permeate flux at each time step',
                         )
        def PFlux_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.PFlux[t] == 
                b.VAGMD_Models_NaCl(b.TEI, 
                                    pyunits.convert(self.feed_props[0].flow_vol_phase["Liq"], 
                                                    to_units = pyunits.L / pyunits.h),
                                    b.TCI, 
                                    b.feed_props[t].conc_mass_phase_comp["Liq","TDS"], 
                                    b.feed_props[t].temperature)[0]
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Permeate flow rate at each time step',
                         )
        def PFR_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return b.permeate_flow_rate[t] == b.PFlux[t] * module_area
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Condenser outlet temperature at each time step',
                         )
        def TCO_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.condenser_out_props[t].temperature == 
                b.VAGMD_Models_NaCl(b.TEI, 
                                    pyunits.convert(self.feed_props[0].flow_vol_phase["Liq"], 
                                                    to_units = pyunits.L / pyunits.h),
                                    b.TCI, 
                                    b.feed_props[t].conc_mass_phase_comp["Liq","TDS"], 
                                    b.feed_props[t].temperature)[1]
                + 273.15
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Evaporator outlet temperature at each time step',
                         )
        def TEO_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.evaporator_out_props[t].temperature == 
                b.VAGMD_Models_NaCl(b.TEI, 
                                    pyunits.convert(self.feed_props[0].flow_vol_phase["Liq"], 
                                                    to_units = pyunits.L / pyunits.h),
                                    b.TCI, 
                                    b.feed_props[t].conc_mass_phase_comp["Liq","TDS"], 
                                    b.feed_props[t].temperature)[2]
                + 273.15
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Log mean temperature difference at each time step',
                         )
        def ATml_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                TCO = b.condenser_out_props[t].temperature - 273.15 * pyunits.K
                TEO = b.evaporator_out_props[t].temperature - 273.15 * pyunits.K
                return (b.log_mean_temp_dif[t] == ((b.TEI - TCO) - (TEO - b.TCI)) / log((b.TEI - TCO)/(TEO - b.TCI + 1e-6))
                )
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Average temperature of the feed flow in the heat source at each time step',
                         )
        def T_F_avg_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return b.avg_feed_props[t].temperature == (b.TEI + 273.15*pyunits.K + b.condenser_out_props[t].temperature) / 2
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Salinity of the feed flow in the heat source at each time step',
                         )
        def S_F_avg_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return (b.avg_feed_props[t].mass_frac_phase_comp["Liq","TDS"] == 
                        b.VAGMD_Models_NaCl(b.TEI, 
                                            pyunits.convert(self.feed_props[0].flow_vol_phase["Liq"], 
                                                            to_units = pyunits.L / pyunits.h),
                                            b.TCI, 
                                            b.feed_props[t].conc_mass_phase_comp["Liq","TDS"], 
                                            b.feed_props[t].temperature)[3] / 1000
                        )

            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Feed flow rate in the heat source at each time step',
                         )
        def V_F_avg_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return b.avg_feed_props[t].flow_vol_phase["Liq"] == b.feed_props[0].flow_vol_phase["Liq"]
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Thermal power at each time step',
                         )
        def ThPower_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                FFR = b.feed_props[t].flow_vol_phase["Liq"]
                CpF = b.avg_feed_props[t].cp_mass_phase["Liq"]
                RhoF = b.avg_feed_props[t].dens_mass_phase["Liq"]
                TCO = b.condenser_out_props[t].temperature - 273.15 * pyunits.K
                return b.thermal_power[t] == pyunits.convert((FFR * CpF * (b.TEI - TCO)) * (RhoF),
                                                              to_units = pyunits.kW)
            else:
                return Constraint.Skip

        @self.Constraint(self.flowsheet().config.time,
                         doc='Thermal energy at each time step',
                         )
        def ThEnergy_i(b, t):
            if t > 0 and t in b.flowsheet().config.time_set:
                return b.thermal_energy[t] == pyunits.convert(b.thermal_power[t] * b.tR,
                                                              to_units = pyunits.kWh)
            else:
                return Constraint.Skip


    '''
    Calculation models
    '''
    def VAGMD_Models_NaCl(self, TEI, FFR, TCI, SgL, Ttank):
    # Model parameters
        PFluxAS26      = [0.798993148477908, 0.314627216640160, 0.559805181621833, -0.146236734128216, -0.659197144919924, 0, 0.185658514024503, -0.107221706014227, 0, 0, -0.187626469717738, 0, 0, 0, 0.128420664686447]
        PFluxAS7_high  = [9.41014468300000, 0, 0, 0, -0.0188989390000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        PFluxAS7_low   = [4.82699491400000, 1.37479276300000, 1.91988177200000, -0.574212905000000, -0.641257664000000, 0.399259954000000, 0, 0, 0, 0, 0, 0, -0.588924321000000, 0, 0]
        TCOAS26        = [65.1084685465240, 9.15474718837607, -0.917460918908258, 0.480070517276181, -1.06168979823129, 0, 0, 0, 0.142552983811052, 0, 0, 0, 0, 0, 0]
        TCOAS7_high    = [67.0068599900000, 0, 0, 0, -0.0145469190000000, 0, 0, 0, 0, 0, 0, 0, 0]
        TCOAS7_low     = [58.8189910000000, 8.79691366400000, -2.06741401400000, 1.63187967600000, -0.914624645000000, -0.536574144000000, -0.249657477000000, 0.398973861000000, -0.153760262000000, 0.102355281000000, 0.696768080000000, -0.300582958000000, -0.557410173000000]
        TEOAS26        = [29.2261439896435, 0.569016152083381, 0.824636694807529, 4.62669502530487, 1.37222105534565, 0, 0, 0.220665657590258, 0, 0, 0, 0, 0, 0, 0]
        TEOAS7_high    = [36.2497021800000, 0, 0, 0, 0.0126951860000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        TEOAS7_low     = [34.4309511000000, 1.55140768500000, 1.85928314600000, 4.52887180500000, 1.10791196800000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.808447211000000]
        Coder          = [[-5.68487382500000, 0.0705622560000000, 0.000152146000000000],
                            [-1.58460599600000, 0.00102338700000000, 1.20000000000000e-06],
                            [-4.27697973100000, 0.175533630000000, -0.000178178000000000],
                            [-1.49331349500000, 0.0146627780000000, 5.62000000000000e-06],
                            [-7, 0.100000000000000, 0],
                            [-2.14285714285714, 0.00285714285714286, 0],
                            [-5, 0.200000000000000, 0],
                            [-1.33333333333333, 0.00950000000000000, 0]]
        a1 = 0.983930048493388
        a2 = -4.8359231959954E-04
        S_c = a1 * SgL + a2 * SgL**2 # [g/kg]

        CoderVars = [ [1, TEI, TEI**2],
                    [1, FFR, FFR**2],
                    [1, TCI, TCI**2],
                    [1, S_c, S_c**2]]

    # Model calculations
        if self.config.module_type == 7:
            if self.config.high_brine_salinity:
                TEI = 0
                FFR = 0
                TCI = 0
                S_r = S_c

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high   
            else:
                TEI = sum(CoderVars[0][j] * Coder[0][j] for j in range(len(Coder[0])))#np.dot(CoderVars[0], Coder[0])
                FFR = sum(CoderVars[1][j] * Coder[1][j] for j in range(len(Coder[1])))
                TCI = sum(CoderVars[2][j] * Coder[2][j] for j in range(len(Coder[2])))#np.dot(CoderVars[2], Coder[2])
                S_r = sum(CoderVars[3][j] * Coder[3][j] for j in range(len(Coder[3])))

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_low, TCOAS7_low, TEOAS7_low  


            VarsAS7  = [1, TEI, FFR, TCI, S_r, FFR*TEI, TCI*TEI, S_r*TEI, FFR*TCI, FFR*S_r, S_r*TCI, TEI**2, FFR**2, TCI**2, S_r**2]
            VarsAS7_TCO = [1, TEI, FFR, TCI, S_r, FFR*TEI, S_r*TEI, FFR*TCI, FFR*S_r, S_r*TCI, FFR**2, S_r**2, FFR**3]

            PFlux = sum(VarsAS7[j] * PFluxAS7[j] for j in range(len(VarsAS7)))
            TCO   = sum(VarsAS7_TCO[j] * TCOAS7[j] for j in range(len(VarsAS7_TCO)))
            TEO   = sum(VarsAS7[j] * TEOAS7[j] for j in range(len(VarsAS7)))       

        else:
            TEI = sum(CoderVars[0][j] * Coder[4][j] for j in range(len(Coder[0])))
            FFR = sum(CoderVars[1][j] * Coder[5][j] for j in range(len(Coder[1])))
            TCI = sum(CoderVars[2][j] * Coder[6][j] for j in range(len(Coder[2])))
            S_r = sum(CoderVars[3][j] * Coder[7][j] for j in range(len(Coder[3])))   

            VarsAS26 = [1, TEI, FFR, TCI, S_r, TCI*TEI, FFR*TEI, S_r*TEI, FFR*TCI, S_r*TCI, FFR*S_r, TEI**2, FFR**2, TCI**2, S_r**2]

            PFlux = sum(VarsAS26[j] * PFluxAS26[j] for j in range(len(VarsAS26)))
            TCO   = sum(VarsAS26[j] * TCOAS26[j] for j in range(len(VarsAS26)))
            TEO   = sum(VarsAS26[j] * TEOAS26[j] for j in range(len(VarsAS26)))

        return [PFlux, TCO, TEO, S_c]
