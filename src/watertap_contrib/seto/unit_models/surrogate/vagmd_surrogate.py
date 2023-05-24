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
    Block,
    Set,
    Var,
    Param,
    value,
    Suffix,
    ConcreteModel,
    PositiveIntegers,
    Reference,
    Constraint,
    units as pyunits,
    exp,
    log,
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
from watertap.core.util.model_diagnostics.infeasible import *

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)


@declare_process_block_class("VAGMDSurrogate")
class VAGMDData(UnitModelBlockData):
    """
    Vacuumed Membrane distillation (air-gapped) - batch operation model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be True",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package_seawater",
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
        "property_package_water",
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
            default="AS7C1.5L",
            domain=In(["AS7C1.5L", "AS26C7.2L"]),
            doc="""Selection of module type (7 for AS7C1.5L and 26 for AS26C7.2L)""",
        ),
    )
    CONFIG.declare(
        "cooling_system_type",
        ConfigValue(
            default="closed",
            domain=In(["open", "closed"]),
            doc="""Selection of cooling system type (open or closed)""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package_seawater.get_metadata().get_derived_units

        """
        Model parameters
        """
        self.gas_constant_1 = Param(
            initialize = 8.314,
            units = pyunits.J / pyunits.K / pyunits.mol,
            doc = "Gas constan in J/K mol"
        )

        self.gas_constant_2 = Param(
            initialize = 0.082,
            units = pyunits.atm * pyunits.L / pyunits.K /pyunits.mol,
            doc = "Gas constant in atm * L /K mol"
        )

        self.pump_efficiency = Param(
            initialize = 0.6,
            units = pyunits.dimensionless,
            doc = "Pump efficiency"
        )

        self.heat_exchanger_area = Param(
            initialize = 1.34,
            units = pyunits.m**2,
            doc = "Effective heat transfer coefficient"
        )

        self.cooling_flow_rate = Param(
            initialize = 1265,
            units = pyunits.L / pyunits.h,
            doc = "Cooling water volumetric flow rate, fixed to 1265L/h, to maintain vacuum pressure inside the MD module around 200 mbar"
        )

        self.thermal_heat_transfer_coeff = Param(
            initialize = 3168,
            units = pyunits.W / pyunits.m**2 / pyunits.K,
            doc = "Overall heat transfer coefficient"
        )

        self.cooling_flow_pressure_drop = Param(
            initialize = 170,
            units = pyunits.mbar,
            doc = "Cooling flow pressure drop"
        )

        """
        MD Module type and corresponding area
        """
        if self.config.module_type == "AS7C1.5L":
            self.module_area = Param(
                initialize=7.2, units=pyunits.m**2, doc="Area of module AS7C1.5L"
            )
        else:  # module_type = "AS26C7.2L"
            self.module_area = Param(
                initialize=25.92, units=pyunits.m**2, doc="Area of module AS26C7.2L"
            )

        """
        Intermediate variables
        """
        self.permeate_flux = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.L / pyunits.h / pyunits.m**2,
            doc="Permeate flux",
        )

        self.feed_flow_pressure_drop = Var(
            initialize = 100,
            bounds = (0,None),
            units = pyunits.mbar,
            doc = "Feed flow pressure drop"
        )

        self.log_mean_temp_dif = Var(
            initialize=0,
            bounds=(0, None),
            units=pyunits.C,
            doc="Log mean temperature difference in the MD module",
        )

        # self.thermal_resistance_hot = Var(
        #     initialize = 1e3,
        #     bounds = (0, None),
        #     units = pyunits.W / pyunits.k,
        #     doc = "Thermal resistance on the hot side"
        # )

        # self.thermal_resistance_cold = Var(
        #     initialize = 1e3,
        #     bounds = (0, None),
        #     units = pyunits.W / pyunits.k,
        #     doc = "Thermal resistance on the cold side"
        # )

        # self.number_transfer_units = Var(
        #     initialize = 10,
        #     bounds = (0, None),
        #     units = pyunits.dimensionless,
        #     doc = "The number of transfer units"
        # )

    

        """
        Output variables
        """
        self.thermal_power = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power",
        )

        self.feed_pump_power_elec = Var(
            initialize = 0.005,
            bounds = (0, None),
            units = pyunits.kW,
            doc = "Electric power for pumping feed water"
        )

        self.cooling_pump_power_elec = Var(
            initialize = 0.005,
            bounds = (0, None),
            units = pyunits.kW,
            doc = "Electric power for pumping cooling water"
        )

        """
        Add block for the batched feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_seawater
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict
        )

        """
        Add block for the permeated water
        """
        tmp_dict["defined_state"] = False

        self.permeate_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of permeated water",
            **tmp_dict
        )
        # salinity in permeate is zero
        self.permeate_props[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)

        # permeate temperature is assumed the same as feed water
        @self.Constraint(doc="distillate temperature")
        def eq_distillate_temp(b):
            return b.permeate_props[0].temperature == b.feed_props[0].temperature

        """
        Add block for water at the evaporator inlet
        """
        tmp_dict["defined_state"] = True

        self.evaporator_in_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of evaporator inlet water",
            **tmp_dict
        )

        """
        Add block for water at the evaporator outlet
        """
        tmp_dict["defined_state"] = False

        self.evaporator_out_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        """
        Add block for water at the condenser inlet
        """
        tmp_dict["defined_state"] = True

        self.condenser_in_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of condenser inlet water",
            **tmp_dict
        )

        """
        Add block for water at the condenser outlet
        """
        tmp_dict["defined_state"] = False

        self.condenser_out_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of condenser outlet water",
            **tmp_dict
        )

        """
        Add block for the average status of the feed flow
        """
        self.avg_feed_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Average properties of the feed water",
            **tmp_dict
        )

        tmp_dict["parameters"] = self.config.property_package_water
        """
        Add block for the average status in the condenser
        """
        self.avg_condenser_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Average properties of water in the condenser",
            **tmp_dict
        )    
        self.avg_condenser_props[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        self.avg_condenser_props[0].pressure.fix(101325)

        """
        Constraint equations
        """
        # Set alias for state properties
        FFR = pyunits.convert(
            self.feed_props[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.h
        )
        S = self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        TEI = self.evaporator_in_props[0].temperature
        TCI = self.condenser_in_props[0].temperature
        Ttank = self.feed_props[0].temperature
        TEO = self.evaporator_out_props[0].temperature
        TCO = self.condenser_out_props[0].temperature

        @self.Constraint(doc="Average temperature of the feed flow in the heat source")
        def eq_avg_temp_heat_tank(b):
            return b.avg_feed_props[0].temperature == (TEI + TCO) / 2

        @self.Constraint(doc="Average temperature of the flow in the condenser")
        def eq_avg_temp_condenser(b):
            return b.avg_condenser_props[0].temperature == (TCI + TCO) / 2

        @self.Constraint(doc="Average salinity of the feed flow")
        def eq_avg_salinity_feed_tank(b):
            return (
                b.avg_feed_props[0].mass_frac_phase_comp["Liq", "TDS"]
                == b._get_membrane_performance(TEI, FFR, TCI, S)[3] / 1000
            )

        @self.Constraint(doc="Flowrate of the average feed flow block")
        def eq_feed_volumetric_flow_rate(b):
            return (
                b.avg_feed_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="Permeate flux")
        def eq_permeate_flux(b):
            return b.permeate_flux == b._get_membrane_performance(TEI, FFR, TCI, S)[0]

        @self.Constraint(doc="Feed flow pressure drop")
        def eq_feed_flow_pressure_drop(b):
            return (
                    b.feed_flow_pressure_drop == 
                    b._get_pressure_drop(pyunits.convert(b.feed_props[0].flow_vol_phase["Liq"], 
                                                         to_units=pyunits.L / pyunits.h), 
                                         b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"])
                    )

        @self.Constraint(doc="Evaporatore outlet temperature")
        def eq_evaporator_outlet_temp(b):
            return (
                b.evaporator_out_props[0].temperature
                == b._get_membrane_performance(TEI, FFR, TCI, S)[2] + 273.15
            )

        @self.Constraint(doc="Condenser outlet temperature")
        def eq_condenser_outlet_temp(b):
            return (
                b.condenser_out_props[0].temperature
                == b._get_membrane_performance(TEI, FFR, TCI, S)[1] + 273.15
            )

        @self.Constraint(doc="Permeate flow rate")
        def eq_permeate_volumetric_flow_rate(b):
            return b.permeate_props[0].flow_vol_phase["Liq"] == pyunits.convert(
                b.permeate_flux * self.module_area,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Constraint(doc="Condenser mass flow rate")
        def eq_condenser_mass_flow_rate(b):
            return b.avg_condenser_props[0].flow_mass_phase_comp["Liq", "H2O"] == b.permeate_props[0].flow_mass_phase_comp["Liq", "H2O"]

        @self.Constraint(doc="initial log mean temperature difference")
        def eq_log_mean_temp_dif(b):
            return b.log_mean_temp_dif == ((TEI - TCO) - (TEO - TCI)) / log(
                (TEI - TCO) / (TEO - TCI + 1e-6)
            )

        @self.Constraint(doc="Thermal power")
        def eq_thermal_power(b):
            CpF = b.avg_feed_props[0].cp_mass_phase["Liq"]
            RhoF = b.avg_feed_props[0].dens_mass_phase["Liq"]
            return b.thermal_power == pyunits.convert(
                (FFR * CpF * (TEI - TCO)) * (RhoF), to_units=pyunits.kW
            )

        @self.Constraint(doc="Electric power for pumping feed water")
        def eq_feed_pump_power_elec(b):
            return (
                b.feed_pump_power_elec == pyunits.convert(
                                        b.gas_constant_1 / b.gas_constant_2 
                                        * b.feed_flow_pressure_drop / b.pump_efficiency
                                        * b.feed_props[0].flow_vol_phase["Liq"],
                                        to_units = pyunits.kW)
            )

        @self.Constraint(doc="Electric power for pumping cooling water")
        def eq_cooling_pump_power_elec(b):
            return (
                b.cooling_pump_power_elec == pyunits.convert(
                                        b.gas_constant_1 / b.gas_constant_2 
                                        * b.cooling_flow_pressure_drop / b.pump_efficiency
                                        * b.cooling_flow_rate,
                                        to_units = pyunits.kW)
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.permeate_flux) is None:
            iscale.set_scaling_factor(self.permeate_flux, 1)

        if iscale.get_scaling_factor(self.log_mean_temp_dif) is None:
            iscale.set_scaling_factor(self.log_mean_temp_dif, 1e-1)

        if iscale.get_scaling_factor(self.thermal_power) is None:
            iscale.set_scaling_factor(self.thermal_power, 1e-1)

        if iscale.get_scaling_factor(self.feed_flow_pressure_drop) is None:
            iscale.set_scaling_factor(self.feed_flow_pressure_drop, 1e-1)

        if iscale.get_scaling_factor(self.feed_pump_power_elec) is None:
            iscale.set_scaling_factor(self.feed_pump_power_elec, 1e3)

        if iscale.get_scaling_factor(self.cooling_pump_power_elec) is None:
            iscale.set_scaling_factor(self.cooling_pump_power_elec, 1e3)


        # Transforming constraint

        sf = iscale.get_scaling_factor(self.feed_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_distillate_temp, sf)

        sf = iscale.get_scaling_factor(self.avg_feed_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_avg_temp_heat_tank, sf)

        sf = iscale.get_scaling_factor(self.avg_condenser_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_avg_temp_condenser, sf)

        sf = iscale.get_scaling_factor(
            self.avg_feed_props[0].mass_frac_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.eq_avg_salinity_feed_tank, sf)

        sf = iscale.get_scaling_factor(
            self.avg_feed_props[0].flow_vol_phase["Liq"]
        )
        iscale.constraint_scaling_transform(self.eq_feed_volumetric_flow_rate, sf)

        sf = iscale.get_scaling_factor(self.permeate_flux)
        iscale.constraint_scaling_transform(self.eq_permeate_flux, sf)

        sf = iscale.get_scaling_factor(
            self.evaporator_out_props[0].temperature
        )
        iscale.constraint_scaling_transform(self.eq_evaporator_outlet_temp, sf)

        sf = iscale.get_scaling_factor(
            self.condenser_out_props[0].temperature
        )
        iscale.constraint_scaling_transform(self.eq_condenser_outlet_temp, sf)

        sf = iscale.get_scaling_factor(
            self.permeate_props[0].flow_vol_phase["Liq"]
        )
        iscale.constraint_scaling_transform(self.eq_permeate_volumetric_flow_rate, sf)

        sf = iscale.get_scaling_factor(
            self.avg_condenser_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )        
        iscale.constraint_scaling_transform(self.eq_condenser_mass_flow_rate, sf)

        sf = iscale.get_scaling_factor(self.log_mean_temp_dif)
        iscale.constraint_scaling_transform(self.eq_log_mean_temp_dif, sf)

        sf = iscale.get_scaling_factor(self.thermal_power)
        iscale.constraint_scaling_transform(self.eq_thermal_power, sf)

        sf = iscale.get_scaling_factor(self.feed_flow_pressure_drop)
        iscale.constraint_scaling_transform(self.eq_feed_flow_pressure_drop, sf)

        sf = iscale.get_scaling_factor(self.feed_pump_power_elec)
        iscale.constraint_scaling_transform(self.eq_feed_pump_power_elec, sf)

        sf = iscale.get_scaling_factor(self.cooling_pump_power_elec)
        iscale.constraint_scaling_transform(self.eq_cooling_pump_power_elec, sf)

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")
        init_log.info_low("Starting initialization...")
        # Iniitialize feed properties

        flags = self.feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )

        init_log.info_low("Initialization for feed Completed.")

        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.feed_props[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value
        init_log.info_low("Starting initialization of evaporator_in_props ")
        self.evaporator_in_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of evaporator_out_props ")
        self.evaporator_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_low("Starting initialization of condenser_in_props ")
        self.condenser_in_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of condenser_out_props ")
        self.condenser_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of avg_feed_props ")
        self.avg_feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        state_args_water = {}
        state_dict_water = self.avg_condenser_props[
            self.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_water.keys():
            if state_dict_water[k].is_indexed():
                state_args_water[k] = {}
                for m in state_dict_water[k].keys():
                    state_args_water[k][m] = state_dict_water[k][m].value
            else:
                state_args_water[k] = state_dict_water[k].value

        for p, j in self.avg_condenser_props.phase_component_set:
            if p == "Vap" and j == "H2O":
                state_args_water["flow_mass_phase_comp"][(p, j)] = 0

        init_log.info_low("Starting initialization of avg_condenser_props ")
        self.avg_condenser_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_water,
        )

        state_args_permeate = deepcopy(state_args)
        for p, j in self.permeate_props.phase_component_set:
            if j == "TDS":
                state_args_permeate["flow_mass_phase_comp"][(p, j)] = 0

        init_log.info_low("Starting initialization of permeate_props ")
        self.permeate_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,
        )

        opt = get_solver(solver, optarg)
        print('unit dof', degrees_of_freedom(self))
        assert degrees_of_freedom(self) == 0        

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            print_close_to_bounds(self)
            print_infeasible_constraints(self) 
        init_log.info("Initialization status {}.".
                      format(idaeslog.condition(res)))

    """
    Equation to calculate pressure drop
    """
    def _get_pressure_drop(self, flow_rate, salinity):
        if self.config.module_type == "AS7C1.5L":
            coefficients = [
                -158.2007422,
                0.39402609,
                0,
                0.000585345,
                8.93618e-5,
                -0.000287828,
            ]
        else:  # self.config.module_type == 'AS26C7.2L'
            coefficients = [
                -72.53793298,
                0.110437201,
                0,
                0.000643495,
                0.000189924,
                -0.001111447,
            ]

        return (
            coefficients[0]
            + coefficients[1] * flow_rate
            + coefficients[2] * salinity
            + coefficients[3] * flow_rate * salinity
            + coefficients[4] * flow_rate**2
            + coefficients[5] * salinity**2
        )

    """
    Surrogate equations for membrane performance
    """

    def _get_membrane_performance(self, TEI, FFR, TCI, SgL):
        # Model parameters
        PFluxAS26 = [
            0.798993148477908,
            0.314627216640160,
            0.559805181621833,
            -0.146236734128216,
            -0.659197144919924,
            0,
            0.185658514024503,
            -0.107221706014227,
            0,
            0,
            -0.187626469717738,
            0,
            0,
            0,
            0.128420664686447,
        ]
        PFluxAS7_high = [
            9.41014468300000,
            0,
            0,
            0,
            -0.0188989390000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        PFluxAS7_low = [
            4.82699491400000,
            1.37479276300000,
            1.91988177200000,
            -0.574212905000000,
            -0.641257664000000,
            0.399259954000000,
            0,
            0,
            0,
            0,
            0,
            0,
            -0.588924321000000,
            0,
            0,
        ]
        TCOAS26 = [
            65.1084685465240,
            9.15474718837607,
            -0.917460918908258,
            0.480070517276181,
            -1.06168979823129,
            0,
            0,
            0,
            0.142552983811052,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_high = [
            67.0068599900000,
            0,
            0,
            0,
            -0.0145469190000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_low = [
            58.8189910000000,
            8.79691366400000,
            -2.06741401400000,
            1.63187967600000,
            -0.914624645000000,
            -0.536574144000000,
            -0.249657477000000,
            0.398973861000000,
            -0.153760262000000,
            0.102355281000000,
            0.696768080000000,
            -0.300582958000000,
            -0.557410173000000,
        ]
        TEOAS26 = [
            29.2261439896435,
            0.569016152083381,
            0.824636694807529,
            4.62669502530487,
            1.37222105534565,
            0,
            0,
            0.220665657590258,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_high = [
            36.2497021800000,
            0,
            0,
            0,
            0.0126951860000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_low = [
            34.4309511000000,
            1.55140768500000,
            1.85928314600000,
            4.52887180500000,
            1.10791196800000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0.808447211000000,
        ]
        Coder = [
            [-5.68487382500000, 0.0705622560000000, 0.000152146000000000],
            [-1.58460599600000, 0.00102338700000000, 1.20000000000000e-06],
            [-4.27697973100000, 0.175533630000000, -0.000178178000000000],
            [-1.49331349500000, 0.0146627780000000, 5.62000000000000e-06],
            [-7, 0.100000000000000, 0],
            [-2.14285714285714, 0.00285714285714286, 0],
            [-5, 0.200000000000000, 0],
            [-1.33333333333333, 0.00950000000000000, 0],
        ]
        a1 = 0.983930048493388
        a2 = -4.8359231959954e-04
        S_c = a1 * SgL + a2 * SgL**2  # [g/kg]

        TEI -= 273.15
        TCI -= 273.15
        
        CoderVars = [
            [1, TEI, TEI**2],
            [1, FFR, FFR**2],
            [1, TCI, TCI**2],
            [1, S_c, S_c**2],
        ]

        # Model calculations
        if self.config.module_type == "AS7C1.5L":
            if self.config.high_brine_salinity:
                TEI = 0
                FFR = 0
                TCI = 0
                S_r = S_c

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high
            else:
                TEI = sum(CoderVars[0][j] * Coder[0][j] for j in range(len(Coder[0])))
                FFR = sum(CoderVars[1][j] * Coder[1][j] for j in range(len(Coder[1])))
                TCI = sum(CoderVars[2][j] * Coder[2][j] for j in range(len(Coder[2])))
                S_r = sum(CoderVars[3][j] * Coder[3][j] for j in range(len(Coder[3])))

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_low, TCOAS7_low, TEOAS7_low

            VarsAS7 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                TCI * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]
            VarsAS7_TCO = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                FFR**2,
                S_r**2,
                FFR**3,
            ]

            PFlux = sum(VarsAS7[j] * PFluxAS7[j] for j in range(len(VarsAS7)))
            TCO = sum(VarsAS7_TCO[j] * TCOAS7[j] for j in range(len(VarsAS7_TCO)))
            TEO = sum(VarsAS7[j] * TEOAS7[j] for j in range(len(VarsAS7)))

        else:
            TEI = sum(CoderVars[0][j] * Coder[4][j] for j in range(len(Coder[0])))
            FFR = sum(CoderVars[1][j] * Coder[5][j] for j in range(len(Coder[1])))
            TCI = sum(CoderVars[2][j] * Coder[6][j] for j in range(len(Coder[2])))
            S_r = sum(CoderVars[3][j] * Coder[7][j] for j in range(len(Coder[3])))

            VarsAS26 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                TCI * TEI,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                S_r * TCI,
                FFR * S_r,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]

            PFlux = sum(VarsAS26[j] * PFluxAS26[j] for j in range(len(VarsAS26)))
            TCO = sum(VarsAS26[j] * TCOAS26[j] for j in range(len(VarsAS26)))
            TEO = sum(VarsAS26[j] * TEOAS26[j] for j in range(len(VarsAS26)))

        return [PFlux, TCO, TEO, S_c]