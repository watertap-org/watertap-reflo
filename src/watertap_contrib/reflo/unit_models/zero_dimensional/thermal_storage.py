#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""
Thermal Energy Storage Tank
"""

# Import Pyomo libraries
from copy import deepcopy
from pyomo.common.config import ConfigBlock, ConfigValue, In

from pyomo.environ import (
    Var,
    Constraint,
    check_optimal_termination,
    Param,
    value,
    log,
    units as pyunits,
    Expression,
    NonNegativeReals,
)

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)

from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

from idaes.core.util.config import is_physical_parameter_block
from watertap.core import InitializationMixin

_log = idaeslog.getLogger(__name__)
__author__ = "Mukta Hardikar"

@declare_process_block_class("ThermalEnergyStorage")
class ThermalEnergyStorageData(UnitModelBlockData):
    """
    Thermal energy storage model Class
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="Battery does not support dynamic models, thus this must be False",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag",
            doc="""Battery does not have defined volume, thus this must be False.""",
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
        "storage material", 
        ConfigValue(default="salt", 
                    doc="Thermal storage material")
    )

    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
            **default** - EnergyBalanceType.useDefault.
            **Valid values:** {
            **EnergyBalanceType.useDefault - refer to property package for default
            balance type
            **EnergyBalanceType.none** - exclude energy balances,
            **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
            **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
            **EnergyBalanceType.energyTotal** - single energy balance for material,
            **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
            **default** - MomentumBalanceType.pressureTotal.
            **Valid values:** {
            **MomentumBalanceType.none** - exclude momentum balances,
            **MomentumBalanceType.pressureTotal** - single pressure balance for material,
            **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
            **MomentumBalanceType.momentumTotal** - single momentum balance for material,
            **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    def build(self):
        super().build()

        # Add hx inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  
        self.hx_inlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of exit stream from the heat source",
            **tmp_dict,
        )

        # Add hx outlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # block is not an inlet
        self.hx_outlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of stream exiting the TES back to the heat source",
            **tmp_dict,
        )

        # Add process inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True 
        self.process_inlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of exit stream from the heat source",
            **tmp_dict,
        )

        # Add process outlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # block is not an inlet
        self.process_outlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of stream exiting the TES back to the heat source",
            **tmp_dict,
        )

        # Add Ports
        self.add_inlet_port(name='tes_hx_inlet',block = self.hx_inlet_block)
        self.add_outlet_port(name='tes_hx_outlet',block = self.hx_outlet_block)
        self.add_inlet_port(name='tes_process_inlet',block = self.process_inlet_block)
        self.add_outlet_port(name='tes_process_outlet',block = self.process_outlet_block)

        # From https://github.com/gmlc-dispatches/dispatches/blob/main/dispatches/properties/solarsalt_properties.py
        # **** Requires Temperature in K
        # Ref: (2015) Chang et al, Energy Procedia 69, 779 - 789
        # Specific Heat Capacity as a function of Temperature, J/kg/K
        # def specific_heat_cp(model):
        #    return model.specific_heat_cp
        #       == 1443 + 0.172 * model.temperature

        self.salt_csp = Param(
            initialize=1443,
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of the salt in the TES",
        )

        self.tes_diameter = Var(
            initialize=10, 
            units=pyunits.m, 
            doc="Diameter of the thermal storage tank"
        )

        self.tes_H_D_ratio = Param(
            initialize = 2,
            doc='Height to diameter ratio of tank'
        )
        # self.tes_height = Var(
        #     initialize=10, 
        #     units=pyunits.m, 
        #     doc="Height of the thermal storage tank"
        # )

        self.tes_volume = Var(
            initialize=100,
            units=pyunits.m**3,
            doc="Volume of the thermal storage tank",
        )

        self.heat_in = Var(
            self.flowsheet().config.time,
            initialize=0,
            units=pyunits.W,
            doc="Thermal energy from a solar power source",
        )

        self.heat_out = Var(
            self.flowsheet().config.time,
            initialize=0,
            units=pyunits.W,
            doc="Thermal energy exiting the thermal storage tank",
        )

        self.tes_initial_temp = Var(
            initialize=30+273.15,
            units=pyunits.K,
            bounds=(25+273.15, 150+273.15),
            doc="Temperature of the thermal storage tank initially from the previous time step",
        )

        self.tes_temp = Var(
            self.flowsheet().config.time,
            initialize=30+273.15,
            units=pyunits.K,
            bounds=(25+273.15, 150+273.15),
            doc="Temperature of the thermal storage tank to track convective heat losses",
        )

        self.dt = Var(
            initialize=1, 
            units=pyunits.s, 
            doc="Time step for multiperiod"
        )

        ## TODO Convective heat loss as a function of tank temperature
        # Can be updated to be a function of temperature
        # self.h_conv = Param(
        #     initialize = 1,
        #     units = pyunits.J/pyunits.m**2/pyunits.K
        # )

        # Excess heat that needs to be dissipated because of limit on tes temp
        # self.heat_dissipation = Var(
        #     self.flowsheet().config.time,
        #     initialize = 10,
        # )

        self.salt_mass = Var(
            initialize=100,
            units=pyunits.kg,
            doc="Mass of salt in the thermal energy storage tank",
        )

        self.salt_packing_density = Param(
            initialize=1,
            units=pyunits.kg / pyunits.m**3,
            doc="Packing density of salt",
        )

        # constraints
        # tank volume- Update to have a ratio between height and diameter so that volume is optimized
        @self.Constraint(doc="Calculate optimal tank size")
        def eq_tes_volume(b):
            return b.tes_volume == b.tes_diameter*b.tes_H_D_ratio * 3.14 * b.tes_diameter**2 / 4

        @self.Constraint(doc="Calculate mass of solar salt")
        def eq_salt_mass(b):
            return b.salt_mass == b.tes_volume * b.salt_packing_density

        # Constraint to calculate the total heat entering
        @self.Constraint(self.flowsheet().config.time)
        def eq_heat_in(b,t):
            return b.heat_in[t] == (
                b.hx_inlet_block[t].enth_flow_phase['Liq']
                + b.process_inlet_block[t].enth_flow_phase['Liq']
            )
                
        # Constraint to calculate the total heat entering
        @self.Constraint(self.flowsheet().config.time)
        def eq_heat_out(b,t):
            return b.heat_out[t] == (
                b.hx_outlet_block[t].enth_flow_phase['Liq']
                + b.process_outlet_block[t].enth_flow_phase['Liq']
            )

        # the temperature of the tank after each time step
        @self.Constraint(self.flowsheet().config.time)
        def eq_tes_temp(b, t):
            return b.tes_temp[t] == b.tes_initial_temp + 1 / (
                b.salt_csp * b.salt_mass
            ) * (b.heat_in[t]*b.dt - b.heat_out[t]*b.dt)

        
    def initialize_build(
            self, 
            state_args=None,
            outlvl=idaeslog.NOTSET, 
            solver=None, 
            optarg=None
    ):

        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Create solver
        opt = get_solver(solver=solver, options=optarg)    
                
        self.hx_inlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.hx_inlet_block[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args_out = deepcopy(state_args)

        self.hx_outlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        self.process_inlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        self.process_outlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        # solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "TES initialization status {}.".format(idaeslog.condition(res))
        )




