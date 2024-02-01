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
Thermal Energy Storage Tank

Molten salt Properties:
From https://github.com/gmlc-dispatches/dispatches/blob/main/dispatches/properties/solarsalt_properties.py
Ref: (2015) Chang et al, Energy Procedia 69, 779 - 789

"""

# Import Pyomo libraries
from copy import deepcopy
from pyomo.common.config import ConfigBlock, ConfigValue, In

from pyomo.environ import (
    Var,
    Param,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)

from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from idaes.core.util.config import is_physical_parameter_block
from watertap.core import InitializationMixin

from watertap_contrib.reflo.costing.solar.thermal_energy_storage import (
    cost_tes,
)

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
        "storage material", ConfigValue(default="salt", doc="Thermal storage material")
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
        self.add_inlet_port(name="tes_hx_inlet", block=self.hx_inlet_block)
        self.add_outlet_port(name="tes_hx_outlet", block=self.hx_outlet_block)
        self.add_inlet_port(name="tes_process_inlet", block=self.process_inlet_block)
        self.add_outlet_port(name="tes_process_outlet", block=self.process_outlet_block)

        # TODO: Add multiple tanks?

        self.salt_csp = Param(
            initialize=1443,
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="Specific heat capacity of the salt in the TES",
        )

        self.hours_storage = Var(
            initialize=8, bounds=(0, 24), units=pyunits.h, doc="Hours of storage"
        )

        self.heat_load = Var(
            initialize=120,
            bounds=(0, None),
            units=pyunits.MW,
            doc="Design thermal output rate",
        )

        self.thermal_energy_capacity = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.MW * pyunits.h,
            doc="Thermal energy storage capacity for hours of storage at design thermal output rate",
        )

        self.tes_volume = Var(
            initialize=100,
            bounds=(0, None),
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

        self.tes_initial_temperature = Var(
            initialize=30 + 273.15,
            units=pyunits.K,
            bounds=(25 + 273.15, 150 + 273.15),
            doc="Temperature of the thermal storage tank initially from the previous time step",
        )

        self.tes_temperature = Var(
            self.flowsheet().config.time,
            initialize=30 + 273.15,
            units=pyunits.K,
            bounds=(25 + 273.15, 150 + 273.15),
            doc="Temperature of the thermal storage tank to track convective heat losses",
        )

        self.dt = Var(initialize=3600, units=pyunits.s, doc="Time step for multiperiod")

        ## TODO Convective heat loss as a function of tank temperature
        # Can be updated to be a function of temperature

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

        self.pump_power = Param(
            initialize=1, units=pyunits.W, mutable=True, doc="Pump power"
        )

        self.pump_eff = Param(
            initialize=1,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Pump efficiency",
        )

        self.electricity = Var(initialize=1, units=pyunits.W, doc="Total electricity")

        self.temperature_design = Param(
            initialize=80 + 273.15, units=pyunits.K, doc="ambient temperature"
        )

        self.temperature_ambient = Param(
            initialize=25 + 273.15, units=pyunits.K, doc="ambient temperature"
        )

        @self.Constraint(doc="Calculate mass of solar salt")
        def eq_salt_mass(b):
            return b.salt_mass == b.tes_volume * b.salt_packing_density

        # Constraint to calculate the total heat entering
        @self.Constraint(self.flowsheet().config.time)
        def eq_heat_in(b, t):
            return b.heat_in[t] == (
                b.hx_inlet_block[t].enth_flow_phase["Liq"]
                + b.process_inlet_block[t].enth_flow_phase["Liq"]
            )

        # Constraint to calculate the total heat entering
        @self.Constraint(self.flowsheet().config.time)
        def eq_heat_out(b, t):
            return b.heat_out[t] == (
                b.hx_outlet_block[t].enth_flow_phase["Liq"]
                + b.process_outlet_block[t].enth_flow_phase["Liq"]
            )

        # the temperature of the tank after each time step
        @self.Constraint(self.flowsheet().config.time)
        def eq_tes_temp(b, t):
            return b.tes_temperature[t] == b.tes_initial_temperature + 1 / (
                b.salt_csp * b.salt_mass
            ) * (b.heat_in[t] * b.dt - b.heat_out[t] * b.dt)

        @self.Constraint()
        def eq_thermal_capacity(b):
            return b.thermal_energy_capacity == b.hours_storage * b.heat_load

        @self.Constraint()
        def eq_tes_volume(b):
            return b.tes_volume == b.thermal_energy_capacity / (
                b.salt_csp
                * b.salt_packing_density
                * (b.temperature_design - b.temperature_ambient)
            )

        @self.Constraint(doc="Pump power")
        def eq_P_pump(b):
            return b.electricity == (b.pump_power / b.pump_eff)

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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
        init_log.info("TES initialization status {}.".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.tes_volume) is None:
            sf = iscale.get_scaling_factor(self.tes_volume, default=1e-1)
            iscale.set_scaling_factor(self.tes_volume, sf)

        if iscale.get_scaling_factor(self.heat_in) is None:
            sf = iscale.get_scaling_factor(self.heat_in, default=1e-3)
            iscale.set_scaling_factor(self.heat_in, sf)

        if iscale.get_scaling_factor(self.heat_out) is None:
            sf = iscale.get_scaling_factor(self.heat_out, default=1e-3)
            iscale.set_scaling_factor(self.heat_out, sf)

        if iscale.get_scaling_factor(self.tes_initial_temperature) is None:
            sf = iscale.get_scaling_factor(self.tes_initial_temperature, default=1e-2)
            iscale.set_scaling_factor(self.tes_initial_temperature, sf)

        if iscale.get_scaling_factor(self.tes_temperature) is None:
            sf = iscale.get_scaling_factor(self.tes_temperature, default=1e-2)
            iscale.set_scaling_factor(self.tes_temperature, sf)

        if iscale.get_scaling_factor(self.dt) is None:
            sf = iscale.get_scaling_factor(self.dt, default=1e-3)
            iscale.set_scaling_factor(self.dt, sf)

        if iscale.get_scaling_factor(self.salt_mass) is None:
            sf = iscale.get_scaling_factor(self.salt_mass, default=1e-2)
            iscale.set_scaling_factor(self.salt_mass, sf)

        if iscale.get_scaling_factor(self.electricity) is None:
            sf = iscale.get_scaling_factor(self.electricity, default=1e-1)
            iscale.set_scaling_factor(self.electricity, sf)

        if iscale.get_scaling_factor(self.heat_load) is None:
            sf = iscale.get_scaling_factor(self.heat_load, default=1e-3)
            iscale.set_scaling_factor(self.heat_load, sf)

        if iscale.get_scaling_factor(self.thermal_energy_capacity) is None:
            sf = iscale.get_scaling_factor(self.thermal_energy_capacity, default=1e-5)
            iscale.set_scaling_factor(self.thermal_energy_capacity, sf)

    @property
    def default_costing_method(self):
        return cost_tes
