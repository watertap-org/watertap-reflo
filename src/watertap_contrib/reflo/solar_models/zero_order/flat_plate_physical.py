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

from copy import deepcopy
from pyomo.environ import Var, Param, Expression, log, exp, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import declare_process_block_class

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap_contrib.reflo.costing.solar.flat_plate import (
    cost_flat_plate,
)

__author__ = "Mukta Hardikar, Matthew Boyd"


@declare_process_block_class("FlatPlatePhysical")
class FlatPlatePhysicalData(SolarEnergyBaseData):
    """
    Physical model for flat plate
    based on equations in Solar Engineering of Thermal Processes, Duffie and Beckman, 4th ed.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

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

        self._tech_type = "flat_plate"
        # Add inlet state block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.inlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of inlet stream",
            **tmp_dict,
        )

        # Add outlet state block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # block is not an inlet
        self.outlet_block = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet stream",
            **tmp_dict,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.inlet_block)
        self.add_outlet_port(name="outlet", block=self.outlet_block)

        # ==========PARAMETERS==========

        self.number_collectors = Param(
            initialize=1,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Number of collectors in array",
        )

        self.FR_ta = Param(
            initialize=0.689,  # optical gain "a" in Hottel-Whillier-Bliss equation [hcoll = a - b*dT]; defaults from SAM
            units=pyunits.dimensionless,
            mutable=True,
            doc="Product of collector heat removal factor (FR), cover transmittance (t), and shortwave absorptivity of absorber (a)",
        )

        self.FR_UL = Param(
            initialize=3.85,  # Thermal loss coeff "b" in Hottel-Whillier-Bliss equation [hcoll = a - b*dT]; defaults from SAM
            units=pyunits.kilowatt / (pyunits.m**2 * pyunits.K),
            mutable=True,
            doc="Product of collector heat removal factor (FR) and overall heat loss coeff. of collector (UL)",
        )

        self.trans_absorb_prod = Param(
            initialize=1,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Effective transmittance-absorption product",
        )

        self.heat_loss_coeff = Param(
            initialize=1,
            units=pyunits.watt / (pyunits.m**2 * pyunits.K),
            mutable=True,
            doc="Overall collector heat loss coefficient",
        )

        self.mdot_test = Param(
            initialize=1,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate of fluid during characterization test",
        )

        self.cp_test = Param(
            initialize=4184,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            doc="Specific heat capacity of fluid during characterization test",
        )

        self.cp_use = Param(
            initialize=4184,
            units=pyunits.J / (pyunits.kg * pyunits.K),
            mutable=True,
            doc="specific heat capacity of fluid during use",
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

        self.max_irradiance = Param(
            initialize=1000,
            units=pyunits.W / pyunits.m**2,
            mutable=True,
            doc="Maximum irradiance at the location",
        )

        self.factor_delta_T = Param(
            initialize=0.03,  # this is a guess as to what the 30 represents in the equation for total collector area in SAM documentation
            units=pyunits.K,
            mutable=True,
            doc="Influent minus ambient temperature",
        )

        self.temperature_ambient = Param(
            initialize=30 + 273.15, units=pyunits.K, doc="Ambient temperature"
        )

        # ==========VARIABLES==========

        self.collector_area_total = Var(
            initialize=1, units=pyunits.m**2, doc="Total collector area"
        )

        self.collector_area = Var(
            initialize=1, units=pyunits.m**2, doc="Area of a single collector"
        )

        self.total_irradiance = Var(
            initialize=900,
            units=pyunits.W / pyunits.m**2,
            doc="Total irradiance",
        )

        self.Fprime_UL = Var(
            initialize=1,
            units=pyunits.W / (pyunits.m**2 * pyunits.K),
            doc="Product of collector efficiency factor and overall heat loss coefficient at test conditions, D&B Eq. 6.20.4",
        )

        self.ratio_FRta = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Ratio of FRta_use to FRta_test, D&B Eq. 6.20.3",
        )

        self.net_heat_gain = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=pyunits.W,
            doc="Useful net heat gain",
        )

        self.heat_load = Var(
            initialize=1, units=pyunits.kW, doc="Rated plant heat capacity in MW"
        )

        self.heat_annual = Var(
            initialize=1, units=pyunits.MWh, doc="Annual heat generated by flat plate"
        )

        self.heat_loss_coeff = Expression(expr=self.FR_UL / self.FR_ta)

        # ==========CONSTRAINTS==========

        @self.Constraint(doc="Total collector area")
        def eq_collector_area_total(b):
            return b.collector_area_total == b.collector_area * b.number_collectors

        @self.Constraint(
            doc="Corrected collector heat loss coefficient, D&B Eq. 6.20.4, calculated at test conditions"
        )
        def eq_Fprime_UL(b):
            return b.Fprime_UL == (
                -b.mdot_test
                * b.cp_test
                / b.collector_area
                * log(
                    1
                    - b.FR_ta
                    * b.heat_loss_coeff
                    * b.collector_area
                    / (b.mdot_test * b.cp_test)
                )
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Ratio of FRta_use to FRta_test, D&B Eq. 6.20.3",
        )
        def eq_ratio_FRta(b, t):
            return b.ratio_FRta == (
                (
                    b.inlet_block[t].flow_mass_phase_comp["Liq", "H2O"]
                    * b.number_collectors
                    * b.cp_use
                    / (b.collector_area * b.number_collectors)
                    * (
                        1
                        - exp(
                            -b.collector_area
                            * b.number_collectors
                            * b.Fprime_UL
                            / (
                                b.inlet_block[t].flow_mass_phase_comp["Liq", "H2O"]
                                * b.number_collectors
                                * b.cp_use
                            )
                        )
                    )
                )
                / (b.FR_ta * b.heat_loss_coeff)
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Useful net heat gain, not accounting for pipe heat losses or a heat exchanger, D&B Eq. 6.8.1",
        )
        def eq_net_heat_gain(b, t):
            return b.net_heat_gain[t] == (
                b.collector_area
                * b.number_collectors
                * b.ratio_FRta
                * (
                    b.FR_ta
                    * b.trans_absorb_prod
                    * b.total_irradiance
                    * b.trans_absorb_prod
                    - b.FR_ta
                    * b.heat_loss_coeff
                    * (b.inlet_block[t].temperature - b.temperature_ambient)
                )
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraint to calculate the outlet temperature from FPC",
        )
        def eq_outlet_temp(b, t):
            return b.outlet_block[t].temperature == (
                b.inlet_block[t].temperature
                + b.net_heat_gain[t]
                / b.inlet_block[t].flow_mass_phase_comp["Liq", "H2O"]
                / b.cp_use
            )

        @self.Constraint(doc="Pump power")
        def eq_P_pump(b):
            return b.electricity == (b.pump_power / b.pump_eff)

        @self.Constraint(self.flowsheet().config.time, doc="Useful heat generated")
        def eq_heat(b, t):
            return b.heat == b.net_heat_gain[t]

        @self.Constraint(doc="Heat generated annually")
        def eq_heat_annual(b):
            return b.heat_annual == b.heat * pyunits.convert(
                1 * pyunits.year, to_units=pyunits.hour
            )

        # Need to check this
        @self.Constraint(doc="Heat load of system")
        def eq_heat_load(b):
            first_term = pyunits.convert(
                b.FR_ta * b.trans_absorb_prod * b.max_irradiance,
                to_units=pyunits.kW / pyunits.m**2,
            )
            return b.heat_load == pyunits.convert(
                b.collector_area_total * (first_term - b.FR_UL * b.factor_delta_T),
                to_units=pyunits.kW,
            )

    def initialize(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):

        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Create solver
        opt = get_solver(solver=solver, options=optarg)

        self.inlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.inlet_block[
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

        self.outlet_block.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        # solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info("FPC initialization status {}.".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.collector_area) is None:
            sf = iscale.get_scaling_factor(self.collector_area, default=1)
            iscale.set_scaling_factor(self.collector_area, sf)

        if iscale.get_scaling_factor(self.collector_area_total) is None:
            sf = iscale.get_scaling_factor(self.collector_area_total, default=1)
            iscale.set_scaling_factor(self.collector_area_total, sf)

        if iscale.get_scaling_factor(self.total_irradiance) is None:
            sf = iscale.get_scaling_factor(self.total_irradiance, default=1e-2)
            iscale.set_scaling_factor(self.total_irradiance, sf)

        if iscale.get_scaling_factor(self.Fprime_UL) is None:
            sf = iscale.get_scaling_factor(self.Fprime_UL, default=1)
            iscale.set_scaling_factor(self.Fprime_UL, sf)

        if iscale.get_scaling_factor(self.ratio_FRta) is None:
            sf = iscale.get_scaling_factor(self.ratio_FRta, default=1)
            iscale.set_scaling_factor(self.ratio_FRta, sf)

        if iscale.get_scaling_factor(self.net_heat_gain) is None:
            sf = iscale.get_scaling_factor(self.net_heat_gain, default=1e-3)
            iscale.set_scaling_factor(self.net_heat_gain, sf)

        if iscale.get_scaling_factor(self.heat_load) is None:
            sf = iscale.get_scaling_factor(self.heat_load, default=1e-3)
            iscale.set_scaling_factor(self.heat_load, sf)

        if iscale.get_scaling_factor(self.heat_annual) is None:
            sf = iscale.get_scaling_factor(self.heat_annual, default=1e-6)
            iscale.set_scaling_factor(self.heat_annual, sf)

        if iscale.get_scaling_factor(self.heat) is None:
            sf = iscale.get_scaling_factor(self.heat, default=1e-3)
            iscale.set_scaling_factor(self.heat, sf)

        if iscale.get_scaling_factor(self.electricity) is None:
            sf = iscale.get_scaling_factor(self.electricity, default=0.1)
            iscale.set_scaling_factor(self.electricity, sf)

    @property
    def default_costing_method(self):
        return cost_flat_plate
