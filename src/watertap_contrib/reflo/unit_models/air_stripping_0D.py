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

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    Block,
    log,
    log10,
    exp,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from enum import Enum, auto

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import InitializationError, ConfigurationError

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin

__author__ = "Kurban Sitterley"


_log = idaeslog.getLogger(__name__)


class PackingMaterial(Enum):
    none = auto()


@declare_process_block_class("AirStripping0D")
class AirStripping0DData(InitializationMixin, UnitModelBlockData):
    """
    Zero dimensional air stripping model
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

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.none,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.none.
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

    CONFIG.declare(
        "target",
        ConfigValue(
            default="X",
            domain=str,
            description="Designates targeted species for removal",
        ),
    )

    def build(self):
        super().build()

        target = self.config.target
        self.target_set = Set(initialize=[target])
        comps = self.config.property_package.component_list
        solutes = self.config.property_package.solute_set
        phase_set = self.config.property_package.phase_list
        self.phase_target_set = phase_set * self.target_set

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.process_flow = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.process_flow.add_state_blocks(has_phase_equilibrium=True)
        self.process_flow.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.process_flow.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=False
        )
        # self.process_flow.add_isothermal_assumption()
        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        prop_in = self.process_flow.properties_in[0]

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.opt_air_water_ratio_param = Param(
            initialize=3.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor multiplied by minimum air-to-water ratio for design",
        )
        self.surface_area_packing_total = Var(
            units=pyunits.m**-1,
            bounds=(0, None),
            doc="Total specific surface area of packing.",
        )
        self.surface_area_packing_wetted = Var(
            units=pyunits.m**-1,
            bounds=(0, None),
            doc="Wetted specific surface area of packing.",
        )

        self.diam_nominal_packing = Var(units=pyunits.m, bounds=(0, None))

        self.surf_tension_packing = Var(
            units=pyunits.kg * pyunits.s**-2,
            bounds=(0, None),
            doc="Surface tension of packing",
        )
        self.surf_tension_water = Var(
            units=pyunits.kg * pyunits.s**-2,
            bounds=(0, None),
            doc="Surface tension of water",
        )

        self.stripping_factor = Var(
            initialize=2, bounds=(1, 20), units=pyunits.dimensionless
        )

        self.min_air_water_ratio = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Minumum air-to-water ratio",
        )

        self.opt_air_water_ratio = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Optimum air-to-water ratio",
        )

        self.height_packed_tower = Var(
            initialize=1,
            bounds=(0, 10),
            units=pyunits.m,
            doc="Height of packed tower",
        )

        self.diam_packed_tower = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Diameter of packed tower",
        )

        self.mass_transfer_coeff = Var(
            self.phase_target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Mass transfer coefficient in tower",
        )

        self.mass_loading_rate = Var(
            self.phase_target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            doc="Mass loading rate in tower",
        )

        self.height_transfer_unit = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Height of one transfer unit",
        )

        self.number_transfer_unit = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.N_Re = Var(
            phase_set,
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Reynolds number",
        )

        self.N_Fr = Var(
            initialize=0.1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Froude number",
        )

        self.N_We = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Weber number",
        )

        self.N_Sc = Var(
            phase_set,
            initialize=700,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.N_We = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Weber number",
        )

        self.wettability_parameter = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Wettability parameter",
        )

        self.packing_efficiency_number = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Packing efficiency number",
        )

        self.build_onda()

    def build_onda(self):

        self.onda = Block()

        self.onda_F = F = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop F term",
        )

        # a0 = -6.6599 + 4.3077*F - 1.3503*F^2 + 0.15931*F^3
        self.onda_a0_param1 = a01 = Param(
            initialize=-6.6599,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a0 term, first parameter",
        )
        self.onda_a0_param2 = a02 = Param(
            initialize=4.3077,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a0 term, second parameter",
        )
        self.onda_a0_param3 = a03 = Param(
            initialize=-1.3503,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a0 term, third parameter",
        )
        self.onda_a0_param4 = a04 = Param(
            initialize=0.15931,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a0 term, fourth parameter",
        )

        self.onda_a0 = a0 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a0 term",
        )

        @self.Constraint(doc="Onda a0 equation")
        def eq_onda_a0(b):
            return a0 == a01 + a02*F + a03*F**2 + a04*F**3


        # a1 = 3.0945 - 4.3512*F + 1.6240*F^2 - 0.20855*F^3
        self.onda_a1_param1 = a11 = Param(
            initialize=3.0945,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a1 term, first parameter",
        )
        self.onda_a1_param2 = a12 = Param(
            initialize=-4.3512,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a1 term, second parameter",
        )
        self.onda_a1_param3 = a13 = Param(
            initialize=1.6240,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a1 term, third parameter",
        )
        self.onda_a1_param4 = a14 = Param(
            initialize=-0.20855,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a1 term, fourth parameter",
        )
        self.onda_a1 = a1 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a1 term",
        )

        @self.Constraint(doc="Onda a1 equation")
        def eq_onda_a1(b):
            return a1 == a11 + a12*F + a13*F**2 + a14*F**3

        # a2 = 1.7611 - 2.3394*F + 0.89914*F^2 - 0.11597*F^3
        self.onda_a2_param1 = a21 = Param(
            initialize=1.7611,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a2 term, first parameter",
        )
        self.onda_a2_param2 = a22 = Param(
            initialize=-2.3394,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a2 term, second parameter",
        )
        self.onda_a2_param3 = a23 = Param(
            initialize=0.89914,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a2 term, third parameter",
        )
        self.onda_a2_param4 = a24 = Param(
            initialize=-0.115971,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a2 term, fourth parameter",
        )
        self.onda_a2 = a2 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Onda correlation: Pressure drop a2 term",
        )

        @self.Constraint(doc="Onda a2 equation")
        def eq_onda_a2(b):
            return a2 == a21 + a22*F + a23*F**2 + a24*F**3
        
    def initialize_build(
        self,
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
        pass
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # flags = self.process_flow.properties_in.initialize(
        #     outlvl=outlvl,
        #     optarg=optarg,
        #     solver=solver,
        #     state_args=state_args,
        #     hold_state=True,
        # )
        init_log.info("Initialization Step 1a Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.process_flow.properties_in[
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

        # for p, j in self.process_flow.properties_out.phase_component_set:
        #     if j in self.target_ion_set:
        #         state_args_out["flow_mol_phase_comp"][(p, j)] = (
        #             state_args["flow_mol_phase_comp"][(p, j)] * 1e-3
        #         )

        # self.process_flow.properties_out.initialize(
        #     outlvl=outlvl,
        #     optarg=optarg,
        #     solver=solver,
        #     state_args=state_args_out,
        # )
        init_log.info("Initialization Step 1b Complete.")

        state_args_regen = deepcopy(state_args)

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        pass

        # if iscale.get_scaling_factor(self.t_breakthru) is None:
        #     iscale.set_scaling_factor(self.t_breakthru, 1e-6)

    def _get_stream_table_contents(self, time_point=0):
        pass

    #     return create_stream_table_dataframe(
    #         {},
    #         time_point=time_point,
    #     )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_ion_exchange
