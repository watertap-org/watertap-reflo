#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
    Var,
    Param,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.core.util.initialization import interval_initializer
from watertap.unit_models.crystallizer import CrystallizationData
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)

from watertap_contrib.reflo.costing.units.multi_effect_crystallizer import (
    cost_multi_effect_crystallizer,
)

__author__ = "Oluwamayowa Amusat, Zhuoran Zhang, Kurban Sitterley"


@declare_process_block_class("CrystallizerEffect")
class CrystallizerEffectData(CrystallizationData):
    """
    Zero dimensional model for crystallizer effect
    """

    CONFIG = CrystallizationData.CONFIG()

    CONFIG.declare(
        "property_package_vapor",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for heating steam properties",
            doc="""Property parameter object used to define steam property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "standalone",
        ConfigValue(
            default=True,
            domain=bool,
            description="Flag to indicate if model is used alone or as part of MultiEffectCrystallizer.",
            doc="""Flag to indicate if model is used alone or as part of MultiEffectCrystallizer unit model,
    **default** - True.""",
        ),
    )

    def build(self):
        super().build()

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.properties_pure_water = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of pure water vapour at outlet",
            **tmp_dict,
        )

        self.add_port(name="pure_water", block=self.properties_pure_water)

        self.properties_pure_water[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        self.properties_pure_water[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0)
        self.properties_pure_water[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        self.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        self.properties_in[0].flow_vol_phase["Liq"]

        self.inlet.temperature.setub(1000)
        self.outlet.temperature.setub(1000)
        self.solids.temperature.setub(1000)
        self.vapor.temperature.setub(1000)
        self.pure_water.temperature.setub(1000)

        self.steam_pressure = Param(
            initialize=3,
            mutable=True,
            units=pyunits.bar,
            doc="Steam pressure (gauge) for crystallizer heating: 3 bar default based on Dutta example",
        )

        self.efficiency_pump = Param(
            initialize=0.7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Crystallizer pump efficiency - assumed",
        )

        self.pump_head_height = Param(
            initialize=1,
            mutable=True,
            units=pyunits.m,
            doc="Crystallizer pump head height -  assumed, unvalidated",
        )

        self.energy_flow_superheated_vapor = Var(
            initialize=1e5,
            bounds=(-5e8, 5e8),
            units=pyunits.watt,
            doc="Energy that could be supplied from vapor",
        )

        self.delta_temperature_in = Var(
            self.flowsheet().time,
            initialize=35,
            bounds=(None, None),
            units=pyunits.degK,
            doc="Temperature difference at the inlet side",
        )

        self.delta_temperature_out = Var(
            self.flowsheet().time,
            initialize=35,
            bounds=(None, None),
            units=pyunits.degK,
            doc="Temperature difference at the outlet side",
        )

        delta_temperature_chen_callback(self)

        self.heat_exchanger_area = Var(
            initialize=1000.0,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Heat exchanger area",
        )

        self.overall_heat_transfer_coefficient = Var(
            initialize=1e3,
            bounds=(0, None),
            units=pyunits.watt / pyunits.m**2 / pyunits.degK,
            doc="Overall heat transfer coefficient for heat exchangers",
        )

        @self.Constraint()
        def eq_pure_water_mass_flow_rate(b):
            return (
                b.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(doc="Thermal energy in the vapor")
        def eq_vapor_energy_constraint(b):
            return b.energy_flow_superheated_vapor == (
                b.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
                * (
                    b.properties_pure_water[
                        0
                    ].dh_vap_mass_solvent  # Latent heat from the vapor
                    + b.properties_vapor[0].enth_mass_solvent["Vap"]
                    - b.properties_pure_water[0].enth_mass_solvent["Vap"]
                )
            )

        self.del_component(self.eq_p_con1)
        self.del_component(self.eq_p_con2)

        @self.Constraint()
        def eq_p_con1(b):
            return b.pressure_operating == b.properties_out[0].pressure

        @self.Constraint()
        def eq_p_con2(b):
            return b.pressure_operating == b.properties_solids[0].pressure

        @self.Constraint()
        def eq_p_con4(b):
            return b.properties_pure_water[0].pressure == self.pressure_operating

        @self.Constraint()
        def eq_p_con5(b):
            return b.properties_pure_water[0].pressure_sat == self.pressure_operating

        if self.config.standalone:
            tmp_dict["parameters"] = self.config.property_package_vapor
            tmp_dict["defined_state"] = False

            self.heating_steam = self.config.property_package_vapor.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of inlet heating steam",
                **tmp_dict,
            )

            self.add_port(name="steam", block=self.heating_steam)
            self.steam.temperature.setub(1000)

            @self.Constraint(doc="Change in temperature at inlet")
            def eq_delta_temperature_inlet(b):
                return (
                    b.delta_temperature_in[0]
                    == b.heating_steam[0].temperature - b.temperature_operating
                )

            @self.Constraint(doc="Change in temperature at outlet")
            def eq_delta_temperature_outlet(b):
                return (
                    b.delta_temperature_out[0]
                    == b.heating_steam[0].temperature - b.properties_in[0].temperature
                )

            @self.Constraint(doc="Heat transfer equation")
            def eq_heat_transfer(b):
                return b.work_mechanical[0] == pyunits.convert(
                    b.overall_heat_transfer_coefficient
                    * b.heat_exchanger_area
                    * b.delta_temperature[0],
                    to_units=pyunits.kJ * pyunits.s**-1,
                )

            @self.Constraint(doc="Calculate mass flow rate of heating steam")
            def eq_heating_steam_flow_rate(b):
                return b.work_mechanical[0] == (
                    pyunits.convert(
                        b.heating_steam[0].dh_vap_mass
                        * b.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"],
                        to_units=pyunits.kJ * pyunits.s**-1,
                    )
                )

        else:
            self.del_component(self.inlet)
            self.del_component(self.outlet)
            self.del_component(self.solids)
            self.del_component(self.vapor)
            self.del_component(self.pure_water)

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure changer initialization routines

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
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")

        if state_args is None:
            state_args = {}
            state_dict = self.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        state_args_solids = deepcopy(state_args)

        for p, j in self.properties_solids.phase_component_set:
            if p == "Sol":
                state_args_solids["flow_mass_phase_comp"][p, j] = state_args[
                    "flow_mass_phase_comp"
                ]["Liq", j]
            elif p in ["Liq", "Vap"]:
                state_args_solids["flow_mass_phase_comp"][p, j] = 1e-8

        self.properties_solids.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_solids,
        )

        state_args_vapor = deepcopy(state_args)

        for p, j in self.properties_vapor.phase_component_set:
            if p == "Vap":
                state_args_vapor["flow_mass_phase_comp"][p, j] = state_args[
                    "flow_mass_phase_comp"
                ]["Liq", j]
            elif p == "Liq" or p == "Sol":
                state_args_vapor["flow_mass_phase_comp"][p, j] = 1e-8

        self.properties_vapor.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_vapor,
        )

        self.properties_pure_water.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_vapor,
        )

        if hasattr(self, "heating_steam"):
            state_args_steam = deepcopy(state_args)

            for p, j in self.properties_vapor.phase_component_set:
                state_args_steam["flow_mass_phase_comp"][p, j] = 1

            self.heating_steam.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_steam,
            )

        init_log.info_high("Initialization Step 2 Complete.")

        interval_initializer(self)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info_high(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        if iscale.get_scaling_factor(self.work_mechanical[0]) is None:
            iscale.set_scaling_factor(
                self.work_mechanical[0],
                iscale.get_scaling_factor(
                    self.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]
                )
                * iscale.get_scaling_factor(
                    self.properties_in[0].enth_mass_solvent["Vap"]
                )
                * 1e3,
            )
        if iscale.get_scaling_factor(self.energy_flow_superheated_vapor) is None:
            iscale.set_scaling_factor(
                self.energy_flow_superheated_vapor,
                iscale.get_scaling_factor(
                    self.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
                )
                * iscale.get_scaling_factor(
                    self.properties_vapor[0].enth_mass_solvent["Vap"]
                ),
            )

        if iscale.get_scaling_factor(self.delta_temperature_in[0]) is None:
            iscale.set_scaling_factor(self.delta_temperature_in[0], 0.1)

        if iscale.get_scaling_factor(self.delta_temperature_out[0]) is None:
            iscale.set_scaling_factor(self.delta_temperature_out[0], 0.1)

        if iscale.get_scaling_factor(self.heat_exchanger_area) is None:
            iscale.set_scaling_factor(self.heat_exchanger_area, 0.1)

        if iscale.get_scaling_factor(self.overall_heat_transfer_coefficient) is None:
            iscale.set_scaling_factor(self.overall_heat_transfer_coefficient, 0.01)

        iscale.set_scaling_factor(self.properties_solids[0].flow_vol_phase["Vap"], 1e5)
        iscale.set_scaling_factor(self.properties_out[0].flow_vol_phase["Vap"], 1e5)

        for _, c in self.eq_p_con4.items():
            sf = iscale.get_scaling_factor(self.properties_pure_water[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        iscale.constraint_scaling_transform(self.eq_vapor_energy_constraint, 1e-6)
        iscale.constraint_scaling_transform(self.eq_operating_pressure_constraint, 1e-5)
        iscale.constraint_scaling_transform(self.eq_p_con5, 1e-5)

    @property
    def default_costing_method(self):
        return cost_multi_effect_crystallizer