#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
    ConcreteModel,
    Var,
    check_optimal_termination,
    Param,
    Constraint,
    Suffix,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowsheetBlock,
)
from watertap.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block

from idaes.core.util.exceptions import InitializationError

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.util.initialization import interval_initializer
from watertap.unit_models.crystallizer import Crystallization, CrystallizationData
from watertap_contrib.reflo.costing.units.crystallizer_watertap import (
    cost_crystallizer_watertap,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)

_log = idaeslog.getLogger(__name__)

__author__ = "Oluwamayowa Amusat, Zhuoran Zhang, Kurban Sitterley"


@declare_process_block_class("CrystallizerEffect")
class CrystallizerEffectData(CrystallizationData):
    """
    Zero-order model for crystallizer effect
    """

    CONFIG = CrystallizationData.CONFIG()

    CONFIG.declare(
        "property_package_vapor",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for heating and motive steam properties",
            doc="""Property parameter object used to define steasm property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )

    def build(self):
        super().build()

        # self.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        # self.properties_in[0].flow_vol_phase["Liq"]

        # self.tsat_constants_dict = {
        # "tsat_c1": (42.6776, pyunits.K),
        # "tsat_c2": (-3892.7, pyunits.K),
        # "tsat_c3": (1000, pyunits.kPa),
        # "tsat_c4": (-9.48654, pyunits.dimensionless),
        # }

        # for c, (v, u) in self.tsat_constants_dict.items():
        #     self.properties_in[0].add_component(c, Param(initialize=v, units=u))

        # self.properties_in[0].add_component("steam_pressure_sat", Param(initialize=150, units=pyunits.kPa, doc="Steam pressure saturated"))

        # self.properties_in[0].add_component("test", Var())

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

        self.properties_pure_water[0].flow_mass_phase_comp["Liq", "H2O"].fix(1e-8)
        self.properties_pure_water[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0)
        self.properties_pure_water[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        self.properties_pure_water[0].mass_frac_phase_comp["Liq", "NaCl"]

        self.heating_steam = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of inlet heating steam",
            **tmp_dict,
        )

        self.add_port(name="steam", block=self.heating_steam)

        # self.heating_steam[0].flow_mass_phase_comp["Liq", "H2O"].fix(1e-8)
        # self.heating_steam[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0)
        self.heating_steam[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        # self.heating_steam[0].mass_frac_phase_comp["Liq", "NaCl"]

        self.inlet.temperature.setub(1000)
        self.outlet.temperature.setub(1000)
        self.solids.temperature.setub(1000)
        self.vapor.temperature.setub(1000)
        self.pure_water.temperature.setub(1000)
        self.steam.temperature.setub(1000)

        self.energy_flow_superheated_vapor = Var(
            initialize=1e5,
            bounds=(-5e6, 5e6),
            units=pyunits.kJ / pyunits.s,
            doc="Energy could be supplied from vapor",
        )

        self.delta_temperature_in = Var(
            self.flowsheet().time,
            initialize=35,
            bounds=(None, None),
            units=pyunits.K,
            doc="Temperature difference at the inlet side",
        )
        self.delta_temperature_out = Var(
            self.flowsheet().time,
            initialize=35,
            bounds=(None, None),
            units=pyunits.K,
            doc="Temperature difference at the outlet side",
        )
        delta_temperature_chen_callback(self)

        self.area = Var(
            initialize=1000.0,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Heat exchange area",
        )

        self.overall_heat_transfer_coefficient = Var(
            # self.flowsheet().time,
            initialize=100.0,
            bounds=(0, None),
            units=pyunits.W / pyunits.m**2 / pyunits.K,
            doc="Overall heat transfer coefficient",
        )

        # self.overall_heat_transfer_coefficient.fix(100)
        # @self.Constraint()
        # def eq_mass_frac_nacl_steam(b):
        #     return (
        #         b.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        #         == b.heating_steam[0].conc_mass_phase_comp["Liq", "NaCl"]
        #     )
        # @self.Constraint(doc="Flow rate of liquid motive steam is zero")
        # def eq_motive_steam_liquid_mass(b):
        #     return b.heating_steam[0].flow_mass_phase_comp["Liq", "H2O"] == 0
        @self.Constraint()
        def eq_pure_vapor_flow_rate(b):
            return (
                b.properties_pure_water[0].flow_mass_phase_comp["Vap", "H2O"]
                == b.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]
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

        @self.Constraint(doc="Change in temperature at inlet")
        def eq_delta_temperature_in(b):
            return (
                b.delta_temperature_in[0]
                == b.heating_steam[0].temperature_sat_solvent - b.temperature_operating
            )

        @self.Constraint(doc="Change in temperature at outlet")
        def eq_delta_temperature_out(b):
            return (
                b.delta_temperature_out[0]
                == b.heating_steam[0].temperature_sat_solvent
                - b.properties_in[0].temperature
            )

        # self.del_component(self.eq_p_con1)
        # self.del_component(self.eq_p_con2)
        # @self.Constraint()
        # def eq_p_con1(b):
        #     return b.pressure_operating == b.properties_out[0].pressure

        # @self.Constraint()
        # def eq_p_con2(b):
        #     return b.pressure_operating == b.properties_solids[0].pressure

        # @self.Constraint()
        # def eq_temp_222(b):
        #     return b.properties_in[0].temperature == b.temperature_operating
        @self.Constraint()
        def eq_p_con4(b):
            return b.properties_pure_water[0].pressure == self.pressure_operating

        @self.Constraint()
        def eq_p_con5(b):
            return b.properties_pure_water[0].pressure_sat == self.pressure_operating

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

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
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
            elif p == "Liq" or p == "Vap":
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

        state_args_steam = deepcopy(state_args)

        for p, j in self.properties_vapor.phase_component_set:
            state_args_steam["flow_mass_phase_comp"][p, j] = 1e-8

        self.heating_steam.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_steam,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        interval_initializer(self)
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

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

        if iscale.get_scaling_factor(self.area) is None:
            iscale.set_scaling_factor(self.area, 0.1)

        if iscale.get_scaling_factor(self.overall_heat_transfer_coefficient) is None:
            iscale.set_scaling_factor(self.overall_heat_transfer_coefficient, 0.1)

        for ind, c in self.eq_p_con4.items():
            sf = iscale.get_scaling_factor(self.properties_pure_water[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

    def _calculate_saturation_temperature(self):
        pass


if __name__ == "__main__":
    import watertap.property_models.unit_specific.cryst_prop_pack as props
    from watertap.property_models.water_prop_pack import WaterParameterBlock

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = props.NaClParameterBlock()
    m.fs.vapor = WaterParameterBlock()

    m.fs.eff = eff = CrystallizerEffect(property_package=m.fs.props)
    # m.fs.eff.display()

    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.15
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    crystallizer_yield = 0.5
    operating_pressure_eff1 = 0.78

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    eff.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    eff.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    eff.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    eff.inlet.pressure[0].fix(feed_pressure)
    eff.inlet.temperature[0].fix(feed_temperature)
    eff.inlet.temperature[0].fix(404.5)
    eff.properties_in[0].pressure_sat

    # eff.steam.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
    #     eps
    # )
    # eff.steam.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
    #     eps
    # )
    # eff.heating_steam[0].temperature.fix(404.5)
    # eff.steam.temperature[0].fix(404.5)
    # eff.heating_steam[0].flow_mass_phase_comp.fix(eps)
    # eff.heating_steam[0].mass_frac_phase_comp.fix(eps)
    # eff.heating_steam[0].pressure.fix()

    eff.heating_steam.calculate_state(
        var_args={
            ("pressure", None): 101325,
            # ("pressure_sat", None): 28100
            ("temperature", None): 393,
            ("mass_frac_phase_comp", ("Liq", "H2O")): 1 - eps,
            ("mass_frac_phase_comp", ("Liq", "NaCl")): eps,
            ("flow_vol_phase", ("Liq")): 10,
        },
        hold_state=True,
    )

    ###################### 
    eff.crystallization_yield["NaCl"].fix(crystallizer_yield)
    eff.crystal_growth_rate.fix()
    eff.souders_brown_constant.fix()
    eff.crystal_median_length.fix()
    eff.overall_heat_transfer_coefficient.fix(100)

    eff.pressure_operating.fix(operating_pressure_eff1 * pyunits.bar)

    print(f"dof = {degrees_of_freedom(m.fs.eff)}")

    eff.initialize()

    print("\nPROPERTIES IN\n")
    eff.properties_in[0].flow_mass_phase_comp.display()
    eff.properties_in[0].mass_frac_phase_comp.display()
    eff.properties_in[0].temperature.display()
    eff.properties_in[0].pressure_sat.display()
    eff.properties_in[0].temperature_sat_solvent.display()

    # print("\nPROPERTIES VAPOR\n")
    # eff.properties_vapor[0].temperature.display()
    # eff.properties_vapor[0].pressure_sat.display()
    # eff.properties_vapor[0].temperature_sat_solvent.display()

    print("\nPROPERTIES HEATING STEAM\n")
    eff.heating_steam[0].flow_mass_phase_comp.display()
    eff.heating_steam[0].mass_frac_phase_comp.display()
    eff.heating_steam[0].temperature.display()
    eff.heating_steam[0].pressure_sat.display()
    eff.heating_steam[0].temperature_sat_solvent.display()

    eff.temperature_operating.display()
    # eff.heating_steam.display()
