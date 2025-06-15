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

from pyomo.environ import (
    check_optimal_termination,
    Var,
    Constraint,
    Expression,
    Suffix,
    RangeSet,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
    FlowsheetBlock,
)
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

from watertap.core import InitializationMixin, ControlVolume0DBlock
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect
from watertap_contrib.reflo.costing.units.multi_effect_crystallizer import (
    cost_multi_effect_crystallizer,
)

_log = idaeslog.getLogger(__name__)

__author__ = "Oluwamayowa Amusat, Zhuoran Zhang, Kurban Sitterley"


@declare_process_block_class("MultiEffectCrystallizer")
class MultiEffectCrystallizerData(InitializationMixin, UnitModelBlockData):
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
        "number_effects",
        ConfigValue(
            default=4,
            domain=PositiveInt,
            description="Number of effects of the multi-effect crystallizer system",
            doc="""Number of effects of the multi-effect crystallizer system.""",
        ),
    )

    def build(self):
        super().build()

        if self.config.number_effects <= 1:
            raise ConfigurationError(
                "The MultiEffectCrystallizer model requires more than 1 effect."
                "To model a crystallizer with one effect, use the CrystallizerEffect model with 'standalone=True'."
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)
        self.control_volume.add_material_balances(
            balance_type=MaterialBalanceType.componentPhase, has_mass_transfer=True
        )
        self.add_inlet_port(name="inlet", block=self.control_volume)
        self.add_outlet_port(name="outlet", block=self.control_volume)

        self.number_effects = self.config.number_effects
        self.Effects = RangeSet(self.config.number_effects)

        self.first_effect = self.Effects.first()
        self.last_effect = self.Effects.last()

        self.effects = FlowsheetBlock(self.Effects, dynamic=False)

        # There is no solid NaCl coming in
        self.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        # All solid NaCl accessed through properties_solids
        self.control_volume.mass_transfer_term[0, "Sol", "NaCl"].fix(0)

        # Local expressions to aggregate material flows for all effects
        total_mass_flow_water_in_expr = 0
        total_mass_flow_salt_in_expr = 0
        total_mass_flow_pure_water_out_expr = 0
        total_mass_flow_water_out_expr = 0
        total_mass_flow_salt_out_expr = 0
        total_flow_vol_in_expr = 0
        total_flow_vol_out_expr = 0

        for n, eff in self.effects.items():
            eff.effect = effect = CrystallizerEffect(
                property_package=self.config.property_package,
                property_package_vapor=self.config.property_package_vapor,
                standalone=False,
            )
            effect.properties_in[0].conc_mass_phase_comp

            total_flow_vol_in_expr += effect.properties_in[0].flow_vol_phase["Liq"]
            total_flow_vol_out_expr += effect.properties_pure_water[0].flow_vol_phase[
                "Liq"
            ]

            total_mass_flow_water_in_expr += effect.properties_in[
                0
            ].flow_mass_phase_comp["Liq", "H2O"]

            total_mass_flow_salt_in_expr += effect.properties_in[
                0
            ].flow_mass_phase_comp["Liq", "NaCl"]

            total_mass_flow_pure_water_out_expr += effect.properties_pure_water[
                0
            ].flow_mass_phase_comp["Liq", "H2O"]

            total_mass_flow_water_out_expr += effect.properties_out[
                0
            ].flow_mass_phase_comp["Liq", "H2O"]

            total_mass_flow_salt_out_expr += (
                effect.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"]
                + effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"]
            )

            if n == self.first_effect:
                tmp_dict = dict(**self.config.property_package_args)
                tmp_dict["has_phase_equilibrium"] = False
                tmp_dict["parameters"] = self.config.property_package_vapor
                tmp_dict["defined_state"] = False

                effect.heating_steam = (
                    self.config.property_package_vapor.state_block_class(
                        self.flowsheet().config.time,
                        doc="Material properties of inlet heating steam",
                        **tmp_dict,
                    )
                )

                @effect.Constraint(
                    doc="Change in temperature at inlet for first effect"
                )
                def eq_delta_temperature_inlet_effect_1(b):
                    return (
                        b.delta_temperature_in[0]
                        == b.heating_steam[0].temperature - b.temperature_operating
                    )

                @effect.Constraint(
                    doc="Change in temperature at outlet for first effect"
                )
                def eq_delta_temperature_outlet_effect_1(b):
                    return (
                        b.delta_temperature_out[0]
                        == b.heating_steam[0].temperature
                        - b.properties_in[0].temperature
                    )

                @effect.Constraint(doc="Heat transfer equation for first effect")
                def eq_heat_transfer_effect_1(b):
                    return b.work_mechanical[0] == pyunits.convert(
                        b.overall_heat_transfer_coefficient
                        * b.heat_exchanger_area
                        * b.delta_temperature[0],
                        to_units=pyunits.kJ * pyunits.s**-1,
                    )

                @effect.Constraint(doc="Calculate mass flow rate of heating steam")
                def eq_heating_steam_flow_rate(b):
                    return b.work_mechanical[0] == (
                        pyunits.convert(
                            b.heating_steam[0].dh_vap_mass
                            * b.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"],
                            to_units=pyunits.kJ * pyunits.s**-1,
                        )
                    )

                self.add_port(name="solids", block=effect.properties_solids)
                self.add_port(name="vapor", block=effect.properties_vapor)
                self.add_port(name="steam", block=effect.heating_steam)
                self.steam.temperature.setub(1000)

            else:
                prev_effect = self.effects[n - 1].effect

                del_temp_in_constr = Constraint(
                    expr=effect.delta_temperature_in[0]
                    == prev_effect.properties_vapor[0].temperature
                    - effect.temperature_operating,
                    doc=f"Change in temperature at inlet for effect {n}",
                )
                self.add_component(
                    f"eq_delta_temperature_inlet_effect_{n}", del_temp_in_constr
                )

                del_temp_out_constr = Constraint(
                    expr=effect.delta_temperature_out[0]
                    == prev_effect.properties_pure_water[0].temperature
                    - effect.properties_in[0].temperature,
                    doc=f"Change in temperature at outlet for effect {n}",
                )
                self.add_component(
                    f"eq_delta_temperature_outlet_effect_{n}", del_temp_out_constr
                )

                hx_constr = Constraint(
                    expr=prev_effect.energy_flow_superheated_vapor
                    == effect.overall_heat_transfer_coefficient
                    * effect.heat_exchanger_area
                    * effect.delta_temperature[0],
                    doc=f"Heat transfer equation for effect {n}",
                )
                self.add_component(f"eq_heat_transfer_effect_{n}", hx_constr)

                energy_flow_constr = Constraint(
                    expr=effect.work_mechanical[0]
                    == pyunits.convert(
                        prev_effect.energy_flow_superheated_vapor,
                        to_units=pyunits.kJ * pyunits.s**-1,
                    ),
                    doc=f"Energy supplied to effect {n}",
                )
                self.add_component(
                    f"eq_energy_for_effect_{n}_from_effect_{n - 1}", energy_flow_constr
                )

        @self.Constraint(doc="Mass transfer term for liquid water")
        def eq_mass_transfer_term_liq_water(b):
            return b.control_volume.mass_transfer_term[0, "Liq", "H2O"] == -1 * (
                total_mass_flow_water_out_expr
            )

        @self.Constraint(doc="Mass transfer term for vapor water")
        def eq_mass_transfer_term_vap_water(b):
            return b.control_volume.mass_transfer_term[0, "Vap", "H2O"] == -1 * (
                b.effects[1].effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(doc="Mass transfer term for salt in liquid phase")
        def eq_mass_transfer_term_liq_salt(b):
            return b.control_volume.mass_transfer_term[0, "Liq", "NaCl"] == -1 * (
                total_mass_flow_salt_out_expr
            )

        @self.Constraint(doc="Steam flow")
        def eq_overall_steam_flow(b):
            return (
                b.control_volume.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]
                == b.effects[1]
                .effect.heating_steam[0]
                .flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(doc="Mass balance of water for all effects")
        def eq_overall_mass_balance_water_in(b):
            return (
                b.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                == total_mass_flow_water_in_expr
            )

        @self.Constraint(doc="Mass balance of salt for all effects")
        def eq_overall_mass_balance_salt_in(b):
            return (
                b.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
                == total_mass_flow_salt_in_expr
            )

        @self.Constraint(doc="Control volume temperature at outlet")
        def eq_temperature_outlet(b):
            return (
                b.control_volume.properties_out[0].temperature
                == b.effects[b.last_effect].effect.properties_pure_water[0].temperature
            )

        @self.Constraint(doc="Control volume pressure at outlet")
        def eq_isobaric(b):
            return (
                b.control_volume.properties_out[0].pressure
                == b.control_volume.properties_in[0].pressure
            )

        self.total_flow_vol_in = Expression(expr=total_flow_vol_in_expr)

        self.recovery_vol_phase = Var(
            ["Liq"],
            initialize=0.75,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Overall unit recovery on volumetric basis",
        )

        @self.Constraint(doc="Volumetric recovery")
        def eq_recovery_vol_phase(b):
            return (
                b.recovery_vol_phase["Liq"]
                == total_flow_vol_out_expr / total_flow_vol_in_expr
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
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

        init_args = dict(
            state_args=state_args,
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        flow_mass_water_initial = (
            self.control_volume.properties_in[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value
            / self.number_effects
        )
        flow_mass_salt_initial = (
            self.control_volume.properties_in[0]
            .flow_mass_phase_comp["Liq", "NaCl"]
            .value
            / self.number_effects
        )

        opt = get_solver(solver, optarg)

        for n, eff in self.effects.items():
            # Each effect is first solved in a vacuum with linking constraints deactivated
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
                flow_mass_water_initial
            )
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
                flow_mass_salt_initial
            )
            if n == 1:
                if not degrees_of_freedom(eff.effect) == 0:
                    raise InitializationError(
                        f"Degrees of freedom in first effect must be zero during initialization, "
                        f"but has {degrees_of_freedom(eff.effect)}. Check inlet conditions and re-initialize."
                    )
                eff.effect.initialize(**init_args)
                inlet_conc = (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .value
                )
                mass_transfer_coeff = eff.effect.overall_heat_transfer_coefficient.value
            else:
                # Deactivate contraint that links energy flow between effects
                linking_constr = getattr(
                    self, f"eq_energy_for_effect_{n}_from_effect_{n - 1}"
                )
                linking_constr.deactivate()
                eff.effect.initialize(**init_args)
                linking_constr.activate()
                # Fix inlet feed concentration to be equal across all effects
                eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                    inlet_conc
                )
                eff.effect.overall_heat_transfer_coefficient.fix(mass_transfer_coeff)

            # Unfix inlet mass flow rates to allow unit model to determine based on energy flows
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

            init_log.info(f"Initialization of Effect {n} Complete.")

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
            iscale.set_scaling_factor(self.recovery_vol_phase, 10)

            iscale.constraint_scaling_transform(self.eq_isobaric, 1e-6)

        if 1 in self.Effects:
            iscale.constraint_scaling_transform(
                self.effects[1].effect.eq_heat_transfer_effect_1, 1e-4
            )
            iscale.constraint_scaling_transform(
                self.effects[1].effect.eq_heating_steam_flow_rate, 1e-4
            )

        if 2 in self.Effects:
            iscale.constraint_scaling_transform(self.eq_heat_transfer_effect_2, 1e-6)
            iscale.constraint_scaling_transform(
                self.eq_energy_for_effect_2_from_effect_1, 1e-4
            )
        if 3 in self.Effects:
            iscale.constraint_scaling_transform(self.eq_heat_transfer_effect_3, 1e-6)
            iscale.constraint_scaling_transform(
                self.eq_energy_for_effect_3_from_effect_2, 1e-4
            )
        if 4 in self.Effects:
            iscale.constraint_scaling_transform(self.eq_heat_transfer_effect_4, 1e-6)
            iscale.constraint_scaling_transform(
                self.eq_energy_for_effect_4_from_effect_3, 1e-4
            )

    @property
    def default_costing_method(self):
        return cost_multi_effect_crystallizer