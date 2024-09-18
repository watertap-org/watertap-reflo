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


# Import Pyomo libraries
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    assert_optimal_termination,
    Constraint,
    Expression,
    Suffix,
    RangeSet,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
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
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
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

    CONFIG.declare(
        "number_effects",
        ConfigValue(
            default=4,
            domain=PositiveInt,
            description="Number of effects of the multi-effect crystallizer system",
            doc="""A ConfigBlock specifying the number of effects, which can only be 4.""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.number_effects = self.config.number_effects
        self.Effects = RangeSet(self.config.number_effects)

        self.first_effect = self.Effects.first()
        self.last_effect = self.Effects.last()

        self.effects = FlowsheetBlock(self.Effects, dynamic=False)

        total_flow_vol_in_expr = 0

        for n, eff in self.effects.items():
            eff.effect = effect = CrystallizerEffect(
                property_package=self.config.property_package,
                property_package_vapor=self.config.property_package_vapor,
                standalone=False,
            )
            effect.properties_in[0].conc_mass_phase_comp
            total_flow_vol_in_expr += effect.properties_in[0].flow_vol_phase["Liq"]

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

                # self.add_port(name="steam", block=effect.heating_steam)

                @effect.Constraint(
                    doc="Change in temperature at inlet for first effect."
                )
                def eq_delta_temperature_inlet_effect_1(b):
                    return (
                        b.delta_temperature_in[0]
                        == b.heating_steam[0].temperature - b.temperature_operating
                    )

                @effect.Constraint(
                    doc="Change in temperature at outlet for first effect."
                )
                def eq_delta_temperature_outlet_effect_1(b):
                    return (
                        b.delta_temperature_out[0]
                        == b.heating_steam[0].temperature
                        - b.properties_in[0].temperature
                    )

                @effect.Constraint(doc="Heat transfer equation for first effect.")
                def eq_heat_transfer_effect_1(b):
                    return b.work_mechanical[0] == (
                        b.overall_heat_transfer_coefficient
                        * b.heat_exchanger_area
                        * b.delta_temperature[0]
                    )

                @effect.Constraint(doc="Calculate mass flow rate of heating steam")
                def eq_heating_steam_flow_rate(b):
                    return b.work_mechanical[0] == (
                        pyunits.convert(
                            b.heating_steam[0].dh_vap_mass
                            * b.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"],
                            to_units=pyunits.kJ / pyunits.s,
                        )
                    )

                self.add_port(name="inlet", block=effect.properties_in)
                self.add_port(name="outlet", block=effect.properties_out)
                self.add_port(name="solids", block=effect.properties_solids)
                self.add_port(name="vapor", block=effect.properties_vapor)
                self.add_port(name="pure_water", block=effect.properties_pure_water)
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
                effect.add_component(
                    f"eq_delta_temperature_inlet_effect_{n}", del_temp_in_constr
                )

                del_temp_out_constr = Constraint(
                    expr=effect.delta_temperature_out[0]
                    == prev_effect.properties_pure_water[0].temperature
                    - effect.properties_in[0].temperature,
                    doc=f"Change in temperature at outlet for effect {n}",
                )
                effect.add_component(
                    f"eq_delta_temperature_outlet_effect_{n}", del_temp_out_constr
                )

                hx_constr = Constraint(
                    expr=prev_effect.energy_flow_superheated_vapor
                    == effect.overall_heat_transfer_coefficient
                    * effect.heat_exchanger_area
                    * effect.delta_temperature[0],
                    doc=f"Heat transfer equation for effect {n}",
                )
                effect.add_component(f"eq_heat_transfer_effect_{n}", hx_constr)

                energy_flow_constr = Constraint(
                    expr=effect.work_mechanical[0]
                    == prev_effect.energy_flow_superheated_vapor,
                    doc=f"Energy supplied to effect {n}",
                )
                effect.add_component(
                    f"eq_energy_for_effect_{n}_from_effect_{n - 1}", energy_flow_constr
                )

        self.total_flow_vol_in = Expression(expr=total_flow_vol_in_expr)

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
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)
        for n, eff in self.effects.items():
            if n == 1:
                assert degrees_of_freedom(eff.effect) == 0
                eff.effect.initialize(**init_args)
                inlet_conc = (
                    eff.effect.properties_in[0]
                    .conc_mass_phase_comp["Liq", "NaCl"]
                    .value
                )
                mass_transfer_coeff = eff.effect.overall_heat_transfer_coefficient.value
            else:
                c = getattr(eff.effect, f"eq_energy_for_effect_{n}_from_effect_{n - 1}")
                c.deactivate()
                eff.effect.initialize(**init_args)
                c.activate()
                eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
                eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
                eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                    inlet_conc
                )
                eff.effect.overall_heat_transfer_coefficient.fix(mass_transfer_coeff)

            init_log.info(f"Initialization of Effect {n} Complete.")

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    @property
    def default_costing_method(self):
        return cost_multi_effect_crystallizer


if __name__ == "__main__":

    from watertap_contrib.reflo.costing import TreatmentCosting
    import watertap.property_models.unit_specific.cryst_prop_pack as props
    from watertap.property_models.water_prop_pack import WaterParameterBlock
    from watertap.core.util.model_diagnostics.infeasible import *
    from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
    from idaes.core.util.scaling import *
    from idaes.core.util.testing import initialization_tester
    from idaes.core.util.scaling import (
        calculate_scaling_factors,
        unscaled_variables_generator,
        badly_scaled_var_generator,
    )

    solver = get_solver()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = props.NaClParameterBlock()
    m.fs.vapor = WaterParameterBlock()

    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Vap", "H2O"))
    m.fs.props.set_default_scaling("flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl"))

    m.fs.mec = mec = MultiEffectCrystallizer(
        property_package=m.fs.props, property_package_vapor=m.fs.vapor
    )
    for n, effects in mec.effects.items():
        print(n, effects)

    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.15
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    crystallizer_yield = 0.5
    operating_pressures = [0.45, 0.25, 0.208, 0.095]
    operating_pressure_eff1 = 0.45  # bar
    operating_pressure_eff2 = 0.25  # bar
    operating_pressure_eff3 = 0.208  # bar
    operating_pressure_eff4 = 0.095  # bar

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff1 = mec.effects[1].effect
    eff2 = mec.effects[2].effect
    eff3 = mec.effects[3].effect
    eff4 = mec.effects[4].effect

    for (n, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        eff.effect.properties_in[0].pressure.fix(feed_pressure)
        eff.effect.properties_in[0].temperature.fix(feed_temperature)

        eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(eps)
        eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(eps)
        eff.effect.properties_in[0].conc_mass_phase_comp[...]
        eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
        # eff.effect.crystal_growth_rate.fix(5e-9)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        # eff.effect.crystal_median_length.fix(0.6e-3)
        eff.effect.crystal_median_length.fix()
        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.set_value(100)
        if n == 1:
            eff.effect.overall_heat_transfer_coefficient.fix(100)
            eff.effect.heating_steam[0].pressure_sat
            eff.effect.heating_steam[0].dh_vap_mass
            eff.effect.heating_steam.calculate_state(
                var_args={
                    ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
                    ("pressure", None): 101325,
                    # ("pressure_sat", None): 28100
                    ("temperature", None): 393,
                },
                hold_state=True,
            )
            eff.effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
            eff.effect.heating_steam[0].flow_vol_phase
        print(f"dof effect {n} = {degrees_of_freedom(eff.effect)}")

    print(f"dof before init = {degrees_of_freedom(m)}")

    calculate_scaling_factors(m)
    try:
        mec.initialize()
    except:
        print_infeasible_constraints(m)

    print(f"DOF before solve  = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    print(f"termination {results.solver.termination_condition}")
    assert_optimal_termination(results)

    # for n, eff in mec.effects.items():
    #     print(f"\nEFFECT {n}\n")
    #     eff.effect.overall_heat_transfer_coefficient.display()
    #     eff.effect.properties_solids[0].flow_mass_phase_comp.display()
    #     eff.effect.temperature_operating.display()

    m.fs.costing = TreatmentCosting()
    m.fs.mec.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_type": "mass_basis"},
    )

    m.fs.costing.nacl_recovered.cost.set_value(-0.024)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(mec.total_flow_vol_in)

    print(f"DOF after costing = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    print(f"termination {results.solver.termination_condition}")
    print(f"LCOW = {m.fs.costing.LCOW()}")

    # m.fs.costing.display()
    # mec.costing.display()
    # eff1.work_mechanical.display()
    # eff1.pressure_operating.display()
    # eff1.heating_steam[0].pressure.display()
    # eff1.heating_steam[0].pressure_sat.display()

    # for n, eff in mec.effects.items():
    #     print(f"\nEFFECT {n}\n")
    #     eff.effect.properties_pure_water[0].temperature.display()
    # eff.effect.height_crystallizer.display()
    # eff.effect.height_slurry.display()
    # eff.effect.volume_suspension.display()
    # eff.effect.t_res.display()
    #     eff.effect.overall_heat_transfer_coefficient.display()
    #     eff.effect.properties_vapor[0].temperature.display()
    #     eff.effect.temperature_operating.display()
    #     eff.effect.properties_in[0].flow_mass_phase_comp.display()
    #     eff.effect.properties_in[0].conc_mass_phase_comp.display()
    #     eff.effect.properties_out[0].flow_mass_phase_comp.display()


##################################################################################################
# "Linking" constraints below that weren't working well but want to keep them here until model is ready for merge
# just in case so I don't have to write them again.

# mass_flow_solid_nacl_constr = Constraint(
#     expr=effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"]
#     == prev_effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"],
#     doc="Mass flow of solid NaCl for effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_mass_flow_sol_nacl_effect_{n}",
#     mass_flow_solid_nacl_constr,
# )

# mass_flow_vap_water_constr = Constraint(
#     expr=effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]
#     == prev_effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"],
#     doc="Mass flow of water vapor for effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_mass_flow_vap_water_effect_{n}",
#     mass_flow_vap_water_constr,
# )

# prop_in_press_constr = Constraint(
#     expr=effect.properties_in[0].pressure
#     == prev_effect.properties_in[0].pressure,
#     doc="Inlet properties pressure for effect {n}",
# )
# effect.add_component(f"eq_equiv_press_effect_{n}", prop_in_press_constr)

# prop_in_press_constr_ub = Constraint(
#     expr=effect.properties_in[0].pressure
#     <= 1.0001 * prev_effect.properties_in[0].pressure,
#     doc="Inlet properties pressure for effect {n}",
# )
# effect.add_component(f"eq_equiv_press_effect_{n}_ub", prop_in_press_constr_ub)

# prop_in_press_constr_lb = Constraint(
#     expr=effect.properties_in[0].pressure
#     >= 0.9999 * prev_effect.properties_in[0].pressure,
#     doc="Inlet properties pressure for effect {n}",
# )
# effect.add_component(f"eq_equiv_press_effect_{n}_lb", prop_in_press_constr_lb)

# prop_in_temp_constr = Constraint(
#     expr=effect.properties_in[0].temperature
#     == prev_effect.properties_in[0].temperature,
#     doc="Inlet properties temperature for effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_temp_effect_{n}", prop_in_temp_constr
# )

# cryst_growth_rate_constr = Constraint(
#     expr=effect.crystal_growth_rate == prev_effect.crystal_growth_rate,
#     doc="Equivalent crystal growth rate effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_crystal_growth_rate_effect_{n}", cryst_growth_rate_constr
# )

# souders_brown_constr = Constraint(
#     expr=effect.souders_brown_constant
#     == prev_effect.souders_brown_constant,
#     doc="Equivalent Sounders Brown constant effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_souders_brown_constant_effect_{n}", souders_brown_constr
# )

# cryst_med_len_constr = Constraint(
#     expr=effect.crystal_median_length
#     == prev_effect.crystal_median_length,
#     doc="Equivalent crystal median length effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_crystal_median_length_effect_{n}", cryst_med_len_constr
# )

# cryst_yield_constr = Constraint(
#     expr=effect.crystallization_yield["NaCl"]
#     == prev_effect.crystallization_yield["NaCl"],
#     doc="Equivalent crystallization yield effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_crystallization_yield_effect_{n}", cryst_yield_constr
# )

# op_press_constr = Constraint(
#     expr=effect.pressure_operating <= prev_effect.pressure_operating,
#     doc=f"Equivalent operating pressure effect {n}",
# )
# effect.add_component(
#     f"eq_equiv_pressure_operating_effect_{n}", op_press_constr
# )

# brine_conc_constr = Constraint(
#     expr=effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
#     >= 0.95* prev_effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
# )
# effect.add_component(f"eq_equiv_brine_conc_effect_{n}_lb", brine_conc_constr)
# brine_conc_constr = Constraint(
#     expr=effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
#     <= 1.05 * prev_effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
# )
# effect.add_component(f"eq_equiv_brine_conc_effect_{n}_ub", brine_conc_constr)
