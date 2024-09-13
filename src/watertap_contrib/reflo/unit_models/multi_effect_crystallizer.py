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

from watertap.core.util.model_diagnostics.infeasible import * 
# Import Pyomo libraries
from pyomo.environ import (
    ConcreteModel,
    Var,
    check_optimal_termination,
    Param,
    Constraint,
    Suffix,
    RangeSet,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
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

from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect


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

        self.current_effect = Var()

        for n, eff in self.effects.items():
            self.current_effect.fix(n)
            eff.effect = effect = CrystallizerEffect(
                property_package=self.config.property_package,
                property_package_vapor=self.config.property_package_vapor,
                standalone=False,
            )
            if n == self.first_effect:

                tmp_dict = dict(**self.config.property_package_args)
                tmp_dict["has_phase_equilibrium"] = False
                tmp_dict["parameters"] = self.config.property_package_vapor
                tmp_dict["defined_state"] = False

                effect.heating_steam = self.config.property_package_vapor.state_block_class(
                    self.flowsheet().config.time,
                    doc="Material properties of inlet heating steam",
                    **tmp_dict,
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
                        * b.area
                        * b.delta_temperature[0]
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
                    * effect.area
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

                mass_flow_solid_nacl_constr = Constraint(
                    expr=effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"]
                    == prev_effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"],
                    doc="Mass flow of solid NaCl for effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_mass_flow_sol_nacl_effect_{n}",
                    mass_flow_solid_nacl_constr,
                )

                mass_flow_vap_water_constr = Constraint(
                    expr=effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"]
                    == prev_effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"],
                    doc="Mass flow of water vapor for effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_mass_flow_vap_water_effect_{n}",
                    mass_flow_vap_water_constr,
                )

                prop_in_press_constr = Constraint(
                    expr=effect.properties_in[0].pressure
                    == prev_effect.properties_in[0].pressure,
                    doc="Inlet properties pressure for effect {n}",
                )
                effect.add_component(f"eq_equiv_press_effect_{n}", prop_in_press_constr)

                prop_in_flow_mass_water_constr = Constraint(
                    expr=effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                    == prev_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"],
                    doc="Inlet properties mass flow water for effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_mass_flow_liq_water_effect_{n}", prop_in_flow_mass_water_constr
                )

                prop_in_flow_mass_nacl_constr = Constraint(
                    expr=effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
                    == prev_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"],
                    doc="Inlet properties mass flow NaCl for effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_mass_flow_liq_nacl_effect_{n}", prop_in_flow_mass_nacl_constr
                )

                prop_in_temp_constr = Constraint(
                    expr=effect.properties_in[0].temperature
                    == prev_effect.properties_in[0].temperature,
                    doc="Inlet properties temperature for effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_temp_effect_{n}", prop_in_temp_constr
                )

                # steam_temp_sat_constr = Constraint(
                #     expr=effect.heating_steam[0].temperature
                #     == prev_effect.heating_steam[0].temperature,
                #     doc="Steam saturation temperature for effect {n}",
                # )
                # effect.add_component(
                #     f"eq_equiv_steam_temp_sat_effect_{n}", steam_temp_sat_constr
                # )

                # steam_press_sat_constr = Constraint(
                #     expr=effect.heating_steam[0].pressure_sat
                #     == prev_effect.heating_steam[0].pressure_sat,
                #     doc="Steam saturation pressure for effect {n}",
                # )
                # effect.add_component(
                #     f"eq_equiv_steam_press_sat_effect_{n}", steam_press_sat_constr
                # )

                cryst_growth_rate_constr = Constraint(
                    expr=effect.crystal_growth_rate == prev_effect.crystal_growth_rate,
                    doc="Equivalent crystal growth rate effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_crystal_growth_rate_effect_{n}", cryst_growth_rate_constr
                )

                souders_brown_constr = Constraint(
                    expr=effect.souders_brown_constant
                    == prev_effect.souders_brown_constant,
                    doc="Equivalent Sounders Brown constant effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_souders_brown_constant_effect_{n}", souders_brown_constr
                )

                cryst_med_len_constr = Constraint(
                    expr=effect.crystal_median_length
                    == prev_effect.crystal_median_length,
                    doc="Equivalent crystal median length effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_crystal_median_length_effect_{n}", cryst_med_len_constr
                )

                cryst_yield_constr = Constraint(
                    expr=effect.crystallization_yield["NaCl"]
                    == prev_effect.crystallization_yield["NaCl"],
                    doc="Equivalent crystallization yield effect {n}",
                )
                effect.add_component(
                    f"eq_equiv_crystallization_yield_effect_{n}", cryst_yield_constr
                )

                # op_press_constr = Constraint(
                #     expr=effect.pressure_operating <= prev_effect.pressure_operating,
                #     doc=f"Equivalent operating pressure effect {n}",
                # )
                # effect.add_component(
                #     f"eq_equiv_pressure_operating_effect_{n}", op_press_constr
                # )


print(list(range(1, 5, 1)))
if __name__ == "__main__":
    import watertap.property_models.unit_specific.cryst_prop_pack as props
    from watertap.property_models.water_prop_pack import WaterParameterBlock

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

    feed_flow_mass = 2
    feed_mass_frac_NaCl = 0.15
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    crystallizer_yield = 0.5
    operating_pressures = [0.78, 0.25, 0.208, 0.095]
    operating_pressure_eff1 = 0.78  # bar
    operating_pressure_eff2 = 0.25  # bar
    operating_pressure_eff3 = 0.208  # bar
    operating_pressure_eff4 = 0.095  # bar

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    eps = 1e-6

    eff1 = mec.effects[1].effect
    eff2 = mec.effects[2].effect
    eff3 = mec.effects[3].effect
    eff4 = mec.effects[4].effect

    mec.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    mec.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    mec.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    mec.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)

    mec.inlet.pressure[0].fix(feed_pressure)
    mec.inlet.temperature[0].fix(feed_temperature)

    eff1.heating_steam[0].pressure_sat

    eff1.heating_steam.calculate_state(
        var_args={
            ("pressure", None): 101325,
            # ("pressure_sat", None): 28100
            ("temperature", None): 393,
        },
        hold_state=True,
    )
    eff1.crystallization_yield["NaCl"].fix(crystallizer_yield)
    eff1.crystal_growth_rate.fix()
    eff1.souders_brown_constant.fix()
    eff1.crystal_median_length.fix()
    eff1.overall_heat_transfer_coefficient.fix(100)

    # eff.pressure_operating.fix(operating_pressure_eff1 * pyunits.bar)
    for (n, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        print(f"dof effect {n} = {degrees_of_freedom(eff.effect)}")

    print(f"dof = {degrees_of_freedom(m)}")
    calculate_scaling_factors(m)
    eff1.initialize()
    # eff1.properties_in[0].display()
    cvc(eff2.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], eff2.eq_equiv_mass_flow_liq_water_effect_2)
    cvc(eff2.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"], eff2.eq_equiv_mass_flow_liq_nacl_effect_2)
    cvc(eff2.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"], eff2.eq_equiv_mass_flow_sol_nacl_effect_2)
    cvc(eff2.properties_in[0].flow_mass_phase_comp["Vap", "H2O"], eff2.eq_equiv_mass_flow_vap_water_effect_2)
    cvc(eff2.properties_in[0].temperature, eff2.eq_equiv_temp_effect_2)
    # eff2.initialize()
    # eff2.properties_in[0].display()
    # mec.current_effect.display()
    # eff2.del_component(eff2.properties_pure_water)
    # eff2.properties_pure_water.display()
    # try:
    flags = eff2.properties_in.initialize(hold_state=True)
    eff2.properties_out.initialize(hold_state=False)
    eff2.properties_solids.initialize(hold_state=False)
    eff2.properties_vapor.initialize(hold_state=False)
    eff2.properties_out[0].display()
    eff2.properties_solids[0].phase_component_set.display()
    # print(flags)
    # except: 
    #     eff2.properties_in[0].display()
    # eff1.properties_in[0].display()
    # eff1.heating_steam[0].display()
    # print(f"dof = {degrees_of_freedom(m)}")
    # for _, eff in mec.effects.items():
    #     try:
    #         eff.effect.initialize()
    #     except:
    #         pass
            # print_infeasible_constraints(m)
            # print_variables_close_to_bounds(m)
            # raise
