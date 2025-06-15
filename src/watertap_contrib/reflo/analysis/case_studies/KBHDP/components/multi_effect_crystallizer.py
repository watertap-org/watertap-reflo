# %%
import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    constraint_scaling_transform,
    get_scaling_factor,
    extreme_jacobian_columns,
    extreme_jacobian_rows,
    badly_scaled_var_generator,
    unscaled_variables_generator,
    extreme_jacobian_entries,
)
from watertap.core.util.model_diagnostics import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver

from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock,
    NaClStateBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock,
    WaterStateBlock,
)
from watertap_contrib.reflo.unit_models.multi_effect_crystallizer import (
    MultiEffectCrystallizer,
)
from watertap_contrib.reflo.unit_models.crystallizer_effect import CrystallizerEffect

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)

__all__ = [
    "build_system",
    "build_mec",
    "set_mec_op_conditions",
    "add_mec_costing",
    "init_mec",
    "unfix_mec",
    "mec_rescaling",
]


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = TreatmentCosting()
    build_mec(m, m.fs)

    return m


def build_mec(m, blk) -> None:
    blk.properties = NaClParameterBlock()
    blk.vapor_properties = WaterParameterBlock()

    blk.unit = MultiEffectCrystallizer(
        property_package=blk.properties, property_package_vapor=blk.vapor_properties
    )


def set_mec_op_conditions(
    m,
    blk,
    operating_pressures=[0.45, 0.25, 0.208, 0.095],
    feed_H2O=153.34422736111105,
    feed_NaCl=38.336056840277756,
    nacl_yield=0.9,
    heat_transfer_coeff=100,
) -> None:
    mec = blk.unit

    # Guessed values for initialization
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    ### TOTAL GOING INTO MEC
    flow_mass_phase_water_total = feed_H2O
    flow_mass_phase_salt_total = feed_NaCl
    ### TOTAL INTO EACH EFFECT INITIAL
    """
    Note: In the initial solve of the system, assume the total feed flow rate is 1 kg/s,
    which is align to the default value, in order to guarantee a solution in the initial solve.
    """
    # flow_mass_phase_water_per = 116.2473764168908 / 100 * pyunits.kg / pyunits.s
    # flow_mass_phase_salt_per = 28.478213652777765 / 100 * pyunits.kg / pyunits.s
    flow_mass_phase_water_per = (
        flow_mass_phase_water_total
        / (flow_mass_phase_water_total + flow_mass_phase_salt_total)
        * pyunits.kg
        / pyunits.s
    )
    flow_mass_phase_salt_per = (
        flow_mass_phase_salt_total
        / (flow_mass_phase_water_total + flow_mass_phase_salt_total)
        * pyunits.kg
        / pyunits.s
    )

    saturated_steam_pressure = 101325 * pyunits.Pa + pyunits.convert(
        3 * pyunits.bar, to_units=pyunits.Pa
    )

    ### FIX UNIT MODEL PARAMETERS
    for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            flow_mass_phase_water_per
        )
        eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            flow_mass_phase_salt_per
        )

        eff.effect.properties_in[0].pressure.fix(feed_pressure)
        eff.effect.properties_in[0].temperature.fix(feed_temperature)

        eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
        eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        eff.effect.properties_in[0].conc_mass_phase_comp[...]

        eff.effect.crystallization_yield["NaCl"].fix(nacl_yield)
        eff.effect.crystal_growth_rate.fix()
        eff.effect.souders_brown_constant.fix()
        eff.effect.crystal_median_length.fix()

        eff.effect.pressure_operating.fix(
            pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
        )
        eff.effect.overall_heat_transfer_coefficient.fix(heat_transfer_coeff)

    first_effect = mec.effects[1].effect

    first_effect.overall_heat_transfer_coefficient.fix(heat_transfer_coeff)
    first_effect.heating_steam[0].pressure_sat
    first_effect.heating_steam[0].dh_vap_mass
    first_effect.heating_steam.calculate_state(
        var_args={
            ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
            ("pressure", None): saturated_steam_pressure,
            ("pressure_sat", None): saturated_steam_pressure,
        },
        hold_state=True,
    )
    first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    ### FIX CV PROPERTIES EXCEPT FOR THE LIQUID FLOW RATES
    mec.control_volume.properties_in[0].pressure.fix(feed_pressure)
    mec.control_volume.properties_in[0].temperature.fix(feed_temperature)
    mec.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)

    """
    Check DOF
    """
    for n, eff in mec.effects.items():
        assert degrees_of_freedom(eff.effect) == 0


def init_mec(blk):
    mec = blk.unit

    ### Scale
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    calculate_scaling_factors(blk)

    # from idaes.core.util.scaling import badly_scaled_var_generator
    # bad_var = badly_scaled_var_generator(blk)
    # print('here')
    # for v,l in bad_var:
    #     print(v.name, l)

    # ej = extreme_jacobian_entries(blk, True)
    # for e, c, v in ej:
    #     print(f"{c.name} --> {v.name}: {e}")

    ### INITIALIZE FOR EACH EFFECT
    for n, eff in mec.effects.items():
        eff.effect.initialize()

    ### UNFIX THE INLET FLOW RATES OF EACH EFFECT

    # mec.feed_effect = Var(
    #     mec.Effects,
    #     initialize=0.25,
    #     bounds=(0, 1),
    #     units=pyunits.dimensionless,
    #     doc="Fraction of feed entering each effect",
    # )

    for n, eff in mec.effects.items():
        if n > 1:
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
            eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
            eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

    solver = get_solver()
    results = solver.solve(blk)
    assert_optimal_termination(results)


def unfix_mec(blk):
    mec = blk.unit

    first_effect = mec.effects[1].effect
    ### Release 1st effect flow rate and fix total flow rate instead
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    mec.inlet.temperature[0].unfix()
    mec.inlet.pressure[0].unfix()


def mec_rescaling(
    blk, flow_mass_phase_water_total=116.247, flow_mass_phase_salt_total=28.478
):
    """
    Note: Rescaling is probably needed for extremely large feed flow,
    """
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
        index=("Liq", "NaCl"),
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Vap", "H2O"),
    )
    blk.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_salt_total),
        index=("Sol", "NaCl"),
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Vap", "H2O"),
    )
    blk.vapor_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(flow_mass_phase_water_total),
        index=("Liq", "H2O"),
    )


def add_mec_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing
    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block,
        costing_method_arguments={"cost_work_as": "heat"},
    )


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    print("DOF before solving: ", degrees_of_freedom(m))
    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results


if __name__ == "__main__":
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = TreatmentCosting()
    build_mec(m, m.fs)

    feed_H2O = 111.14146762116363
    feed_NaCl = 28.47821365277779
    set_mec_op_conditions(
        m,
        m.fs,
        operating_pressures=[0.45, 0.25, 0.208, 0.095],
        feed_H2O=feed_H2O,
        feed_NaCl=feed_NaCl,
        nacl_yield=0.9,
        heat_transfer_coeff=1300,
    )

    init_mec(m.fs)
    unfix_mec(m.fs)

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_H2O)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(feed_NaCl)

    m.fs.unit.inlet.temperature[0].fix(273.15 + 30.51)
    m.fs.unit.inlet.pressure[0].fix(101325)
    mec_rescaling(
        m.fs, flow_mass_phase_water_total=feed_H2O, flow_mass_phase_salt_total=feed_NaCl
    )
    add_mec_costing(m, m.fs)

    # ej = extreme_jacobian_entries(m, True)
    # for e, c, v in ej:
    #     print(f"{c.name} --> {v.name}: {e}")

    m.fs.unit.inlet.display()
    print("dof", degrees_of_freedom(m))
    solve(m,tee=False)
    print("feed conc", feed_NaCl / (feed_H2O + feed_NaCl))

    m.fs.costing._find_flow_unit('electricity')

#     # %%
#     conc = []
#     fail = []
#     for i in [110, 120, 130, 140, 150]:
#         for j in [20, 25, 30, 35, 40]:
#             try:
#                 m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(i)
#                 m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(j)
#                 solve(m)
#                 conc.append((i, j, j / (i + j)))
#             except:
#                 fail.append((i, j, j / (i + j)))


# if __name__ == "__main__":
#     feed_pressure = 101325
#     feed_temperature = 273.15 + 20
#     eps = 1e-12
#     m = ConcreteModel()
#     m.fs = FlowsheetBlock(dynamic=False)

#     m.fs.properties = NaClParameterBlock()
#     m.fs.vapor_properties = WaterParameterBlock()

#     m.fs.unit = mec = MultiEffectCrystallizer(
#         property_package=m.fs.properties, property_package_vapor=m.fs.vapor_properties
#     )

#     operating_pressures = [0.45, 0.25, 0.208, 0.095]

#     feed_flow_mass = 1
#     feed_mass_frac_NaCl = 0.15
#     crystallizer_yield = 0.5
#     feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

#     total_feed_flow_mass = 4

#     atm_pressure = 101325 * pyunits.Pa
#     saturated_steam_pressure_gage = 3 * pyunits.bar
#     saturated_steam_pressure = atm_pressure + pyunits.convert(
#         saturated_steam_pressure_gage, to_units=pyunits.Pa
#     )

#     for (_, eff), op_pressure in zip(mec.effects.items(), operating_pressures):
#         eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
#             feed_flow_mass * feed_mass_frac_H2O
#         )
#         eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
#             feed_flow_mass * feed_mass_frac_NaCl
#         )
#         eff.effect.properties_in[0].pressure.fix(feed_pressure)
#         eff.effect.properties_in[0].temperature.fix(feed_temperature)

#         eff.effect.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(eps)
#         eff.effect.properties_in[0].flow_mass_phase_comp["Vap", "H2O"].fix(eps)
#         eff.effect.properties_in[0].conc_mass_phase_comp[...]
#         eff.effect.crystallization_yield["NaCl"].fix(crystallizer_yield)
#         eff.effect.crystal_growth_rate.fix()
#         eff.effect.souders_brown_constant.fix()
#         eff.effect.crystal_median_length.fix()
#         eff.effect.pressure_operating.fix(
#             pyunits.convert(op_pressure * pyunits.bar, to_units=pyunits.Pa)
#         )
#         eff.effect.overall_heat_transfer_coefficient.fix(1300)

#     first_effect = m.fs.unit.effects[1].effect

#     first_effect.overall_heat_transfer_coefficient.fix(1300)
#     first_effect.heating_steam[0].pressure_sat
#     first_effect.heating_steam[0].dh_vap_mass
#     first_effect.heating_steam.calculate_state(
#         var_args={
#             ("flow_mass_phase_comp", ("Liq", "H2O")): 0,
#             ("pressure", None): saturated_steam_pressure,
#             ("pressure_sat", None): saturated_steam_pressure,
#         },
#         hold_state=True,
#     )
#     first_effect.heating_steam[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

#     # m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
#     #     total_feed_flow_mass * feed_mass_frac_H2O
#     # )
#     # m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
#     #     total_feed_flow_mass * feed_mass_frac_NaCl
#     # )

#     m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
#     m.fs.unit.control_volume.properties_in[0].pressure.fix(feed_pressure)
#     m.fs.unit.control_volume.properties_in[0].temperature.fix(feed_temperature)

#     print("dof after setup", degrees_of_freedom(m.fs.unit))

#     for n, eff in mec.effects.items():
#         eff.effect.initialize()

#     for n, eff in mec.effects.items():
#         if n > 1:
#             eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
#             eff.effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
#             eff.effect.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()

#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
#     )
#     m.fs.properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
#     )
#     m.fs.vapor_properties.set_default_scaling(
#         "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
#     )
#     m.fs.vapor_properties.set_default_scaling(
#         "flow_mass_phase_comp", 1, index=("Liq", "H2O")
#     )

#     calculate_scaling_factors(m)

#     first_effect.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
#     first_effect.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

#     mec.inlet.temperature[0].unfix()
#     mec.inlet.pressure[0].unfix()

#     feed_H2O = 0.8 * 200
#     feed_NaCl = 0.2 * 200

#     m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_H2O)
#     m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(feed_NaCl)

#     m.fs.unit.inlet.temperature[0].fix(273.15 + 30.51)
#     m.fs.unit.inlet.pressure[0].fix(101325)
#     # m.fs.unit.initialize()

#     mec_rescaling(
#         m.fs, flow_mass_phase_water_total=feed_H2O, flow_mass_phase_salt_total=feed_NaCl
#     )

#     print("dof before solving: ", degrees_of_freedom(m))

#     # Deactivate unnecessary constraints
#     # for n, eff in m.fs.unit.effects.items():
#     #     eff.effect.properties_vapor[0].eq_flow_vol_phase["Liq"].deactivate()
#     #     eff.effect.properties_vapor[0].eq_flow_vol_phase["Sol"].deactivate()
#     #     eff.effect.properties_solids[0].eq_flow_vol_phase["Liq"].deactivate()
#     #     eff.effect.properties_solids[0].eq_flow_vol_phase["Vap"].deactivate()
#     #     eff.effect.properties_out[0].eq_flow_vol_phase["Vap"].deactivate()
#     #     eff.effect.properties_out[0].eq_flow_vol_phase["Sol"].deactivate()

#     #     eff.effect.eq_enthalpy_balance.deactivate()

#     ej = extreme_jacobian_entries(m, True)
#     for e, c, v in ej:
#         print(f"{c.name} --> {v.name}: {e}")

#     # m.dummy_obj = Objective(expr=0)
#     # from idaes.core.util import DiagnosticsToolbox

#     # dt = DiagnosticsToolbox(m)
#     # dt.report_structural_issues()
#     # dt.display_variables_with_extreme_values()
#     # dt.report_numerical_issues()
#     # dt.display_variables_with_extreme_jacobians()
#     # dt.display_constraints_with_extreme_jacobians()

#     solver = get_solver()
#     results = solver.solve(m)
#     assert_optimal_termination(results)
# # %%