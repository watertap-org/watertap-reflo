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

# Pyomo imports
from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
<<<<<<< Updated upstream
=======
    Expression,
    Reference,
>>>>>>> Stashed changes
    Var,
    value,
    assert_optimal_termination,
)
<<<<<<< Updated upstream
=======
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
    )
>>>>>>> Stashed changes
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# IDAES imports
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.scaling import (
    calculate_scaling_factors,
<<<<<<< Updated upstream
)
from idaes.core import FlowsheetBlock, MaterialBalanceType
=======
    list_badly_scaled_variables,
)
from idaes.core import FlowsheetBlock, MaterialBalanceType, EnergyBalanceType
>>>>>>> Stashed changes
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models import (
    Mixer,
    Separator,
<<<<<<< Updated upstream
    HeatExchanger
)
# from idaes.models.unit_models.heat_exchanger import HX0DInitializer
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
=======
    HeatExchanger,
    Heater,
)
# from idaes.models.unit_models.heat_exchanger import HX0DInitializer
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback, delta_temperature_underwood_callback
>>>>>>> Stashed changes

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import FODrawSolutionParameterBlock
from watertap_contrib.reflo.unit_models.zero_order.forward_osmosis_zo import ForwardOsmosisZO

<<<<<<< Updated upstream
def build_fo_treviflowsheet(
=======
from idaes.core.util import DiagnosticsToolbox

def build_fo_trevi_flowsheet(
>>>>>>> Stashed changes
    m=None,

):
    """
    This function builds a flowsheet as a representative of Trevi's FO system configuration

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.seawater_properties = SeawaterParameterBlock()
    m.fs.draw_solution_properties = FODrawSolutionParameterBlock()
<<<<<<< Updated upstream

    add_fo(m.fs)

    add_HX(m.fs)

    return m

def add_HX(fs):
    fs.HX1A = HeatExchanger(
        delta_temperature_callback=delta_temperature_lmtd_callback,
        hot_side_name="product_water",
        cold_side_name="weak_draw",
        product_water={"property_package": fs.draw_solution_properties},
        weak_draw={"property_package": fs.draw_solution_properties}
    )
=======
    mfs = m.fs

    RO_recovery = 0.9
    NF_recovery = 0.8


    # Add FO module
    add_fo(mfs)

    def separation_heat(t):
        return 105000
    
    # Add heat exchangers
    mfs.HX1A = HeatExchanger(
        delta_temperature_callback=delta_temperature_underwood_callback,
        hot_side_name="product_water",
        cold_side_name="weak_draw",
        product_water={"property_package": m.fs.draw_solution_properties},
        weak_draw={"property_package": m.fs.draw_solution_properties},
    ) 
    
    # mfs.HX1A.cold_side.properties_in[0].liquid_separation = 0
    # mfs.HX1A.cold_side.properties_out[0].liquid_separation = 0

    mfs.HX2A = HeatExchanger(
        delta_temperature_callback=delta_temperature_underwood_callback,
        hot_side_name="reg_draw",
        cold_side_name="weak_draw",
        reg_draw={"property_package": m.fs.draw_solution_properties},
        weak_draw={"property_package": m.fs.draw_solution_properties}
    )

    
    @mfs.Constraint(doc = 'Same outlet temperature from HX 1A and 2A')
    def outlet_temp_HX1(b):
        return (b.HX1A.weak_draw_outlet.temperature[0] == b.HX2A.weak_draw_outlet.temperature[0])
    
    # Add a separator to represent RO and NF fed with the product water to remove remained
    # draw solution (NF) and boron (RO), if existing
    m.fs.S2 = Separator(
        property_package=m.fs.draw_solution_properties,
        mixed_state_block=None,
        outlet_list=["RO_reject", "NF_reject", "fresh_water"],
        split_basis = 3, # Component flow
    )

    @mfs.Constraint(mfs.draw_solution_properties.component_list,
                    doc = 'NF permeate is sent to RO for further treatment')
    def S2_RO_reject(b,j):
        permeate_coeff = {"H2O": NF_recovery * (1 - RO_recovery),
                          "DrawSolution": 0}
        return (b.S2.RO_reject.flow_mass_phase_comp[0, "Liq", j]
                == b.S2.inlet.flow_mass_phase_comp[0, "Liq", j]
                 * permeate_coeff[j]
            )

    @mfs.Constraint(mfs.draw_solution_properties.component_list,
                    doc = 'NF removes all draw solution and the reject is recirculated')
    def S2_NF_reject(b,j):
        permeate_coeff = {"H2O": 1 - NF_recovery,
                          "DrawSolution": 1}
        return (b.S2.NF_reject.flow_mass_phase_comp[0, "Liq", j]
                == b.S2.inlet.flow_mass_phase_comp[0, "Liq", j]
                 * permeate_coeff[j]
                )

    # Add mixer to mix NF reject that contains draw solution and weak draw
    mfs.M1 = Mixer(
        property_package=m.fs.draw_solution_properties,
        material_balance_type=MaterialBalanceType.componentPhase,
        # momentum_mixing_type = 2,
        energy_mixing_type=1,
        inlet_list=["NF_reject", "weak_draw"],
    )

    # Add a separator to diverge the weak draw solution from FO module for heat recovery
    mfs.S1 = Separator(
        property_package=m.fs.draw_solution_properties,
        mixed_state_block=m.fs.M1.mixed_state,
        outlet_list=["to_HX1A", "to_HX2A"],
        split_basis = 1, # Total flow
    )

    # Add mixer to combine pre-heated weak draw solution
    mfs.M2 = Mixer(
        property_package=m.fs.draw_solution_properties,
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=1,
        inlet_list=["HX1A", "HX2A"],
    )

    # Add connections
    mfs.HX1A_cold_inlet = Arc(
        source=m.fs.S1.to_HX1A, destination=m.fs.HX1A.cold_side_inlet
    )
    mfs.HX1A_hot_inlet = Arc(
        source=m.fs.fo.product, destination=m.fs.HX1A.hot_side_inlet
    )
    mfs.HX2A_cold_inlet = Arc(
        source=m.fs.S1.to_HX2A, destination=m.fs.HX2A.cold_side_inlet
    )
    mfs.HX2A_hot_inlet = Arc(
        source=m.fs.fo.reg_draw, destination=m.fs.HX2A.hot_side_inlet
    )
    mfs.S2_inlet = Arc(
        source=m.fs.HX1A.hot_side_outlet, destination=m.fs.S2.inlet
    )
    mfs.NF_reject_to_M1 = Arc(
        source=m.fs.S2.NF_reject, destination=m.fs.M1.NF_reject
    )
    mfs.weak_draw_to_M1 = Arc(
        source=m.fs.fo.weak_draw, destination=m.fs.M1.weak_draw
    )
    mfs.HX1A_to_M2 = Arc(
        source=m.fs.HX1A.cold_side_outlet, destination=m.fs.M2.HX1A
    )
    mfs.HX2A_to_M2 = Arc(
        source=m.fs.HX2A.cold_side_outlet, destination=m.fs.M2.HX2A
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


    mfs.HX1A.area.fix(50000)    
    mfs.HX1A.overall_heat_transfer_coefficient[0].fix(100)
    mfs.HX2A.area.fix(50000)    
    mfs.HX2A.overall_heat_transfer_coefficient[0].fix(100)
    # mfs.HX1A.heat_transfer_equation.deactivate()
    # mfs.HX2A.heat_transfer_equation.deactivate()

    iscale.set_scaling_factor(mfs.HX1A.area, 1e-5)
    iscale.set_scaling_factor(mfs.HX1A.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(mfs.HX1A.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX1A.cold_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX2A.area, 1e-5)
    iscale.set_scaling_factor(mfs.HX2A.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(mfs.HX2A.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX2A.cold_side.heat, 1e-7)

    return m


>>>>>>> Stashed changes

def add_fo(fs):
    """
    This function adds a FO module to the current flowsheet

    Returns:
        object: A FO zero-order model
    """
    fs.fo = ForwardOsmosisZO(
        property_package_water = fs.seawater_properties,
        property_package_draw_solution = fs.draw_solution_properties,
    )

    fo = fs.fo

    # System specifications
    recovery_ratio = 0.3  # Assumed FO recovery ratio
<<<<<<< Updated upstream
=======
    nanofiltration_recovery_ratio = 0.8 # Nanofiltration recovery ratio
>>>>>>> Stashed changes
    dp_brine = 0  # Required pressure over brine osmotic pressure (Pa)
    heat_mixing = 105  # Heat of mixing in the membrane (MJ/m3 product)
    reneration_temp = 90  # Separation temperature of the draw solution (C)
    separator_temp_loss = 1  # Temperature loss in the separator (K)
    feed_temperature = 13  # Feed water temperature (C)
    feed_vol_flow = 3.704  # Feed water volumetric flow rate (m3/s)
    feed_TDS_mass = 0.035  # TDS mass fraction of feed
    strong_draw_temp = 20  # Strong draw solution inlet temperature (C)
    strong_draw_mass = 0.8  # Strong draw solution mass fraction
    product_draw_mass = 0.01  # Mass fraction of draw in the product water


    fo.recovery_ratio.fix(recovery_ratio)
<<<<<<< Updated upstream
=======
    fo.nanofiltration_recovery_ratio.fix(nanofiltration_recovery_ratio)
>>>>>>> Stashed changes
    fo.dp_brine.fix(dp_brine)
    fo.heat_mixing.fix(heat_mixing)
    fo.regeneration_temp.fix(reneration_temp + 273.15)
    fo.separator_temp_loss.fix(separator_temp_loss)

    # Specifyf strong draw solution properties
    fo.strong_draw_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "DrawSolution")): strong_draw_mass,
            ("temperature", None): strong_draw_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    fo.strong_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].unfix()

    # Specifyf product water properties
    fo.product_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1,
            ("mass_frac_phase_comp", ("Liq", "DrawSolution")): product_draw_mass,
            ("temperature", None): reneration_temp - separator_temp_loss + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    fo.product_props[0].temperature.unfix()

    # Specify feed properties
    fo.feed_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): feed_vol_flow,
            ("mass_frac_phase_comp", ("Liq", "TDS")): feed_TDS_mass,
            ("temperature", None): feed_temperature + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

<<<<<<< Updated upstream
    # Set scaling factors for mass flow rates
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "DrawSolution")
=======
    fo.weak_draw_props[0].pressure.fix(101325)
    fo.reg_draw_props[0].pressure.fix(101325)

    # Set scaling factors for mass flow rates
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "TDS")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Liq", "DrawSolution")
>>>>>>> Stashed changes
    )


def fix_dof_and_initialize(m, outlvl=idaeslog.WARNING):
    calculate_scaling_factors(m)

    m.fs.fo.initialize()
    # Unfix the state variables and fix mass fractrion of two state blocks
    strong_draw_mass = 0.8  # Strong draw solution mass fraction
    product_draw_mass = 0.01  # Mass fraction of draw in the product water
    m.fs.fo.unfix_and_fix_freedom(strong_draw_mass, product_draw_mass)


<<<<<<< Updated upstream
if __name__=="__main__":
    m = build_fo_treviflowsheet()
    fix_dof_and_initialize(m)
=======
    m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value  
    m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value  
    m.fs.S1.to_HX1A.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value / 2 
    m.fs.S1.to_HX1A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value / 2 
    m.fs.S1.to_HX2A.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value / 2 
    m.fs.S1.to_HX2A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value / 2
    m.fs.S1.to_HX1A_state[0].flow_vol_phase["Liq"]
    m.fs.S1.to_HX2A_state[0].flow_vol_phase["Liq"]
    m.fs.S1.to_HX1A_state[0].cp_mass_phase["Liq"]    
    m.fs.S1.initialize()

    state_args_1A_hot = {"flow_mass_phase_comp":{("Liq", "H2O"): m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value,
                                            ("Liq", "DrawSolution"): m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value},
                    "temperature": m.fs.fo.product_props[0].temperature.value,
                    "pressure": m.fs.fo.product_props[0].pressure.value}
    state_args_1A_cold = {"flow_mass_phase_comp":{("Liq", "H2O"): m.fs.S1.to_HX1A.flow_mass_phase_comp[0, "Liq", "H2O"].value,
                                            ("Liq", "DrawSolution"): m.fs.S1.to_HX1A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value},
                    "temperature": m.fs.S1.to_HX1A.temperature[0].value,
                    "pressure": m.fs.S1.to_HX1A.pressure[0].value}
    m.fs.HX1A.initialize(state_args_1 = state_args_1A_hot, state_args_2 = state_args_1A_cold)
    m.fs.HX1A.cold_side.properties_out[0].liquid_separation = 1

    state_args_2A_hot = {"flow_mass_phase_comp":{("Liq", "H2O"): m.fs.fo.reg_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value,
                                            ("Liq", "DrawSolution"): m.fs.fo.reg_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value},
                    "temperature": m.fs.fo.reg_draw_props[0].temperature.value,
                    "pressure": m.fs.fo.reg_draw_props[0].pressure.value}
    state_args_2A_cold = {"flow_mass_phase_comp":{("Liq", "H2O"): m.fs.S1.to_HX2A.flow_mass_phase_comp[0, "Liq", "H2O"].value,
                                            ("Liq", "DrawSolution"): m.fs.S1.to_HX2A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value},
                    "temperature": m.fs.S1.to_HX2A.temperature[0].value,
                    "pressure": m.fs.S1.to_HX2A.pressure[0].value}
    m.fs.HX2A.initialize(state_args_1 = state_args_2A_hot, state_args_2 = state_args_2A_cold)
    m.fs.HX2A.cold_side.properties_out[0].liquid_separation = 1

    m.fs.S2.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value  
    m.fs.S2.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value 
    m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value  
    m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value 
    m.fs.S2.RO_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value * 1e-1
    m.fs.S2.RO_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value  * 1e-1
    m.fs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    m.fs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value 
    m.fs.S2.initialize()

    m.fs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    m.fs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    m.fs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value
    m.fs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value
    m.fs.M1.outlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    m.fs.M1.outlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    m.fs.M1.mixed_state[0].flow_vol_phase["Liq"]
    
    m.fs.M1.initialize(outlvl = idaeslog.DEBUG)



    m.fs.M2.HX1A.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    m.fs.M2.HX1A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    m.fs.M2.HX2A.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value
    m.fs.M2.HX2A.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value
    m.fs.M2.outlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    m.fs.M2.outlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    m.fs.M2.initialize()

    print('flowsheet initalization completed')

if __name__=="__main__":
    # from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
    # from pyomo.environ import Objective
    # m.fs.o = Objective(expr = 0)
    # nlp = PyomoNLP(m)
    # jac = nlp.evaluate_jacobian()
    # var_list = nlp.primals_names()
    # con_list = nlp.constraint_names()
    # jac_array = jac.toarray()
    # import pandas as pd
    # df = pd.DataFrame(jac_array, index = con_list, columns=var_list)
    # # if has_touched_vars:
    # #     msg= "_with_touched_vars"
    # # else:
    # #     msg = ""    
    # df.to_csv(f'jacobian.csv')

    m = build_fo_trevi_flowsheet()
    fix_dof_and_initialize(m)
    
    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()

    # try:
    #     fix_dof_and_initialize(m)
    # except:
    #     pass

    # dt.report_numerical_issues()

    # display_constraints_with_large_residuals()
    # display_variables_at_or_outside_bounds()
    # display_variables_with_extreme_jacobians()
    # display_constraints_with_extreme_jacobians()




    # print(m.fs.HX1A._get_stream_table_contents())
    # print(m.fs.HX1A.heat_duty[0].value)

    print('no. of var', number_variables(m))
    print('no. of const', 
        number_total_constraints(m),)
    print('no of unused var', number_unused_variables(m))

>>>>>>> Stashed changes

    from watertap.core.util.initialization import check_dof
    check_dof(m, fail_flag=True)

    solver = get_solver()

<<<<<<< Updated upstream
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(m.fs.fo.strong_draw_props[0].flow_vol_phase["Liq"].value)

=======
    # results = solver.solve(m)
    # assert_optimal_termination(results)

    m.fs.HX1A.overall_heat_transfer_coefficient[0].unfix()
    m.fs.HX2A.overall_heat_transfer_coefficient[0].unfix()
    m.fs.HX1A.weak_draw_outlet.temperature.fix(70 + 273.15)
    m.fs.HX1A.product_water_outlet.temperature.fix(23 + 273.15)
    # m.fs.HX2A.reg_draw_outlet.temperature.fix(23 + 273.15)


    check_dof(m, fail_flag=True)
    results = solver.solve(m)
    assert_optimal_termination(results)

    # def separation_heat(t):
    #     return 0
    
    # mfs = m.fs
    # cold_side = m.fs.HX1A.cold_side
    # # cold_side.del_component(cold_side.heat)
    # # mfs.HX1A.del_component(mfs.HX1A.heat_duty)
    # cold_side.del_component(cold_side.enthalpy_balances)
    
    # cold_side.add_total_enthalpy_balances(
    #     # balance_type=EnergyBalanceType.enthalpyTotal, 
    #     has_heat_transfer=True,
    #     custom_term = separation_heat,
    # )
    # mfs.HX1A.heat_duty = Reference(cold_side.heat)

    # mfs = m.fs
    # mfs.HX1A.cold_side.del_component(mfs.HX1A.cold_side.enthalpy_balances)
    
    # mfs.HX1A.cold_side.add_energy_balances(
    #     balance_type=EnergyBalanceType.enthalpyTotal, 
    #     has_heat_transfer=True,
    #     custom_term = separation_heat,
    # )
    # mfs.HX1A.heat_duty = Reference(cold_side.heat)

    # cold_side = mfs.HX1A.cold_side

    # results = solver.solve(m)
    # assert_optimal_termination(results)



    print('1A cold in temp', m.fs.HX1A.weak_draw_inlet.temperature[0].value - 273.15)
    print('1A cold out temp', m.fs.HX1A.weak_draw_outlet.temperature[0].value - 273.15)
    print('1A hot in temp', m.fs.HX1A.product_water_inlet.temperature[0].value - 273.15)
    print('1A hot out temp', m.fs.HX1A.product_water_outlet.temperature[0].value - 273.15)
    print('')
    print('2A cold in temp', m.fs.HX2A.weak_draw_inlet.temperature[0].value - 273.15)
    print('2A cold out temp', m.fs.HX2A.weak_draw_outlet.temperature[0].value - 273.15)
    print('2A hot in temp', m.fs.HX2A.reg_draw_inlet.temperature[0].value - 273.15)
    print('2A hot out temp', m.fs.HX2A.reg_draw_outlet.temperature[0].value - 273.15)
    print('')
    print('FO product vol', m.fs.fo.product_props[0].flow_vol_phase["Liq"].value)
    print('FO product cp', m.fs.fo.product_props[0].cp_mass_phase["Liq"].value)
    print('FO reg draw vol', m.fs.fo.reg_draw_props[0].flow_vol_phase["Liq"].value)
    print('FO reg draw cp', m.fs.fo.reg_draw_props[0].cp_mass_phase["Liq"].value)
    print('FO weak draw vol', m.fs.fo.weak_draw_props[0].flow_vol_phase["Liq"].value)
    print('FO weak draw cp', m.fs.fo.weak_draw_props[0].cp_mass_phase["Liq"].value)
    print('')
    print('FO weak H2O', m.fs.fo.weak_draw.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('FO weak draw', m.fs.fo.weak_draw.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)
    print('FO weak temp', m.fs.fo.weak_draw.temperature[0].value - 273.15)
    print('')
    print('FO mixed H2O', m.fs.S1.inlet.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('FO mixed draw', m.fs.S1.inlet.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)
    print('FO mixed temp', m.fs.S1.inlet.temperature[0].value - 273.15)

    print('')
    print('1A cold in H2O', m.fs.HX1A.cold_side_inlet.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('1A cold in Draw', m.fs.HX1A.cold_side_inlet.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)
    print('1A cold vol',   m.fs.S1.to_HX1A_state[0].flow_vol_phase["Liq"].value)
    print('1A hot in H2O', m.fs.HX1A.hot_side_inlet.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('1A hot in Draw', m.fs.HX1A.hot_side_inlet.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)
    print('')
    print('2A cold in H2O', m.fs.HX2A.cold_side_inlet.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('2A cold in Draw', m.fs.HX2A.cold_side_inlet.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)
    print('2A cold vol',   m.fs.S1.to_HX2A_state[0].flow_vol_phase["Liq"].value)
    print('2A hot in H2O', m.fs.HX2A.hot_side_inlet.flow_mass_phase_comp[0,"Liq","H2O"].value)
    print('2A hot in Draw', m.fs.HX2A.hot_side_inlet.flow_mass_phase_comp[0,"Liq","DrawSolution"].value)

    print('')
    print('HX1A hot side heat load (MJ)', value(m.fs.HX1A.hot_side.heat[0]) / 1e6)
    print('HX1A cold side heat load (MJ)', value(m.fs.HX1A.cold_side.heat[0]) / 1e6)
    print('HX2A hot side heat load (MJ)', value(m.fs.HX2A.hot_side.heat[0]) / 1e6)
    print('HX2A cold side heat load (MJ)', value(m.fs.HX2A.cold_side.heat[0]) / 1e6)
    print('HX1A approach T', value(m.fs.HX1A.delta_temperature[0]))

    print('')
    print('HX1A area', m.fs.HX1A.area.value)
    print('HX2A area', m.fs.HX2A.area.value)
    print('HX1A coef', m.fs.HX2A.overall_heat_transfer_coefficient[0].value)
    print('HX2A coef', m.fs.HX2A.overall_heat_transfer_coefficient[0].value)

    print('')
    print('Fresh water mass', m.fs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    print('RO reject mass', m.fs.S2.RO_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    print('NF reject H2O mass', m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    print('NF reject draw mass', m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value)

    print('')
    print('NF reject temp', m.fs.S2.NF_reject.temperature[0].value - 273.15)
    print('weak draw temp', m.fs.fo.weak_draw_props[0].temperature.value - 273.15)
    print('M1 temperature', m.fs.M1.mixed_state[0].temperature.value -273.15)
    print('M1 H2O mass', m.fs.M1.mixed_state[0].flow_mass_phase_comp['Liq', 'H2O'].value)
    print('M1 H2O mass', m.fs.M1.mixed_state[0].flow_mass_phase_comp['Liq', 'DrawSolution'].value)

    print('')
    print('M2 mixed H2O mass', m.fs.M2.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"].value)
    print('M2 mixed draw mass', m.fs.M2.mixed_state[0].flow_mass_phase_comp["Liq", "DrawSolution"].value)
    print('M2 mixed temp', m.fs.M2.mixed_state[0].temperature.value - 273.15)



    # m = ConcreteModel()
    # m.fs = FlowsheetBlock(dynamic=False)
    # m.fs.seawater_properties = SeawaterParameterBlock()
    # m.fs.draw_solution_properties = FODrawSolutionParameterBlock()
    # mfs = m.fs
    # # Add mixer to mix NF reject that contains draw solution and weak draw
    # mfs.M1 = Mixer(
    #     property_package=m.fs.draw_solution_properties,
    #     material_balance_type=MaterialBalanceType.componentTotal,
    #     energy_mixing_type=1,
    #     inlet_list=["NF_reject", "weak_draw"],
    # )

    # mfs.M1.NF_reject.temperature[0].fix(47.8+273.15)
    # mfs.M1.NF_reject.pressure[0].fix(101325)
    # mfs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].fix(223.2)
    # mfs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].fix(0.01)

    # mfs.M1.weak_draw.temperature[0].fix(20.81+273.15)
    # mfs.M1.weak_draw.pressure[0].fix(101325)
    # mfs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1479.96)
    # mfs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "DrawSolution"].fix(1516.85)

    # # Set scaling factors for mass flow rates
    # mfs.seawater_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
    # )
    # mfs.seawater_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-3, index=("Liq", "TDS")
    # )

    # mfs.draw_solution_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
    # )
    # mfs.draw_solution_properties.set_default_scaling(
    #     "flow_mass_phase_comp", 1e-1, index=("Liq", "DrawSolution")
    # )

    # from idaes.core.util.scaling import (
    #     calculate_scaling_factors,
    #     constraint_scaling_transform,
    #     unscaled_variables_generator,
    #     unscaled_constraints_generator,
    #     badly_scaled_var_generator,
    #     list_badly_scaled_variables,
    # )
    # badly_scaled_var_lst = list(badly_scaled_var_generator(m))
    # list_badly_scaled_variables(m)

    # # @mfs.Constraint(
    # #     doc="Temperature difference")
    # # def temp_bound2(b):
    # #     return b.M1.weak_draw_state[0].mass_frac_phase_comp["Liq", "DrawSolution"] <= 0.6
    
    # # @mfs.Constraint(
    # #     doc="Temperature difference")
    # # def temp_bound3(b):
    # #     return b.M1.mixed_state[0].mass_frac_phase_comp["Liq", "DrawSolution"] <= 0.6



    # mfs.M1.weak_draw.temperature[0].value = 273.15 + 23
    # mfs.M1.weak_draw.pressure[0].value = 101325
    # mfs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "H2O"].value = 1700
    # mfs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = 1528
    # mfs.M1.initialize()

    # from watertap.core.util.initialization import check_dof
    # check_dof(m, fail_flag=True)

    # solver = get_solver()

    # results = solver.solve(m)
    # assert_optimal_termination(results)   

    # print(mfs.M1.mixed_state[0].temperature.value - 273.15)
    # print(mfs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"].value)
    # print(mfs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "DrawSolution"].value)
    # print(mfs.M1.NF_reject_state[0].mass_frac_phase_comp["Liq", "DrawSolution"].value)
    # print(mfs.M1.weak_draw_state[0].mass_frac_phase_comp["Liq", "DrawSolution"].value)
    # print(mfs.M1.mixed_state[0].mass_frac_phase_comp["Liq", "DrawSolution"].value)
>>>>>>> Stashed changes
