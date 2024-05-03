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
    Var,
    value,
    assert_optimal_termination,
)
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# IDAES imports
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models import (
    Mixer,
    Separator,
    HeatExchanger
)
# from idaes.models.unit_models.heat_exchanger import HX0DInitializer
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import FODrawSolutionParameterBlock
from watertap_contrib.reflo.unit_models.zero_order.forward_osmosis_zo import ForwardOsmosisZO

def build_fo_treviflowsheet(
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
    )


def fix_dof_and_initialize(m, outlvl=idaeslog.WARNING):
    calculate_scaling_factors(m)

    m.fs.fo.initialize()
    # Unfix the state variables and fix mass fractrion of two state blocks
    strong_draw_mass = 0.8  # Strong draw solution mass fraction
    product_draw_mass = 0.01  # Mass fraction of draw in the product water
    m.fs.fo.unfix_and_fix_freedom(strong_draw_mass, product_draw_mass)


if __name__=="__main__":
    m = build_fo_treviflowsheet()
    fix_dof_and_initialize(m)

    from watertap.core.util.initialization import check_dof
    check_dof(m, fail_flag=True)

    solver = get_solver()

    results = solver.solve(m)
    assert_optimal_termination(results)

    print(m.fs.fo.strong_draw_props[0].flow_vol_phase["Liq"].value)

