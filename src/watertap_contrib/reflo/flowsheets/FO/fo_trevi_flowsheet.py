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

# Pyomo imports
from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
    value,
)
from pyomo.network import Arc

# IDAES imports
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import (
    Mixer,
    Separator,
    HeatExchanger,
    Heater,
)
from idaes.models.unit_models.heat_exchanger import delta_temperature_underwood_callback

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

# WaterTAP REFLO imports
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import (
    FODrawSolutionParameterBlock,
)

from watertap_contrib.reflo.unit_models.zero_order.forward_osmosis_zo import (
    ForwardOsmosisZO,
)


def build_fo_trevi_flowsheet(
    m=None,
    recovery_ratio=0.3,  # Assumed FO recovery ratio
    RO_recovery_ratio=0.9,  # RO recovery ratio
    NF_recovery_ratio=0.8,  # Nanofiltration recovery ratio
    dp_brine=0,  # Required pressure over brine osmotic pressure (Pa)
    heat_mixing=105,  # Heat of mixing in the membrane (MJ/m3 product)
    separation_temp=90,  # Separation temperature of the draw solution (C)
    separator_temp_loss=1,  # Temperature loss in the separator (K)
    feed_temperature=13,  # Feed water temperature (C)
    feed_vol_flow=0.022,  # Feed water volumetric flow rate (m3/s)
    feed_TDS_mass=0.035,  # TDS mass fraction of feed
    strong_draw_temp=20,  # Strong draw solution inlet temperature (C)
    strong_draw_mass=0.8,  # Strong draw solution mass fraction
    product_draw_mass=0.01,  # Mass fraction of draw in the product water
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
    mfs = m.fs

    # Add FO module
    add_fo(
        mfs,
        recovery_ratio,
        NF_recovery_ratio,
        dp_brine,
        heat_mixing,
        separation_temp,
        separator_temp_loss,
        feed_temperature,
        feed_vol_flow,
        feed_TDS_mass,
        strong_draw_temp,
        strong_draw_mass,
        product_draw_mass,
    )

    # Add heat exchangers
    mfs.HX1 = HeatExchanger(
        delta_temperature_callback=delta_temperature_underwood_callback,
        hot_side_name="product_water",
        cold_side_name="weak_draw",
        product_water={"property_package": m.fs.draw_solution_properties},
        weak_draw={"property_package": m.fs.draw_solution_properties},
    )

    mfs.HX2 = HeatExchanger(
        delta_temperature_callback=delta_temperature_underwood_callback,
        hot_side_name="reg_draw",
        cold_side_name="weak_draw",
        reg_draw={"property_package": m.fs.draw_solution_properties},
        weak_draw={"property_package": m.fs.draw_solution_properties},
    )

    @mfs.Constraint(doc="Same outlet temperature from HX1 and HX2")
    def outlet_temp_HX1(b):
        return (
            b.HX1.weak_draw_outlet.temperature[0]
            == b.HX2.weak_draw_outlet.temperature[0]
        )

    # Add a separator to represent RO and NF fed with the product water to remove remained
    # draw solution (NF) and boron (RO), if existing
    m.fs.S2 = Separator(
        property_package=m.fs.draw_solution_properties,
        mixed_state_block=None,
        outlet_list=["RO_reject", "NF_reject", "fresh_water"],
        split_basis=3,  # Component flow
    )

    @mfs.Constraint(
        mfs.draw_solution_properties.component_list,
        doc="NF permeate is sent to RO for further treatment",
    )
    def S2_RO_reject(b, j):
        permeate_coeff = {
            "H2O": NF_recovery_ratio * (1 - RO_recovery_ratio),
            "DrawSolution": 0,
        }
        return (
            b.S2.RO_reject.flow_mass_phase_comp[0, "Liq", j]
            == b.S2.inlet.flow_mass_phase_comp[0, "Liq", j] * permeate_coeff[j]
        )

    @mfs.Constraint(
        mfs.draw_solution_properties.component_list,
        doc="NF removes all draw solution and the reject is recirculated",
    )
    def S2_NF_reject(b, j):
        permeate_coeff = {"H2O": 1 - NF_recovery_ratio, "DrawSolution": 1}
        return (
            b.S2.NF_reject.flow_mass_phase_comp[0, "Liq", j]
            == b.S2.inlet.flow_mass_phase_comp[0, "Liq", j] * permeate_coeff[j]
        )

    # Add mixer to mix NF reject that contains draw solution and weak draw
    mfs.M1 = Mixer(
        property_package=m.fs.draw_solution_properties,
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=1,
        inlet_list=["NF_reject", "weak_draw"],
    )

    # Add a separator to diverge the weak draw solution from FO module for heat recovery
    mfs.S1 = Separator(
        property_package=m.fs.draw_solution_properties,
        mixed_state_block=m.fs.M1.mixed_state,
        outlet_list=["to_HX1", "to_HX2"],
        split_basis=1,  # Total flow
    )

    # Add a heater to heat up weak draw solution from HX1 and HX2
    mfs.H1 = Heater(
        property_package=m.fs.draw_solution_properties,
    )
    mfs.H2 = Heater(
        property_package=m.fs.draw_solution_properties,
    )

    @mfs.Constraint(doc="Temperature at the separator")
    def H1_out_temp(b):
        return b.H1.outlet.temperature[0] == b.fo.regeneration_temp

    @mfs.Constraint(doc="Temperature at the separator")
    def H2_out_temp(b):
        return b.H2.outlet.temperature[0] == b.fo.regeneration_temp

    # Add a cooler that cools regenerated draw solution using supplemental source water
    mfs.Cooler = Heater(
        property_package=m.fs.draw_solution_properties,
    )

    @mfs.Constraint(doc="Strong draw solution temperature entering FO module")
    def Cooler_out_temp(b):
        return b.Cooler.outlet.temperature[0] == b.fo.feed_props[0].temperature

    # Add connections
    mfs.HX1_cold_inlet = Arc(
        source=m.fs.S1.to_HX1, destination=m.fs.HX1.cold_side_inlet
    )
    mfs.HX1_hot_inlet = Arc(source=m.fs.fo.product, destination=m.fs.HX1.hot_side_inlet)
    mfs.HX2_cold_inlet = Arc(
        source=m.fs.S1.to_HX2, destination=m.fs.HX2.cold_side_inlet
    )
    mfs.HX2_hot_inlet = Arc(
        source=m.fs.fo.reg_draw, destination=m.fs.HX2.hot_side_inlet
    )
    mfs.S2_inlet = Arc(source=m.fs.HX1.hot_side_outlet, destination=m.fs.S2.inlet)
    mfs.NF_reject_to_M1 = Arc(source=m.fs.S2.NF_reject, destination=m.fs.M1.NF_reject)
    mfs.weak_draw_to_M1 = Arc(source=m.fs.fo.weak_draw, destination=m.fs.M1.weak_draw)
    mfs.HX1_to_H1 = Arc(source=m.fs.HX1.cold_side_outlet, destination=m.fs.H1.inlet)
    mfs.HX2_to_H2 = Arc(source=m.fs.HX2.cold_side_outlet, destination=m.fs.H2.inlet)
    mfs.HX2_to_Cooler = Arc(
        source=m.fs.HX2.hot_side_outlet, destination=m.fs.Cooler.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Calculate system capacity (m3/day)
    fo = mfs.fo

    @mfs.Expression(doc="Calculate system capacity (m3/day)")
    def system_capacity(b):
        return pyunits.convert(
            mfs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "H2O"]
            / (1000 * pyunits.kg / pyunits.m**3),  # Fresh water density
            to_units=pyunits.m**3 / pyunits.day,
        )

    # Calculate specific thermal energy consumption of the system (k)
    @mfs.Expression(doc="Calculate specific thermal energy consumption (kWh/m3)")
    def specific_energy_consumption_thermal(b):
        return pyunits.convert(
            (mfs.H1.heat_duty[0] + mfs.H2.heat_duty[0]) / mfs.system_capacity,
            to_units=pyunits.kWh / pyunits.m**3,
        )

    # Specify heat exchanger parameters
    mfs.HX1.area.fix(50000)
    mfs.HX1.overall_heat_transfer_coefficient[0].fix(100)
    mfs.HX2.area.fix(50000)
    mfs.HX2.overall_heat_transfer_coefficient[0].fix(100)

    iscale.set_scaling_factor(mfs.HX1.area, 1e-5)
    iscale.set_scaling_factor(mfs.HX1.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(mfs.HX1.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX1.cold_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX2.area, 1e-5)
    iscale.set_scaling_factor(mfs.HX2.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(mfs.HX2.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(mfs.HX2.cold_side.heat, 1e-7)

    return m


def add_fo(
    fs,
    recovery_ratio=0.3,  # Assumed FO recovery ratio
    NF_recovery_ratio=0.8,  # Nanofiltration recovery ratio
    dp_brine=0,  # Required pressure over brine osmotic pressure (Pa)
    heat_mixing=105,  # Heat of mixing in the membrane (MJ/m3 product)
    separation_temp=90,  # Separation temperature of the draw solution (C)
    separator_temp_loss=1,  # Temperature loss in the separator (K)
    feed_temperature=13,  # Feed water temperature (C)
    feed_vol_flow=3.704,  # Feed water volumetric flow rate (m3/s)
    feed_TDS_mass=0.035,  # TDS mass fraction of feed
    strong_draw_temp=20,  # Strong draw solution inlet temperature (C)
    strong_draw_mass=0.8,  # Strong draw solution mass fraction
    product_draw_mass=0.01,  # Mass fraction of draw in the product water
):
    """
    This function adds a FO module to the current flowsheet

    Returns:
        object: A FO zero-order model
    """
    fs.fo = ForwardOsmosisZO(
        property_package_water=fs.seawater_properties,
        property_package_draw_solution=fs.draw_solution_properties,
    )

    fo = fs.fo

    fo.recovery_ratio.fix(recovery_ratio)
    fo.nanofiltration_recovery_ratio.fix(NF_recovery_ratio)
    fo.dp_brine.fix(dp_brine)
    fo.heat_mixing.fix(heat_mixing)
    fo.regeneration_temp.fix(separation_temp + 273.15)
    fo.separator_temp_loss.fix(separator_temp_loss)

    # Specify strong draw solution properties
    fo.strong_draw_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): feed_vol_flow,
            ("mass_frac_phase_comp", ("Liq", "DrawSolution")): strong_draw_mass,
            ("temperature", None): strong_draw_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    fo.strong_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].unfix()

    # Specify product water properties
    fo.product_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): feed_vol_flow,
            ("mass_frac_phase_comp", ("Liq", "DrawSolution")): product_draw_mass,
            ("temperature", None): separation_temp - separator_temp_loss + 273.15,
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

    fo.weak_draw_props[0].pressure.fix(101325)
    fo.reg_draw_props[0].pressure.fix(101325)

    # Set scaling factors for mass flow rates
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / feed_vol_flow / 100, index=("Liq", "H2O")
    )
    fs.seawater_properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / feed_vol_flow / 10, index=("Liq", "TDS")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / feed_vol_flow / 100, index=("Liq", "H2O")
    )
    fs.draw_solution_properties.set_default_scaling(
        "flow_mass_phase_comp", 1 / feed_vol_flow / 10, index=("Liq", "DrawSolution")
    )


def fix_dof_and_initialize(
    m,
    strong_draw_mass_frac=0.8,
    product_draw_mass_frac=0.01,
    NF_recovery_ratio=0.8,
    RO_recovery_ratio=0.9,
):

    iscale.calculate_scaling_factors(m)

    m.fs.fo.initialize()
    # Unfix the state variables and fix mass fraction of two state blocks
    m.fs.fo.unfix_and_fix_freedom(strong_draw_mass_frac, product_draw_mass_frac)

    m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    )
    m.fs.S1.to_HX1.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value / 2
    )
    m.fs.S1.to_HX1.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value / 2
    )
    m.fs.S1.to_HX2.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value / 2
    )
    m.fs.S1.to_HX2.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.S1.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value / 2
    )
    m.fs.S1.to_HX1_state[0].flow_vol_phase["Liq"]
    m.fs.S1.to_HX2_state[0].flow_vol_phase["Liq"]
    m.fs.S1.to_HX1_state[0].cp_mass_phase["Liq"]
    m.fs.S1.initialize()

    state_args_HX1_hot = {
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): m.fs.fo.product_props[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            ("Liq", "DrawSolution"): m.fs.fo.product_props[0]
            .flow_mass_phase_comp["Liq", "DrawSolution"]
            .value,
        },
        "temperature": m.fs.fo.product_props[0].temperature.value,
        "pressure": m.fs.fo.product_props[0].pressure.value,
    }
    state_args_HX1_cold = {
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): m.fs.S1.to_HX1.flow_mass_phase_comp[0, "Liq", "H2O"].value,
            ("Liq", "DrawSolution"): m.fs.S1.to_HX1.flow_mass_phase_comp[
                0, "Liq", "DrawSolution"
            ].value,
        },
        "temperature": m.fs.S1.to_HX1.temperature[0].value,
        "pressure": m.fs.S1.to_HX1.pressure[0].value,
    }
    m.fs.HX1.initialize(
        state_args_1=state_args_HX1_hot, state_args_2=state_args_HX1_cold
    )
    # Cold side has liquid separation
    m.fs.HX1.cold_side.properties_out[0].liquid_separation.fix(1)
    m.fs.HX1.cold_side.properties_out[0].mass_frac_after_separation = (
        strong_draw_mass_frac
    )
    iscale.set_scaling_factor(
        m.fs.HX1.cold_side.properties_out[0].heat_separation_phase,
        1e-5,
    )

    state_args_HX2_hot = {
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): m.fs.fo.reg_draw_props[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .value,
            ("Liq", "DrawSolution"): m.fs.fo.reg_draw_props[0]
            .flow_mass_phase_comp["Liq", "DrawSolution"]
            .value,
        },
        "temperature": m.fs.fo.reg_draw_props[0].temperature.value,
        "pressure": m.fs.fo.reg_draw_props[0].pressure.value,
    }
    state_args_HX2_cold = {
        "flow_mass_phase_comp": {
            ("Liq", "H2O"): m.fs.S1.to_HX2.flow_mass_phase_comp[0, "Liq", "H2O"].value,
            ("Liq", "DrawSolution"): m.fs.S1.to_HX2.flow_mass_phase_comp[
                0, "Liq", "DrawSolution"
            ].value,
        },
        "temperature": m.fs.S1.to_HX2.temperature[0].value,
        "pressure": m.fs.S1.to_HX2.pressure[0].value,
    }
    m.fs.HX2.initialize(
        state_args_1=state_args_HX2_hot, state_args_2=state_args_HX2_cold
    )
    # Cold side has liquid separation
    m.fs.HX2.cold_side.properties_out[0].liquid_separation.fix(1)
    m.fs.HX2.cold_side.properties_out[0].mass_frac_after_separation = (
        strong_draw_mass_frac
    )

    # Initialize separator S2
    m.fs.S2.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    m.fs.S2.inlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    )
    m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
        * (1 - NF_recovery_ratio)
    )
    m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
        * (1 - NF_recovery_ratio)
    )
    m.fs.S2.RO_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
        * NF_recovery_ratio
        * (1 - RO_recovery_ratio)
    )
    m.fs.S2.RO_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
        * NF_recovery_ratio
        * (1 - RO_recovery_ratio)
    )
    m.fs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "H2O"].value
        * NF_recovery_ratio
        * RO_recovery_ratio
    )
    m.fs.S2.fresh_water.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
        * NF_recovery_ratio
        * RO_recovery_ratio
    )
    m.fs.S2.initialize()

    # Initialize mixer M1
    m.fs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    m.fs.M1.weak_draw.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
    )
    m.fs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value
    )
    m.fs.M1.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value
    )
    m.fs.M1.outlet.flow_mass_phase_comp[0, "Liq", "H2O"].value = (
        m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].value
        + m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "H2O"].value
    )
    m.fs.M1.outlet.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value = (
        m.fs.fo.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].value
        + m.fs.S2.NF_reject.flow_mass_phase_comp[0, "Liq", "DrawSolution"].value
    )

    m.fs.M1.mixed_state[0].flow_vol_phase["Liq"]
    m.fs.M1.initialize()

    # Initialize heaters
    m.fs.H1.initialize()
    m.fs.H2.initialize()
    m.fs.Cooler.initialize()


def get_flowsheet_performance(m):
    overall_performance = {
        "Production capacity (m3/day)": value(m.fs.system_capacity),
        "Specific thermal energy consumption (kWh/m3)": value(
            m.fs.specific_energy_consumption_thermal
        ),
        "Thermal power requirement (kW)": value(
            m.fs.specific_energy_consumption_thermal
        )
        * value(m.fs.system_capacity)
        / 24,
        "LCOW ($/m3)": value(m.fs.costing.LCOW),
    }

    operational_parameters = {
        # Heat exchanger HX1 temperatures
        "HX1 cold in temp": m.fs.HX1.weak_draw_inlet.temperature[0].value - 273.15,
        "HX1 cold out temp": m.fs.HX1.weak_draw_outlet.temperature[0].value - 273.15,
        "HX1 hot in temp": m.fs.HX1.product_water_inlet.temperature[0].value - 273.15,
        "HX1 hot out temp": m.fs.HX1.product_water_outlet.temperature[0].value - 273.15,
        # Heat exchanger HX2 temperatures
        "HX2 cold in temp": m.fs.HX2.weak_draw_inlet.temperature[0].value - 273.15,
        "HX2 cold out temp": m.fs.HX2.weak_draw_outlet.temperature[0].value - 273.15,
        "HX2 hot in temp": m.fs.HX2.reg_draw_inlet.temperature[0].value - 273.15,
        "HX2 hot out temp": m.fs.HX2.reg_draw_outlet.temperature[0].value - 273.15,
        # Product water flow rates
        "FO product vol (m3/s)": m.fs.fo.product_props[0].flow_vol_phase["Liq"].value,
        "FO product H2O mass (kg/s)": m.fs.fo.product_props[0]
        .flow_mass_phase_comp["Liq", "H2O"]
        .value,
        "FO product draw mass (kg/s)": m.fs.fo.product_props[0]
        .flow_mass_phase_comp["Liq", "DrawSolution"]
        .value,
        # Regenerated draw solution flow rates
        "FO reg draw vol (m3/s)": m.fs.fo.reg_draw_props[0].flow_vol_phase["Liq"].value,
        "FO reg H2O mass (kg/s)": m.fs.fo.reg_draw_props[0]
        .flow_mass_phase_comp["Liq", "H2O"]
        .value,
        "FO reg draw mass (kg/s)": m.fs.fo.reg_draw_props[0]
        .flow_mass_phase_comp["Liq", "DrawSolution"]
        .value,
        # Weak draw solution properties at FO module outlet
        "FO weak H2O mass(kg/s)": m.fs.fo.weak_draw.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "FO weak draw mass(kg/s)": m.fs.fo.weak_draw.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        "FO weak temp": m.fs.fo.weak_draw.temperature[0].value - 273.15,
        # Weak draw solution properties after mixing with NF reject
        "Mixed weak H2O mass(kg/s)": m.fs.S1.inlet.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "Mixed weak draw mass (kg/s)": m.fs.S1.inlet.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        "Mixed weak temp": m.fs.S1.inlet.temperature[0].value - 273.15,
        # Flows going into HX1
        "HX1 cold inlet H2O(kg/s)": m.fs.HX1.cold_side_inlet.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "HX1 cold inlet Draw(kg/s)": m.fs.HX1.cold_side_inlet.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        "HX1 hot inlet H2O(kg/s)": m.fs.HX1.hot_side_inlet.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "HX1 hot inlet Draw(kg/s)": m.fs.HX1.hot_side_inlet.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        # Flows going into HX2
        "HX2 cold in H2O(kg/s)": m.fs.HX2.cold_side_inlet.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "HX2 cold in Draw(kg/s)": m.fs.HX2.cold_side_inlet.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        "HX2 hot in H2O(kg/s)": m.fs.HX2.hot_side_inlet.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "HX2 hot in Draw(kg/s)": m.fs.HX2.hot_side_inlet.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        # Energy balances in heat exchangers
        "HX1 hot side heat load (MJ)": value(m.fs.HX1.hot_side.heat[0]) / 1e6,
        "HX1 cold side heat load (MJ)": value(m.fs.HX1.cold_side.heat[0]) / 1e6,
        "HX2 hot side heat load (MJ)": value(m.fs.HX2.hot_side.heat[0]) / 1e6,
        "HX2 cold side heat load (MJ)": value(m.fs.HX2.cold_side.heat[0]) / 1e6,
        "HX1 approach temp": value(m.fs.HX1.delta_temperature[0]),
        "HX2 approach temp": value(m.fs.HX2.delta_temperature[0]),
        # Separated flow rates after S2
        "Fresh water H2O mass (kg/s)": m.fs.S2.fresh_water.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "RO reject H2O mass (kg/s)": m.fs.S2.RO_reject.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "NF reject H2O mass (kg/s)": m.fs.S2.NF_reject.flow_mass_phase_comp[
            0, "Liq", "H2O"
        ].value,
        "NF reject draw mass (kg/s)": m.fs.S2.NF_reject.flow_mass_phase_comp[
            0, "Liq", "DrawSolution"
        ].value,
        # Mixed flow properties after M1
        "M1 temperature": m.fs.M1.mixed_state[0].temperature.value - 273.15,
        "M1 H2O mass (kg/s)": m.fs.M1.mixed_state[0]
        .flow_mass_phase_comp["Liq", "H2O"]
        .value,
        "M1 draw mass (kg/s)": m.fs.M1.mixed_state[0]
        .flow_mass_phase_comp["Liq", "DrawSolution"]
        .value,
        # Heat loads
        "Heater 1 heat load": m.fs.H1.heat_duty[0].value,
        "Heater 2 heat load": m.fs.H2.heat_duty[0].value,
        "Cooler heat load": m.fs.Cooler.heat_duty[0].value,
    }

    return overall_performance, operational_parameters
