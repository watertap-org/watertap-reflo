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

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
    Var,
)
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import (
    Mixer,
    Separator,
)
from idaes.core.util.scaling import (
    set_scaling_factor,
    get_scaling_factor,
    calculate_scaling_factors,
)
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

from watertap_contrib.reflo.unit_models.surrogate import (
    VAGMDSurrogateBase,
    LTMEDSurrogate,
)

solver = get_solver()


def build_med_md_flowsheet(
    m=None,
    phase="processing",
    med_feed_salinity=30,
    med_feed_temp=25,
    med_steam_temp=80,
    med_capacity=1,
    med_recovery_ratio=0.5,
    batch_volume=50,
    md_feed_flow_rate=600,
    md_evap_inlet_temp=80,
    md_cond_inlet_temp=25,
    md_module_type="AS26C7.2L",
    md_cooling_system_type="closed",
    md_cooling_inlet_temp=25,
    md_high_brine_salinity=False,
    dt=None,
):
    if m is None:
        m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.liquid_prop = SeawaterParameterBlock()
    m.fs.vapor_prop = WaterParameterBlock()
    med_inputs = {
        "feed_salinity": med_feed_salinity,
        "feed_temperature": med_feed_temp,
        "steam_temperature": med_steam_temp,
        "sys_capacity": med_capacity,
        "recovery_ratio": med_recovery_ratio,
    }
    add_med(m.fs, med_inputs)  # add m.fs.med compnonent

    # Add a separator that splits MED brine into overflow and mixer input
    m.fs.S1 = Separator(
        property_package=m.fs.liquid_prop,
        mixed_state_block=m.fs.med.brine_props,
        outlet_list=["separator_to_mixer", "overflow"],
    )

    # Add a mixer that takes in brine from MED and MD
    m.fs.M1 = Mixer(
        property_package=m.fs.liquid_prop,
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=1,
        inlet_list=["MED_brine", "remained_liquid", "MD_brine"],
    )

    vagmd_inputs = {
        "feed_flow_rate": md_feed_flow_rate,
        "evap_inlet_temp": md_evap_inlet_temp,
        "cond_inlet_temp": md_cond_inlet_temp,
        "feed_temp": 25,
        "feed_salinity": med_feed_salinity / med_recovery_ratio,
        "md_high_brine_salinity": md_high_brine_salinity,
        "module_type": md_module_type,
        "cooling_system_type": md_cooling_system_type,
        "cooling_inlet_temp": md_cooling_inlet_temp,
    }

    add_vagmd(m.fs, vagmd_inputs)

    # Add a separator the splits the mixed flow into MD feed and remained liquid
    m.fs.S2 = Separator(
        property_package=m.fs.liquid_prop,
        mixed_state_block=m.fs.M1.mixed_state,
        outlet_list=["remained_liquid", "MD_feed"],
    )

    # Add connections
    m.fs.separator_to_mixer = Arc(
        source=m.fs.S1.separator_to_mixer, destination=m.fs.M1.MED_brine
    )

    m.fs.vagmd_to_mixer = Arc(source=m.fs.vagmd.brine, destination=m.fs.M1.MD_brine)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Add system constraints
    add_processing_phase_constraint(m.fs, phase, batch_volume, dt)

    # Initiate the properties that needs to be calculated
    property_initial_value(m.fs)

    return m


def property_initial_value(mfs):
    """
    Touch the properties that need to be calculated
    """
    mfs.S1.overflow_state[0].conc_mass_phase_comp
    mfs.S1.overflow_state[0].flow_vol_phase
    mfs.S1.separator_to_mixer_state[0].conc_mass_phase_comp
    mfs.S1.separator_to_mixer_state[0].flow_vol_phase
    mfs.M1.mixed_state[0].conc_mass_phase_comp
    mfs.M1.mixed_state[0].flow_vol_phase
    mfs.M1.MED_brine_state[0].flow_vol_phase
    mfs.M1.MD_brine_state[0].flow_vol_phase
    mfs.S2.remained_liquid_state[0].conc_mass_phase_comp
    mfs.S2.remained_liquid_state[0].flow_vol_phase
    mfs.M1.remained_liquid_state[0].conc_mass_phase_comp
    mfs.M1.remained_liquid_state[0].flow_vol_phase


def add_processing_phase_constraint(
    mfs,
    phase,
    batch_volume,
    dt=None,
):
    # Add status from the previous step
    mfs.volume_in_tank = Var(
        initialize=1,
        bounds=(0, None),
        units=pyunits.m**3,
        doc="Overflow volume in tank",
    )

    mfs.volume_in_tank_previous = Var(
        initialize=1,
        bounds=(0, None),
        units=pyunits.m**3,
        doc="Overflow volume in tank from previous time step",
    )

    # Define time interval
    mfs.dt = Var(
        initialize=30,
        bounds=(0, None),
        units=pyunits.s,
        doc="Time step of the simulation (s)",
    )

    @mfs.Constraint(doc="Calculate time interval of each period")
    def eq_dt(b):
        if dt is not None:
            return b.dt == dt

        elif mfs.vagmd.config.module_type == "AS7C1.5L":
            return b.dt == 20352.55 / pyunits.convert(
                mfs.vagmd.feed_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.L / pyunits.h,
            )
        else:  # module_type == "AS26C7.2L"
            return b.dt == 73269.19 / pyunits.convert(
                mfs.vagmd.feed_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.L / pyunits.h,
            )

    @mfs.Constraint(doc="remained volume is the batch volume")
    def eq_S2_remained_volume(b):
        return b.S2.remained_liquid_state[0].flow_vol_phase[
            "Liq"
        ] * b.dt == pyunits.convert(batch_volume * pyunits.L, to_units=pyunits.m**3)

    @mfs.Constraint(doc="Split mixed flow to MD feed")
    def eq_S2_to_vagmd(b):
        return (
            b.S2.MD_feed_state[0].flow_vol_phase["Liq"]
            == b.vagmd.feed_props[0].flow_vol_phase["Liq"]
        )

    @mfs.Constraint(doc="Current volume in the mixer")
    def eq_tank_volume(b):
        return b.volume_in_tank == b.volume_in_tank_previous + pyunits.convert(
            b.S1.overflow_state[0].flow_vol_phase["Liq"] * b.dt,
            to_units=pyunits.m**3,
        )

    @mfs.Expression(
        doc="Calculate specific thermal energy consumption in VAGMD (kWh/m3)"
    )
    def vagmd_specific_energy_consumption_thermal(b):
        return mfs.vagmd.thermal_power / pyunits.convert(
            mfs.vagmd.permeate_flux * mfs.vagmd.module_area,
            to_units=pyunits.m**3 / pyunits.h,
        )

    @mfs.Expression(
        doc="Calculate specific electric energy consumption in VAGMD (kWh/m3)"
    )
    def vagmd_specific_energy_consumption_electric(b):
        return (
            mfs.vagmd.feed_pump_power_elec + mfs.vagmd.cooling_pump_power_elec
        ) / pyunits.convert(
            mfs.vagmd.permeate_flux * mfs.vagmd.module_area,
            to_units=pyunits.m**3 / pyunits.h,
        )


def add_med(fs, inputs):
    """Method to add an MED component to an exisitng flowsheet
    Args:
        fs: exisitng flowsheet
        inputs: a dictionary depicting the MED configurations
    """
    fs.med = LTMEDSurrogate(
        property_package_liquid=fs.liquid_prop,
        property_package_vapor=fs.vapor_prop,
        number_effects=12,
    )
    steam = fs.med.steam_props[0]

    # System specification
    # Input variable 1: Feed salinity (30-60 g/L = kg/m3)
    feed_salinity = inputs["feed_salinity"] * pyunits.kg / pyunits.m**3  # g/L = kg/m3

    # Input variable 2: Feed temperature (15-35 deg C)
    feed_temperature = inputs["feed_temperature"]  # degC

    # Input variable 3: Heating steam temperature (60-85 deg C)
    steam_temperature = inputs["steam_temperature"]  # degC

    # Input variable 4: System capacity (> 2000 m3/day)
    sys_capacity = inputs["sys_capacity"] * pyunits.m**3 / pyunits.day  # m3/day

    # Input variable 5: Recovery ratio (30%- 50%)
    recovery_ratio = 0.5 * pyunits.dimensionless  # dimensionless

    feed_flow = pyunits.convert(
        (sys_capacity / recovery_ratio), to_units=pyunits.m**3 / pyunits.s
    )

    fs.med.feed_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): feed_flow,
            ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
            ("temperature", None): feed_temperature + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )
    # Fix input steam temperature
    steam.temperature.fix(steam_temperature + 273.15)

    sf = (
        1
        / 1000
        / pyunits.convert_value(
            inputs["sys_capacity"],
            from_units=pyunits.m**3 / pyunits.day,
            to_units=pyunits.m**3 / pyunits.s,
        )
    )

    # Fix target recovery rate
    fs.med.recovery_vol_phase[0, "Liq"].fix(recovery_ratio)
    fs.liquid_prop.set_default_scaling("flow_mass_phase_comp", sf, index=("Liq", "H2O"))
    fs.liquid_prop.set_default_scaling(
        "flow_mass_phase_comp", 10 * sf, index=("Liq", "TDS")
    )
    fs.vapor_prop.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    fs.vapor_prop.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Vap", "H2O"))


def add_vagmd(fs, inputs):
    """Method to add a VAGMD component to an existing flowsheet
    Args:
        fs: existing flowsheet
        inputs: a dictionary depicting the MD configurations
    """

    # System specification (Input variables)
    feed_flow_rate = inputs["feed_flow_rate"]  # 400 - 1100 L/h
    evap_inlet_temp = inputs["evap_inlet_temp"]  # 60 - 80 deg C
    cond_inlet_temp = inputs["cond_inlet_temp"]  # 20 - 30 deg C
    feed_temp = inputs["feed_temp"]  # 20 - 30 deg C
    feed_salinity = inputs["feed_salinity"]  # 35 - 292 g/L
    md_high_brine_salinity = inputs["md_high_brine_salinity"]  # -
    module_type = inputs["module_type"]
    cooling_system_type = inputs["cooling_system_type"]
    cooling_inlet_temp = inputs["cooling_inlet_temp"]
    # deg C, not required when cooling system type is "closed"

    # Identify if the final brine salinity is larger than 175.3 g/L for module "AS7C1.5L"
    # If yes, then operational parameters need to be fixed at a certain value,
    # and coolying circuit is closed to maintain condenser inlet temperature constant
    if module_type == "AS7C1.5L" and md_high_brine_salinity:
        cooling_system_type = "closed"
        feed_flow_rate = 1100  # L/h
        evap_inlet_temp = 80  # deg C
        cond_inlet_temp = 25  # deg C

    fs.vagmd = VAGMDSurrogateBase(
        property_package_seawater=fs.liquid_prop,
        property_package_water=fs.vapor_prop,
        module_type=module_type,
        high_brine_salinity=md_high_brine_salinity,
        cooling_system_type=cooling_system_type,
    )

    # Specify evaporator inlet temperature
    fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

    # Identify cooling system type
    # Closed circuit, in which TCI is forced to be constant and the cooling water temperature can be adjusted.
    if cooling_system_type == "closed":
        fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
    # Open circuit, in which cooling is available at a constant water temperature and condenser inlet temperature varies.
    else:  # "open"
        fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)


def fix_dof_and_initialize(
    m,
    outlvl=idaeslog.WARNING,
):
    """Fix degrees of freedom and initialize the flowsheet
    This function fixes the degrees of freedom of each unit and initializes the entire flowsheet.
    Args:
        m: Pyomo `Block` or `ConcreteModel` containing the flowsheet
        outlvl: Logger (default: idaeslog.WARNING)
    """

    calculate_scaling_factors(m)
    # Add scaling factors for variable and constraints constructed in the flowsheet
    if get_scaling_factor(m.fs.volume_in_tank) is None:
        set_scaling_factor(m.fs.volume_in_tank, 1e2)
        set_scaling_factor(m.fs.volume_in_tank_previous, 1e2)

    if get_scaling_factor(m.fs.dt) is None:
        set_scaling_factor(m.fs.dt, 1e-1)

    m.fs.med.initialize()
    m.fs.S1.initialize()
    # Specify feed flow state properties
    m.fs.vagmd.feed_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): pyunits.convert(
                600 * pyunits.L / pyunits.h,
                to_units=pyunits.m**3 / pyunits.s,
            ),
            ("conc_mass_phase_comp", ("Liq", "TDS")): 70,
            ("temperature", None): 35 + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )
    m.fs.vagmd.initialize()

    calculate_variable_from_constraint(m.fs.dt, m.fs.eq_dt)

    m.fs.M1.remained_liquid_state.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 50 / 1000 / m.fs.dt.value,  # m3/s
            ("conc_mass_phase_comp", ("Liq", "TDS")): 70 * pyunits.kg / pyunits.m**3,
            ("temperature", None): 35.7 + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    m.fs.M1.initialize()
    m.fs.S2.initialize()

    m.fs.volume_in_tank_previous.fix(0)

    return
