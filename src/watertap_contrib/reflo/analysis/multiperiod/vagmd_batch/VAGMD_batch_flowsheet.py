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
    Var,
    units as pyunits,
)

# IDAES imports
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap_contrib.reflo.unit_models.surrogate import VAGMDSurrogate


def build_vagmd_flowsheet(
    m=None,
    dt=None,
    system_capacity=2000,  # m3/day
    feed_flow_rate=600,
    evap_inlet_temp=80,
    cond_inlet_temp=25,
    feed_temp=25,
    feed_salinity=35,
    recovery_ratio=0.5,
    initial_batch_volume=50,
    module_type="AS7C1.5L",
    cooling_system_type="closed",
    high_brine_salinity=False,  # True if brine salinity > 175.3 g/L
    cooling_inlet_temp=25,  # Not required when cooling system type is "closed"
):
    """
    This function builds a unit model for a certain time period

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.seawater_properties = SeawaterParameterBlock()
    m.fs.water_properties = WaterParameterBlock()

    m.fs.vagmd = VAGMDSurrogate(
        property_package_seawater=m.fs.seawater_properties,
        property_package_water=m.fs.water_properties,
        module_type=module_type,
        high_brine_salinity=high_brine_salinity,
        cooling_system_type=cooling_system_type,
    )

    # Run helper function to determine the salinity mode
    (
        feed_flow_rate,
        evap_inlet_temp,
        cond_inlet_temp,
        cooling_system_type,
    ) = m.fs.vagmd._determine_salinity_mode(
        feed_flow_rate,
        evap_inlet_temp,
        cond_inlet_temp,
        module_type,
        high_brine_salinity,
        cooling_system_type,
    )

    # Specify system capacity
    m.fs.vagmd.system_capacity.fix(system_capacity)
    # Specify feed flow state properties
    m.fs.vagmd.feed_props.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): pyunits.convert(
                feed_flow_rate * pyunits.L / pyunits.h,
                to_units=pyunits.m**3 / pyunits.s,
            ),
            ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
            ("temperature", None): feed_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    # Fix the model inputs
    m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

    if cooling_system_type == "closed":
        m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
    else:  # "open"
        m.fs.vagmd.cooling_in_props[0].temperature.fix(cooling_inlet_temp + 273.15)

    # Define time interval
    m.fs.dt = Var(
        initialize=30,
        bounds=(0, None),
        units=pyunits.s,
        doc="Time step of the simulation (s)",
    )

    """
    Variables inherited from previous period
    """
    m.fs.pre_feed_salinity = Var(
        initialize=35,
        bounds=(0, 292),
        units=pyunits.g / pyunits.L,
        doc="Feed salinity from previous time step",
    )
    m.fs.pre_feed_temperature = Var(
        initialize=298.15,
        bounds=(273.15, 373.15),
        units=pyunits.K,
        doc="Feed temperature from previous time step",
    )
    m.fs.pre_evap_out_temp = Var(
        initialize=308.15,
        bounds=(273.15, 373.15),
        units=pyunits.K,
        doc="Feed temperature from previous time step",
    )
    m.fs.pre_permeate_flow_rate = Var(
        initialize=1e-5,
        units=pyunits.m**3 / pyunits.s,
        doc="Permeate flow rate from previous time step",
    )
    m.fs.pre_acc_distillate_volume = Var(
        initialize=0,
        units=pyunits.L,
        doc="Accumulated volume of distillate from previous time step",
    )

    m.fs.pre_acc_thermal_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated thermal energy consumption from previous time step",
    )

    m.fs.pre_acc_electric_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated electric energy consumption from previous time step",
    )

    m.fs.pre_acc_cooling_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated cooling energy consumption from previous time step",
    )

    m.fs.pre_thermal_power = Var(
        initialize=10,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Thermal power consumption from previous time step",
    )

    m.fs.pre_cooling_power = Var(
        initialize=10,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Cooling power consumption from previous time step",
    )

    m.fs.pre_cooling_pump_power_elec = Var(
        initialize=10,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Electric power for pumping cooling water from previous time step",
    )

    m.fs.pre_feed_pump_power_elec = Var(
        initialize=10,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Electric power for pumping feed water from previous time step",
    )

    """
    Variables aggrated through periods
    """
    m.fs.acc_distillate_volume = Var(
        initialize=0,
        units=pyunits.L,
        doc="Accumulated volume of distillate",
    )

    m.fs.acc_recovery_ratio = Var(
        initialize=0,
        units=pyunits.dimensionless,
        doc="Accumulated recovery rate",
    )

    m.fs.acc_thermal_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated thermal energy consumption",
    )

    m.fs.acc_cooling_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated cooling energy consumption",
    )

    m.fs.acc_electric_energy = Var(
        initialize=0,
        units=pyunits.kWh,
        doc="Accumulated electric energy consumption",
    )

    m.fs.specific_energy_consumption_thermal = Var(
        initialize=100,
        units=pyunits.kWh / pyunits.m**3,
        doc="Specific thermal power consumption (kWh/m3)",
    )

    m.fs.specific_energy_consumption_electric = Var(
        initialize=1,
        units=pyunits.kWh / pyunits.m**3,
        doc="Specific electric power consumption (kWh/m3)",
    )

    m.fs.gain_output_ratio = Var(
        initialize=10,
        units=pyunits.dimensionless,
        doc="Gain output ratio",
    )

    """
    Constraint equations
    """
    mfs = m.fs

    @mfs.Constraint(doc="Calculate time interval of each period")
    def eq_dt(b):
        if dt:  # If dt has been specified
            return b.dt == dt
        elif module_type == "AS7C1.5L":
            return b.dt == 20352.55 / feed_flow_rate
        else:  # module_type == "AS26C7.2L"
            return b.dt == 73269.19 / feed_flow_rate

    @mfs.Constraint(doc="Calculate current feed salinity")
    def eq_feed_salinity(b):
        return b.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"] * (
            initial_batch_volume - b.acc_distillate_volume
        ) == b.pre_feed_salinity * (initial_batch_volume - b.pre_acc_distillate_volume)

    @mfs.Constraint(doc="Calculate current feed temperature")
    def eq_feed_temp(b):
        return b.vagmd.feed_props[0].temperature == (
            feed_flow_rate
            * pyunits.convert(b.dt, to_units=pyunits.h)
            * b.pre_evap_out_temp
            + (initial_batch_volume - b.pre_acc_distillate_volume)
            * b.pre_feed_temperature
        ) / (
            feed_flow_rate * pyunits.convert(b.dt, to_units=pyunits.h)
            + initial_batch_volume
            - b.pre_acc_distillate_volume
        )

    @mfs.Constraint(doc="Calculate accumulated distillate volume")
    def eq_acc_distillate_volume(b):
        return b.acc_distillate_volume == b.pre_acc_distillate_volume + pyunits.convert(
            b.pre_permeate_flow_rate * b.dt, to_units=pyunits.L
        )

    @mfs.Constraint(doc="Calculate accmulated recovery ratio")
    def eq_acc_recovery_ratio(b):
        return (
            b.acc_recovery_ratio
            == 1
            - feed_salinity / b.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )

    @mfs.Constraint(doc="Calculate accmulated thermal energy consumption")
    def eq_acc_thermal_energy(b):
        return b.acc_thermal_energy == b.pre_acc_thermal_energy + pyunits.convert(
            b.pre_thermal_power * b.dt, to_units=pyunits.kWh
        )

    @mfs.Constraint(doc="Calculate accmulated thermal energy consumption")
    def eq_acc_cooling_energy(b):
        return b.acc_cooling_energy == b.pre_acc_cooling_energy + pyunits.convert(
            b.pre_cooling_power * b.dt, to_units=pyunits.kWh
        )

    @mfs.Constraint(doc="Calculate specific thermal energy consumption")
    def eq_specific_energy_consumption_thermal(b):
        return b.specific_energy_consumption_thermal == pyunits.convert(
            b.acc_thermal_energy / (b.acc_distillate_volume + 1e-8 * pyunits.L),
            to_units=pyunits.kWh / pyunits.m**3,
        )

    @mfs.Constraint(doc="Calculate accmulated electric energy consumption")
    def eq_acc_electric_energy(b):
        return b.acc_electric_energy == b.pre_acc_electric_energy + pyunits.convert(
            (b.pre_feed_pump_power_elec + b.pre_cooling_pump_power_elec) * b.dt,
            to_units=pyunits.kWh,
        )

    @mfs.Constraint(doc="Calculate specific electric energy consumption")
    def eq_specific_energy_consumption_electric(b):
        return b.specific_energy_consumption_electric == pyunits.convert(
            b.acc_electric_energy / (b.acc_distillate_volume + 1e-8 * pyunits.L),
            to_units=pyunits.kWh / pyunits.m**3,
        )

    @mfs.Constraint(doc="Calculate gain output ratio")
    def eq_gain_output_ratio(b):
        return b.gain_output_ratio == pyunits.convert(
            b.acc_distillate_volume
            * b.vagmd.avg_condenser_props[0].dens_mass_phase["Liq"]
            * b.vagmd.avg_condenser_props[0].dh_vap_mass
            / (b.acc_thermal_energy + 1e-8 * pyunits.kWh),
            to_units=pyunits.dimensionless,
        )

    return m


def fix_dof_and_initialize(
    m,
    feed_flow_rate=600,
    feed_salinity=35,
    feed_temp=25,
    outlvl=idaeslog.WARNING,
):

    # Initialize the flowsheet to the beginning of the batch operation (t = 0)
    m.fs.pre_feed_temperature.fix(feed_temp + 273.15)
    m.fs.pre_permeate_flow_rate.fix(0)
    m.fs.acc_distillate_volume.fix(0)
    m.fs.acc_thermal_energy.fix(0)
    m.fs.acc_cooling_energy.fix(0)
    m.fs.acc_electric_energy.fix(0)
    m.fs.pre_thermal_power.fix(0)
    m.fs.pre_cooling_power.fix(0)
    m.fs.pre_feed_pump_power_elec.fix(0)
    m.fs.pre_cooling_pump_power_elec.fix(0)

    iscale.calculate_scaling_factors(m.fs.vagmd)
    m.fs.vagmd.initialize_build(outlvl=outlvl)

    if iscale.get_scaling_factor(m.fs.dt) is None:
        iscale.set_scaling_factor(m.fs.dt, 1e-1)

    if iscale.get_scaling_factor(m.fs.pre_feed_salinity) is None:
        iscale.set_scaling_factor(m.fs.pre_feed_salinity, 1e-2)

    if iscale.get_scaling_factor(m.fs.pre_feed_temperature) is None:
        iscale.set_scaling_factor(m.fs.pre_feed_temperature, 1e-2)

    if iscale.get_scaling_factor(m.fs.pre_evap_out_temp) is None:
        iscale.set_scaling_factor(m.fs.pre_evap_out_temp, 1e-2)

    if iscale.get_scaling_factor(m.fs.pre_permeate_flow_rate) is None:
        iscale.set_scaling_factor(m.fs.pre_permeate_flow_rate, 1e5)

    if iscale.get_scaling_factor(m.fs.pre_acc_distillate_volume) is None:
        iscale.set_scaling_factor(m.fs.pre_acc_distillate_volume, 1e-1)

    if iscale.get_scaling_factor(m.fs.pre_acc_thermal_energy) is None:
        iscale.set_scaling_factor(m.fs.pre_acc_thermal_energy, 1e0)

    if iscale.get_scaling_factor(m.fs.pre_acc_electric_energy) is None:
        iscale.set_scaling_factor(m.fs.pre_acc_electric_energy, 1e2)

    if iscale.get_scaling_factor(m.fs.pre_acc_cooling_energy) is None:
        iscale.set_scaling_factor(m.fs.pre_acc_cooling_energy, 1e0)

    if iscale.get_scaling_factor(m.fs.acc_distillate_volume) is None:
        iscale.set_scaling_factor(m.fs.acc_distillate_volume, 1e-1)

    if iscale.get_scaling_factor(m.fs.acc_recovery_ratio) is None:
        iscale.set_scaling_factor(m.fs.acc_recovery_ratio, 1e2)

    if iscale.get_scaling_factor(m.fs.acc_thermal_energy) is None:
        iscale.set_scaling_factor(m.fs.acc_thermal_energy, 1e0)

    if iscale.get_scaling_factor(m.fs.acc_electric_energy) is None:
        iscale.set_scaling_factor(m.fs.acc_electric_energy, 1e2)

    if iscale.get_scaling_factor(m.fs.acc_cooling_energy) is None:
        iscale.set_scaling_factor(m.fs.acc_cooling_energy, 1e0)

    if iscale.get_scaling_factor(m.fs.pre_thermal_power) is None:
        iscale.set_scaling_factor(m.fs.pre_thermal_power, 1e-1)

    if iscale.get_scaling_factor(m.fs.pre_cooling_power) is None:
        iscale.set_scaling_factor(m.fs.pre_cooling_power, 1e-1)

    if iscale.get_scaling_factor(m.fs.pre_feed_pump_power_elec) is None:
        iscale.set_scaling_factor(m.fs.pre_feed_pump_power_elec, 1e3)

    if iscale.get_scaling_factor(m.fs.pre_cooling_pump_power_elec) is None:
        iscale.set_scaling_factor(m.fs.pre_cooling_pump_power_elec, 1e3)

    if iscale.get_scaling_factor(m.fs.specific_energy_consumption_thermal) is None:
        iscale.set_scaling_factor(m.fs.specific_energy_consumption_thermal, 1e-2)

    if iscale.get_scaling_factor(m.fs.specific_energy_consumption_electric) is None:
        iscale.set_scaling_factor(m.fs.specific_energy_consumption_electric, 1)

    if iscale.get_scaling_factor(m.fs.gain_output_ratio) is None:
        iscale.set_scaling_factor(m.fs.gain_output_ratio, 1e-1)
