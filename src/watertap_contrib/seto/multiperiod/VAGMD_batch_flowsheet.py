# Pyomo imports
from pyomo.environ import (
    ConcreteModel,
    Var,
    Objective,
    Param,
    Expression,
    Constraint,
    Block,
    log10,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)

# WaterTAP imports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.unit_models.surrogate import VAGMDSurrogate


def build_vagmd_flowsheet(
    m=None,
    feed_flow_rate=600,
    evap_inlet_temp=80,
    cond_inlet_temp=25,
    feed_temp=25,
    feed_salinity=35,
    initial_batch_volume=50,
    module_type="AS7C1.5L",
    high_brine_salinity=False,
    cooling_system_type="closed",
):
    """
    This function builds a unit model for a certain time period

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.vagmd = VAGMDSurrogate(
        property_package=m.fs.properties,
        module_type=module_type,
        high_brine_salinity=high_brine_salinity,
        cooling_system_type=cooling_system_type,
    )

    # Fix the model inputs
    m.fs.vagmd.feed_props[0].flow_vol_phase["Liq"].fix(
        pyunits.convert(
            feed_flow_rate * pyunits.L / pyunits.h, to_units=pyunits.m**3 / pyunits.s
        )
    )
    m.fs.vagmd.evaporator_in_props[0].temperature.fix(evap_inlet_temp + 273.15)

    if cooling_system_type == "closed":  # TODO: update closed cooling
        m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)
    else:  # "open"
        m.fs.vagmd.condenser_in_props[0].temperature.fix(cond_inlet_temp + 273.15)

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

    m.fs.specific_energy_consumption_thermal = Var(
        initialize=100,
        units=pyunits.kWh / pyunits.m**3,
        doc="Specific thermal power consumption (kWh/m3)",
    )

    """
    Constriant equatiosn
    """
    mfs = m.fs

    @mfs.Constraint(doc="Calculate time interval of each period")
    def eq_dt(b):
        if module_type == "AS7C1.5L":
            return b.vagmd.dt == 20352.55 / feed_flow_rate
        else:  # module_type == "AS26C7.2L"
            return b.vagmd.dt == 73269.19 / feed_flow_rate

    @mfs.Constraint(doc="Calculate current feed salinity")
    def eq_feed_salinity(b):
        return b.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"] * (
            initial_batch_volume - b.acc_distillate_volume
        ) == b.pre_feed_salinity * (initial_batch_volume - b.pre_acc_distillate_volume)

    @mfs.Constraint(doc="Calculate current feed temperature")
    def eq_feed_temp(b):
        return b.vagmd.feed_props[0].temperature == (
            feed_flow_rate
            * pyunits.convert(b.vagmd.dt, to_units=pyunits.h)
            * b.pre_evap_out_temp
            + (initial_batch_volume - b.pre_acc_distillate_volume)
            * b.pre_feed_temperature
        ) / (
            feed_flow_rate * pyunits.convert(b.vagmd.dt, to_units=pyunits.h)
            + initial_batch_volume
            - b.pre_acc_distillate_volume
        )

    @mfs.Constraint(doc="Calculate accumulated distillate volume")
    def eq_acc_distillate_volume(b):
        return b.acc_distillate_volume == b.pre_acc_distillate_volume + pyunits.convert(
            b.pre_permeate_flow_rate * b.vagmd.dt, to_units=pyunits.L
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
        return b.acc_thermal_energy == b.pre_acc_thermal_energy + b.vagmd.thermal_energy

    # @m.Constraint(doc="Calculate specific thermal energy consumption")
    # def eq_specific_thermal_energy_consumption(b):
    #     return (b.fs.specific_energy_consumption_thermal == pyunits.convert(b.fs.acc_thermal_energy / (b.fs.acc_distillate_volume+1e-8*pyunits.L),
    #                                                                         to_units= pyunits.kWh / pyunits.m**3))

    return m


def fix_dof_and_initialize(
    m,
    outlvl=idaeslog.WARNING,
):

    m.fs.vagmd.calculate_scaling_factors()

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

    if iscale.get_scaling_factor(m.fs.acc_distillate_volume) is None:
        iscale.set_scaling_factor(m.fs.acc_distillate_volume, 1e-1)

    if iscale.get_scaling_factor(m.fs.acc_recovery_ratio) is None:
        iscale.set_scaling_factor(m.fs.acc_recovery_ratio, 1e2)

    if iscale.get_scaling_factor(m.fs.acc_thermal_energy) is None:
        iscale.set_scaling_factor(m.fs.acc_thermal_energy, 1e0)

    if iscale.get_scaling_factor(m.fs.specific_energy_consumption_thermal) is None:
        iscale.set_scaling_factor(m.fs.specific_energy_consumption_thermal, 1e-2)

    # Transforming constraint
    sf = iscale.get_scaling_factor(m.fs.vagmd.dt)
    iscale.constraint_scaling_transform(m.fs.eq_dt, sf)

    # sf = (iscale.get_scaling_factor(m.fs.vagmd.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]) *
    #     iscale.get_scaling_factor(m.fs.acc_distillate_volume))
    iscale.constraint_scaling_transform(m.fs.eq_feed_salinity, 1e4)

    # sf = iscale.get_scaling_factor(m.fs.vagmd.feed_props[0].temperature )
    iscale.constraint_scaling_transform(m.fs.eq_feed_temp, 1e-2)

    # sf = iscale.get_scaling_factor(m.fs.eq_acc_distillate_volume )
    iscale.constraint_scaling_transform(m.fs.eq_acc_distillate_volume, 1e-1)

    # sf = iscale.get_scaling_factor(m.fs.acc_recovery_ratio )
    iscale.constraint_scaling_transform(m.fs.eq_acc_recovery_ratio, 1e2)

    # sf = iscale.get_scaling_factor(m.fs.acc_thermal_energy )
    iscale.constraint_scaling_transform(m.fs.eq_acc_thermal_energy, 1e0)
