import os
import numpy as np

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Param,
    Constraint,
    Var,
    Block,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import *
from idaes.core import UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)
from watertap.core.util.infeasible import *

from watertap_contrib.seto.analysis.net_metering.util import *
from watertap_contrib.seto.costing import (
    TreatmentCosting,
    EnergyCosting,
    SETOSystemCosting,
)
from watertap_contrib.seto.solar_models.zero_order import PhotovoltaicZO
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP


solver = get_solver()

absolute_path = os.path.dirname(__file__)
print(absolute_path)

tech_config_file = "/pysam_data/pvsamv1.json"
tech_config_file = absolute_path + tech_config_file
grid_config_file = "/pysam_data/grid.json"
grid_config_file = absolute_path + grid_config_file
rate_config_file = "/pysam_data/utilityrate5.json"
rate_config_file = absolute_path + rate_config_file
cash_config_file = "/pysam_data/singleowner.json"
cash_config_file = absolute_path + cash_config_file
weather_file = "/pysam_data/phoenix_az_33.450495_-111.983688_psmv3_60_tmy.csv"
weather_file = absolute_path + weather_file


def build_ro_pv():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = SETODatabase()
    m.pysam = PySAMWaterTAP(
        pysam_model="pv",
        tech_config_file=tech_config_file,
        grid_config_file=grid_config_file,
        rate_config_file=rate_config_file,
        cash_config_file=cash_config_file,
        weather_file=weather_file,
    )
    m.fs.properties = NaClParameterBlock()

    treatment = m.fs.treatment = Block()
    energy = m.fs.energy = Block()

    energy.pv = PhotovoltaicZO(property_package=m.fs.properties, database=m.db)
    treatment.feed = Feed(property_package=m.fs.properties)
    treatment.product = Product(property_package=m.fs.properties)
    treatment.disposal = Product(property_package=m.fs.properties)

    treatment.p1 = Pump(property_package=m.fs.properties)

    treatment.ro = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )

    treatment.erd = EnergyRecoveryDevice(property_package=m.fs.properties)

    treatment.a1 = Arc(source=treatment.feed.outlet, destination=treatment.p1.inlet)
    treatment.a2 = Arc(source=treatment.p1.outlet, destination=treatment.ro.inlet)
    treatment.a3 = Arc(
        source=treatment.ro.permeate, destination=treatment.product.inlet
    )
    treatment.a4 = Arc(source=treatment.ro.retentate, destination=treatment.erd.inlet)
    treatment.a5 = Arc(
        source=treatment.erd.outlet, destination=treatment.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(treatment)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    set_scaling_factor(treatment.p1.control_volume.work, 1e-3)
    set_scaling_factor(treatment.ro.area, 1e-2)
    treatment.feed.properties[0].flow_vol_phase["Liq"]
    treatment.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    energy.pv.properties[0].flow_vol_phase["Liq"]
    energy.pv.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
    set_scaling_factor(treatment.erd.control_volume.work, 1e-3)
    calculate_scaling_factors(m)

    return m


def set_operating_conditions(
    flow_in=1e-3, tds=35, ro_area_guess=50, water_recovery=0.5
):
    pv = m.fs.energy.pv
    p1 = m.fs.treatment.p1
    feed = m.fs.treatment.feed
    feed.properties[0].pressure.fix(101325)
    feed.properties[0].temperature.fix(298.15)
    feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): tds * 1e-3,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )
    pv.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): tds * 1e-3,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )
    p1.efficiency_pump.fix(0.8)
    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.treatment.feed.properties[0],
        solver=solver,
        over_pressure=0.25,
        water_recovery=water_recovery,
        NaCl_passage=0.01,
    )
    operating_pressure_psi = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.psi
    )()
    operating_pressure_bar = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.bar
    )()
    print(
        f"\nOperating Pressure Estimate = {round(operating_pressure_bar, 2)} bar = {round(operating_pressure_psi, 2)} psi\n"
    )
    p1.control_volume.properties_out[0].pressure.fix(operating_pressure)
    m.db.get_unit_operation_parameters("solar_energy")
    pv.load_parameters_from_database(use_default_removal=True)


def initialize_treatment(ro_area_guess=50, water_recovery=0.5):
    ro = m.fs.treatment.ro
    p1 = m.fs.treatment.p1
    feed = m.fs.treatment.feed
    erd = m.fs.treatment.erd
    p1 = m.fs.treatment.p1

    ro.A_comp.fix(4.2e-12)
    ro.B_comp.fix(3.5e-8)
    ro.feed_side.channel_height.fix(1e-3)
    ro.feed_side.spacer_porosity.fix(0.97)
    ro.permeate.pressure[0].fix(101325)
    ro.width.fix(5)

    ro.feed_side.properties_in[0].flow_mass_phase_comp[
        "Liq", "H2O"
    ] = p1.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]()
    ro.feed_side.properties_in[0].flow_mass_phase_comp[
        "Liq", "NaCl"
    ] = p1.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"]()
    ro.feed_side.properties_in[0].temperature = feed.properties[0].temperature()
    ro.feed_side.properties_in[0].pressure = p1.control_volume.properties_out[
        0
    ].pressure()
    ro.area.fix(ro_area_guess)
    ro.initialize()
    ro.area.unfix()
    ro.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)
    propagate_state(m.fs.treatment.a4)

    erd.efficiency_pump.fix(0.95)
    erd.control_volume.properties_out[0].pressure.fix(101325)
    erd.initialize()

    propagate_state(m.fs.treatment.a1)
    p1.initialize()
    propagate_state(m.fs.treatment.a2)


def initialize_energy():
    m.fs.energy.pv.initialize()


def initialize_sys(ro_area_guess=50, water_recovery=0.5):
    optarg = solver.options
    m.fs.treatment.feed.initialize(optarg=optarg)
    initialize_treatment(ro_area_guess=ro_area_guess, water_recovery=water_recovery)
    initialize_energy()


def optimize_setup(
    press_lb=125,
    press_ub=1200,
    area_lb=1,
    area_ub=150,
    prod_salinity=500e-6,
    min_flux=2.5e-4,
):
    ro = m.fs.treatment.ro
    p1 = m.fs.treatment.p1

    press_lb_Pa = pyunits.convert(press_lb * pyunits.psi, to_units=pyunits.Pa)()
    press_ub_Pa = pyunits.convert(press_ub * pyunits.psi, to_units=pyunits.Pa)()

    p1.control_volume.properties_out[0].pressure.unfix()
    p1.control_volume.properties_out[0].pressure.setlb(press_lb_Pa)
    p1.control_volume.properties_out[0].pressure.setub(press_ub_Pa)
    p1.deltaP.setlb(0)

    ro.area.unfix()
    ro.area.setlb(area_lb)
    ro.area.setub(area_ub)

    m.fs.treatment.prod_salinity = Param(initialize=prod_salinity, mutable=True)
    m.fs.treatment.min_flux = Param(initialize=min_flux, mutable=True)

    m.fs.treatment.eq_prod_salinity = Constraint(
        expr=m.fs.treatment.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.treatment.prod_salinity
    )

    constraint_scaling_transform(m.fs.treatment.eq_prod_salinity, 1e4)

    m.fs.treatment.eq_min_flux = Constraint(
        expr=m.fs.treatment.ro.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.treatment.min_flux
    )

    assert_degrees_of_freedom(m.fs.treatment, 1)


def add_costing():
    treatment = m.fs.treatment
    energy = m.fs.energy
    treatment.costing = TreatmentCosting()
    energy.costing = EnergyCosting()
    m.db.get_unit_operation_parameters("photovoltaic")
    energy.pv.load_parameters_from_database(use_default_removal=True)

    energy.pv.costing = UnitModelCostingBlock(flowsheet_costing_block=energy.costing)
    treatment.ro.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )
    treatment.erd.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )
    treatment.p1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=treatment.costing
    )

    treatment.costing.cost_process()
    energy.costing.cost_process()

    m.fs.sys_costing = SETOSystemCosting()
    flow_out = treatment.product.properties[0].flow_vol
    m.fs.sys_costing.add_LCOW(flow_out)
    m.fs.sys_costing.add_specific_electric_energy_consumption(flow_out)

    treatment.costing.initialize()
    energy.costing.initialize()


def solve_it(solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")
    return results


def fix_pv_costing():
    pvc = m.fs.energy.pv.costing
    pvc.system_capacity.fix(0)
    pvc.annual_generation.fix(0)
    pvc.land_area.fix(0)


def fix_treatment_global_params():
    tc = m.fs.treatment.costing
    tc.factor_total_investment.fix(1)
    tc.factor_maintenance_labor_chemical.fix(0)


def fix_pysam_costing():
    tech_model = m.pysam.tech_model
    cash_model = m.pysam.cash_model

    avg_gen = np.mean(m.pysam.hourly_energy)
    m.fs.energy.pv.electricity.fix(-1 * avg_gen)

    annual_gen = m.pysam.annual_energy * pyunits.kWh
    annual_gen = pyunits.convert(annual_gen, to_units=pyunits.MWh)()
    land_area = cash_model.LandLease.land_area * pyunits.acres

    pvc = m.fs.energy.pv.costing
    pvc.land_area.fix(land_area)
    pv_params = m.fs.energy.costing.photovoltaic
    pv_params.fixed_operating_by_capacity.fix(cash_model.SystemCosts.om_capacity[0])
    m.fs.energy.pv.costing.system_capacity.fix(m.pysam.nameplate_dc * 1000)
    pv_params.variable_operating_by_generation.fix(
        cash_model.SystemCosts.om_production[0]
    )
    pvc.annual_generation.fix(annual_gen)


flow_in = 4.38e-3  # m3/s
tds = 50  # g/L
ro_area_guess = 50
water_recovery = 0.5
press_lb = 125  # psi
press_ub = 5500  # psi
area_lb = 1  # m2
area_ub = 1200  # m2
prod_salinity = 200e-4
min_flux = 0.5e-6

oversize_factor = 1

m = build_ro_pv()
set_operating_conditions(
    flow_in=flow_in,
    tds=tds,
    ro_area_guess=ro_area_guess,
    water_recovery=water_recovery,
)
initialize_sys()
optimize_setup(
    press_lb=press_lb,
    press_ub=press_ub,
    area_lb=area_lb,
    area_ub=area_ub,
    prod_salinity=prod_salinity,
    min_flux=min_flux,
)
add_costing()

m.fs.obj = Objective(expr=m.fs.sys_costing.LCOW)

m.fs.energy.pv.oversize_factor = oversize_factor

fix_pv_costing()
fix_treatment_global_params()
m.results = solve_it()
display_ro_pv_results(m)
desired_pv_size = (
    m.fs.treatment.costing.aggregate_flow_electricity() * m.fs.energy.pv.oversize_factor
)

cash_model_kwargs = {"om_fixed": 1e4, "om_production": 20}
m.pysam.run_pv_single_owner(
    desired_size=desired_pv_size, cash_model_kwargs=cash_model_kwargs
)

m.fs.sys_costing.add_LCOE()

fix_pysam_costing()

m.fs.treatment.ro.area.fix()
m.fs.treatment.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix()
m.fs.treatment.p1.control_volume.properties_out[0].pressure.unfix()
# m.fs.treatment.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].unfix()
m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
# m.fs.treatment.ro.recovery_mass_phase_comp[0, 'Liq', 'H2O'].fix(water_recovery)

m.results = solve_it(check_termination=False)
display_ro_pv_results(m)
