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

import os
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    Param,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    log10,
    Block,
    value,
    Objective,
    Constraint,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    set_scaling_factor,
    calculate_scaling_factors,
    constraint_scaling_transform,
)
from idaes.core.util.initialization import propagate_state

from watertap.core.solvers import get_solver
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

from watertap_contrib.reflo.analysis.net_metering.util import (
    display_ro_pv_results,
    display_pv_results,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.solar_models.zero_order import Photovoltaic
from watertap_contrib.reflo.core import PySAMWaterTAP


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
    """Builds the structure of the PV-RO system

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
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

    energy.pv = Photovoltaic()
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

    treatment.feed.properties[0].flow_vol_phase["Liq"]
    treatment.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    set_scaling_factor(treatment.erd.control_volume.work, 1e-3)

    return m


def set_operating_conditions(m, flow_in=1e-2, conc_in=30, water_recovery=0.5):
    """Sets operating condition for the PV-RO system

    Args:
        m (obj): Pyomo model
        flow_in (float, optional): feed volumetric flow rate [m3/s]. Defaults to 1e-2.
        conc_in (int, optional): solute concentration [g/L]. Defaults to 30.
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
    m.fs.treatment.feed.properties[0].pressure.fix(101325)
    m.fs.treatment.feed.properties[0].temperature.fix(298.15)
    m.fs.treatment.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): conc_in * 1e-3,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    operating_pressure = calculate_operating_pressure(
        feed_state_block=m.fs.treatment.feed.properties[0],
        solver=solver,
        over_pressure=0.0,
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

    p1_work_scale_est = np.floor(
        log10(flow_in * operating_pressure)
    )  # Estimating pump work to predict p1 scaling factor
    set_scaling_factor(
        m.fs.treatment.p1.control_volume.work, 10 ** (-p1_work_scale_est)
    )

    m.fs.treatment.p1.control_volume.properties_out[0].pressure.fix(operating_pressure)
    m.fs.treatment.p1.efficiency_pump.fix(0.8)


def initialize_treatment(m, water_recovery=0.5):
    """

    Args:
        m (obj): Pyomo model
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
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
    ro.feed_side.velocity[0, 0].fix(0.55)  # crossflow velocity (m/s)
    # ro.feed_side.velocity[0, 0].unfix() # crossflow velocity (m/s)

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

    est_flux = 30.0  # [LMH] Use to help initialization. If initial guess for mem area is too far off, model will fail
    ro_area_guess = (
        water_recovery
        * 1000
        * (
            value(
                pyunits.convert(
                    m.fs.treatment.feed.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.m**3 / pyunits.hr,
                )
            )
        )
    ) / est_flux
    ro_area_scale = np.floor(
        log10(ro_area_guess)
    )  # Estimating membrane area to predict scaling factor

    print(
        f"\nRO Membrane Area Estimate = {round(ro_area_guess, 2)} m^2 assuming a {round(est_flux, 2)} LMH Water Flux\n"
    )
    ro.area.fix(10**ro_area_scale)
    set_scaling_factor(ro.area, 10 ** (-ro_area_scale))
    calculate_scaling_factors(m)

    # solve feed
    solver.solve(feed)

    # # initialize pump
    propagate_state(m.fs.treatment.a1)
    p1.initialize()

    propagate_state(m.fs.treatment.a2)
    ro.initialize()
    ro.area.unfix()
    ro.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)

    propagate_state(m.fs.treatment.a4)
    erd.efficiency_pump.fix(0.95)
    erd.control_volume.properties_out[0].pressure.fix(101325)
    erd.initialize()

    propagate_state(m.fs.treatment.a3)
    propagate_state(m.fs.treatment.a5)

    print(
        f"\nFeed Flowrate = {pyunits.convert(ro.feed_side.properties_in[0].flow_vol, to_units=pyunits.Mgallons/pyunits.day)():<2.3f} MGD ; Feed Side Velocity = {value(ro.feed_side.velocity[0, 0])} m/s ; Solute Conc. In = {pyunits.convert(ro.feed_side.properties_in[0].conc_mass_phase_comp['Liq', 'NaCl'], to_units=pyunits.g/pyunits.L)():<3.1f} g/L \n"
    )


def initialize_energy(m):
    m.fs.energy.pv.initialize()
    m.fs.energy.pv.oversize_factor.set_value(1)


def initialize_system(m, water_recovery=0.5):
    optarg = solver.options
    m.fs.treatment.feed.initialize(optarg=optarg)
    initialize_treatment(m, water_recovery=water_recovery)
    initialize_energy(m)


def optimize_setup(
    m,
    opt_target,
    press_lb=100,
    press_ub=4500,
    area_lb=1,
    area_ub=20000,
    prod_salinity=0.5,
    min_flux=2.5e-4,
):
    m.fs.obj = Objective(expr=opt_target)
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

    ro.feed_side.velocity[0, 0].unfix()
    ro.feed_side.velocity[0, 0].setlb(0.01)
    ro.feed_side.velocity[0, 0].setub(1)

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


def add_costing(m):
    treatment = m.fs.treatment
    energy = m.fs.energy
    treatment.costing = TreatmentCosting()
    energy.costing = EnergyCosting()

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

    m.fs.sys_costing = REFLOCosting()
    m.fs.sys_costing.add_LCOW(treatment.product.properties[0].flow_vol)
    m.fs.sys_costing.add_specific_electric_energy_consumption(
        treatment.product.properties[0].flow_vol
    )

    treatment.costing.initialize()
    energy.costing.initialize()


def fix_pv_costing(m):
    m.fs.energy.pv.costing.system_capacity.fix(0)
    m.fs.energy.pv.costing.annual_generation.fix(0)
    m.fs.energy.pv.costing.land_area.fix(0)


def fix_treatment_global_params(m):
    m.fs.treatment.costing.total_investment_factor.fix(1)
    m.fs.treatment.costing.maintenance_labor_chemical_factor.fix(0)


def size_pv(m):
    desired_pv_size = (
        m.fs.treatment.costing.aggregate_flow_electricity()
        * m.fs.energy.pv.oversize_factor()
    )
    cash_model_kwargs = {"om_fixed": 1e4, "om_production": 20}
    m.pysam.run_pv_single_owner(
        desired_size=desired_pv_size, cash_model_kwargs=cash_model_kwargs
    )
    m.fs.sys_costing.add_LCOE()


def fix_pysam_costing(m):
    tech_model = m.pysam.tech_model
    cash_model = m.pysam.cash_model

    avg_gen = np.mean(m.pysam.hourly_energy)
    m.fs.energy.pv.electricity.fix(-1 * avg_gen)

    annual_gen = pyunits.convert(
        (m.pysam.annual_energy * pyunits.kWh), to_units=pyunits.MWh
    )()
    land_area = cash_model.LandLease.land_area * pyunits.acres

    m.fs.energy.pv.costing.land_area.fix(land_area)
    m.fs.energy.costing.photovoltaic.fixed_operating_by_capacity.fix(
        cash_model.SystemCosts.om_capacity[0]
    )
    m.fs.energy.pv.costing.system_capacity.fix(m.pysam.nameplate_dc * 1000)
    m.fs.energy.costing.photovoltaic.variable_operating_by_generation.fix(
        cash_model.SystemCosts.om_production[0]
    )
    m.fs.energy.pv.costing.annual_generation.fix(annual_gen)
    m.fs.energy.costing.maintenance_labor_chemical_factor.fix(0)


def solve(m, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=tee)
    size_pv(m)
    fix_pysam_costing(m)
    results = solver.solve(m, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")

    return results


def model_setup(Q, conc, recovery):
    m = build_ro_pv()
    set_operating_conditions(m, flow_in=Q, conc_in=conc, water_recovery=recovery)
    initialize_system(m)
    add_costing(m)
    fix_pv_costing(m)
    fix_treatment_global_params(m)
    optimize_setup(m, m.fs.sys_costing.LCOW)

    return m


def run(m):
    results = solve(m)
    assert_optimal_termination(results)
    display_ro_pv_results(m)
    display_pv_results(m)

    return m, results


def main():
    m = model_setup(6.375e-2, 75, 0.5)
    m, results = run(m)

    return m, results


if __name__ == "__main__":
    m, results = main()
