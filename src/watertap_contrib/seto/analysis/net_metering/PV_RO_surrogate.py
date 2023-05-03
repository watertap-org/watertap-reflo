import os
import numpy as np

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Param,
    Constraint,
    Block,
    Var,
    TransformationFactory,
    assert_optimal_termination,
    check_optimal_termination,
    value,
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
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
from watertap_contrib.seto.solar_models.surrogate.pv import PVSurrogate

solver = get_solver()

absolute_path = os.path.dirname(__file__)


def build_ro_pv():
    """Builds the structure of the PV-RO system

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    print(f'\n{"=======> BUILDING FLOWSHEET <=======":^60}\n')

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    treatment = m.fs.treatment = Block()
    energy = m.fs.energy = Block()

    energy.pv = PVSurrogate()
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
    set_scaling_factor(treatment.p1.control_volume.work, 1e-3)   
    calculate_scaling_factors(m)

    return m


def set_operating_conditions(m, flow_in=1e-2, conc_in=30, water_recovery=0.5):
    """Sets operating condition for the PV-RO system

    Args:
        m (obj): Pyomo model
        flow_in (float, optional): feed volumetric flow rate [m3/s]. Defaults to 1e-2.
        conc_in (int, optional): solute concentration [g/L]. Defaults to 30.
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    m.fs.treatment.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): conc_in * 1e-3,
            ("pressure", None): 101325,
            ("temperature", None): 298.15
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
    ro.feed_side.velocity[0, 0].fix(0.4)  # crossflow velocity (m/s)

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
    ro.area.setub(10**(ro_area_scale+1))
    ro.width.setub(10**(ro_area_scale))
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


def initialize_energy(m):
    m.fs.energy.pv.load_surrogate()


def initialize_sys(m, water_recovery=0.5):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    optarg = solver.options
    m.fs.treatment.feed.initialize(optarg=optarg)
    initialize_treatment(m, water_recovery=water_recovery)
    initialize_energy(m)


def optimize_setup(
    m,
    opt_target,
    press_lb=10,
    press_ub=1500,
    area_lb=1,
    area_ub=6000000,
    prod_salinity=0.5,
    min_flux=1.0e-5,
    min_flux_mass_phase_comp = 1.0e-5,
    min_eq_mass_frac_permeate = 1e-6
):
    
    print(f'\n{"=======> OPTIMIZING <=======":^60}\n')

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
    ro.feed_side.velocity[0, 0].setlb(0.0)
    ro.feed_side.velocity[0, 0].setub(0.5)

    m.fs.treatment.prod_salinity = Param(initialize=prod_salinity, mutable=True)
    m.fs.treatment.min_flux = Param(initialize=min_flux, mutable=True)
    m.fs.treatment.flux_mass_phase_comp = Param(initialize=min_flux_mass_phase_comp, mutable=True)

    m.fs.treatment.eq_prod_salinity = Constraint(
        expr=m.fs.treatment.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"]
        <= m.fs.treatment.prod_salinity
    )

    constraint_scaling_transform(m.fs.treatment.eq_prod_salinity, 1e4)

    m.fs.treatment.eq_min_flux_constraint = Constraint(
        expr=m.fs.treatment.ro.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        >= m.fs.treatment.min_flux
    )


def add_costing(m, cap_max = None):
    treatment = m.fs.treatment
    energy = m.fs.energy
    treatment.costing = TreatmentCosting()
    energy.costing = EnergyCosting()

    energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing
        )
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
        
    # m.fs.energy.pv_design_constraint = Constraint(
    #     expr=m.fs.energy.pv.design_size == m.fs.treatment.costing.aggregate_flow_electricity
    # )

    m.fs.energy.pv_electricity_constraint = Constraint(
        expr=m.fs.energy.pv.electricity <= 0
    )
    
    m.fs.energy.pv.costing.land_constraint = Constraint(
        expr=m.fs.energy.pv.costing.land_area == m.fs.energy.pv.land_req
    )
    
    energy.costing.cost_process()
    m.fs.sys_costing = SETOSystemCosting()
    
    m.fs.sys_costing.add_LCOW(treatment.product.properties[0].flow_vol)
    m.fs.sys_costing.add_specific_electric_energy_consumption(
        treatment.product.properties[0].flow_vol
    )
    m.fs.sys_costing.add_LCOE(e_model="surrogate")
    if cap_max != None:
        m.fs.sys_costing.total_capital_cost.setlb(0)
        m.fs.sys_costing.total_capital_cost.setub(cap_max)
    
    treatment.costing.initialize()
    energy.costing.initialize()
    
def fix_treatment_global_params(m):
    m.fs.treatment.costing.factor_total_investment.fix(1)
    m.fs.energy.costing.factor_total_investment.fix(1)
    m.fs.treatment.costing.factor_maintenance_labor_chemical.fix(0)

def automate_rescale_variables(self, rescale_factor=1, default=1):
        if rescale_factor is None:
            rescale_factor = 1
        for var, sv in badly_scaled_var_generator(self):
            sf = get_scaling_factor(var)
            if get_scaling_factor(var) is None:
                print(f"{var} is missing a scaling factor")
                sf = default
                set_scaling_factor(var, sf, data_objects=False)

            set_scaling_factor(var, sf / sv * rescale_factor)
            calculate_scaling_factors(self)

def debug(m, solver=None, automate_rescale=False, resolve=False):
    if solver is None:
        solver = get_solver()

    print(f'\n{"=======> DEBUGGING <=======":^60}\n')
    print(f'\n{"=======> BADLY SCALED VARIABLES <=======":^60}\n')
    badly_scaled_vars = list(badly_scaled_var_generator(m))
    print([print(i[0], i[1]) for i in badly_scaled_vars])
    print(f'\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
    print_infeasible_bounds(m)
    print(f'\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
    print_infeasible_constraints(m)
    print(f'\n{"=======> CONSTRAINTS CLOSE TO BOUNDS <=======":^60}\n')
    print_close_to_bounds(m)

    if automate_rescale:
        print(
            f"\n{len(badly_scaled_vars)} poorly scaled "
            f"variable(s) will be rescaled so that each scaled variable value = 1\n"
        )
        automate_rescale_variables(m)
        badly_scaled_vars = list(badly_scaled_var_generator(m))

    print(
            f"\nNow {len(badly_scaled_vars)} poorly scaled\n"
        )

    if resolve:
        results = solver.solve(m, tee=True)
        debug(m)
        return results


def solve(m, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()

    try:
        results = solver.solve(m, tee=tee)
    except:
        results = debug(m, automate_rescale=True, resolve=True)
    if check_termination:
        termination = check_optimal_termination(results)
        if termination != True:
            results = debug(m, automate_rescale=True, resolve=True)
        assert_optimal_termination(results)
    print(f"\nDOF = {degrees_of_freedom(m)}")
    print(f"MODEL SOLVE = {results.solver.termination_condition.swapcase()}")

    return results

def model_setup(Q, conc, recovery):
    m = build_ro_pv()
    set_operating_conditions(m, flow_in=Q, conc_in=conc, water_recovery=recovery)
    initialize_sys(m, water_recovery=recovery)
    add_costing(m, cap_max = 9.0E6)
    fix_treatment_global_params(m)
    optimize_setup(m, m.fs.sys_costing.LCOW)
    return m

def run(m):
    results = solve(m)
    assert_optimal_termination(results)
    display_ro_pv_results(m)
    return m, results

def mgd_to_m3s(Q):
    return 0.0438*Q

def main():
    m = model_setup(mgd_to_m3s(10), 35, 0.5)
    m, results = run(m)

    return m, results

if __name__ == "__main__":
    m, results = main()
