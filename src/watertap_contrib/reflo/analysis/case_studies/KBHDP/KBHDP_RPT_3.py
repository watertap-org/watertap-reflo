import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    TransformationFactory,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
    REFLOSystemCosting,
)

from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.MD import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.FPC import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.deep_well_injection import *


def propagate_state(arc):
    _prop_state(arc)


def build_system(inlet_cond, n_time_points):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.treatment = Block()
    m.fs.energy = Block()

    m.fs.treatment.costing = TreatmentCosting()
    m.fs.energy.costing = EnergyCosting()

    # Property package
    m.fs.params = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.params)
    m.fs.product = Product(property_package=m.fs.params)
    # m.fs.disposal = Product(property_package=m.fs.params)

    # Create MD unit model at flowsheet level
    m.fs.treatment.md = FlowsheetBlock(dynamic=False)
    model_options, n_time_points = build_md(
        m, m.fs.treatment.md, inlet_cond, n_time_points
    )
    m.fs.treatment.dwi = FlowsheetBlock(dynamic=False)
    build_DWI(m, m.fs.treatment.dwi, m.fs.params)

    m.fs.energy.fpc = FlowsheetBlock(dynamic=False)
    build_fpc(m.fs.energy.fpc)

    add_connections(m)
    return m, model_options, n_time_points


def add_connections(m):

    m.fs.feed_to_md = Arc(
        source=m.fs.feed.outlet, destination=m.fs.treatment.md.feed.inlet
    )

    m.fs.md_to_product = Arc(
        source=m.fs.treatment.md.permeate.outlet, destination=m.fs.product.inlet
    )

    m.fs.md_to_dwi = Arc(
        source=m.fs.treatment.md.concentrate.outlet,
        destination=m.fs.treatment.dwi.unit.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


# TODO: Should the recovery be a constraint of a variable thats fixed?
def add_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.feed_flow_vol = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.L / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.product.properties[0].conc_mass_phase_comp
    # m.fs.disposal.properties[0].conc_mass_phase_comp


# TODO: Add fraction of heat as solar
def add_energy_constraints(m):

    @m.Constraint()
    def eq_thermal_req(b):
        return (
            b.fs.energy.costing.aggregate_flow_heat
            == -1 * b.fs.treatment.costing.aggregate_flow_heat * 0.5
        )


def add_costing(m, treatment_costing_block, energy_costing_block):
    add_fpc_costing(m.fs.energy.fpc, energy_costing_block)
    add_md_costing(m.fs.treatment.md.unit, treatment_costing_block)
    add_DWI_costing(m.fs.treatment, m.fs.treatment.dwi, treatment_costing_block)


def calc_costing(m):

    # Treatment costing
    m.fs.treatment.costing.capital_recovery_factor.fix(0.08764)
    m.fs.treatment.costing.wacc.unfix()
    m.fs.treatment.costing.cost_process()

    m.fs.treatment.costing.initialize()

    m.fs.treatment.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol
    )
    m.fs.treatment.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    # Energy costing
    m.fs.energy.costing.cost_process()

    m.fs.energy.costing.initialize()
    m.fs.energy.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.energy.costing.add_LCOW(m.fs.product.properties[0].flow_vol)


def add_system_costing(m):
    # System costing
    m.fs.sys_costing = REFLOSystemCosting()
    # m.fs.sys_costing.wacc.unfix()
    m.fs.sys_costing.frac_from_grid.fix(0.01)
    m.fs.sys_costing.heat_cost_sell.set_value(0)
    m.fs.sys_costing.cost_process()

    m.fs.sys_costing.initialize()

    m.fs.sys_costing.add_LCOW(m.fs.product.properties[0].flow_vol)


def set_inlet_conditions(blk, model_options):

    print(f'\n{"=======> SETTING FEED CONDITIONS <=======":^60}\n')
    feed_flow_rate = model_options["feed_flow_rate"]
    feed_salinity = model_options["feed_salinity"]
    feed_temp = model_options["feed_temp"]

    blk.feed.properties.calculate_state(
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


def set_operating_conditions(m):

    set_md_op_conditions(m.fs.treatment.md)
    set_fpc_op_conditions(m.fs.energy.fpc, hours_storage=8, temperature_hot=80)


def init_system(m, blk, model_options, n_time_points, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_md)

    init_md(blk.treatment.md, model_options, n_time_points)

    propagate_state(m.fs.md_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.md_to_dwi)
    # m.fs.disposal.initialize()

    init_DWI(m, blk.treatment.dwi, verbose=True, solver=None)

    init_fpc(blk.energy.fpc)


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print_infeasible_bounds(m)
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print(msg)
        return results


def report_sys_costing(blk):

    print(f"\n\n-------------------- System Costing Report --------------------\n")
    print("\n")

    print(f'{"LCOW":<30s}{value(blk.LCOW):<20,.2f}{pyunits.get_units(blk.LCOW)}')

    print(
        f'{"Capital Cost":<30s}{value(blk.total_capital_cost):<20,.2f}{pyunits.get_units(blk.total_capital_cost)}'
    )

    print(
        f'{"Total Operating Cost":<30s}{value(blk.total_operating_cost):<20,.2f}{pyunits.get_units(blk.total_operating_cost)}'
    )

    print(
        f'{"Agg Fixed Operating Cost":<30s}{value(blk.aggregate_fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.aggregate_fixed_operating_cost)}'
    )

    print(
        f'{"Agg Variable Operating Cost":<30s}{value(blk.aggregate_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.aggregate_variable_operating_cost)}'
    )

    print(
        f'{"Heat flow":<30s}{value(blk.aggregate_flow_heat):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat)}'
    )

    # print(
    #     f'{"Total heat cost":<30s}{value(blk.total_heat_operating_cost):<20,.2f}{pyunits.get_units(blk.total_heat_operating_cost)}'
    # )

    print(
        f'{"Heat purchased":<30s}{value(blk.aggregate_flow_heat_purchased):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat_purchased)}'
    )

    print(
        f'{"Heat sold":<30s}{value(blk.aggregate_flow_heat_sold):<20,.2f}{pyunits.get_units(blk.aggregate_flow_heat_sold)}'
    )

    print(
        f'{"Elec Flow":<30s}{value(blk.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity)}'
    )

    # print(
    #     f'{"Total elec cost":<30s}{value(blk.total_electric_operating_cost):<20,.2f}{pyunits.get_units(blk.total_electric_operating_cost)}'
    # )

    print(
        f'{"Elec purchased":<30s}{value(blk.aggregate_flow_electricity_purchased):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity_purchased)}'
    )

    print(
        f'{"Elec sold":<30s}{value(blk.aggregate_flow_electricity_sold):<20,.2f}{pyunits.get_units(blk.aggregate_flow_electricity_sold)}'
    )


if __name__ == "__main__":

    solver = get_solver()
    solver = SolverFactory("ipopt")

    Qin = 4 * pyunits.Mgal / pyunits.day  # KBHDP flow rate
    Qin = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.day)

    feed_salinity = 12 * pyunits.g / pyunits.L

    overall_recovery = 0.5

    inlet_cond = {
        "inlet_flow_rate": Qin,
        "inlet_salinity": feed_salinity,
        "recovery": overall_recovery,
    }

    n_time_points = None

    # Build  MD, DWI and FPC
    m, model_options, n_time_points = build_system(
        inlet_cond, n_time_points=n_time_points
    )
    set_inlet_conditions(m.fs, model_options)

    init_system(m, m.fs, model_options, n_time_points)

    set_operating_conditions(m)

    print(f"\nBefore Costing System Degrees of Freedom: {degrees_of_freedom(m)}")

    results = solver.solve(m)

    # Touching variables to solve for volumetric flow rate
    m.fs.product.properties[0].flow_vol_phase
    # m.fs.disposal.properties[0].flow_vol_phase

    add_costing(
        m,
        treatment_costing_block=m.fs.treatment.costing,
        energy_costing_block=m.fs.energy.costing,
    )

    calc_costing(m)
    add_system_costing(m)
    iscale.set_scaling_factor(
        m.fs.sys_costing.aggregate_flow_electricity_purchased, 1e4
    )
    # add_energy_constraints(m)

    print(f"\nAfter Costing System Degrees of Freedom: {degrees_of_freedom(m)}")
    print("Treatment block:", degrees_of_freedom(m.fs.treatment))
    print("Treatment costing block:", degrees_of_freedom(m.fs.treatment.costing))
    print("Energy block:", degrees_of_freedom(m.fs.energy))
    print("Energy costing block:", degrees_of_freedom(m.fs.energy.costing))

    results = solver.solve(m)

    print_infeasible_constraints(m)

    print(f"\nAfter Solve System Degrees of Freedom: {degrees_of_freedom(m)}")
    print("Treatment block:", degrees_of_freedom(m.fs.treatment))
    print("Treatment costing block:", degrees_of_freedom(m.fs.treatment.costing))
    print("Energy block:", degrees_of_freedom(m.fs.energy))
    print("Energy costing block:", degrees_of_freedom(m.fs.energy.costing))

    m.fs.sys_costing.frac_from_grid.unfix()

    m.fs.obj = Objective(expr=m.fs.sys_costing.LCOW)
    results = solver.solve(m)

    print(f"\nAfter Optimization System Degrees of Freedom: {degrees_of_freedom(m)}")

    print("\n")
    print(
        f'{"Treatment LCOW":<30s}{value(m.fs.treatment.costing.LCOW):<10.2f}{pyunits.get_units(m.fs.treatment.costing.LCOW)}'
    )

    print("\n")
    print(
        f'{"Energy LCOW":<30s}{value(m.fs.energy.costing.LCOW):<10.2f}{pyunits.get_units(m.fs.energy.costing.LCOW)}'
    )

    print("\n")
    print(
        f'{"System LCOW":<30s}{value(m.fs.sys_costing.LCOW):<10.2f}{pyunits.get_units(m.fs.sys_costing.LCOW)}'
    )

    report_MD(m, m.fs.treatment.md)
    report_md_costing(m, m.fs.treatment)

    print_DWI_costing_breakdown(m.fs.treatment, m.fs.treatment.dwi)

    report_fpc(m, m.fs.energy.fpc.unit)
    report_fpc_costing(m, m.fs.energy)
    report_sys_costing(m.fs.sys_costing)

    # print("\nTreatment")
    # m.fs.treatment.costing.aggregate_fixed_operating_cost.display()
    # m.fs.treatment.costing.aggregate_variable_operating_cost.display()
    # m.fs.treatment.costing.total_operating_cost.display()
    # m.fs.treatment.costing.total_capital_cost.display()

    # print("\nEnergy")
    # m.fs.energy.costing.aggregate_fixed_operating_cost.display()
    # m.fs.energy.costing.aggregate_variable_operating_cost.display()
    # m.fs.energy.costing.total_capital_cost.display()
    # m.fs.energy.costing.total_operating_cost.display()

    # m.fs.sys_costing.display()
    # m.fs.energy.costing.display()
    # m.fs.energy.display()
