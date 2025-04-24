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
    SolverFactory,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state

# import idaes.core.util.scaling as iscale
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
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
import pandas as pd

import pathlib

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/kbhdp_case_study.yaml"

__all__ = [
    "build_system",
    "add_connections",
    "add_costing",
    "add_system_costing",
    "set_inlet_conditions",
    "set_operating_conditions",
    "init_system",
    "print_results_summary",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(__location__, "swh-kbhdp.json")


def propagate_state(arc):
    _prop_state(arc)


def build_system(Qin=4, Cin=12, water_recovery=0.5):
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.treatment = Block()
    m.fs.energy = Block()

    m.inlet_flow_rate = pyunits.convert(
        Qin * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = water_recovery

    m.fs.treatment.costing = TreatmentCosting()
    m.fs.energy.costing = EnergyCosting()

    # Property package
    m.fs.properties = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.treatment.feed = Feed(property_package=m.fs.properties)
    m.fs.treatment.product = Product(property_package=m.fs.properties)
    # m.fs.disposal = Product(property_package=m.fs.properties)

    # Create MD unit model at flowsheet level
    m.fs.treatment.md = FlowsheetBlock()

    build_md(
        m,
        m.fs.treatment.md,
    )
    m.fs.treatment.dwi = FlowsheetBlock()
    build_DWI(m, m.fs.treatment.dwi, m.fs.properties)

    m.fs.energy.FPC = FlowsheetBlock()
    build_fpc(m)

    return m


def add_connections(m):
    treatment = m.fs.treatment

    treatment.feed_to_md = Arc(
        source=treatment.feed.outlet, destination=treatment.md.feed.inlet
    )

    treatment.md_to_product = Arc(
        source=treatment.md.permeate.outlet, destination=treatment.product.inlet
    )

    treatment.md_to_dwi = Arc(
        source=treatment.md.concentrate.outlet,
        destination=treatment.dwi.unit.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_costing(m, treatment_costing_block=None, energy_costing_block=None):
    # Solving the system before adding costing
    # solver = SolverFactory("ipopt")
    if treatment_costing_block is None:
        treatment_costing_block = m.fs.treatment.costing
    if energy_costing_block is None:
        energy_costing_block = m.fs.energy.costing
    solver = get_solver()
    solve(m, solver=solver, tee=False)
    add_fpc_costing(m, costing_block=energy_costing_block)
    add_md_costing(m.fs.treatment.md.mp, treatment_costing_block)
    add_DWI_costing(
        m.fs.treatment, m.fs.treatment.dwi, costing_blk=treatment_costing_block
    )
    # System costing
    treatment_costing_block.cost_process()
    energy_costing_block.cost_process()
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()

    print("\n--------- INITIALIZING SYSTEM COSTING ---------\n")

    treatment_costing_block.initialize()
    energy_costing_block.initialize()
    m.fs.costing.initialize()
    m.fs.costing.add_annual_water_production(
        m.fs.treatment.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)
    m.fs.costing.add_LCOH()


def calc_costing(m, heat_price=0.01, electricity_price=0.07):
    # Touching variables to solve for volumetric flow rate
    m.fs.product.properties[0].flow_vol_phase

    # Treatment costing
    # Overwriting values in yaml
    m.fs.treatment.costing.heat_cost.fix(heat_price)
    m.fs.treatment.costing.electricity_cost.fix(electricity_price)
    m.fs.treatment.costing.cost_process()

    m.fs.treatment.costing.initialize()

    m.fs.treatment.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol
    )
    m.fs.treatment.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    # Energy costing
    m.fs.energy.costing.electricity_cost.fix(electricity_price)
    m.fs.energy.costing.cost_process()

    m.fs.energy.costing.initialize()
    m.fs.energy.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.energy.costing.add_LCOH()


def add_constraints(m):
    treatment = m.fs.treatment

    m.fs.water_recovery = Var(
        initialize=m.water_recovery,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=treatment.feed.properties[0].flow_vol * m.fs.water_recovery
        == treatment.product.properties[0].flow_vol
    )


def set_scaling(m):
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "TDS"))

    calculate_scaling_factors(m)


def set_inlet_conditions(m):
    print(f'\n{"=======> SETTING FEED CONDITIONS <=======":^60}\n')

    m.fs.treatment.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.inlet_flow_rate,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.inlet_salinity,
            ("temperature", None): 298.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def set_operating_conditions(m, hours_storage=8):
    set_inlet_conditions(m)
    set_fpc_op_conditions(m, hours_storage=hours_storage, temperature_hot=80)


def init_system(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    treatment.feed.initialize()

    init_md(treatment.md)

    propagate_state(treatment.md_to_product)
    treatment.product.initialize()

    propagate_state(treatment.md_to_dwi)
    # m.fs.disposal.initialize()

    init_DWI(m, blk.treatment.dwi, verbose=True, solver=None)

    init_fpc(m.fs.energy)


def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()
    # solver.options["max_iter"] = 5000
    # solver.options["halt_on_ampl_error"] = "yes"

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


def optimize(m):
    m.fs.costing.frac_heat_from_grid.unfix()
    m.fs.obj = Objective(expr=m.fs.costing.LCOT)


def report_costing(blk):
    print(f"\n\n-------------------- System Costing Report --------------------\n")
    print("\n")

    print(f'{"LCOT":<30s}{value(blk.LCOT):<20,.2f}{pyunits.get_units(blk.LCOT)}')

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


def main(
    water_recovery=0.5,
    heat_price=0.07,
    electricity_price=0.07,
    frac_heat_from_grid=0.01,
    hours_storage=8,
    run_optimization=True,
):
    # Build  MD, DWI and FPC
    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    set_scaling(m)
    # check_jac(m)
    init_system(m, m.fs)
    print(f"dof = {degrees_of_freedom(m)}")
    results = solve(m.fs.treatment.md.mp)
    # results = solve(m.fs.treatment.md)
    # results = solve(m.fs.treatment)
    # results = solve(m.fs.energy)
    results = solve(m)
    print(f"termination {results.solver.termination_condition}")
    add_costing(m)

    print(f"dof = {degrees_of_freedom(m)}")
    results = solve(m)
    print(f"termination 2 {results.solver.termination_condition}")


def print_results_summary(m):
    print(f"\nAfter Optimization System Degrees of Freedom: {degrees_of_freedom(m)}")

    print("\n")
    print(
        f'{"Treatment LCOW":<30s}{value(m.fs.treatment.costing.LCOW):<10.2f}{pyunits.get_units(m.fs.treatment.costing.LCOW)}'
    )

    print("\n")
    print(
        f'{"Energy LCOH":<30s}{value(m.fs.energy.costing.LCOH):<10.2f}{pyunits.get_units(m.fs.energy.costing.LCOH)}'
    )

    print("\n")
    print(
        f'{"System LCOT":<30s}{value(m.fs.costing.LCOT) :<10.2f}{pyunits.get_units(m.fs.costing.LCOT)}'
    )

    print("\n")
    print(
        f'{"Percent from the grid":<30s}{value(m.fs.costing.frac_heat_from_grid):<10.2f}{pyunits.get_units(m.fs.costing.frac_heat_from_grid)}'
    )

    report_MD(m, m.fs.treatment.md)
    report_md_costing(m, m.fs.treatment)

    print_DWI_costing_breakdown(m.fs.treatment, m.fs.treatment.dwi)

    report_fpc(m, m.fs.energy.fpc.unit)
    report_fpc_costing(m, m.fs.energy)
    report_costing(m.fs.costing)


def save_results(m):
    results_df = pd.DataFrame(
        columns=[
            "water_recovery",
            "heat_price",
            "LCOH",
            "hours_storage",
            "frac_heat_from_grid",
            "product_annual_production",
            "utilization_factor",
            "capital_recovery_factor",
            "unit",
            "cost_component",
            "cost",
            "norm_cost_component",
        ]
    )

    capex_output = {
        "FPC": value(m.fs.energy.fpc.unit.costing.capital_cost)
        * value(m.fs.costing.capital_recovery_factor),
        "MD": value(
            m.fs.treatment.md.unit.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.capital_cost
        )
        * value(m.fs.costing.capital_recovery_factor),
        "DWI": 0,
        "Heat": 0,
        "Electricity": 0,
    }

    fixed_opex_output = {
        "FPC": value(m.fs.energy.fpc.unit.costing.fixed_operating_cost)
        + value(m.fs.energy.fpc.unit.costing.capital_cost)
        * value(m.fs.energy.costing.maintenance_labor_chemical_factor),
        "MD": value(
            m.fs.treatment.md.unit.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.fixed_operating_cost
        )
        + value(
            m.fs.treatment.md.unit.get_active_process_blocks()[
                -1
            ].fs.vagmd.costing.capital_cost
        )
        * value(m.fs.treatment.costing.maintenance_labor_chemical_factor),
        "DWI": 0,
        "Heat": 0,
        "Electricity": 0,
    }
    variable_opex_output = {
        "FPC": 0,
        "MD": 0,
        "DWI": value(m.fs.treatment.dwi.unit.costing.variable_operating_cost),
        "Heat": value(m.fs.costing.total_heat_operating_cost),
        "Electricity": value(m.fs.costing.total_electric_operating_cost),
    }

    for unit in ["FPC", "MD", "DWI", "Heat", "Electricity"]:
        # Add fixed_opex
        temp = {
            "water_recovery": value(m.fs.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.fpc.unit.hours_storage),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "product_annual_production": value(m.fs.costing.annual_water_production),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "fixed_opex",
            "cost": fixed_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)
        # Add variable opex
        temp = {
            "water_recovery": value(m.fs.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.fpc.unit.hours_storage),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "product_annual_production": value(m.fs.costing.annual_water_production),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "variable_opex",
            "cost": variable_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

        # Add opex
        temp = {
            "water_recovery": value(m.fs.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.fpc.unit.hours_storage),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "product_annual_production": value(m.fs.costing.annual_water_production),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "opex",
            "cost": variable_opex_output[unit] + fixed_opex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

        # Add capex
        temp = {
            "water_recovery": value(m.fs.water_recovery),
            "heat_price": value(m.fs.costing.heat_cost_buy),
            "LCOH": value(m.fs.energy.costing.LCOH),
            "hours_storage": value(m.fs.energy.fpc.unit.hours_storage),
            "frac_heat_from_grid": value(m.fs.costing.frac_heat_from_grid),
            "product_annual_production": value(m.fs.costing.annual_water_production),
            "utilization_factor": value(m.fs.costing.utilization_factor),
            "capital_recovery_factor": value(m.fs.costing.capital_recovery_factor),
            "unit": unit,
            "cost_component": "capex",
            "cost": capex_output[unit],
        }
        results_df = results_df.append(temp, ignore_index=True)

    results_df["norm_cost_component"] = (
        results_df["cost"]
        / results_df["product_annual_production"]
        / results_df["utilization_factor"]
    )

    file_name = (
        "RPT3_water_recovery_"
        + str(value(m.fs.water_recovery))
        + "_heat_price_"
        + str(value(m.fs.costing.heat_cost_buy))
        + "_hours_storage_"
        + str(value(m.fs.energy.fpc.unit.hours_storage))
    )

    # results_df.to_csv(
    #     r"C:\Users\mhardika\Documents\SETO\Case Studies\RPT3\RPT3_results\\"
    #     + file_name
    #     + ".csv"
    # )
    # Flow cost


if __name__ == "__main__":
    main(
        water_recovery=0.8,
        heat_price=0.08,
        electricity_price=0.07,
        frac_heat_from_grid=0.5,
        hours_storage=6,
        run_optimization=False,
    )
