from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
    Block,
    Constraint,
    SolverFactory,
    check_optimal_termination,
)
import os

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver

from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.scaling import *
from idaes.core.util.model_statistics import *
from watertap_contrib.reflo.solar_models.surrogate.flat_plate.flat_plate_surrogate import (
    FlatPlateSurrogate,
)

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from watertap_contrib.reflo.costing import (
    EnergyCosting,
)
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    check_jac,
    calc_scale,
)

__all__ = [
    "build_fpc",
    "init_fpc",
    "set_fpc_op_conditions",
    "add_fpc_costing",
    "add_FPC_scaling",
    "add_fpc_costing_scaling",
    "report_fpc",
    "print_FPC_costing_breakdown",
]
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.abspath(os.path.join(__location__, ".."))
weather_file = os.path.join(parent_dir, "el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(parent_dir, 'data/fpc', "solar_water_heating-kbhdp.json")
dataset_filename = os.path.join(parent_dir, 'data/fpc', "FPC_KBHDP_el_paso.pkl")
surrogate_filename = os.path.join(parent_dir, 'data/fpc', "FPC_KBHDP_el_paso.json")


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()
    m.fs.energy = Block()

    m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)

    m.fs.fpc = FlowsheetBlock(dynamic=False)

    return m


def build_fpc(m):
    energy = m.fs.energy

    print(f'\n{"=======> BUILDING FPC SYSTEM <=======":^60}\n')

    input_bounds = dict(heat_load=[0.5, 200])
    input_units = dict(heat_load="MW")
    input_variables = {
        "labels": ["heat_load"],
        "bounds": input_bounds,
        "units": input_units,
    }

    output_units = dict(heat_annual="kWh", electricity_annual="kWh")
    output_variables = {
        "labels": ["heat_annual", "electricity_annual"],
        "units": output_units,
    }


    energy.FPC = FlatPlateSurrogate(
        surrogate_model_file=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=False,
    )


def init_fpc(blk):
    blk.FPC.initialize()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()


def set_fpc_op_conditions(m, hours_storage=12, temperature_hot=80):
    energy = m.fs.energy
    # energy.FPC.load_surrogate()

    # energy.FPC.hours_storage.fix(hours_storage)
    # # Assumes the hot temperature to the inlet of a 'MD HX'
    # energy.FPC.temperature_hot.fix(temperature_hot)
    # Assumes the cold temperature from the outlet temperature of a 'MD HX'
    # energy.FPC.temperature_cold.set_value(20)
    energy.FPC.heat_load.fix(50)


def add_fpc_costing(m, costing_block=None):
    energy = m.fs.energy
    if costing_block is None:
        energy.costing = EnergyCosting()

    energy.FPC.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
    )


def add_FPC_scaling(m, blk):
    # set_scaling_factor(blk.heat_annual_scaled, 1e-8)
    # set_scaling_factor(blk.electricity_annual_scaled, 1e-2)
    set_scaling_factor(blk.heat_load, 1e-6)
    set_scaling_factor(blk.hours_storage, 1 / 10)

    # constraint_scaling_transform(blk.heat_constraint, 1e-3)
    # constraint_scaling_transform(blk.electricity_constraint, 1e-4)


def add_fpc_costing_scaling(m, blk):
    set_scaling_factor(blk.yearly_heat_production, 1e-9)
    set_scaling_factor(blk.lifetime_heat_production, 1e-9)
    set_scaling_factor(blk.aggregate_flow_heat, 1e3)


def breakdown_dof(blk):
    equalities = [c for c in activated_equalities_generator(blk)]
    active_vars = variables_in_activated_equalities_set(blk)
    fixed_active_vars = fixed_variables_in_activated_equalities_set(blk)
    unfixed_active_vars = unfixed_variables_in_activated_equalities_set(blk)
    print("\n ===============DOF Breakdown================\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(blk)}")
    print(f"Activated Variables: ({len(active_vars)})")
    for v in active_vars:
        print(f"   {v}")
    print(f"Activated Equalities: ({len(equalities)})")
    for c in equalities:
        print(f"   {c}")

    print(f"Fixed Active Vars: ({len(fixed_active_vars)})")
    for v in fixed_active_vars:
        print(f"   {v}")

    print(f"Unfixed Active Vars: ({len(unfixed_active_vars)})")
    for v in unfixed_active_vars:
        print(f"   {v}")
    print("\n")
    print(f" {f' Active Vars':<30s}{len(active_vars)}")
    print(f"{'-'}{f' Fixed Active Vars':<30s}{len(fixed_active_vars)}")
    print(f"{'-'}{f' Activated Equalities':<30s}{len(equalities)}")
    print(f"{'='}{f' Degrees of Freedom':<30s}{degrees_of_freedom(blk)}")
    print("\nSuggested Variables to Fix:")

    if degrees_of_freedom != 0:
        unfixed_vars_without_constraint = [
            v for v in active_vars if v not in unfixed_active_vars
        ]
        for v in unfixed_vars_without_constraint:
            if v.fixed is False:
                print(f"   {v}")


def report_fpc(m):
    blk = m.fs.energy.FPC

    print(f"\n\n-------------------- FPC Report --------------------\n")
    print(
        f'{"Number of collectors":<30s}{value(blk.number_collectors):<20,.2f}{pyunits.get_units(blk.number_collectors)}'
    )

    print(
        f'{"Collector area":<30s}{value(blk.collector_area_total):<20,.2f}{pyunits.get_units(blk.collector_area_total)}'
    )

    print(
        f'{"Storage volume":<30s}{value(blk.storage_volume):<20,.2f}{pyunits.get_units(blk.storage_volume)}'
    )
    print(
        f'{"Hours of Storage":<30s}{value(blk.hours_storage):<20,.2f}{pyunits.get_units(blk.hours_storage)}'
    )
    print(
        f'{"Heat load":<30s}{value(blk.heat_load):<20,.2f}{pyunits.get_units(blk.heat_load)}'
    )

    print(
        f'{"Heat annual":<30s}{value(blk.heat_annual):<20,.2f}{pyunits.get_units(blk.heat_annual)}'
    )

    print(f'{"Heat":<30s}{value(blk.heat):<20,.2f}{pyunits.get_units(blk.heat)}')

    print(
        f'{"Electricity annual":<30s}{value(blk.electricity_annual):<20,.2f}{pyunits.get_units(blk.electricity_annual)}'
    )

    print(
        f'{"Electricity":<30s}{value(blk.electricity):<20,.2f}{pyunits.get_units(blk.electricity)}'
    )


def print_FPC_costing_breakdown(m, blk):
    print(f"\n\n-------------------- FPC Costing Report --------------------\n")
    energy = m.fs.energy

    print(
        f'{"Capital Cost":<30s}{value(blk.costing.capital_cost):<20,.2f}{pyunits.get_units(blk.costing.capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(blk.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.fixed_operating_cost)}'
    )
    print("")
    print(
        f'{"Fixed Op. Cost by Capacity":<30s}{value(energy.costing.flat_plate.fixed_operating_by_capacity):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.fixed_operating_by_capacity)}'
    )
    print(
        f'{"Cost Per Collector Area":<30s}{value(energy.costing.flat_plate.cost_per_area_collector):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.cost_per_area_collector)}'
    )
    print(
        f'{"Cost Per Volume Storage":<30s}{value(energy.costing.flat_plate.cost_per_volume_storage):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.cost_per_volume_storage)}'
    )
    print(
        f'{"Cost Per Land":<30s}{value(energy.costing.flat_plate.land_cost_per_acre):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.land_cost_per_acre)}'
    )
    print(
        f'{"Land Req.":<30s}{value(blk.costing.land_area):<20,.2f}{pyunits.get_units(blk.costing.land_area)}'
    )
    print("")

    # print(
    #     f'{"Aggregated Variable Operating Cost":<30s}{value(blk.costing.aggregate_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.aggregate_variable_operating_cost)}'
    # )

    # print(
    #     f'{"Heat flow":<30s}{value(blk.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_heat)}'
    # )

    # print(
    #     f'{"Heat Cost":<30s}{value(blk.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["heat"])}'
    # )

    # print(
    #     f'{"Elec Flow":<30s}{value(blk.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_electricity)}'
    # )

    # print(
    #     f'{"Elec Cost":<30s}{value(blk.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["electricity"])}'
    # )


def solve(m, solver=None, tee=True, raise_on_failure=True, debug=False):
    # ---solving---
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            print("\n--------- TREATMENT ---------\n")
            print("\n--------- ENERGY ---------\n")
            check_jac(m.fs.energy)

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(m)
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print('\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print('\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print('\n{"=======> CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print("\n--------- FAILED SOLVE!!! ---------\n")
        print(msg)
        assert False


if __name__ == "__main__":

    solver = get_solver()
    # solver = SolverFactory("ipopt")

    m = build_system()

    build_fpc(m)
    set_fpc_op_conditions(m)
    add_FPC_scaling(m, m.fs.energy.FPC)
    init_fpc(m.fs.energy)

    add_fpc_costing(m)
    results = solve(m)
    report_fpc(m)
    print("\n")
    print(
        f'{"Heat Load FPC":<40s}{value(m.fs.energy.FPC.heat_load):<5.1f}{pyunits.get_units(m.fs.energy.FPC.heat_load)}'
    )
    print(
        f'{"Heat Flow FPC":<40s}{value(pyunits.convert(m.fs.energy.FPC.heat, to_units=pyunits.MW)):<5.1f}{pyunits.get_units(pyunits.convert(m.fs.energy.FPC.heat, to_units=pyunits.MW))}'
    )
    print(
        f'{"Load Utilization %":<40s}{(100*value(pyunits.convert(m.fs.energy.FPC.heat, to_units=pyunits.MW)))/(value(pyunits.convert(m.fs.energy.FPC.heat_load, to_units=pyunits.MW))):<5.1f}{"%"}'
    )
