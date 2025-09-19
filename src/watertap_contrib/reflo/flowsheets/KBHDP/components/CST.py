import os
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.solar_models import TroughSurrogate
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve


__all__ = [
    "build_CST",
    "init_CST",
    "set_CST_op_conditions",
    "add_CST_costing",
    "add_CST_costing_scaling",
    "report_CST",
    "print_CST_costing_breakdown",
]


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
par_dir = os.path.dirname(__location__)
dataset_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate_data.pkl"
surrogate_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate.json"
# dataset_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate_data-unscaled.pkl"
# surrogate_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate-unscaled.json"


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.costing = EnergyCosting()

    m.fs.cst = FlowsheetBlock()

    return m


def build_CST(blk):

    print(f'\n{"=======> BUILDING CST SYSTEM <=======":^60}\n')

    # Create input_variables for configuring surrogate model
    input_units = dict(system_capacity="MW")
    input_variables = {
        "labels": ["system_capacity"],
        "units": input_units,
    }

    # Create output_variables for configuring surrogate model
    output_units = dict(
        heat_annual="kWh/year",
        electricity_annual="kWh/year",
        total_aperture_area="m**2",
    )
    output_variables = {
        "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
        "units": output_units,
    }

    # Create surrogate model configuration dictionary
    trough_dict = dict(
        # surrogate_filename_save=surrogate_filename,
        surrogate_model_file=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    blk.unit = TroughSurrogate(**trough_dict)


def init_CST(blk):
    # Fix input variables for initialization
    blk.unit.initialize()


def set_CST_op_conditions(blk, system_capacity=50):

    blk.unit.system_capacity.fix(system_capacity)
    # 24 hours storage and 300C temperature loop used to create surrogate, so this must be fixed
    blk.unit.hours_storage.set_value(24)
    blk.unit.temperature_loop.set_value(300)


def add_CST_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def add_CST_costing_scaling(blk):
    unit = blk.unit
    iscale.constraint_scaling_transform(unit.costing.direct_cost_constraint, 1e-8)
    iscale.constraint_scaling_transform(unit.costing.indirect_cost_constraint, 1e-6)
    iscale.constraint_scaling_transform(unit.costing.capital_cost_constraint, 1e-8)


def report_CST(blk, w=30):
    unit = blk.unit
    title = "CST Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")

    print(
        f'{"System Capacity":<{w}s}{value(unit.system_capacity):<{w},.2f}{pyunits.get_units(unit.system_capacity)}'
    )

    print(
        f'{"Heat annual":<{w}s}{value(unit.heat_annual):<{w},.2f}{pyunits.get_units(unit.heat_annual)}'
    )

    print(
        f'{"Thermal Power Produced":<{w}s}{value(unit.heat):<{w},.2f}{pyunits.get_units(unit.heat)}'
    )

    print(
        f'{"Electricity annual":<{w}s}{value(unit.electricity_annual):<{w},.2f}{pyunits.get_units(unit.electricity_annual)}'
    )

    print(
        f'{"Electric Power Consumed":<{w}s}{value(unit.electricity):<{w},.2f}{pyunits.get_units(unit.electricity)}'
    )


def print_CST_costing_breakdown(blk, w=30):
    unit = blk.unit
    title = "CST Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")

    print(
        f'{"Capital Cost":<{w}s}{value(unit.costing.capital_cost):<{w},.2f}{pyunits.get_units(unit.costing.capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<{w}s}{value(unit.costing.fixed_operating_cost):<{w},.2f}{pyunits.get_units(unit.costing.fixed_operating_cost)}'
    )

    print(
        f'{"Variable Operating Cost":<{w}s}{value(unit.costing.variable_operating_cost):<{w},.2f}{pyunits.get_units(unit.costing.variable_operating_cost)}'
    )
    print("\n\n")


def main():

    m = build_system()
    build_CST(m.fs.cst)
    set_CST_op_conditions(m.fs.cst)
    assert degrees_of_freedom(m) == 0
    init_CST(m.fs.cst)

    results = solve(m)
    assert_optimal_termination(results)

    add_CST_costing(m.fs.cst)
    m.fs.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.add_LCOH()
    m.fs.costing.initialize()

    from idaes.core.solvers import get_solver

    solver = get_solver()
    results = solve(m, solver=solver)
    assert_optimal_termination(results)

    report_CST(m.fs.cst)
    print_CST_costing_breakdown(m.fs.cst)

    return m


if __name__ == "__main__":
    m = main()
