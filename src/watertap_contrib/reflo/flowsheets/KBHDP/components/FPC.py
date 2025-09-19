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
from watertap_contrib.reflo.solar_models import FlatPlateSurrogate
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_FPC",
    "init_FPC",
    "set_FPC_op_conditions",
    "add_FPC_costing",
    "add_FPC_scaling",
    "report_FPC",
    "print_FPC_costing_breakdown",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
par_dir = os.path.dirname(__location__)
dataset_filename = f"{par_dir}/data/fpc/kbhdp_fpc_surrogate_data.pkl"
surrogate_filename = f"{par_dir}/data/fpc/kbhdp_fpc_surrogate.json"


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.costing = EnergyCosting()

    m.fs.fpc = FlowsheetBlock()

    return m


def build_FPC(blk):

    print(f'\n{"=======> BUILDING FPC SYSTEM <=======":^60}\n')
    input_bounds = dict(system_capacity=[0.5, 200])
    input_units = dict(system_capacity="MW")
    input_variables = {
        "labels": ["system_capacity"],
        "bounds": input_bounds,
        "units": input_units,
    }

    output_units = dict(heat_annual="kWh/year", electricity_annual="kWh/year")
    output_variables = {
        "labels": ["heat_annual", "electricity_annual"],
        "units": output_units,
    }
    fpc_dict = {
        "dataset_filename": dataset_filename,
        "input_variables": input_variables,
        "output_variables": output_variables,
        "scale_training_data": True,
        "surrogate_model_file": surrogate_filename,
        # "surrogate_filename_save": surrogate_filename,
    }

    blk.unit = FlatPlateSurrogate(**fpc_dict)


def init_FPC(blk):
    blk.unit.initialize()


def set_FPC_op_conditions(blk, system_capacity=50):

    blk.unit.system_capacity.fix(system_capacity)
    # 24 hours storage and 80C used to create surrogate, so these must be set
    blk.unit.hours_storage.set_value(24)
    blk.unit.temperature_hot.set_value(80)


def add_FPC_costing(blk, costing_block=None):

    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
    )


def add_FPC_scaling(blk):
    iscale.set_scaling_factor(blk.unit.system_capacity, 1e-6)


def report_FPC(blk, w=25):
    title = "FPC Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    unit = blk.unit

    print(f'\n{f"Parameter":<{w}}{f"Value":<{w}}{f"Units":<{w}}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"Number of collectors":<{w}s}{value(unit.number_collectors):<{w},.2f}{pyunits.get_units(unit.number_collectors)}'
    )

    print(
        f'{"Collector area":<{w}s}{value(unit.collector_area_total):<{w},.2f}{pyunits.get_units(unit.collector_area_total)}'
    )

    print(
        f'{"Storage volume":<{w}s}{value(unit.storage_volume):<{w},.2f}{pyunits.get_units(unit.storage_volume)}'
    )
    print(
        f'{"Hours of Storage":<{w}s}{value(unit.hours_storage):<{w},.2f}{pyunits.get_units(unit.hours_storage)}'
    )
    print(
        f'{"System Capacity":<{w}s}{value(unit.system_capacity):<{w},.2f}{pyunits.get_units(unit.system_capacity)}'
    )

    print(
        f'{"Heat Annual":<{w}s}{value(unit.heat_annual):<{w},.2f}{pyunits.get_units(unit.heat_annual)}'
    )

    print(
        f'{"Thermal Power":<{w}s}{value(unit.heat):<{w},.2f}{pyunits.get_units(unit.heat)}'
    )

    print(
        f'{"Electricity Annual":<{w}s}{value(unit.electricity_annual):<{w},.2f}{pyunits.get_units(unit.electricity_annual)}'
    )

    print(
        f'{"Electric Power Reqd":<{w}s}{value(unit.electricity):<{w},.2f}{pyunits.get_units(unit.electricity)}'
    )


def print_FPC_costing_breakdown(blk, w=25):
    unit = blk.unit
    title = "FPC Costing Breakdown"
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
        f'{"Land Req.":<{w}s}{value(unit.costing.land_area):<{w},.2f}{pyunits.get_units(unit.costing.land_area)}'
    )
    print("\n\n")


def main():

    m = build_system()

    build_FPC(m.fs.fpc)
    set_FPC_op_conditions(m.fs.fpc)
    assert degrees_of_freedom(m) == 0
    add_FPC_scaling(m.fs.fpc)
    init_FPC(m.fs.fpc)

    add_FPC_costing(m.fs.fpc)
    m.fs.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.add_LCOH()
    m.fs.costing.initialize()

    results = solve(m)
    assert_optimal_termination(results)

    report_FPC(m.fs.fpc)

    return m


if __name__ == "__main__":
    m = main()
