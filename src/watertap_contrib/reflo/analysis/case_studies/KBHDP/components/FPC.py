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

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.solar_models import FlatPlateSurrogate

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
par_dir = os.path.dirname(__location__)
dataset_filename = f"{par_dir}/data/fpc/kbhdp_fpc_surrogate_data.pkl"
surrogate_filename = f"{par_dir}/data/fpc/kbhdp_fpc_surrogate.json"


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.costing = EnergyCosting()

    m.fs.fpc = FlowsheetBlock()

    return m


def build_fpc(blk):
    # energy = m.fs.energy

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


def init_fpc(blk):
    blk.unit.initialize()


def set_fpc_op_conditions(blk, system_capacity=50):

    blk.unit.system_capacity.fix(system_capacity)
    # 24 hours storage and 80C used to create surrogate, so these must be set
    blk.unit.hours_storage.set_value(24)
    blk.unit.temperature_hot.set_value(80)


def add_fpc_costing(blk, costing_block=None):

    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
    )
    m.fs.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.initialize()


def add_FPC_scaling(blk):
    iscale.set_scaling_factor(blk.unit.system_capacity, 1e-6)


def report_fpc(blk):
    m = blk.model()
    unit = blk.unit

    print(f"\n\n-------------------- FPC Report --------------------\n")
    print(
        f'{"Number of collectors":<30s}{value(unit.number_collectors):<20,.2f}{pyunits.get_units(unit.number_collectors)}'
    )

    print(
        f'{"Collector area":<30s}{value(unit.collector_area_total):<20,.2f}{pyunits.get_units(unit.collector_area_total)}'
    )

    print(
        f'{"Storage volume":<30s}{value(unit.storage_volume):<20,.2f}{pyunits.get_units(unit.storage_volume)}'
    )
    print(
        f'{"Hours of Storage":<30s}{value(unit.hours_storage):<20,.2f}{pyunits.get_units(unit.hours_storage)}'
    )
    print(
        f'{"System Capacity":<30s}{value(unit.system_capacity):<20,.2f}{pyunits.get_units(unit.system_capacity)}'
    )

    print(
        f'{"Heat annual":<30s}{value(unit.heat_annual):<20,.2f}{pyunits.get_units(unit.heat_annual)}'
    )

    print(f'{"Heat":<30s}{value(unit.heat):<20,.2f}{pyunits.get_units(unit.heat)}')

    print(
        f'{"Electricity annual":<30s}{value(unit.electricity_annual):<20,.2f}{pyunits.get_units(unit.electricity_annual)}'
    )

    print(
        f'{"Electricity":<30s}{value(unit.electricity):<20,.2f}{pyunits.get_units(unit.electricity)}'
    )

    if hasattr(unit, "costing"):
        print(f"\n\n-------------------- FPC Costing Report --------------------\n")

        print(
            f'{"Capital Cost":<30s}{value(unit.costing.capital_cost):<20,.2f}{pyunits.get_units(unit.costing.capital_cost)}'
        )

        print(
            f'{"Fixed Operating Cost":<30s}{value(unit.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(unit.costing.fixed_operating_cost)}'
        )
        print("")
        print(
            f'{"Fixed Op. Cost by Capacity":<30s}{value(m.fs.costing.flat_plate.fixed_operating_by_capacity):<20,.2f}{pyunits.get_units(m.fs.costing.flat_plate.fixed_operating_by_capacity)}'
        )
        print(
            f'{"Cost Per Collector Area":<30s}{value(m.fs.costing.flat_plate.cost_per_area_collector):<20,.2f}{pyunits.get_units(m.fs.costing.flat_plate.cost_per_area_collector)}'
        )
        print(
            f'{"Cost Per Volume Storage":<30s}{value(m.fs.costing.flat_plate.cost_per_volume_storage):<20,.2f}{pyunits.get_units(m.fs.costing.flat_plate.cost_per_volume_storage)}'
        )
        print(
            f'{"Cost Per Land":<30s}{value(m.fs.costing.land_cost):<20,.2f}{pyunits.get_units(m.fs.costing.land_cost)}'
        )
        print(
            f'{"Land Req.":<30s}{value(unit.costing.land_area):<20,.2f}{pyunits.get_units(unit.costing.land_area)}'
        )
        print("")

        print(
            f'{"Thermal Power":<30s}{value(m.fs.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_heat)}'
        )

        print(
            f'{"Electric Power Reqd":<30s}{value(m.fs.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_electricity)}'
        )

        print(
            f'{"Electricity Cost":<30s}{value(m.fs.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_costs["electricity"])}'
        )


def main():

    solver = get_solver()
    m = build_system()

    build_fpc(m.fs.fpc)
    set_fpc_op_conditions(m.fs.fpc)
    assert degrees_of_freedom(m) == 0
    add_FPC_scaling(m.fs.fpc)
    init_fpc(m.fs.fpc)

    add_fpc_costing(m.fs.fpc)
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_fpc(m.fs.fpc)
    m.fs.fpc.unit.heat_annual_scaling.display()
    print(
        pyunits.get_units(
            m.fs.fpc.unit.heat_annual_scaled / m.fs.fpc.unit.heat_annual_scaling
        )
    )
    print(pyunits.get_units(m.fs.fpc.unit.heat_annual))


if __name__ == "__main__":
    main()
