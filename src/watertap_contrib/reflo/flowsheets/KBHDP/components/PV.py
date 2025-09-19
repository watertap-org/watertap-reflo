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
from watertap_contrib.reflo.solar_models import PVSurrogate
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve

__all__ = [
    "build_pv",
    "add_pv_scaling",
    "set_pv_op_conditions",
    "init_pv",
    "add_pv_costing",
    "add_pv_costing_scaling",
    "print_PV_costing_breakdown",
    "report_PV",
]
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
par_dir = os.path.dirname(__location__)
dataset_filename = f"{par_dir}/data/pv/kbhdp_pv_surrogate_data.pkl"
surrogate_filename = f"{par_dir}/data/pv/kbhdp_pv_surrogate.json"


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.costing = EnergyCosting()

    m.fs.pv = FlowsheetBlock()

    build_pv(m.fs.pv)

    return m


def build_pv(blk):
    pv_dict = {
        "input_variables": {
            "labels": ["system_capacity"],
            "units": {"system_capacity": "kW"},
        },
        "output_variables": {
            "labels": ["electricity_annual", "land_req"],
            "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
        },
        "scale_training_data": False,
        "dataset_filename": dataset_filename,
        "surrogate_model_file": surrogate_filename,
    }
    blk.unit = PVSurrogate(**pv_dict)


def set_pv_op_conditions(blk, system_capacity=1000):
    blk.unit.system_capacity.fix(system_capacity)


def add_pv_costing(blk, costing_block=None):
    if costing_block is None:
        m = blk.model()
        costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block,
    )


def add_pv_scaling(blk):

    iscale.set_scaling_factor(blk.unit.system_capacity, 1)
    iscale.set_scaling_factor(blk.unit.electricity, 1000)
    iscale.set_scaling_factor(blk.unit.land_req, 100)


def add_pv_costing_scaling(blk):
    iscale.set_scaling_factor(blk.unit.costing.system_capacity, 1e-5)


def print_PV_costing_breakdown(blk, w=30):
    title = "PV Costing Breakdown"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"PV Capital Cost":<35s}{f"${value(blk.unit.costing.capital_cost):<{w},.0f}"}{pyunits.get_units(blk.unit.costing.capital_cost)}'
    )
    print(
        f'{"PV Operating Cost":<35s}{f"${value(blk.unit.costing.fixed_operating_cost):<{w},.0f}"}{pyunits.get_units(blk.unit.costing.fixed_operating_cost)}'
    )
    print("\n\n")


def report_PV(blk, w=30):
    title = "PV System Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"PV System Capacity":<{w}s}{value(blk.unit.system_capacity):<{w}.1f}{pyunits.get_units(blk.unit.system_capacity)}'
    )
    print(
        f'{"PV Annual Energy":<{w}s}{value(blk.unit.electricity_annual):<{w},.0f}{pyunits.get_units(blk.unit.electricity_annual)}'
    )
    print(f'{"PV Power Generation":<{w}s}{f"{blk.unit.electricity():<{w},.0f}"}{"kW"}')
    print(
        f'{"Land Requirement":<{w}s}{value(blk.unit.land_req):<{w}.1f}{pyunits.get_units(blk.unit.land_req)}'
    )


def init_pv(blk):
    blk.unit.initialize()


def main():
    m = build_system()
    set_pv_op_conditions(m.fs.pv)
    assert degrees_of_freedom(m) == 0
    add_pv_scaling(m.fs.pv)
    iscale.calculate_scaling_factors(m)
    init_pv(m.fs.pv)
    add_pv_costing(m.fs.pv)
    m.fs.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.add_LCOE()
    m.fs.costing.initialize()
    results = solve(m)
    assert_optimal_termination(results)
    report_PV(m.fs.pv)
    print_PV_costing_breakdown(m.fs.pv)

    return m


if __name__ == "__main__":
    m = main()
