from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
    Block,
    Constraint,
    SolverFactory,
)
import os

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver

from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.scaling import *

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from watertap_contrib.reflo.costing import (
    EnergyCosting,
)
from watertap_contrib.reflo.solar_models.zero_order import TroughPySAM

__all__ = [
    "build_cst",
    "init_cst",
    "set_cst_op_conditions",
    "add_cst_costing",
    "report_cst",
    "report_cst_costing",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.abspath(os.path.join(__location__, ".."))
dataset_filename = os.path.join(
    parent_dir,
    "data/cst/trough_data_heat_load_1_100_hours_storage_0_24.pkl",
)
surrogate_filename = dataset_filename.replace(".pkl", ".json")
def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()
    m.fs.energy = Block()

    # m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)

    # m.fs.cst = FlowsheetBlock(dynamic=False)

    return m


def build_cst_pysam(m, energy_blk=None):

    if energy_blk is None:
        energy_blk = m.fs.energy

    # energy = m.fs.energy

    print(f'\n{"=======> BUILDING CST SYSTEM <=======":^60}\n')

    # input_bounds = dict(heat_load=[1, 100], hours_storage=[0, 24])

    energy_blk.trough = TroughPySAM(solar_model_type="pysam")


def init_cst(blk, hours_storage=24, heat_load=10):
    # Fix input variables for initialization
    blk.trough.hours_storage.fix(hours_storage)
    blk.trough.heat_load.fix(heat_load)
    blk.trough.initialize()

    # blk.trough.heat_load.unfix()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()


def set_cst_op_conditions(blk, hours_storage=6):
    blk.trough.hours_storage.fix(hours_storage)


def add_cst_costing(m, costing_block=None):
    if costing_block is None:
        costing_block = m.fs.costing
    m.fs.energy.trough.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def calc_costing(m, blk):
    m.fs.costing.electricity_cost.fix(0.07)
    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    # TODO: Connect to the treatment volume
    # m.fs.costing.add_LCOH()


def report_cst(m, blk):
    # blk = m.fs.cst
    print(f"\n\n-------------------- CST Report --------------------\n")
    print("\n")

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


def report_cst_costing(m, blk):
    print(f"\n\n-------------------- CST Costing Report --------------------\n")
    print("\n")

    # print(
    #     f'{"LCOH":<30s}{value(blk.costing.LCOH):<20,.5f}{pyunits.get_units(blk.costing.LCOH)}'
    # )

    print(
        f'{"Capital Cost":<30s}{value(blk.costing.total_capital_cost):<20,.2f}{pyunits.get_units(blk.costing.total_capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(blk.costing.total_fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.total_fixed_operating_cost)}'
    )

    print(
        f'{"Variable Operating Cost":<30s}{value(blk.costing.total_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.total_variable_operating_cost)}'
    )

    print(
        f'{"Total Operating Cost":<30s}{value(blk.costing.total_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.total_operating_cost)}'
    )

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


if __name__ == "__main__":

    solver = get_solver()
    solver = SolverFactory("ipopt")

    m = build_system()

    build_cst(m)
    set_cst_op_conditions(m.fs.energy)

    init_cst(m.fs.energy, heat_load=50)


    add_cst_costing(m, costing_block=m.fs.costing)
    calc_costing(m, m.fs)
    # m.fs.costing.aggregate_flow_heat.fix(-70000)
    # m.fs.energy.trough.heat_load.unfix()
    # m.fs.energy.trough.heat_annual.fix(0.52714)
    print(f" dof = {degrees_of_freedom(m)}")
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(degrees_of_freedom(m))
    report_cst(m, m.fs.energy.trough)
    report_cst_costing(m, m.fs)
    # m.fs.
    # # m.fs.costing.display()
    m.fs.energy.trough.display()
    
    # m.fs.costing.display()
    # m.fs.costing.used_flows.display()
