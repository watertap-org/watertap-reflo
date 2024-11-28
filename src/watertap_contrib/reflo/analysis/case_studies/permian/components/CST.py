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

from watertap_contrib.reflo.solar_models.surrogate.trough.trough_surrogate import (
    TroughSurrogate,
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

__all__ = [
    "build_cst",
    "init_cst",
    "set_cst_op_conditions",
    "add_cst_costing",
    "report_cst",
    "report_cst_costing",
]


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()

    m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)

    m.fs.cst = FlowsheetBlock(dynamic=False)

    return m


def build_cst(blk, __file__=None):

    print(f'\n{"=======> BUILDING CST SYSTEM <=======":^60}\n')

    if __file__ == None:
        cwd = os.getcwd()
        __file__ = cwd + r"\src\watertap_contrib\reflo\solar_models\surrogate\trough\\"

    dataset_filename = os.path.join(
        os.path.dirname(__file__), r"data\test_trough_data.pkl"
    )
    surrogate_filename = os.path.join(
        os.path.dirname(__file__),
        r"data\test_trough_data_heat_load_100_500_hours_storage_0_26.json",
    )

    input_bounds = dict(heat_load=[100, 500], hours_storage=[0, 26])
    input_units = dict(heat_load="MW", hours_storage="hour")
    input_variables = {
        "labels": ["heat_load", "hours_storage"],
        "bounds": input_bounds,
        "units": input_units,
    }

    output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
    output_variables = {
        "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
        "units": output_units,
    }

    blk.unit = TroughSurrogate(
        surrogate_model_file=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )


def init_cst(blk):
    # Fix input variables for initialization
    blk.unit.hours_storage.fix()
    blk.unit.heat_load.fix()
    blk.unit.initialize()

    blk.unit.heat_load.unfix()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()


def set_cst_op_conditions(blk, hours_storage=6):
    blk.unit.hours_storage.fix(hours_storage)


def add_cst_costing(blk, costing_block):
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def calc_costing(m, blk):
    blk.costing.heat_cost.set_value(0)
    blk.costing.cost_process()
    blk.costing.initialize()

    # TODO: Connect to the treatment volume
    blk.costing.add_LCOW(m.fs.system_capacity)


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

    print(
        f'{"LCOW":<30s}{value(blk.costing.LCOW):<20,.2f}{pyunits.get_units(blk.costing.LCOW)}'
    )

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

    build_cst(m.fs.cst)

    init_cst(m.fs.cst)

    set_cst_op_conditions(m.fs.cst)

    add_cst_costing(m.fs.cst, costing_block=m.fs.costing)
    calc_costing(m, m.fs)
    m.fs.costing.aggregate_flow_heat.fix(-70000)
    results = solver.solve(m)

    print(degrees_of_freedom(m))
    report_cst(m, m.fs.cst.unit)
    report_cst_costing(m, m.fs)

    # m.fs.costing.used_flows.display()
