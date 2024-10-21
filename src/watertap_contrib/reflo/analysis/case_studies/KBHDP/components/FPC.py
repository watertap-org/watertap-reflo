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


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()

    m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)

    m.fs.fpc = FlowsheetBlock(dynamic=False)

    return m


def build_fpc(blk, __file__=None):

    print(f'\n{"=======> BUILDING FPC SYSTEM <=======":^60}\n')

    if __file__ == None:
        __file__ = r"\Users\mhardika\Documents\watertap_seto\watertap-seto\src\watertap_contrib\reflo\solar_models\surrogate\flat_plate\\"

    dataset_filename = os.path.join(
        os.path.dirname(__file__), r"data\flat_plate_data_heat_load_5_200.pkl"
    )
    surrogate_filename = os.path.join(
        os.path.dirname(__file__), r"flat_plate_surrogate.json"
    )

    input_bounds = dict(
        heat_load=[5, 200], hours_storage=[0, 26], temperature_hot=[50, 100]
    )
    input_units = dict(heat_load="MW", hours_storage="hour", temperature_hot="degK")
    input_variables = {
        "labels": ["heat_load", "hours_storage", "temperature_hot"],
        "bounds": input_bounds,
        "units": input_units,
    }

    output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
    output_variables = {
        "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
        "units": output_units,
    }
    fpc_dict = dict(
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    blk.unit = FlatPlateSurrogate(**fpc_dict)


def init_fpc(blk):
    blk.unit.initialize()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()


def set_fpc_op_conditions(blk, hours_storage=4, temperature_hot=80):

    blk.unit.hours_storage.fix(hours_storage)
    # Assumes the hot temperature to the inlet of a 'MD HX'
    blk.unit.temperature_hot.fix(temperature_hot)
    # Assumes the cold temperature from the outlet temperature of a 'MD HX'
    blk.unit.temperature_cold.set_value(20)


def add_fpc_costing(blk, costing_block):
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def calc_costing(m, blk):
    blk.costing.maintenance_labor_chemical_factor.fix(0)
    blk.costing.total_investment_factor.fix(1)
    blk.costing.cost_process()
    blk.costing.initialize()

    # TODO: Connect to the treatment volume
    blk.costing.add_LCOW(m.fs.system_capacity)


def report_fpc(m, blk):
    # blk = m.fs.fpc
    print(f"\n\n-------------------- FPC Report --------------------\n")
    print("\n")

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
        f'{"Heat load":<30s}{value(blk.heat_load):<20,.2f}{pyunits.get_units(blk.heat_load)}'
    )

    print(
        f'{"Heat annual":<30s}{value(blk.heat_annual):<20,.2f}{pyunits.get_units(blk.heat_annual)}'
    )

    print(
        f'{"Electricity annual":<30s}{value(blk.electricity_annual):<20,.2f}{pyunits.get_units(blk.electricity_annual)}'
    )


def report_fpc_costing(m, blk):
    print(f"\n\n-------------------- FPC Costing Report --------------------\n")
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

    print(
        f'{"Aggregated Variable Operating Cost":<30s}{value(blk.costing.aggregate_variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.aggregate_variable_operating_cost)}'
    )

    print(
        f'{"Heat flow":<30s}{value(blk.costing.aggregate_flow_heat):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_heat)}'
    )

    print(
        f'{"Heat Cost":<30s}{value(blk.costing.aggregate_flow_costs["heat"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["heat"])}'
    )

    print(
        f'{"Elec Flow":<30s}{value(blk.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_electricity)}'
    )

    print(
        f'{"Elec Cost":<30s}{value(blk.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(blk.costing.aggregate_flow_costs["electricity"])}'
    )


if __name__ == "__main__":

    solver = get_solver()
    solver = SolverFactory("ipopt")

    m = build_system()

    build_fpc(m.fs.fpc)
    init_fpc(m.fs.fpc)
    set_fpc_op_conditions(m.fs.fpc)

    # TODO: Connect to overall_thermal_requirement
    m.fs.fpc.unit.heat_load.fix(10)

    print("Degrees of Freedom:", degrees_of_freedom(m))

    results = solver.solve(m)
    add_fpc_costing(m.fs.fpc, costing_block=m.fs.costing)
    calc_costing(m, m.fs)

    results = solver.solve(m)

    # print("LCOW:", m.fs.costing.LCOW())

    report_fpc(m, m.fs.fpc.unit)
    report_fpc_costing(m, m.fs)
