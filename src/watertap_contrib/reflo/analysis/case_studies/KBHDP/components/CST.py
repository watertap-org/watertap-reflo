from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
    Block,
    Constraint,
    SolverFactory,
    Param,
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

import pickle

__all__ = [
    "build_cst",
    "init_cst",
    "set_cst_op_conditions",
    "add_cst_costing",
    "add_cst_costing_scaling",
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
        __file__ = (
            cwd + "/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/data/cst/"
        )

    dataset_filename = os.path.join(
        os.path.dirname(__file__),
        r"trough_kbhdp_heat_load_1_100_hours_storage_24_T_loop_out_300.pkl",
    )

    # Updating pickle file output column names
    with open(dataset_filename, "rb") as f:
        df = pickle.load(f)

    # Rename the columns
    df.rename(columns={"annual_energy": "heat_annual"}, inplace=True)
    df.rename(columns={"electrical_load": "electricity_annual"}, inplace=True)

    # Save the modified DataFrame back as a pickle
    with open(dataset_filename, "wb") as f:
        pickle.dump(df, f)

    surrogate_filename = os.path.join(
        os.path.dirname(__file__),
        r"trough_kbhdp_heat_load_1_100_hours_storage_24_T_loop_out_300.json",
    )

    input_bounds = dict(heat_load=[1, 50])  # , hours_storage=[23, 24])
    input_units = dict(heat_load="MW")  # , hours_storage="hour")
    input_variables = {
        "labels": ["heat_load"],  # "hours_storage"],
        "bounds": input_bounds,
        "units": input_units,
    }

    output_units = dict(
        heat_annual_scaled="kWh",
        electricity_annual_scaled="kWh",
        total_aperture_area_scaled="m**2",
    )
    output_variables = {
        "labels": [
            "heat_annual_scaled",
            "electricity_annual_scaled",
            "total_aperture_area_scaled",
        ],
        "units": output_units,
    }

    blk.unit = TroughSurrogate(
        # surrogate_model_file=surrogate_filename,
        surrogate_filename_save=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    if hasattr(blk.unit, "hours_storage"):
        print("Hours of storage is already any input parameter")
    else:
        print("Creating hours of storage parameter")
        blk.unit.hours_storage = Param(initialize=24, units=pyunits.h, mutable=True)


def init_cst(blk):
    # Fix input variables for initialization
    blk.unit.initialize()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()


def set_cst_op_conditions(blk, heat_load=10, hours_storage=6):

    if isinstance(m.fs.cst.unit.hours_storage, Param):
        blk.unit.hours_storage.set_value(hours_storage)

    if isinstance(m.fs.cst.unit.hours_storage, Var):
        blk.unit.hours_storage.fix(hours_storage)
    blk.unit.heat_load.fix(heat_load)


def add_cst_costing(blk, costing_block):
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)


def calc_costing(m, blk):
    blk.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    blk.costing.maintenance_labor_chemical_factor.fix(0)
    blk.costing.initialize()

def add_cst_costing_scaling(m,blk):
    constraint_scaling_transform(blk.costing.direct_cost_constraint, 1e-8)
    constraint_scaling_transform(blk.costing.indirect_cost_constraint, 1e-6)
    constraint_scaling_transform(blk.costing.capital_cost_constraint, 1e-8)


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
    #     f'{"LCOW":<30s}{value(blk.costing.LCOW):<20,.2f}{pyunits.get_units(blk.costing.LCOW)}'
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


if __name__ == "__main__":

    solver = get_solver()
    solver = SolverFactory("ipopt")

    m = build_system()

    build_cst(m.fs.cst)

    set_cst_op_conditions(m.fs.cst, heat_load=50, hours_storage=24)
    init_cst(m.fs.cst)

    results = solver.solve(m)
    report_cst(m, m.fs.cst.unit)

    add_cst_costing(m.fs.cst, costing_block=m.fs.costing)
    calc_costing(m, m.fs)

    add_cst_costing_scaling(m, m.fs.cst.unit)

    try:
        results = solver.solve(m)
    except:
        print_infeasible_constraints(m)

    print(degrees_of_freedom(m))

    report_cst(m, m.fs.cst.unit)
    report_cst_costing(m, m.fs)

    print(
        f'{"Elec Flow":<30s}{value(m.fs.costing.aggregate_flow_electricity):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_electricity)}'
    )

    print(
        f'{"Elec Cost":<30s}{value(m.fs.costing.aggregate_flow_costs["electricity"]):<20,.2f}{pyunits.get_units(m.fs.costing.aggregate_flow_costs["electricity"])}'
    )

    m.fs.costing.add_LCOH()
    print("LCOH:", m.fs.costing.LCOH())
    print("Hours of storage:", m.fs.cst.unit.hours_storage())
    print("Aperture area:", m.fs.cst.unit.total_aperture_area())

    print("CST fixed cost:", m.fs.cst.unit.costing.fixed_operating_cost())

    # Calcualating LCOH like SAM
    cost = m.fs.costing
    lcoh = (
        cost.total_capital_cost * cost.capital_recovery_factor
        + cost.total_operating_cost
    ) / m.fs.cst.unit.heat_annual

    print("\nManual LCOH check\n")
    print("CRF:", cost.capital_recovery_factor())
    print(
        "Numerator:",
        (
            cost.total_capital_cost * cost.capital_recovery_factor
            + cost.total_operating_cost
        )(),
    )
    print("Denominator:", m.fs.cst.unit.heat_annual())
    print("Calculated LCOH:", lcoh())
