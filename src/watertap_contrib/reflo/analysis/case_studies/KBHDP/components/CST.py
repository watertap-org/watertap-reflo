import os
from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    assert_optimal_termination,
    units as pyunits,
)

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.solar_models import TroughSurrogate


__all__ = [
    "build_cst",
    "init_cst",
    "set_cst_op_conditions",
    "add_cst_costing",
    "add_cst_costing_scaling",
    "report_cst",
    "report_cst_costing",
]


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
par_dir = os.path.dirname(__location__)
dataset_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate_data.pkl"
surrogate_filename = f"{par_dir}/data/cst/kbhdp_cst_surrogate.json"


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()

    m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)

    m.fs.cst = FlowsheetBlock(dynamic=False)

    return m


def build_cst(blk):

    print(f'\n{"=======> BUILDING CST SYSTEM <=======":^60}\n')

    # Create input_variables for configuring surrogate model
    input_units = dict(system_capacity="MW")
    input_bounds = dict(system_capacity=[1, 50])
    input_variables = {
        "labels": ["system_capacity"],
        "units": input_units,
        "bounds": input_bounds,
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
        # surrogate_filename_save=dataset_filename,
        surrogate_model_file=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    blk.unit = TroughSurrogate(**trough_dict)


def init_cst(blk):
    # Fix input variables for initialization
    blk.unit.initialize()


def set_cst_op_conditions(blk, heat_load=10):

    blk.unit.system_capacity.fix(heat_load)
    # 24 hours storage used to create surrogate, so this must be fixed
    blk.unit.hours_storage.set_value(24)


def add_cst_costing(blk, costing_block):
    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_block)
    m = blk.model()
    m.fs.costing.cost_process()
    # Updated to be 0 because this factor is not included in SAM
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.initialize()


def add_cst_costing_scaling(m, blk):
    iscale.constraint_scaling_transform(blk.costing.direct_cost_constraint, 1e-8)
    iscale.constraint_scaling_transform(blk.costing.indirect_cost_constraint, 1e-6)
    iscale.constraint_scaling_transform(blk.costing.capital_cost_constraint, 1e-8)


def report_cst(m, blk):
    
    print(f"\n\n-------------------- CST Report --------------------\n\n")

    print(
        f'{"System Capacity":<30s}{value(blk.system_capacity):<20,.2f}{pyunits.get_units(blk.system_capacity)}'
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
    print(f"\n\n-------------------- CST Costing Report --------------------\n\n")

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


def main():

    solver = get_solver()

    m = build_system()

    build_cst(m.fs.cst)

    set_cst_op_conditions(m.fs.cst, heat_load=50)

    assert degrees_of_freedom(m) == 0

    init_cst(m.fs.cst)

    results = solver.solve(m)
    assert_optimal_termination(results)

    report_cst(m, m.fs.cst.unit)

    add_cst_costing(m.fs.cst, costing_block=m.fs.costing)

    add_cst_costing_scaling(m, m.fs.cst.unit)
    
    results = solver.solve(m)
    assert_optimal_termination(results)

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



if __name__ == "__main__":
    main()