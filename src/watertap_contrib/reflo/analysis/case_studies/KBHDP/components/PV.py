import os
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import *
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap_contrib.reflo.solar_models.surrogate.pv import PVSurrogate
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
    REFLOSystemCosting,
)
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    check_jac,
    calc_scale,
)

__all__ = [
    "build_pv",
    "train_pv_surrogate",
    "set_pv_constraints",
    "add_pv_scaling",
    "add_pv_costing_scaling",
    "print_PV_costing_breakdown",
    "report_PV",
]


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    energy = m.fs.energy = Block()

    print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
    return m


def build_pv(m):
    energy = m.fs.energy

    parent_dir = os.path.abspath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "..", "..")
    )

    surrogate_dir = os.path.join(
        parent_dir,
        "solar_models",
        "surrogate",
        "pv",
    )

    dataset_filename = os.path.join(surrogate_dir, "data", "dataset.pkl")

    surrogate_filename = os.path.join(
        surrogate_dir,
        "pv_surrogate.json",
    )

    energy.pv = PVSurrogate(
        surrogate_model_file=surrogate_filename,
        dataset_filename=dataset_filename,
        input_variables={
            "labels": ["design_size"],
            "bounds": {"design_size": [1, 200000]},
            "units": {"design_size": "kW"},
        },
        output_variables={
            "labels": ["annual_energy", "land_req"],
            "units": {"annual_energy": "kWh", "land_req": "acre"},
        },
        scale_training_data=False,
    )


def train_pv_surrogate(m):
    energy = m.fs.energy

    energy.pv.create_rbf_surrogate()

    assert False


def set_pv_constraints(m, focus="Size"):
    energy = m.fs.energy
    m.fs.energy.pv.load_surrogate()

    m.fs.energy.pv.heat.fix(0)

    if focus == "Size":
        m.fs.energy.pv.design_size.fix(500)
    elif focus == "Energy":
        m.fs.energy.pv.annual_energy.fix(10000000)

    m.fs.energy.pv.initialize()


def add_pv_costing(m, blk):
    energy = m.fs.energy
    energy.costing = EnergyCosting()

    energy.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
        costing_method_arguments={
            "cost_method": "simple"
        }
    )


def add_pv_scaling(m, blk):
    pv = blk

    # print(calc_scale(value(pv.annual_energy)))

    iscale.set_scaling_factor(pv.design_size, 1000)
    # iscale.set_scaling_factor(pv.annual_energy, 1)
    iscale.set_scaling_factor(pv.electricity, 1000)
    iscale.set_scaling_factor(pv.land_req, 100)


def add_pv_costing_scaling(m, blk):

    # iscale.set_scaling_factor(blk.system_capacity, 1e-5)
    iscale.set_scaling_factor(blk.yearly_electricity_production, 1e7)
    iscale.set_scaling_factor(blk.lifetime_electricity_production, 1e8)
    iscale.set_scaling_factor(blk.aggregate_flow_electricity, 1e3)
    # iscale.constraint_scaling_transform(blk.lifetime_electricity_production_constraint, 1e2)


def print_PV_costing_breakdown(pv):
    print(f"\n\n-------------------- PV Costing Breakdown --------------------\n")
    print(f'{"PV Capital Cost":<35s}{f"${value(pv.costing.capital_cost):<25,.0f}"}')
    print(
        f'{"PV Operating Cost":<35s}{f"${value(pv.costing.fixed_operating_cost):<25,.0f}"}'
    )


def report_PV(m):
    elec = "electricity"
    print(f"\n\n-------------------- PHOTOVOLTAIC SYSTEM --------------------\n\n")
    print(
        f'{"Land Requirement":<30s}{value(m.fs.energy.pv.land_req):<10.1f}{pyunits.get_units(m.fs.energy.pv.land_req)}'
    )
    print(
        f'{"System Agg. Flow Electricity":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}'
    )
    print(
        f'{"PV Agg. Flow Elec.":<30s}{value(m.fs.energy.pv.design_size):<10.1f}{pyunits.get_units(m.fs.energy.pv.design_size)}'
    )
    print(
        f'{"Treatment Agg. Flow Elec.":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}'
    )

    print(
        f'{"PV Annual Energy":<30s}{value(m.fs.energy.pv.annual_energy):<10,.0f}{pyunits.get_units(m.fs.energy.pv.annual_energy)}'
    )
    # print(
    #     f'{"Treatment Annual Energy":<30s}{value(m.fs.annual_treatment_energy):<10,.0f}{"kWh/yr"}'
    # )
    print("\n")
    print(
        f'{"PV Annual Generation":<25s}{f"{pyunits.convert(m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}'
    )
    print(
        f'{"Treatment Annual Demand":<25s}{f"{pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}'
    )
    print(
        f'{"Grid Electricity Frac":<25s}{f"{100*value(m.fs.costing.frac_elec_from_grid):<25,.3f} %"}'
    )
    print(
        f'{"Treatment Elec Cost":<25s}{f"${value(m.fs.treatment.costing.aggregate_flow_costs[elec]):<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"Energy Elec Cost":<25s}{f"${value(m.fs.energy.costing.aggregate_flow_costs[elec]):<25,.0f}"}{"$/yr":<10s}'
    )
    print("\nEnergy Balance")
    print(
        f'{"Treatment Agg. Flow Elec.":<30s}{value(m.fs.treatment.costing.aggregate_flow_electricity):<10.1f}{"kW"}'
    )
    print(
        f'{"PV Agg. Flow Elec.":<30s}{value(m.fs.energy.costing.aggregate_flow_electricity):<10.1f}{"kW"}'
    )
    print(
        f'{"Electricity Buy":<30s}{f"{value(m.fs.costing.aggregate_flow_electricity_purchased):<10,.0f}"}{"kW":<10s}'
    )
    # print(
    #     f'{"Electricity Sold":<30s}{f"{value(m.fs.costing.aggregate_flow_electricity_sold):<10,.0f}"}{"kW":<10s}'
    # )
    print(
        f'{"Electricity Cost":<29s}{f"${value(m.fs.costing.total_electric_operating_cost):<10,.0f}"}{"$/yr":<10s}'
    )

    # print(m.fs.energy.pv.annual_energy.display())
    # # print(m.fs.energy.pv.costing.annual_generation.display())
    # print(m.fs.costing.total_electric_operating_cost.display())


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


def initialize(m):
    energy = m.fs.energy

    energy.costing.cost_process()
    energy.costing.initialize()


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
            check_jac(m)

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
    m = build_system()
    build_pv(m)
    set_pv_constraints(m, focus="Energy")
    solve(m, debug=True)
    add_pv_costing(m, m.fs.energy.pv)
    add_pv_scaling(m, m.fs.energy.pv)
    iscale.calculate_scaling_factors(m)
    initialize(m)
    solve(m, debug=True)
    print(
        f"{f'Design Size (W):':<30s}{value(pyunits.convert(m.fs.energy.pv.design_size, to_units=pyunits.watt)):<10,.1f}"
    )
    if m.fs.energy.pv.costing.cost_method == "simple":
        print(
            f"{f'Cost Per Watt ($/W):':<30s}{value(m.fs.energy.costing.pv_surrogate.cost_per_watt_installed):<10,.1f}"
        )
    else:
        print(
            f"{f'Cost Per Watt Module ($/W):':<30s}{value(m.fs.energy.costing.pv_surrogate.cost_per_watt_module):<10,.1f}"
        )
    # print(f"{f'Direct Cost Should Be ($):':<30s}{value(pyunits.convert(m.fs.energy.pv.design_size, to_units=pyunits.watt))*value(m.fs.energy.costing.pv_surrogate.cost_per_watt_module):<10,.1f}")
    print(
        f"{f'Direct Cost Currently Is ($):':<30s}{value(m.fs.energy.pv.costing.capital_cost):<10,.1f}"
    )

    # print(m.fs.energy.pv.costing.direct_capital_cost_constraint.pprint())
    # print(m.fs.energy.pv.design_size())
    # print(m.fs.energy.costing.pv_surrogate.cost_per_watt_module())