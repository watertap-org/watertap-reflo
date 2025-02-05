import os
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
    Block,
    Constraint,
    SolverFactory,
)


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
from watertap_contrib.reflo.analysis.case_studies.KBHDP.data.run_pysam_kbhdp_trough import *

__all__ = [
    "build_trough_pysam",
    "add_pysam_trough_model",
    "set_trough_pysam_op_conditions",
    "add_pysam_trough_costing",
    "run_pysam_trough",
    "get_trough_heat_load",
    "report_trough",
    "report_trough_costing",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.abspath(os.path.join(__location__, ".."))
weather_file = os.path.join(parent_dir, "data/el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(
    parent_dir, "data/cst/trough_physical_process_heat-kbhdp.json"
)


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()
    m.fs.energy = Block()

    m.fs.system_capacity = Var(initialize=6000, units=pyunits.m**3 / pyunits.day)
    m.fs.aggregate_flow_heat_treatment = Var(
        initialize=18000,  # about what is required for RPT2
        bounds=(0, None),
        units=pyunits.kilowatt,
        doc="Represents steady-state heat requirement of treatment system",
    )

    return m


def build_trough_pysam(m, energy_blk=None):

    if energy_blk is None:
        energy_blk = m.fs.energy

    print(f'\n{"=======> BUILDING CST SYSTEM <=======":^60}\n')
    energy_blk.trough = TroughPySAM(solar_model_type="pysam")
    add_pysam_trough_model(m)


def add_pysam_trough_model(m):
    config_data = read_trough_module_datafile(param_file)

    m.tech_model = setup_pysam_trough_model(
        weather_file=weather_file, config_data=config_data
    )


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()
    m.fs.aggregate_flow_heat_treatment.fix()


def set_trough_pysam_op_conditions(m, blk, hours_storage=24):
    # These are just initial values
    blk.trough.heat_load.fix(10)
    blk.trough.electricity_annual.fix(1e5)
    blk.trough.heat_annual.fix(1e5)
    
    m.hours_storage = hours_storage
    blk.trough.hours_storage.fix(hours_storage)


def add_pysam_trough_costing(m, costing_block=None):
    if costing_block is None:
        costing_block = m.fs.costing
    m.fs.energy.trough.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_block
    )
    costing_block.cost_process()

def run_pysam_trough(m, heat_load=None):

    results = run_pysam_trough_model(
        m.tech_model,
        heat_load=heat_load,
        hours_storage=m.hours_storage
    )
    return results


def get_trough_heat_load(m, heat_annual_desired, heat_load_start=5, increment_heat_load=2):
    pysam_results = run_pysam_trough(m, heat_load=heat_load_start)
    heat_annual_trough = pysam_results["annual_energy"]
    heat_load = heat_load_start
    while heat_annual_trough < heat_annual_desired:
        heat_load += increment_heat_load
        print(f"Trying Trough heat load = {heat_load:.2f} MW...")
        pysam_results = run_pysam_trough(m, heat_load=heat_load)
        heat_annual_trough = pysam_results["annual_energy"]
        if np.isnan(heat_annual_trough):
            continue
        print(f"Heat annual desired = {heat_annual_desired:.2f} kWh")
        print(f"Heat annual calculated = {heat_annual_trough:.2f} kWh")
        print(f"Calc/Desired = {heat_annual_trough/heat_annual_desired:.2f} kWh\n")
    heat_load_required = heat_load
    return heat_load_required, pysam_results

def report_trough(m, blk):
    # blk = m.fs.trough
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


def report_trough_costing(m, blk):

    print(f"\n\n-------------------- CST Costing Report --------------------\n")
    print("\n")

    # print(
    #     f'{"LCOH":<30s}{value(blk.costing.LCOH):<20,.5f}{pyunits.get_units(blk.costing.LCOH)}'
    # )

    print(
        f'{"Capital Cost":<30s}{value(blk.costing.capital_cost):<20,.2f}{pyunits.get_units(blk.costing.capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(blk.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.fixed_operating_cost)}'
    )

    print(
        f'{"Variable Operating Cost":<30s}{value(blk.costing.variable_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.variable_operating_cost)}'
    )

    print(
        f'{"Total Operating Cost":<30s}{value(m.fs.costing.total_operating_cost):<20,.2f}{pyunits.get_units(m.fs.costing.total_operating_cost)}'
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

    build_trough_pysam(m)

    set_trough_pysam_op_conditions(m)
    add_pysam_trough_costing(m, costing_block=m.fs.costing)
    set_system_op_conditions(m)
    print(f"dof = {degrees_of_freedom(m)}")

    heat_load = value(
        pyunits.convert(m.fs.aggregate_flow_heat_treatment, to_units=pyunits.MW)
    )
    
    pysam_results = run_pysam_trough(m, heat_load=heat_load)
    m.fs.energy.trough.heat_annual.fix(pysam_results["annual_energy"])
    m.fs.energy.trough.electricity_annual.fix(pysam_results["electrical_load"])
    m.fs.energy.trough.heat_load.fix(heat_load)

    heat_annual_required = value(
        pyunits.convert(
            m.fs.aggregate_flow_heat_treatment, to_units=pyunits.kWh * pyunits.year**-1
        )
    )
    heat_load_required, pysam_results = get_trough_heat_load(
        m, heat_annual_required, heat_load_start=29, increment_heat_load=0.25
    )
    m.fs.energy.trough.heat_annual.fix(pysam_results["annual_energy"])
    m.fs.energy.trough.electricity_annual.fix(pysam_results["electrical_load"])
    m.fs.energy.trough.heat_load.fix(heat_load_required)
    results = solver.solve(m)

    report_trough_costing(m, m.fs.energy.trough)
    m.fs.energy.trough.costing.display()

