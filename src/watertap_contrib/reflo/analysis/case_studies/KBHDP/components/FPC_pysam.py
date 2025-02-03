import os
import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    Block,
    Constraint,
    SolverFactory,
    check_optimal_termination,
    units as pyunits,
)

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver

from idaes.core.util.scaling import *
from idaes.core.util.model_statistics import *

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from watertap.core.util.model_diagnostics.infeasible import *
from watertap_contrib.reflo.costing import (
    EnergyCosting,
)
from watertap_contrib.reflo.solar_models.zero_order import FlatPlatePySAM
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    check_jac,
    calc_scale,
)

from watertap_contrib.reflo.analysis.case_studies.KBHDP.data import *

__all__ = [
    "build_fpc_pysam",
    "add_pysam_fpc_model",
    "run_pysam_fpc",
    "get_fpc_heat_load",
    "init_fpc",
    "set_fpc_pysam_op_conditions",
    "add_fpc_pysam_costing",
    "report_fpc",
    "print_FPC_costing_breakdown",
]
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
parent_dir = os.path.abspath(os.path.join(__location__, ".."))
weather_file = os.path.join(parent_dir, "data/el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(parent_dir, "data/fpc/solar_water_heating-kbhdp.json")


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

    build_fpc_pysam(m)

    return m


def build_fpc_pysam(m):
    m.fs.energy.FPC = FlatPlatePySAM(solar_model_type="pysam")
    add_pysam_fpc_model(m)
    # run_pysam_fpc_model(tech_model)


def add_pysam_fpc_model(m):

    temperatures = {
        "T_cold": 20,
        "T_hot": 70,  # this will be overwritten by temperature_hot value
        "T_amb": 18,
    }

    config_data = read_module_datafile(param_file)

    if "solar_resource_file" in config_data:
        del config_data["solar_resource_file"]

    m.tech_model = setup_pysam_fpc_model(
        temperatures=temperatures,
        weather_file=weather_file,
        config_data=config_data,
    )


def run_pysam_fpc(m, heat_load=60):

    results = run_pysam_fpc_model(
        m.tech_model,
        heat_load_mwt=heat_load,
        hours_storage=m.hours_storage,
        temperature_hot=m.temperature_hot,
    )

    return results


def get_fpc_heat_load(m, heat_annual_desired, heat_load_start=5, increment_heat_load=1):
    """
    Use to get the heat_load input required to
    get the desired annual heat output (heat_annual)
    """
    pysam_results = run_pysam_fpc(m, heat_load=heat_load_start)
    heat_annual_fpc = pysam_results["heat_annual"]
    heat_load = heat_load_start
    while heat_annual_fpc < heat_annual_desired:
        heat_load += increment_heat_load
        print(f"Trying FPC heat load = {heat_load:.2f} MW...")
        pysam_results = run_pysam_fpc(m, heat_load=heat_load)
        heat_annual_fpc = pysam_results["heat_annual"]
        print(f"Heat annual desired = {heat_annual_desired:.2f} kWh")
        print(f"Heat annual calculated = {heat_annual_fpc:.2f} kWh")
        print(f"Calc/Desired = {heat_annual_fpc/heat_annual_desired:.2f} kWh\n")
    heat_load_required = heat_load
    return heat_load_required, pysam_results


def init_fpc(blk):
    blk.FPC.initialize()


def set_system_op_conditions(m):
    m.fs.system_capacity.fix()
    m.fs.aggregate_flow_heat_treatment.fix()


def set_fpc_pysam_op_conditions(m, hours_storage=12, temperature_hot=80):
    energy = m.fs.energy
    energy.FPC.heat_load.fix(10)
    energy.FPC.electricity_annual.fix(1e5)
    energy.FPC.heat_annual.fix(1e5)
    # Assumes the cold temperature from the outlet temperature of a 'MD HX'
    energy.FPC.temperature_cold.set_value(20)
    # Assumes the hot temperature to the inlet of a 'MD HX'
    energy.FPC.temperature_hot.set_value(temperature_hot)
    energy.FPC.hours_storage.set_value(hours_storage)
    # NOTE: if you try to pass evaluated Pyomo objects
    # into PySAM (e.g., value(m.fs.energy.FPC.temperature_hot)) you will get errors.
    # Since we need these parameter values to both calculate costing and run PySAM,
    # they are set here in both ways
    m.hours_storage = hours_storage
    m.temperature_hot = temperature_hot


def add_fpc_pysam_costing(m, costing_block=None):
    energy = m.fs.energy
    if costing_block is None:
        energy.costing = EnergyCosting()

    energy.FPC.costing = UnitModelCostingBlock(
        flowsheet_costing_block=energy.costing,
    )


def report_fpc(m):
    blk = m.fs.energy.FPC

    print(f"\n\n-------------------- FPC Report --------------------\n")
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

    print(f'{"Heat":<30s}{value(blk.heat):<20,.2f}{pyunits.get_units(blk.heat)}')

    print(
        f'{"Electricity annual":<30s}{value(blk.electricity_annual):<20,.2f}{pyunits.get_units(blk.electricity_annual)}'
    )

    print(
        f'{"Electricity":<30s}{value(blk.electricity):<20,.2f}{pyunits.get_units(blk.electricity)}'
    )


def print_FPC_costing_breakdown(m, blk):
    print(f"\n\n-------------------- FPC Costing Report --------------------\n")
    energy = m.fs.energy

    print(
        f'{"Capital Cost":<30s}{value(blk.costing.capital_cost):<20,.2f}{pyunits.get_units(blk.costing.capital_cost)}'
    )

    print(
        f'{"Fixed Operating Cost":<30s}{value(blk.costing.fixed_operating_cost):<20,.2f}{pyunits.get_units(blk.costing.fixed_operating_cost)}'
    )
    print("")
    print(
        f'{"Fixed Op. Cost by Capacity":<30s}{value(energy.costing.flat_plate.fixed_operating_by_capacity):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.fixed_operating_by_capacity)}'
    )
    print(
        f'{"Cost Per Collector Area":<30s}{value(energy.costing.flat_plate.cost_per_area_collector):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.cost_per_area_collector)}'
    )
    print(
        f'{"Cost Per Volume Storage":<30s}{value(energy.costing.flat_plate.cost_per_volume_storage):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.cost_per_volume_storage)}'
    )
    print(
        f'{"Cost Per Land":<30s}{value(energy.costing.flat_plate.land_cost_per_acre):<20,.2f}{pyunits.get_units(energy.costing.flat_plate.land_cost_per_acre)}'
    )
    print(
        f'{"Land Req.":<30s}{value(blk.costing.land_area):<20,.2f}{pyunits.get_units(blk.costing.land_area)}'
    )
    print("")

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
            print("\n--------- TREATMENT ---------\n")
            print("\n--------- ENERGY ---------\n")
            check_jac(m.fs.energy)

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

    solver = get_solver()

    m = build_system()

    set_fpc_pysam_op_conditions(m)
    add_fpc_pysam_costing(m)
    set_system_op_conditions(m)
    print(f"dof = {degrees_of_freedom(m)}")

    heat_load = value(
        pyunits.convert(m.fs.aggregate_flow_heat_treatment, to_units=pyunits.MW)
    )

    pysam_results = run_pysam_fpc(m, heat_load=heat_load)
    m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
    m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
    m.fs.energy.FPC.heat_load.fix(heat_load)
    print(f"dof = {degrees_of_freedom(m)}")
    results = solve(m)
    print_FPC_costing_breakdown(m, m.fs.energy.FPC)

    heat_annual_required = value(
        pyunits.convert(
            m.fs.aggregate_flow_heat_treatment, to_units=pyunits.kWh * pyunits.year**-1
        )
    )

    heat_load_required, pysam_results = get_fpc_heat_load(
        m, heat_annual_required, heat_load_start=100, increment_heat_load=2
    )
    print(heat_annual_required, pysam_results["heat_annual"], heat_load_required)
    m.fs.energy.FPC.heat_annual.fix(pysam_results["heat_annual"])
    m.fs.energy.FPC.electricity_annual.fix(pysam_results["electricity_annual"])
    m.fs.energy.FPC.heat_load.fix(heat_load_required)
    results = solve(m)
    print_FPC_costing_breakdown(m, m.fs.energy.FPC)
