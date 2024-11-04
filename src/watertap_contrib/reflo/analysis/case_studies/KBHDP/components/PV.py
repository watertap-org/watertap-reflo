from pyomo.environ import (
    value,
    Constraint,
    units as pyunits,
)

from idaes.core.util.model_statistics import *
import idaes.core.util.scaling as iscale
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap_contrib.reflo.solar_models.surrogate.pv import PVSurrogate

__all__ = [
    "build_pv",
    "train_pv_surrogate",
    "set_pv_constraints",
    "add_pv_scaling",
    "add_pv_costing_scaling",
    "print_PV_costing_breakdown",
    "report_PV",
]


def build_pv(m):
    energy = m.fs.energy

    energy.pv = PVSurrogate(
        surrogate_model_file="/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/solar_models/surrogate/pv/pv_surrogate.json",
        dataset_filename="/Users/zbinger/watertap-reflo/src/watertap_contrib/reflo/solar_models/surrogate/pv/data/dataset.pkl",
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

    if focus == "Size":
        energy.pv_design_constraint = Constraint(
            expr=m.fs.energy.pv.design_size
            == m.fs.treatment.costing.aggregate_flow_electricity
        )
    elif focus == "Energy":
        # energy.pv_design_constraint = Constraint(
        #     expr= m.fs.energy.pv.annual_energy == 43000000
        # )
        m.fs.energy.pv.annual_energy.fix(10000000)

        # == pyunits.convert(
        #         m.fs.treatment.costing.aggregate_flow_electricity,
        #         to_units=pyunits.kWh / pyunits.year,
        #     )

    m.fs.energy.pv.costing.land_constraint = Constraint(
        expr=m.fs.energy.pv.costing.land_area == m.fs.energy.pv.land_req
    )

    m.fs.energy.pv.costing.annual_generation = Expression(
        expr=m.fs.energy.pv.annual_energy
    )

    m.fs.energy.pv.costing.system_capacity_constraint = Constraint(
        expr= m.fs.energy.pv.costing.system_capacity == m.fs.energy.pv.design_size * 1000
    )

    m.fs.energy_balance = Expression(
        expr= 100 * (m.fs.energy.pv.annual_energy) / (m.fs.annual_treatment_energy)
    )


def add_pv_scaling(m, blk):
    pv = blk

    iscale.set_scaling_factor(pv.design_size, 1e-4)
    iscale.set_scaling_factor(pv.annual_energy, 1e-8)
    iscale.set_scaling_factor(pv.electricity, 1e-6)

def add_pv_costing_scaling(m, blk):
    pv = blk

    iscale.set_scaling_factor(m.fs.energy.pv.costing.system_capacity, 1e-5)


def print_PV_costing_breakdown(pv):
    print(f'{"PV Capital Cost":<35s}{f"${value(pv.costing.capital_cost):<25,.0f}"}')
    print(
        f'{"PV Operating Cost":<35s}{f"${value(pv.costing.fixed_operating_cost):<25,.0f}"}'
    )


def report_PV(m):
    elec = "electricity"
    print(f"\n\n-------------------- PHOTOVOLTAIC SYSTEM --------------------\n\n")
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
        f'{"Land Requirement":<30s}{value(m.fs.energy.pv.land_req):<10.1f}{pyunits.get_units(m.fs.energy.pv.land_req)}'
    )
    print(
        f'{"PV Annual Energy":<30s}{value(m.fs.energy.pv.annual_energy):<10,.0f}{pyunits.get_units(m.fs.energy.pv.annual_energy)}'
    )
    print(
        f'{"Treatment Annual Energy":<30s}{value(m.fs.annual_treatment_energy):<10,.0f}{"kWh/yr"}'
    )
    print("\n")
    print(
        f'{"PV Annual Generation":<25s}{f"{pyunits.convert(-1*m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}'
    )
    print(
        f'{"Treatment Annual Demand":<25s}{f"{pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}'
    )
    print(
        f'{"Energy Balance":<25s}{f"{value(m.fs.energy_balance):<25,.2f}"}'
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
    print(
        f'{"Electricity Sold":<30s}{f"{value(m.fs.costing.aggregate_flow_electricity_sold):<10,.0f}"}{"kW":<10s}'
    )
    print(
        f'{"Electricity Cost":<29s}{f"${value(m.fs.costing.total_electric_operating_cost):<10,.0f}"}{"$/yr":<10s}'
    )

    print(m.fs.energy.pv.annual_energy.display())
    print(m.fs.energy.pv.costing.annual_generation.display())
    print(m.fs.costing.total_electric_operating_cost.display())
