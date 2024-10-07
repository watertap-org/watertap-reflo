import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    TransformationFactory,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import *

# Flowsheet function imports
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet import (
    build_vagmd_flowsheet,
    fix_dof_and_initialize,
)

from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_design_model import (
    get_n_time_points,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.MD import *


# TODO: Add translator block after feed block
# TODO: Add opex calculator at each time step and capex on the last time step

__all__ = [
    "build_system",
    "add_connections",
    "add_constraints",
    "add_costing",
    "relax_constraints",
    "set_inlet_conditions",
    "set_operating_conditions",
    "report_MCAS_stream_conc",
    "display_system_stream_table",
    "display_system_build",
    "init_system",
]


def propagate_state(arc):
    _prop_state(arc)


def main():

    m, model_options, n_time_points = build_system()
    add_connections(m)

    # add_constraints(m)

    set_operating_conditions(m)
    init_system(m)
    add_costing(m)
    solve(m)


    display_costing_breakdown(m)


def build_system():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # m.fs.costing = REFLOCosting()
    # m.fs.costing.base_currency = pyunits.USD_2020

    # Property package
    m.fs.params = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.params)
    m.fs.product = Product(property_package=m.fs.params)
    m.fs.disposal = Product(property_package=m.fs.params)

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    model_options, n_time_points= build_md(m, m.fs.md)
    # m.fs.dwi = FlowsheetBlock(dynamic=False)

    return m, model_options, n_time_points


def add_connections(m):

    m.fs.feed_to_md = Arc(
        source = m.fs.feed.outlet,
        destination = m.fs.md.feed.inlet
    )


    m.fs.md_to_product = Arc(
        source = m.fs.md.permeate.outlet,
        destination = m.fs.product.inlet
    )

    m.fs.md_to_disposal = Arc(
        source = m.fs.md.concentrate.outlet,
        destination = m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def add_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.feed_flow_vol = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.L / pyunits.s,
        doc="System Feed Flowrate",
    )

    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.product.properties[0].conc_mass_phase_comp
    m.fs.disposal.properties[0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_in[0.0].conc_mass_phase_comp
    m.fs.MCAS_to_NaCl_translator.properties_out[0.0].conc_mass_phase_comp


# def add_costing(m):

#     m.fs.costing.cost_process()
#     m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
#     m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)


def set_inlet_conditions(
    blk, model_options
):

    print(f'\n{"=======> SETTING FEED CONDITIONS <=======":^60}\n')
    feed_flow_rate = model_options["feed_flow_rate"]
    feed_salinity = model_options["feed_salinity"]
    feed_temp = model_options["feed_temp"]

    blk.feed.properties.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    feed_flow_rate* pyunits.L / pyunits.h,
                    to_units=pyunits.m**3 / pyunits.s,
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
                ("temperature", None): feed_temp + 273.15,
                ("pressure", None): 101325,
            },
            hold_state = True
        )


def set_operating_conditions(m, model_options):
    
    set_inlet_conditions(m.fs, model_options)
    init_md(m.fs.md,  model_options,n_time_points, verbose=True, solver=None)
    set_md_op_conditions(m.fs.md)


def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize()






def solve(m, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print_infeasible_bounds(m)
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print(msg)
        return results



if __name__ == "__main__":
    main()