import pandas as pd
import numpy as np
import time
import pathlib
from collections import defaultdict
from functools import partial
from pyomo.environ import (
    ConcreteModel,
    value,
    Any,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core import MaterialFlowBasis
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import (
    Product,
    Feed,
    StateJunction,
    Separator,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *
from idaes.core.util.initialization import propagate_state

from watertap.core.solvers import get_solver
from watertap.core import Database
from watertap_contrib.reflo.core.wt_reflo_database import REFLODatabase
from watertap.core.zero_order_properties import WaterParameterBlock as ZO
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock as MCAS,
)

# from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOCosting,
)
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
from watertap_contrib.reflo.analysis.case_studies.permian import *
from watertap_contrib.reflo.unit_models.deep_well_injection import DeepWellInjection
from watertap_contrib.reflo.analysis.case_studies.permian.utils.results_dict import (
    build_results_dict,
    results_dict_append,
)

skips = [
    "diffus_phase",
    "diffus_param",
    "dens_mass_param",
    "dh_vap_w_param",
    "cp_phase_param",
    "pressure_sat_param_psatw",
    "enth_mass_param",
    "osm_coeff_param",
    "visc_d_param",
    "therm_cond_phase_param",
    "pressure_sat_param",
    "bpe_",
    "TIC",
    "TPEC",
    "blocks[",
    "yearly_heat_production",
    "yearly_electricity_production",
    "cp_param_NaCl_liq",
    "_translator",
    "permeate_side",
    "properties_interface",
    "material_flow_dx",
    "._flow_terms",
    "pressure_dx",
    "MCAS_properties",
    "cp_param_NaCl_solid",
    "cp_vap_param",
    "temp_sat_solvent",
    "m_mec.fs.properties",
    "cp_mass_phase",
]

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1125 * pyunits.kg / pyunits.m**3
rho_water = 995 * pyunits.kg / pyunits.m**3
electricity_cost_base = 0.0434618999  # USD_2018/kWh equivalent to 0.0575 USD_2023/kWh
heat_cost_base = 0.00894

solver = get_solver()

__all__ = ["build_and_run_permian_SOA", "solve_permian_SOA"]


def build_and_run_permian_SOA(
    pretreatment_recovery=0.98,
    recovery=0.5,
    Qin=5,
    tds=130,
    electricity_cost=electricity_cost_base,
    heat_cost=heat_cost_base,
    **kwargs,
):
    m_pre = build_and_run_permian_pretreatment(
        Qin=Qin,
        tds=tds,
        electricity_cost=electricity_cost,
        heat_cost=heat_cost,
        **kwargs,
    )


    flow_to_mvc = Qin * value(m_pre.fs.treatment.EC.unit.recovery_frac_mass_H2O[0])
    flow_to_mvc = flow_to_mvc * value(
        m_pre.fs.treatment.cart_filt.unit.recovery_frac_mass_H2O[0]
    )
    tds_to_mvc = value(
        pyunits.convert(
            m_pre.fs.treatment.product.properties[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.gram / pyunits.liter,
        )
    )
    # print(f"Flow to MVC: {flow_to_mvc} MGD")
    # print(f"TDS to MVC: {tds_to_mvc} g/L")

    m_mvc = build_and_run_mvc(
        recovery=recovery,
        Qin=Qin * pretreatment_recovery,
        tds=tds_to_mvc,
        electricity_cost=electricity_cost,
        heat_cost=heat_cost,
        **kwargs,
    )

    m_mvc.fs.disposal.properties[0].flow_vol_phase
    m_mvc.fs.disposal.properties[0].conc_mass_phase_comp
    m_mvc.fs.disposal.initialize()

    flow_to_dwi = value(
        pyunits.convert(
            m_mvc.fs.disposal.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.Mgallons / pyunits.day,
        )
    )
    tds_to_dwi = value(
        pyunits.convert(
            m_mvc.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.gram / pyunits.liter,
        )
    )

    m_dwi = build_and_run_dwi(
        Qin=flow_to_dwi,
        tds=tds_to_dwi,
        electricity_cost=electricity_cost,
        heat_cost=heat_cost,
        **kwargs,
    )

    m = build_agg_model(m_pre, m_mvc, m_dwi)

    assert degrees_of_freedom(m) == 0
    results = solve_permian_SOA(m)

    # print(f"System CAPEX = {m.fs.costing.total_capital_cost()}")
    # print(f"System OPEX = {m.fs.costing.total_operating_cost()}")
    # print(
    #     f"System Aggregate Flow Electricity = {m.fs.costing.aggregate_flow_electricity()}"
    # )
    # print(f"numerator = {numerator()}")
    # print(f"product_flow = {product_flow}")
    print(f"System LCOW = {m.fs.costing.LCOW()}")

    print(
        f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    )
    print(
        f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINISHED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    )
    print(
        f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% {Qin}, {tds}, {recovery} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    )
    print(
        f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    )
    return m, m_pre, m_mvc, m_dwi


def build_agg_costing_blk(
    b,
    models=None,
    base_currency=None,
    base_period=pyunits.year,
    inlet_flow=None,
    product_flow=None,
    disposal_flow=None,
):
    if base_currency is None:
        base_currency = pyunits.USD_2023

    b.models = models
    b.costing_blks = []
    for m in b.models:
        print(f"\n\n\n{m.name}\n\n\n")
        if m.name == "pretreatment":
            # continue
            costing_blk = m.fs.treatment.costing
        else:
            costing_blk = m.fs.costing
        b.costing_blks.append(costing_blk)
    b.base_currency = base_currency
    b.base_period = base_period
    b.capital_recovery_factor = Param(
        initialize=(value(b.costing_blks[0].capital_recovery_factor)),
        units=b.base_period**-1,
    )
    b.total_investment_factor = Param(
        initialize=(value(b.costing_blks[0].total_investment_factor)),
        units=b.base_period**-1,
    )
    b.utilization_factor = Param(initialize=value(b.costing_blks[0].utilization_factor))
    b.electricity_cost = Param(initialize=value(b.costing_blks[0].electricity_cost))
    b.heat_cost = Param(initialize=value(b.costing_blks[0].heat_cost))
    b.maintenance_labor_chemical_factor = Param(
        initialize=value(b.costing_blks[0].maintenance_labor_chemical_factor)
    )

    for costing_blk in b.costing_blks:

        assert value(costing_blk.capital_recovery_factor) == value(
            b.capital_recovery_factor
        )
        assert value(costing_blk.electricity_cost) == value(b.electricity_cost)
        assert value(costing_blk.heat_cost) == value(b.heat_cost)
        assert value(costing_blk.total_investment_factor) == value(
            b.total_investment_factor
        )
        assert value(costing_blk.maintenance_labor_chemical_factor) == value(
            b.maintenance_labor_chemical_factor
        )

    b.registered_unit_costing = Set()
    b.flow_types = Set()
    b.used_flows = Set()
    b.registered_flows = defaultdict(list)
    b.registered_flow_costs = defaultdict(list)
    b.used_flow_costs = dict()
    b.all_used_flow_costs = dict()

    for costing_blk in b.costing_blks:

        agg_flow_cost = getattr(costing_blk, "aggregate_flow_costs")

        for f in costing_blk.flow_types:
            if f not in b.flow_types:
                b.flow_types.add(f)

        b.all_used_flow_costs[m.name] = defaultdict(list)

        for f in costing_blk.used_flows:

            if f not in b.used_flows:
                b.used_flows.add(f)
                used_flow_cost = getattr(costing_blk, f"{f}_cost")
                b.used_flow_costs[f] = used_flow_cost

            if f in costing_blk._registered_flows.keys():
                b.registered_flows[f].extend(costing_blk._registered_flows[f])
                b.registered_flow_costs[f].append(value(agg_flow_cost[f]))

            used_flow_cost = getattr(costing_blk, f"{f}_cost")
            b.all_used_flow_costs[m.name][f].append(used_flow_cost)

        for ruc in costing_blk._registered_unit_costing:
            b.registered_unit_costing.add(ruc)

    b.total_capital_cost = Var(
        initialize=1e6, units=b.base_currency, bounds=(None, None)
    )

    @b.Constraint()
    def eq_total_capital_cost(blk):
        return blk.total_capital_cost == sum(
            value(pyunits.convert(cb.total_capital_cost, to_units=b.base_currency))
            for cb in blk.costing_blks
        )

    b.total_operating_cost = Var(
        initialize=1e4, units=b.base_currency / b.base_period, bounds=(None, None)
    )

    @b.Constraint()
    def eq_total_operating_cost(blk):
        return blk.total_operating_cost == sum(
            value(
                pyunits.convert(
                    cb.total_operating_cost,
                    to_units=b.base_currency / b.base_period,
                )
            )
            for cb in blk.costing_blks
        )

    b.aggregate_flow_costs = Var(b.used_flows, bounds=(None, None))

    @b.Constraint(b.used_flows)
    def eq_aggregate_flow_cost(blk, f):
        return blk.aggregate_flow_costs[f] == value(sum(blk.registered_flow_costs[f]))

    def agg_flow_rule(blk, f, funits):
        e = 0
        for flow in blk.registered_flows[f]:
            e += value(pyunits.convert(flow, to_units=funits))
        agg_flow = getattr(blk, f"aggregate_flow_{f}")
        return agg_flow == e

    for k, f in b.registered_flows.items():
        funits = pyunits.get_units(f[0])
        agg_var = Var(units=funits, doc=f"Aggregate flow for {k}")
        b.add_component(f"aggregate_flow_{k}", agg_var)
        agg_const = Constraint(rule=partial(agg_flow_rule, f=k, funits=funits))
        b.add_component(f"aggregate_flow_{k}_constraint", agg_const)

    @b.Expression()
    def LCOW(blk):
        numerator = (
            blk.total_capital_cost * blk.capital_recovery_factor
            + blk.total_operating_cost
        )
        denominator = pyunits.convert(
            product_flow * pyunits.m**3 / pyunits.s,
            to_units=pyunits.m**3 / b.base_period,
        )
        return pyunits.convert(
            numerator / denominator, to_units=b.base_currency / pyunits.m**3
        )

    @b.Expression()
    def LCOT(blk):
        numerator = (
            blk.total_capital_cost * blk.capital_recovery_factor
            + blk.total_operating_cost
        )
        denominator = pyunits.convert(
            product_flow * pyunits.m**3 / pyunits.s,
            to_units=pyunits.m**3 / b.base_period,
        )
        return pyunits.convert(
            numerator / denominator, to_units=b.base_currency / pyunits.m**3
        )

    flow_rate = product_flow * pyunits.m**3 / pyunits.s

    denominator = (
        pyunits.convert(flow_rate, to_units=pyunits.m**3 / base_period)
        * b.utilization_factor
    )
    direct_capex_lcows = Expression(
        Any,
        doc=f"Levelized Cost of Water - direct CAPEX",
        initialize=0 * base_currency / pyunits.m**3,
    )
    indirect_capex_lcows = Expression(
        Any,
        doc=f"Levelized Cost of Water - indirect CAPEX",
        initialize=0 * base_currency / pyunits.m**3,
    )
    fixed_opex_lcows = Expression(
        Any,
        doc=f"Levelized Cost of Water - fixed OPEX",
        initialize=0 * base_currency / pyunits.m**3,
    )
    variable_opex_lcows = Expression(
        Any,
        doc=f"Levelized Cost of Water - variable OPEX",
        initialize=0 * base_currency / pyunits.m**3,
    )
    flow_lcows = Expression(
        Any,
        doc=f"Levelized Cost of Water - chem/material/energy flows",
        initialize=0 * base_currency / pyunits.m**3,
    )
    flow_lcows_check = Expression(
        Any,
        doc=f"Levelized Cost of Water - chem/material/energy flows",
        initialize=0 * base_currency / pyunits.m**3,
    )
    unit_capex_lcows = Expression(
        Any,
        doc=f"Unit CAPEX-based LCOWs",
        initialize=0 * base_currency / pyunits.m**3,
    )
    unit_opex_lcows = Expression(
        Any,
        doc=f"Unit OPEX-based LCOWs",
        initialize=0 * base_currency / pyunits.m**3,
    )
    b.add_component("LCOW_component_direct_capex", direct_capex_lcows)
    b.add_component("LCOW_component_indirect_capex", indirect_capex_lcows)
    b.add_component("LCOW_component_fixed_opex", fixed_opex_lcows)
    b.add_component("LCOW_component_variable_opex", variable_opex_lcows)
    b.add_component("LCOW_flow", flow_lcows)
    b.add_component("LCOW_flow_check", flow_lcows_check)
    b.add_component("LCOW_capex", unit_capex_lcows)
    b.add_component("LCOW_opex", unit_opex_lcows)

    for u in b.registered_unit_costing:
        # print(u.name)
        direct_capex_numerator = 0 * base_currency / base_period
        indirect_capex_numerator = 0 * base_currency / base_period
        fixed_opex_numerator = 0 * base_currency / base_period
        variable_opex_numerator = 0 * base_currency / base_period
        total_capex_numerator = 0 * base_currency / base_period
        total_opex_numerator = 0 * base_currency / base_period
        if hasattr(u, "capital_cost"):
            direct_capital_cost = pyunits.convert(
                u.direct_capital_cost, to_units=base_currency
            )

            capital_cost = pyunits.convert(u.capital_cost, to_units=base_currency)

            direct_capex_numerator += b.capital_recovery_factor * direct_capital_cost
            capex_numerator = b.capital_recovery_factor * (
                b.total_investment_factor * capital_cost
            )
            indirect_capex_numerator += capex_numerator - direct_capex_numerator

            total_capex_numerator += capex_numerator
            # total_capex_numerator += indirect_capex_numerator

            fixed_opex_numerator += b.maintenance_labor_chemical_factor * capital_cost
        if hasattr(u, "fixed_operating_cost"):
            fixed_opex_numerator += pyunits.convert(
                u.fixed_operating_cost, to_units=base_currency / base_period
            )

        total_opex_numerator += fixed_opex_numerator

        if hasattr(u, "variable_operating_cost"):
            variable_opex_numerator += pyunits.convert(
                u.variable_operating_cost, to_units=base_currency / base_period
            )
            total_opex_numerator += variable_opex_numerator

        if ".chem_addition" in u.unit_model.name:
            unit_model_name = "chem_addition"
        elif ".EC." in u.unit_model.name:
            unit_model_name = "EC"
        elif ".cart_filt." in u.unit_model.name:
            unit_model_name = "cart_filt"
        elif ".MVC." in u.unit_model.name:
            unit_model_name = "MVC"
        elif ".DWI." in u.unit_model.name:
            unit_model_name = "DWI"
        else:
            raise ValueError()

        direct_capex_lcows[unit_model_name] = direct_capex_numerator / denominator
        indirect_capex_lcows[unit_model_name] = indirect_capex_numerator / denominator
        fixed_opex_lcows[unit_model_name] = fixed_opex_numerator / denominator
        variable_opex_lcows[unit_model_name] = variable_opex_numerator / denominator
        unit_capex_lcows[unit_model_name] = total_capex_numerator / denominator
        unit_opex_lcows[unit_model_name] = total_opex_numerator / denominator

    for f, fc in b.aggregate_flow_costs.items():
        i = fc / denominator
        flow_lcows[f] = i

    for ftype, flows in b.registered_flows.items():
        # part of total_variable_operating_cost
        # print()
        # print(ftype)
        # print()

        cost_var = b.used_flow_costs[ftype]
        for flow_expr in flows:
            if any(
                s in flow_expr.to_string()
                for s in ["fs.pv.electricity", "fs.trough.heat"]
            ):
                continue
            # print(flow_expr.to_string())
            flow_cost = pyunits.convert(
                flow_expr * cost_var, to_units=base_currency / base_period
            )
            i = (flow_cost * b.utilization_factor) / denominator
            flow_lcows_check[ftype] += i


def build_agg_model(m_pre, m_mvc, m_dwi):

    m_pre.name = "pretreatment"
    m_mvc.name = "MVC"
    m_dwi.name = "DWI"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.feed_flows = [m_pre.fs.treatment.feed.properties[0].flow_vol]
    m.fs.product_flows = [m_mvc.fs.product.properties[0].flow_vol_phase["Liq"]]

    m.fs.costing = Block()

    m.fs.flow_in = Expression(expr=sum(value(f) for f in m.fs.feed_flows))
    m.fs.flow_treated = Expression(expr=sum(value(f) for f in m.fs.product_flows))
    m.fs.flow_waste = Expression(expr=m.fs.flow_in - m.fs.flow_treated)
    m.fs.system_recovery = Expression(expr=m.fs.flow_treated / m.fs.flow_in)

    build_agg_costing_blk(
        m.fs.costing,
        models=[m_pre, m_mvc, m_dwi],
        # inlet_flows=m.fs.inlet_flows,
        product_flow=m.fs.flow_treated,
    )

    return m


def solve_permian_SOA(m):
    solver = get_solver()
    results = solver.solve(m)
    if not check_optimal_termination(results):
        results = solver.solve(m)
        if not check_optimal_termination(results):
            m.fs.optimal_solve_system.fix(0)
    return results


def run_permian_SOA_recovery_sweep(num_pts=5):

    recoveries = [0.4, 0.45, 0.5, 0.55, 0.6]

    m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA()

    rd, rd_pre, rd_mvc, rd_dwi = build_rds(m, m_pre, m_mvc, m_dwi, merge_cols=["recovery"])

    for r in recoveries:
        m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(recovery=r)
        rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"recovery": r})
    
    col_replace = ["m_agg.fs", "m_pre.fs", "m_mvc.fs", "m_dwi.fs"]
    for i, (r, cr) in enumerate(zip([rd, rd_pre, rd_mvc, rd_dwi], col_replace)):

        if i == 0:
            df = pd.DataFrame.from_dict(r)
            
            df.rename(
                columns={c: c.replace("fs", cr) for c in df.columns}, inplace=True
            )
            print(cr , len(df))
        else:
            df1 = pd.DataFrame.from_dict(r)
            df1.rename(
                columns={c: c.replace("fs", cr) for c in df1.columns}, inplace=True
            )
            df = pd.merge(df, df1, on=["recovery"])
            print(cr , len(df), len(df1))
    
    df.to_csv("permian_SOA_mvc_recovery_sweep.csv", index=False)


def run_permian_SOA_electricity_sweep(num_pts=5):
    # base case treatment energy component sweep
    min_rel = 0.5
    max_rel = 1.25

    ### electricity cost
    elec_cost = np.linspace(
        electricity_cost_base * min_rel, electricity_cost_base * max_rel, num_pts
    ).tolist()
    elec_cost.append(electricity_cost_base)

    ### heat cost
    heat_cost = np.linspace(
        heat_cost_base * min_rel, heat_cost_base * max_rel, num_pts
    ).tolist()
    heat_cost.append(heat_cost_base)

    m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA()

    rd, rd_pre, rd_mvc, rd_dwi = build_rds(m, m_pre, m_mvc, m_dwi, merge_cols=["electricity_cost"])

    for ec in elec_cost:
        m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(electricity_cost=ec)
        rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"electricity_cost": ec})
    
    col_replace = ["m_agg.fs", "m_pre.fs", "m_mvc.fs", "m_dwi.fs"]
    for i, (r, cr) in enumerate(zip([rd, rd_pre, rd_mvc, rd_dwi], col_replace)):

        if i == 0:
            df = pd.DataFrame.from_dict(r)
            
            df.rename(
                columns={c: c.replace("fs", cr) for c in df.columns}, inplace=True
            )
            print(cr , len(df))
        else:
            df1 = pd.DataFrame.from_dict(r)
            df1.rename(
                columns={c: c.replace("fs", cr) for c in df1.columns}, inplace=True
            )
            df = pd.merge(df, df1, on=["electricity_cost"])
            print(cr , len(df), len(df1))
    
    df.to_csv("permian_SOA_electricity_cost_sweep.csv", index=False)


def build_rds(m, m_pre, m_mvc, m_dwi, merge_cols=[]):

    rd = build_results_dict(m, skips=skips)
    rd_pre = build_results_dict(m_pre, skips=skips)
    rd_mvc = build_results_dict(m_mvc, skips=skips)
    rd_dwi = build_results_dict(m_dwi, skips=skips)

    for r in [rd, rd_pre, rd_mvc, rd_dwi]:
        for mc in merge_cols:
            r[mc] = list()
    
    return rd, rd_pre, rd_mvc, rd_dwi
    

def append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={}):
    
    rd = results_dict_append(m, rd)
    rd_pre = results_dict_append(m_pre, rd_pre)
    rd_mvc = results_dict_append(m_mvc, rd_mvc)
    rd_dwi = results_dict_append(m_dwi, rd_dwi)

    for r in [rd, rd_pre, rd_mvc, rd_dwi]:
        for c, v in merge_col_dict.items():
            r[c].append(v)
    

    return rd, rd_pre, rd_mvc, rd_dwi



if __name__ == "__main__":
    
    # run_permian_SOA_electricity_sweep()
    run_permian_SOA_recovery_sweep()

    # m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(Qin=5, tds=130)

    # rd, rd_pre, rd_mvc, rd_dwi = build_rds(m, m_pre, m_mvc, m_dwi, merge_cols=["electricity_cost"])
    # print(len(rd["electricity_cost"]))
    # rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"electricity_cost": electricity_cost_base})
    # print(len(rd["electricity_cost"]))
    # rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"electricity_cost": electricity_cost_base * 0.1})
    # print(len(rd["electricity_cost"]))
    # rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"electricity_cost": electricity_cost_base * 0.2})
    # print(len(rd["electricity_cost"]))
    # rd, rd_pre, rd_mvc, rd_dwi = append_rds(m, m_pre, m_mvc, m_dwi, rd, rd_pre, rd_mvc, rd_dwi, merge_col_dict={"electricity_cost": electricity_cost_base * 0.3})
    # print(len(rd["electricity_cost"]))

    # col_replace = ["m_agg.fs", "m_pre.fs", "m_mvc.fs", "m_dwi.fs"]

    # for i, (r, cr) in enumerate(zip([rd, rd_pre, rd_mvc, rd_dwi], col_replace)):
    #     if i == 0:
    #         df = pd.DataFrame.from_dict(r)
            
    #         df.rename(
    #             columns={c: c.replace("fs", cr) for c in df.columns}, inplace=True
    #         )
    #         print(cr , len(df))
    #     else:
    #         df1 = pd.DataFrame.from_dict(r)
    #         df1.rename(
    #             columns={c: c.replace("fs", cr) for c in df1.columns}, inplace=True
    #         )
    #         df = pd.merge(df, df1, on=["electricity_cost"])
    #         print(cr , len(df), len(df1))
    
    # df.to_csv("test_permian_SOA_agg.csv", index=False)
    # assert False

    # results_dict = build_results_dict(m)
    # results_dict["flow_mgd"] = list()
    # results_dict["tds"] = list()

    # recovery = [0.5, 0.55]
    # qs = [5]
    # salt = [90]

    # recovery = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65]
    # # recovery = [0.45]
    # qs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # salt = [75, 90, 105, 120, 130, 135, 150, 165]

    # results_df = pd.DataFrame()
    # results_df_pre = pd.DataFrame()
    # results_df_mvc = pd.DataFrame()
    # results_df_dwi = pd.DataFrame()

    # for r in recovery:

    #     tmp_results_dict = build_results_dict(m)
    #     tmp_results_dict["flow_mgd"] = list()
    #     tmp_results_dict["tds"] = list()
    #     tmp_results_dict["rerun"] = list()

    #     tmp_results_dict_pre = build_results_dict(m_pre)
    #     tmp_results_dict_pre["flow_mgd"] = list()
    #     tmp_results_dict_pre["tds"] = list()
    #     tmp_results_dict_pre["rerun"] = list()

    #     tmp_results_dict_mvc = build_results_dict(m_mvc)
    #     tmp_results_dict_mvc["flow_mgd"] = list()
    #     tmp_results_dict_mvc["tds"] = list()
    #     tmp_results_dict_mvc["rerun"] = list()

    #     tmp_results_dict_dwi = build_results_dict(m_dwi)
    #     tmp_results_dict_dwi["flow_mgd"] = list()
    #     tmp_results_dict_dwi["tds"] = list()
    #     tmp_results_dict_dwi["rerun"] = list()
    #     for q in qs:
    #         for s in salt:
    #             rerun = 0
    #             tds = s
    #             try:
    #                 m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(
    #                     Qin=q, tds=tds, recovery=r
    #                 )
    #             except:
    #                 rerun = 1
    #                 # tds = s * 1.01
    #                 m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(
    #                     Qin=q + 0.01, tds=s * 1.01, recovery=r
    #                 )
    #             tmp_results_dict = results_dict_append(m, tmp_results_dict)
    #             tmp_results_dict["flow_mgd"].append(q)
    #             tmp_results_dict["tds"].append(tds)
    #             tmp_results_dict["rerun"].append(rerun)

    #             tmp_results_dict_pre = results_dict_append(m_pre, tmp_results_dict_pre)
    #             tmp_results_dict_pre["flow_mgd"].append(q)
    #             tmp_results_dict_pre["tds"].append(tds)
    #             tmp_results_dict_pre["rerun"].append(rerun)

    #             tmp_results_dict_mvc = results_dict_append(m_mvc, tmp_results_dict_mvc)
    #             tmp_results_dict_mvc["flow_mgd"].append(q)
    #             tmp_results_dict_mvc["tds"].append(tds)
    #             tmp_results_dict_mvc["rerun"].append(rerun)

    #             tmp_results_dict_dwi = results_dict_append(m_dwi, tmp_results_dict_dwi)
    #             tmp_results_dict_dwi["flow_mgd"].append(q)
    #             tmp_results_dict_dwi["tds"].append(tds)
    #             tmp_results_dict_dwi["rerun"].append(rerun)

    #     results_df = pd.concat([results_df, pd.DataFrame.from_dict(tmp_results_dict)])
    #     results_df_pre = pd.concat(
    #         [results_df_pre, pd.DataFrame.from_dict(tmp_results_dict_pre)]
    #     )
    #     results_df_mvc = pd.concat(
    #         [results_df_mvc, pd.DataFrame.from_dict(tmp_results_dict_mvc)]
    #     )
    #     results_df_dwi = pd.concat(
    #         [results_df_dwi, pd.DataFrame.from_dict(tmp_results_dict_dwi)]
    #     )
    #     # results_dict = results_dict_append(m, results_dict, tmp_results_dict=tmp_results_dict)
    # save_dir = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/permian/SOA_results"
    # # pprint.pprint(results_dict)
    # timestr = time.strftime("%Y%m%d-%H%M%S")
    # # results_df.to_csv(f"permian_soa_results_{timestr}.csv", index=False)
    # # results_df_pre.to_csv(f"permian_soa_results_pre_{timestr}.csv", index=False)
    # # results_df_mvc.to_csv(f"permian_soa_results_mvc_{timestr}.csv", index=False)
    # # results_df_dwi.to_csv(f"permian_soa_results_dwi_{timestr}.csv", index=False)
    # # # m = build_and_run_permian_SOA(Qin=5, tds=105)
    # mvc_col_dict = dict()
    # for c in results_df_mvc.columns:
    #     mvc_col_dict[c] = c.replace("fs.MVC.", "MVC.").replace(
    #         "fs.costing", "MVC_costing"
    #     )
    # results_df_mvc.rename(columns=mvc_col_dict, inplace=True)

    # dwi_col_dict = dict()
    # for c in results_df_dwi.columns:
    #     dwi_col_dict[c] = c.replace("fs.DWI.", "DWI.").replace(
    #         "fs.costing", "DWI_costing"
    #     )
    # results_df_dwi.rename(columns=dwi_col_dict, inplace=True)
    # #     results_df.rename(
    # #     columns={
    # #         "fs.costing.LCOW": "system_LCOW",
    # #         "fs.costing.total_capital_cost": "system_total_capital_cost",
    # #         "fs.costing.total_operating_cost": "system_total_operating_cost",
    # #     }, inplace=True
    # # )
    # results_merged = pd.merge(
    #     results_df, results_df_pre, on=["flow_mgd", "tds", "rerun"]
    # )
    # results_merged = pd.merge(
    #     results_merged, results_df_mvc, on=["flow_mgd", "tds", "rerun"]
    # )
    # results_merged = pd.merge(
    #     results_merged, results_df_dwi, on=["flow_mgd", "tds", "rerun"]
    # )

    # results_merged.to_csv(
    #     f"{save_dir}/permian_soa_results_MERGED_{timestr}.csv", index=False
    # )
