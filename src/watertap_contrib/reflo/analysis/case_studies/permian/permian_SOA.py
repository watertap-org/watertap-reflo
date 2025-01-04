import pandas as pd
import numpy as np
import time
import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
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
# from watertap_contrib.reflo.kurby import *
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

reflo_dir = pathlib.Path(__file__).resolve().parents[3]
case_study_yaml = f"{reflo_dir}/data/technoeconomic/permian_case_study.yaml"
rho = 1125 * pyunits.kg / pyunits.m**3
rho_water = 995 * pyunits.kg / pyunits.m**3

solver = get_solver()

__all__ = ["build_and_run_permian_SOA", "solve_permian_SOA"]


def build_and_run_permian_SOA(
    pretreatment_recovery=0.98, recovery=0.5, Qin=5, tds=130, **kwargs
):
    m_pre = build_and_run_permian_pretreatment(Qin=Qin, tds=tds, **kwargs)

    # results = solve_permian_SOA(m_pre)
    # m.fs.treatment.product.display()
    flow_to_mvc = value(
        pyunits.convert(
            m_pre.fs.treatment.product.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.Mgallons / pyunits.day,
        )
    )
    flow_to_mvc = Qin * pretreatment_recovery
    tds_to_mvc = value(
        pyunits.convert(
            m_pre.fs.treatment.product.properties[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.gram / pyunits.liter,
        )
    )
    # print(f"Flow to MVC: {flow_to_mvc} MGD")
    # print(f"TDS to MVC: {tds_to_mvc} g/L")
    # assert False
    # assert False
    m_mvc = build_and_run_mvc(
        recovery=recovery, Qin=Qin * pretreatment_recovery, tds=tds_to_mvc, **kwargs
    )
    # results = solve_permian_SOA(m_mvc)
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
    # m_mvc.fs.disposal.properties[0].display()
    # assert False
    m_dwi = build_and_run_dwi(Qin=flow_to_dwi, tds=tds_to_dwi, **kwargs)
    # results = solve_permian_SOA(m_dwi)

    # print(f"Flow to pretreatment: {Qin} MGD")
    # print(f"TDS to pretreatment: {tds} g/L")

    # print(f"Flow to MVC: {flow_to_mvc} MGD")
    # print(f"TDS to MVC: {tds_to_mvc} g/L")

    # print(f"Flow to DWI: {flow_to_dwi} MGD")
    # print(f"TDS to DWI: {tds_to_dwi} g/L")

    product_flow = pyunits.convert(
        m_mvc.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.m**3 / pyunits.year,
    )()
    product_flow_mgd = pyunits.convert(
        m_mvc.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )()
    system_recovery = product_flow_mgd / Qin
    flow_to_product = product_flow_mgd

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = Block()
    m.fs.optimal_solve_system = Var(initialize=1)
    m.fs.optimal_solve_system.fix(1)

    m.fs.recovery = Var(initialize=recovery)
    m.fs.recovery.fix()
    m.fs.system_recovery = Var(initialize=system_recovery)
    m.fs.system_recovery.fix()
    m.fs.flow_mgd = Var(initialize=Qin)
    m.fs.flow_mgd.fix()
    m.fs.flow_to_mvc = Var(initialize=flow_to_mvc)
    m.fs.flow_to_mvc.fix()
    m.fs.flow_to_dwi = Var(initialize=flow_to_dwi)
    m.fs.flow_to_dwi.fix()
    m.fs.flow_to_product = Var(initialize=flow_to_product)
    m.fs.flow_to_product.fix()

    m.fs.tds = Var(initialize=tds)
    m.fs.tds.fix()
    m.fs.tds_to_mvc = Var(initialize=tds_to_mvc)
    m.fs.tds_to_mvc.fix()
    m.fs.tds_to_dwi = Var(initialize=tds_to_dwi)
    m.fs.tds_to_dwi.fix()

    m.fs.costing.total_capital_cost = Var(
        initialize=1e6, bounds=(0, None), units=pyunits.USD_2023
    )
    m.fs.costing.pretreatment_capital_cost = Var(
        initialize=m_pre.fs.treatment.costing.total_capital_cost(), bounds=(0, None)
    )
    m.fs.costing.MVC_capital_cost = Var(
        initialize=m_mvc.fs.costing.total_capital_cost(), bounds=(0, None)
    )
    m.fs.costing.DWI_capital_cost = Var(
        initialize=m_dwi.fs.costing.total_capital_cost(), bounds=(0, None)
    )
    m.fs.costing.LCOW = Var(initialize=5, bounds=(0, None))

    m.fs.costing.pretreatment_capital_cost.fix()
    m.fs.costing.MVC_capital_cost.fix()
    m.fs.costing.DWI_capital_cost.fix()

    m.fs.costing.aggregate_flow_electricity = Expression(
        expr=m_pre.fs.treatment.costing.aggregate_flow_electricity()
        + m_mvc.fs.costing.aggregate_flow_electricity()
        # + m_dwi.fs.costing.aggregate_flow_electricity()
    )
    m.fs.costing.total_operating_cost = Expression(
        expr=m_pre.fs.treatment.costing.total_operating_cost()
        + m_mvc.fs.costing.total_operating_cost()
        + m_dwi.fs.costing.total_operating_cost(),
    )
    m.fs.costing.annualized_capital_cost = Expression(
        expr=m.fs.costing.total_capital_cost
        * m_dwi.fs.costing.capital_recovery_factor()
    )

    m.fs.costing.total_capital_cost_constr = Constraint(
        expr=m.fs.costing.total_capital_cost
        == m.fs.costing.pretreatment_capital_cost
        + m.fs.costing.MVC_capital_cost
        + m.fs.costing.DWI_capital_cost
    )

    numerator = m.fs.costing.annualized_capital_cost + m.fs.costing.total_operating_cost

    m.fs.costing.LCOW_constr = Constraint(
        expr=m.fs.costing.LCOW
        == (m.fs.costing.annualized_capital_cost + m.fs.costing.total_operating_cost)
        / product_flow
    )

    # m_mvc.fs.costing.LCOW.display()
    # results = solver.solve(m)
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


def solve_permian_SOA(m):
    solver = get_solver()
    results = solver.solve(m)
    if not check_optimal_termination(results):
        results = solver.solve(m)
        if not check_optimal_termination(results):
            m.fs.optimal_solve_system.fix(0)
    return results


if __name__ == "__main__":

    from watertap_contrib.reflo.kurby import *

    m, m_pre, m_mvc, m_dwi = build_and_run_permian_SOA(Qin=5, tds=130)

    # results_dict = build_results_dict(m)
    # results_dict["flow_mgd"] = list()
    # results_dict["tds"] = list()

    # recovery = [0.4, 0.45]
    # qs = [5]
    # salt = [90, 130]

    # recovery = np.arange(0.4, 0.61, 0.01)
    # qs = np.arange(1, 10 + 0.5, 0.5)
    # salt = np.arange(75, 180, 5)

    # recovery = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65]
    # recovery = [0.45]
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

    # results_merged = pd.merge(
    #     results_df, results_df_pre, on=["flow_mgd", "tds", "rerun"], suffixes=("_system", "_pre")
    # )
    # results_merged = pd.merge(
    #     results_merged, results_df_mvc, on=["flow_mgd", "tds", "rerun"], suffixes=("_pre", "_mvc")
    # )
    # results_merged = pd.merge(
    #     results_merged, results_df_dwi, on=["flow_mgd", "tds", "rerun"], suffixes=("_mvc", "_dwi")
    # )

    # results_merged.to_csv(f"{save_dir}/permian_soa_results_MERGED_{timestr}.csv", index=False)
