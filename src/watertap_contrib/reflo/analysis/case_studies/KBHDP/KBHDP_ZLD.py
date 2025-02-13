import os
import multiprocessing
import itertools
import time
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from functools import partial
from collections import defaultdict
from pyomo.environ import (
    ConcreteModel,
    value,
    Any,
    Set,
    Var,
    Param,
    SolverFactory,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc

import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.scaling import *
from idaes.core.util.model_statistics import *

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.model_diagnostics.infeasible import *

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock as ZOParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock as CrystallizerParameterBlock,
)

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.analysis.case_studies.KBHDP import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.FPC import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import (
    build_results_dict,
    results_dict_append,
)

_log = idaeslog.getLogger(__name__)

__all__ = [
    "build_md_system",
    "build_and_run_kbhdp_zld",
    "build_and_run_md_system",
    "build_and_run_mec_system",
    "set_md_operating_conditions",
    "init_md_system",
    "build_and_run_primary_treatment",
    "build_primary_system",
    "build_primary_treatment",
    "add_primary_connections",
    "add_primary_system_scaling",
    "set_primary_operating_conditions",
    "init_primary_treatment",
    "run_kbhdp_base_case",
]

electricity_cost_base = 0.04988670259431006  # equivalent to 0.066 in USD2023
heat_cost_base = 0.0166

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

merge_cols = [
    "ro_recovery",
    "md_recovery",
    "crystallization_yield",
    "nacl_recovered_cost",
    "frac_elec_from_grid",
    "frac_heat_from_grid",
    "electricity_cost",
    "heat_cost",
]


save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/ZLD_results"
import idaes.logger as idaeslogger
import warnings

warnings.filterwarnings("ignore")

idaeslogger.getLogger("ideas.core").setLevel("CRITICAL")
idaeslogger.getLogger("ideas.core.util.scaling").setLevel("CRITICAL")
idaeslogger.getLogger("idaes.init").setLevel("CRITICAL")
idaeslogger.getLogger("idaes.watertap.core.util.initialization").setLevel("CRITICAL")
idaeslogger.getLogger("idaes.watertap.property_models.seawater_prop_pack").setLevel(
    "CRITICAL"
)
idaeslogger.getLogger("idaes.apps.grid_integration.multiperiod.multiperiod").setLevel(
    "CRITICAL"
)
idaeslogger.getLogger("idaes.core.util.scaling").setLevel("CRITICAL")
idaeslogger.getLogger("idaes.core.base.costing_base").setLevel("CRITICAL")
idaeslogger.getLogger(
    "idaes.watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_multiperiod_unit_model"
).setLevel("CRITICAL")
idaeslogger.getLogger("idaes.watertap.property_models.water_prop_pack").setLevel(
    "CRITICAL"
)
idaeslogger.getLogger("idaes.core.base.property_meta").setLevel("CRITICAL")


############################################################
############################################################
################### PRIMARY TREATMENT SYSTEM ###############
############################################################
############################################################


def build_and_run_primary_treatment(
    Qin=4,
    water_recovery=0.8,
    fixed_pressure=None,
    ro_mem_area=None,
    objective=None,
    electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
    heat_cost=heat_cost_base,
    aluminum_cost=2.23,
    **kwargs,
):

    m = build_primary_system()
    add_primary_treatment_costing(m)
    m.fs.costing.electricity_cost.fix(electricity_cost)
    m.fs.costing.heat_cost.fix(heat_cost)
    m.fs.costing.aluminum_cost.fix(aluminum_cost)
    add_primary_connections(m)
    set_primary_operating_conditions(m, Qin=Qin)
    add_primary_system_scaling(m)
    init_primary_system(m)
    optimize(
        m,
        water_recovery=water_recovery,
        fixed_pressure=fixed_pressure,
        ro_mem_area=ro_mem_area,
        objective=objective,
    )
    m.fs.brine.properties[0].conc_mass_phase_comp
    m.fs.brine.properties[0].flow_vol_phase
    m.fs.brine.flow_mgd = Expression(
        expr=pyunits.convert(
            m.fs.brine.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.Mgallons / pyunits.day,
        )
    )
    results = solve(m)

    return m


def build_primary_system():
    """
    Build primary treatment system flowsheet.
    """
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.MCAS_properties = MCASParameterBlock(
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Cl_-",
            "Mg_2+",
            "K_+",
            "SiO2",
            "Na_+",
            "SO2_-4+",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )

    m.fs.RO_properties = NaClParameterBlock()
    m.fs.ZO_properties = ZOParameterBlock(solute_list=["tds", "tss"])
    m.fs.MD_properties = SeawaterParameterBlock()

    build_primary_treatment(m)

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )

    return m


def build_primary_treatment(m):
    """
    Build treatment train through RO.
    """
    # treatment = m.fs.treatment = Block()

    m.fs.feed = Feed(property_package=m.fs.MCAS_properties)
    m.fs.product = Product(property_package=m.fs.RO_properties)
    m.fs.sludge = Product(property_package=m.fs.ZO_properties)
    m.fs.UF_waste = Product(property_package=m.fs.ZO_properties)
    m.fs.brine = Product(property_package=m.fs.MD_properties)

    m.fs.EC = FlowsheetBlock(dynamic=False)
    m.fs.UF = FlowsheetBlock(dynamic=False)
    m.fs.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.RO = FlowsheetBlock(dynamic=False)

    m.fs.MCAS_to_TDS_translator = Translator_MCAS_to_TDS(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.ZO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    m.fs.TDS_to_NaCl_translator = Translator_TDS_to_NACL(
        inlet_property_package=m.fs.ZO_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.RO_brine_to_MD_translator = Translator_NaCl_to_TDS(
        inlet_property_package=m.fs.RO_properties,
        outlet_property_package=m.fs.MD_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    build_ec(m, m.fs.EC, prop_package=m.fs.ZO_properties)
    build_UF(m, m.fs.UF, prop_package=m.fs.ZO_properties)
    build_ro(m, m.fs.RO, prop_package=m.fs.RO_properties)


def add_primary_connections(m):

    m.fs.feed_to_translator = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.MCAS_to_TDS_translator.inlet,
    )

    m.fs.translator_to_EC = Arc(
        source=m.fs.MCAS_to_TDS_translator.outlet,
        destination=m.fs.EC.feed.inlet,
    )

    m.fs.EC_to_UF = Arc(
        source=m.fs.EC.product.outlet,
        destination=m.fs.UF.feed.inlet,
    )

    m.fs.EC_to_sludge = Arc(
        source=m.fs.EC.disposal.outlet,
        destination=m.fs.sludge.inlet,
    )

    m.fs.UF_to_translator3 = Arc(
        source=m.fs.UF.product.outlet,
        destination=m.fs.TDS_to_NaCl_translator.inlet,
    )

    m.fs.UF_to_waste = Arc(
        source=m.fs.UF.disposal.outlet,
        destination=m.fs.UF_waste.inlet,
    )

    m.fs.translator_to_pump = Arc(
        source=m.fs.TDS_to_NaCl_translator.outlet,
        destination=m.fs.pump.inlet,
    )

    m.fs.pump_to_ro = Arc(
        source=m.fs.pump.outlet,
        destination=m.fs.RO.feed.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.RO.product.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.ro_to_brine_translator = Arc(
        source=m.fs.RO.disposal.outlet,
        destination=m.fs.RO_brine_to_MD_translator.inlet,
    )

    m.fs.brine_translator_to_brine = Arc(
        source=m.fs.RO_brine_to_MD_translator.outlet,
        destination=m.fs.brine.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_treatment_scaling(m):

    # set default scaling for MCAS
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    # set default scaling for ZO
    m.fs.ZO_properties.set_default_scaling("flow_mass_comp", 1e-2, index=("H2O"))
    m.fs.ZO_properties.set_default_scaling("flow_mass_comp", 1, index=("tds"))
    m.fs.ZO_properties.set_default_scaling("flow_mass_comp", 1e5, index=("tss"))

    # set default scaling for SW
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def add_primary_treatment_costing(m):

    m.fs.costing = TreatmentCosting()
    m.fs.costing.electricity_cost.fix(electricity_cost_base)
    m.fs.costing.heat_cost.fix(heat_cost_base)

    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    add_ec_costing(m, m.fs.EC, m.fs.costing)
    add_UF_costing(m, m.fs.UF, m.fs.costing)
    add_ro_costing(m, m.fs.RO, m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])


def add_primary_system_scaling(m):
    set_treatment_scaling(m)
    add_ec_scaling(m, m.fs.EC)
    add_UF_scaling(m.fs.UF)
    add_ro_scaling(m, m.fs.RO)
    calculate_scaling_factors(m)


def set_inlet_conditions(
    m,
    Qin=None,
    supply_pressure=101325,
):

    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    # treatment = m.fs.treatment

    # Convert Q_in from MGD to kg/s
    Qin = pyunits.convert(
        Qin * pyunits.Mgallon * pyunits.day**-1, to_units=pyunits.m**3 / pyunits.s
    )
    feed_density = 1000 * pyunits.kg / pyunits.m**3
    print("\n=======> SETTING FEED CONDITIONS <=======\n")
    print(f"Flow Rate: {value(Qin):<10.2f}{pyunits.get_units(Qin)}")

    if Qin is None:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    else:
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            Qin * feed_density
        )

    inlet_dict = {
        "Ca_2+": 0.61 * pyunits.kg / pyunits.m**3,
        "Mg_2+": 0.161 * pyunits.kg / pyunits.m**3,
        "Alkalinity_2-": 0.0821 * pyunits.kg / pyunits.m**3,
        "SiO2": 0.13 * pyunits.kg / pyunits.m**3,
        "Cl_-": 5.5 * pyunits.kg / pyunits.m**3,
        "Na_+": 5.5 * pyunits.kg / pyunits.m**3,
        "K_+": 0.016 * pyunits.kg / pyunits.m**3,
        "SO2_-4+": 0.23 * pyunits.kg / pyunits.m**3,
    }

    for solute, solute_conc in inlet_dict.items():
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            pyunits.convert(
                (
                    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
                    / (1000 * pyunits.kg / pyunits.m**3)
                )
                * solute_conc,
                to_units=pyunits.kg / pyunits.s,
            )
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", solute]),
            index=("Liq", solute),
        )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )

    feed_temperature = 273.15 + 20

    # # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)


def set_primary_operating_conditions(m, Qin=4, RO_pressure=20e5, **kwargs):
    # treatment = m.fs.treatment
    pump_efi = 0.8  # pump efficiency [-]
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=Qin)
    set_ec_operating_conditions(m, m.fs.EC)
    set_UF_op_conditions(m.fs.UF)
    m.fs.pump.efficiency_pump.fix(pump_efi)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(m, m.fs.RO, mem_area=10000)


def init_primary_treatment(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    # treatment = m.fs.treatment

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # assert_no_degrees_of_freedom(m)
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_translator)

    m.fs.MCAS_to_TDS_translator.initialize(optarg=optarg)
    propagate_state(m.fs.translator_to_EC)

    init_ec(m, m.fs.EC)
    propagate_state(m.fs.EC_to_UF)

    init_UF(m, m.fs.UF)
    propagate_state(m.fs.UF_to_translator3)
    propagate_state(m.fs.UF_to_waste)

    m.fs.TDS_to_NaCl_translator.initialize(optarg=optarg)
    propagate_state(m.fs.translator_to_pump)

    m.fs.pump.initialize(optarg=optarg)

    propagate_state(m.fs.pump_to_ro)

    init_ro_system(m, m.fs.RO)
    propagate_state(m.fs.ro_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.ro_to_brine_translator)

    m.fs.RO_brine_to_MD_translator.initialize()

    propagate_state(m.fs.brine_translator_to_brine)

    m.fs.brine.initialize()
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    # display_system_stream_table(m)


def init_primary_system(m, verbose=True, solver=None):
    print(f'\n{"=======> SYSTEM INITIALIZATION <=======":^60}\n')
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    init_primary_treatment(m)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")


############################################################
############################################################
################### MD  SYSTEM #############################
############################################################
############################################################


def build_and_run_md_system(
    Qin=1,
    Cin=60,
    water_recovery=0.7,
    electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
    heat_cost=heat_cost_base,
):

    m = build_md_system(Qin=Qin, Cin=Cin, water_recovery=water_recovery)
    m.fs.costing.electricity_cost.fix(electricity_cost)
    m.fs.costing.heat_cost.fix(heat_cost)
    set_md_operating_conditions(m)
    init_md_system(m)

    results = solve(m, tee=False)

    m.fs.feed.properties[0].flow_vol_phase
    m.fs.md.feed.properties[0].flow_vol_phase
    m.fs.disposal.properties[0].flow_vol_phase
    m.fs.disposal.properties[0].conc_mass_phase_comp
    m.fs.disposal.flow_mgd = Expression(
        expr=pyunits.convert(
            m.fs.disposal.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.Mgallons / pyunits.day,
        )
    )

    m.fs.md.unit.add_costing_module(m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])

    print("\nMD System Degrees of Freedom:", degrees_of_freedom(m), "\n")

    assert degrees_of_freedom(m) == 0

    results = solve(m)

    # For reporting purposes
    m.fs.costing.md_capital_cost = Param(
        initialize=value(m.fs.md.unit.costing.capital_cost), mutable=True
    )
    m.fs.costing.md_fixed_operating_cost = Param(
        initialize=value(m.fs.md.unit.costing.fixed_operating_cost), mutable=True
    )
    m.fs.costing.md_module_cost = Param(
        initialize=value(m.fs.md.unit.costing.module_cost), mutable=True
    )
    m.fs.costing.md_other_capital_cost = Param(
        initialize=value(m.fs.md.unit.costing.other_capital_cost), mutable=True
    )

    return m


def build_md_system(Qin=4, Cin=12, water_recovery=0.5):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.fs.costing.electricity_cost.fix(electricity_cost_base)
    m.fs.costing.heat_cost.fix(heat_cost_base)

    m.inlet_flow_rate = pyunits.convert(
        Qin * pyunits.Mgallons / pyunits.day, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = water_recovery

    # Property package
    m.fs.properties = SeawaterParameterBlock()

    # Create feed, product and concentrate state blocks
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Create MD unit model at flowsheet level
    m.fs.md = FlowsheetBlock(dynamic=False)
    build_md(m, m.fs.md, prop_package=m.fs.properties)

    m.fs.feed_to_md = Arc(source=m.fs.feed.outlet, destination=m.fs.md.feed.inlet)

    m.fs.md_to_product = Arc(
        source=m.fs.md.permeate.outlet, destination=m.fs.product.inlet
    )

    m.fs.md_to_disposal = Arc(
        source=m.fs.md.concentrate.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_md_operating_conditions(m):
    feed_flow_rate = m.fs.md.model_input["feed_flow_rate"]
    feed_salinity = m.fs.md.model_input["feed_salinity"]
    feed_temp = m.fs.md.model_input["feed_temp"]

    # m.fs.feed.properties.calculate_state(
    #     var_args={
    #         ("flow_vol_phase", "Liq"): pyunits.convert(
    #             feed_flow_rate * pyunits.L / pyunits.h,
    #             to_units=pyunits.m**3 / pyunits.s,
    #         ),
    #         ("conc_mass_phase_comp", ("Liq", "TDS")): feed_salinity,
    #         ("temperature", None): feed_temp + 273.15,
    #         ("pressure", None): 101325,
    #     },
    #     hold_state=True,
    # )

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): m.inlet_flow_rate,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.inlet_salinity,
            ("temperature", None): feed_temp + 273.15,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def init_md_system(m, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n-------------------- INITIALIZING MEMBRANE DISTILLATION --------------------\n\n"
    )
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    print("\n\n")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_md)

    init_md(m, m.fs.md, verbose=True, solver=None)

    propagate_state(m.fs.md_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.md_to_disposal)
    m.fs.disposal.initialize()


############################################################
############################################################
################### MEC SYSTEM #############################
############################################################
############################################################


def build_and_run_mec_system(
    Qin=None,
    Cin=None,
    number_effects=4,
    electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
    heat_cost=heat_cost_base,
    mec_kwargs=dict(),
    **kwargs,
):
    m = build_mec_system(number_effects=number_effects, **mec_kwargs)

    m.fs.costing.electricity_cost.fix(electricity_cost)
    m.fs.costing.heat_cost.fix(heat_cost)
    set_mec_operating_conditions(m, Qin=Qin, Cin=Cin, **mec_kwargs)
    init_mec_system(m)
    results = solve(m)

    # display_mec_streams(m, m.fs.MEC)
    return m


def build_mec_system(number_effects=4, nacl_recovered_cost=None, **kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = TreatmentCosting()

    m.fs.costing.electricity_cost.fix(electricity_cost_base)
    m.fs.costing.heat_cost.fix(heat_cost_base)
    m.nacl_recovered_cost = nacl_recovered_cost

    m.fs.properties = CrystallizerParameterBlock()
    m.fs.vapor_properties = WaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.solids = Product(property_package=m.fs.properties)

    m.fs.MEC = mec = FlowsheetBlock(dynamic=False)

    build_mec(m, m.fs.MEC, number_effects=number_effects)

    m.fs.feed_to_unit = Arc(source=m.fs.feed.outlet, destination=mec.unit.inlet)

    m.fs.mec_to_product = Arc(source=mec.product.outlet, destination=m.fs.product.inlet)

    m.fs.mec_to_solids = Arc(source=mec.solids.outlet, destination=m.fs.solids.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_mec_operating_conditions(
    m,
    Qin=None,  # MGD
    Cin=None,  # g/L
    flow_mass_water=None,
    flow_mass_tds=None,
    operating_pressures=[0.45, 0.25, 0.208, 0.095],
    crystallization_yield=0.5,
    saturated_steam_pressure_gage=3,
    heat_transfer_coefficient=0.1,
    **kwargs,
):
    atm_pressure = 101325 * pyunits.Pa
    rho = 1000 * pyunits.kg / pyunits.m**3

    saturated_steam_pressure = atm_pressure + pyunits.convert(
        saturated_steam_pressure_gage * pyunits.bar, to_units=pyunits.Pa
    )

    m.operating_pressures = operating_pressures
    m.crystallization_yield = crystallization_yield
    m.heat_transfer_coefficient = heat_transfer_coefficient
    m.saturated_steam_pressure = saturated_steam_pressure
    m.saturated_steam_pressure_gage = saturated_steam_pressure_gage

    if flow_mass_water is None:
        Qin = Qin * pyunits.Mgallons / pyunits.day
        m.flow_mass_water = pyunits.convert(Qin * rho, to_units=pyunits.kg / pyunits.s)
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            m.flow_mass_water
        )
    else:
        m.flow_mass_water = flow_mass_water
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_water)
    if flow_mass_tds is None:
        Cin = Cin * pyunits.gram / pyunits.liter
        m.flow_mass_tds = pyunits.convert(Qin * Cin, to_units=pyunits.kg / pyunits.s)
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(m.flow_mass_tds)
    else:
        m.flow_mass_tds = flow_mass_tds
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(flow_mass_tds)

    m.fs.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].flow_vol_phase["Liq"]


def init_mec_system(m):
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)

    init_mec(m, m.fs.MEC)

    propagate_state(m.fs.MEC.unit_to_product)
    m.fs.MEC.product.initialize()

    propagate_state(m.fs.mec_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.MEC.unit_to_solids)
    m.fs.MEC.solids.initialize()

    propagate_state(m.fs.mec_to_solids)
    m.fs.solids.initialize()

    add_mec_costing(m, m.fs.MEC)

    if m.nacl_recovered_cost is not None:
        m.fs.costing.nacl_recovered.cost.set_value(m.nacl_recovered_cost)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])

    m.fs.costing.initialize()


############################################################
############################################################
############################################################
############################################################


def solve(m, solver=None, tee=False, raise_on_failure=True, debug=False):
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
            check_jac(m.fs.treatment)
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


def optimize(
    m,
    water_recovery=0.5,
    fixed_pressure=None,
    ro_mem_area=None,
    objective="LCOW",
):
    # treatment = m.fs.treatment
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if objective == "LCOW":
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    else:
        m.fs.membrane_area_objective = Objective(expr=m.fs.RO.stage[1].module.area)

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.01
        upper_bound = 0.99
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        m.fs.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        m.fs.pump.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        m.fs.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        for idx, stage in m.fs.RO.stage.items():
            stage.module.area.fix(ro_mem_area)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for idx, stage in m.fs.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)


############################################################
############################################################
#################### ENERGY ################################
############################################################
############################################################


def build_and_run_zld_energy_model(
    m1,
    m_md,
    m_mec,
    frac_heat_from_grid=1,
    frac_elec_from_grid=1,
    design_size_start=1000,
    design_size_end=10000,
    increment_design_size=10,
    hours_storage=24,
    pv_cost_per_watt_installed=1.6,
    cst_cost_per_capacity_capital=560,
    cst_cost_per_storage_capital=62,
    **kwargs,
):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.frac_heat_from_grid = Param(initialize=frac_heat_from_grid)
    m.fs.frac_elec_from_grid = Param(initialize=frac_elec_from_grid)

    m.fs.costing = EnergyCosting()

    build_pv(m, energy_blk=m.fs)
    m.fs.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    build_trough_pysam(m, energy_blk=m.fs)
    m.fs.trough.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    set_trough_pysam_op_conditions(m, m.fs, hours_storage=hours_storage)

    m.fs.costing.pv_surrogate.cost_per_watt_installed.fix(pv_cost_per_watt_installed)
    m.fs.costing.trough_surrogate.cost_per_capacity_capital.fix(
        cst_cost_per_capacity_capital
    )
    m.fs.costing.trough_surrogate.cost_per_storage_capital.fix(
        cst_cost_per_storage_capital
    )

    # Find heat load for trough
    if frac_heat_from_grid != 1:
        heat_annual_required = value(
            pyunits.convert(
                (
                    m_md.fs.costing.aggregate_flow_heat
                    + m_mec.fs.costing.aggregate_flow_heat
                )
                * (1 - frac_heat_from_grid),
                to_units=pyunits.kWh * pyunits.year**-1,
            )
        )

        heat_load_required, pysam_results = get_zld_trough_heat_load(
            m, heat_annual_required, **kwargs
        )
        m.fs.trough.slope_of_fit = Param(initialize=pysam_results["slope"])
        m.fs.trough.heat_annual.fix(pysam_results["annual_energy"])
        m.fs.trough.electricity_annual.fix(pysam_results["electrical_load"])
        m.fs.trough.heat_load.fix(heat_load_required)
    else:
        pysam_results = None
        m.fs.trough.slope_of_fit = Param(initialize=0)
        m.fs.trough.heat_annual.fix(0)
        m.fs.trough.electricity_annual.fix(0)
        m.fs.trough.heat_load.fix(0)

    # presolve the trough model to get the electricity requirement
    assert degrees_of_freedom(m.fs.trough) == 0
    results = solve(m.fs.trough)
    assert_optimal_termination(results)
    if frac_elec_from_grid != 1:

        design_sizes = np.arange(
            design_size_start,
            design_size_end + increment_design_size,
            increment_design_size,
        )

        init_data = pd.DataFrame({"design_size": design_sizes})
        surr_results = m.fs.pv.surrogate.evaluate_surrogate(init_data)
        surr_results["design_size"] = design_sizes

        electricity_annual_required = value(
            pyunits.convert(
                (
                    m1.fs.costing.aggregate_flow_electricity
                    + m_md.fs.costing.aggregate_flow_electricity
                    + m_mec.fs.costing.aggregate_flow_electricity
                    + m.fs.trough.electricity
                )
                * (1 - frac_elec_from_grid),
                to_units=pyunits.kWh * pyunits.year**-1,
            )
        )

        i = 0
        annual_energy = surr_results.loc[i].annual_energy
        while annual_energy < electricity_annual_required:
            i += 1
            if i > max(surr_results.index):
                raise ValueError()
            annual_energy = surr_results.loc[i].annual_energy
            design_size = surr_results.loc[i].design_size

        m.fs.pv.design_size.fix(design_size)
    else:
        m.fs.pv.design_size.setlb(None)
        m.fs.pv.design_size.fix(0)
    assert degrees_of_freedom(m) == 0

    m.fs.costing.cost_process()

    results = solve(m)
    assert_optimal_termination(results)

    return m, pysam_results


def get_zld_trough_heat_load(
    m,
    heat_annual_required,
    heat_load_to_fit=50,
    fit_guess=False,
    slope=None,
    lb=0.98,
    ub=1.02,
    tries=None,
    max_tries=3,
    **kwargs,
):
    def _fit_guess():
        pysam_results = run_pysam_trough(m, heat_load=heat_load_to_fit)

        x = [0, heat_load_to_fit]  # input heat_load in MW
        y = [0, pysam_results["annual_energy"]]  # PySAM heat_annual in kWh

        popt, _ = curve_fit(heat_annual_linear, x, y)
        slope = popt[0]
        return slope

    def heat_annual_linear(heat_load, slope):
        return heat_load * slope

    if tries is None:
        tries = 1
    else:
        tries += 1

    if slope is None:
        fit_guess = True

    if fit_guess:
        slope = _fit_guess()

    heat_load_guess = heat_annual_required / slope
    pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
    heat_annual_guess = pysam_results["annual_energy"]

    if (
        not lb * heat_annual_required
        < pysam_results["annual_energy"]
        < heat_annual_required * ub
    ):
        ratio = heat_annual_required / heat_annual_guess
        heat_load_guess *= ratio
        pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
        heat_annual_guess = pysam_results["annual_energy"]

        # if (
        #     not lb * heat_annual_required
        #     < pysam_results["annual_energy"]
        #     < heat_annual_required * ub
        # ):
        #     ratio = heat_annual_required / heat_annual_guess
        #     heat_load_guess *= ratio
        #     pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
        #     heat_annual_guess = pysam_results["annual_energy"]

        #     if (
        #         not lb * heat_annual_required
        #         < pysam_results["annual_energy"]
        #         < heat_annual_required * ub
        #     ):
        #         ratio = heat_annual_required / heat_annual_guess
        #         heat_load_guess *= ratio
        #         pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
        #         heat_annual_guess = pysam_results["annual_energy"]

    # if np.isnan(heat_annual_guess):
    #     heat_load_guess *= 0.98
    #     pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
    #     heat_annual_guess = pysam_results["annual_energy"]

    # if np.isnan(heat_annual_guess):
    #     heat_load_guess *= 0.98
    #     pysam_results = run_pysam_trough(m, heat_load=heat_load_guess)
    #     heat_annual_guess = pysam_results["annual_energy"]
    pysam_results["slope"] = slope

    if (
        not lb * heat_annual_required
        < pysam_results["annual_energy"]
        < heat_annual_required * ub
    ):
        msg = f"The calculated heat_annual is not within the acceptable bounds after {tries} tries:"
        msg += f"\n\t{lb * heat_annual_required:.2f} < {heat_annual_guess:.2f} < {heat_annual_required * ub:2f}"

        if tries > max_tries:
            msg += f"\nYou have reached the max tries. Returning initial guess based on slope."
            # raise ValueError(msg)
            print(msg)
            heat_load_required = heat_annual_required / slope
            pysam_results["annual_energy"] = heat_annual_required
            return heat_load_required, pysam_results

        else:
            msg += "\nTrying again."
            print(msg)
            heat_load_required, pysam_results = get_zld_trough_heat_load(
                m,
                heat_annual_required,
                heat_load_to_fit=heat_load_guess * ub,
                slope=None,
                fit_guess=True,
                lb=lb,
                ub=ub,
                tries=tries,
                **kwargs,
            )
    heat_annual_guess = pysam_results["annual_energy"]
    msg = f"Found acceptable heat_load after {tries} tries!!"
    msg += f"\n\t{lb * heat_annual_required:.2f} < {heat_annual_guess:.2f} < {heat_annual_required * ub:2f}"
    print(msg)

    heat_load_required = heat_load_guess

    return heat_load_required, pysam_results


############################################################
############################################################
#################### FULL ZLD TRAIN ########################
############################################################
############################################################


def build_and_run_kbhdp_zld(
    primary_kwargs=dict(),
    md_kwargs=dict(),
    mec_op_kwargs=dict(),
    en_kwargs=dict(),
):

    m1 = build_and_run_primary_treatment(**primary_kwargs)

    flow_to_md = value(m1.fs.brine.flow_mgd)
    tds_to_md = value(m1.fs.brine.properties[0].conc_mass_phase_comp["Liq", "TDS"])

    m_md = build_and_run_md_system(Qin=flow_to_md, Cin=tds_to_md, **md_kwargs)

    m_md.fs.md_recovery = Param(initialize=m_md.water_recovery)

    flow_to_mec = value(m_md.fs.disposal.flow_mgd)
    tds_to_mec = value(
        m_md.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "TDS"]
    )

    flow_mass_water = value(
        m_md.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    flow_mass_tds = value(
        m_md.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    mec_kwargs = dict(flow_mass_water=flow_mass_water, flow_mass_tds=flow_mass_tds)
    mec_kwargs.update(mec_op_kwargs)

    m_mec = build_and_run_mec_system(mec_kwargs=mec_kwargs, **mec_op_kwargs)

    m_mec.fs.crystallization_yield = Param(
        initialize=value(
            m_mec.fs.MEC.unit.effects[1].effect.crystallization_yield["NaCl"]
        )
    )
    m_mec.fs.nacl_recovered_cost = Param(
        initialize=value(m_mec.fs.costing.nacl_recovered.cost)
    )

    # m_mec.fs.costing.display()
    m_en, pysam_results = build_and_run_zld_energy_model(
        m1,
        m_md,
        m_mec,
        slope=5495407.843363707,  # fit with 50 MW heat_load, 24 hr storage
        **en_kwargs,
    )

    m_agg = build_agg_model(m1, m_md, m_mec, m_en)

    print(f" dof = {degrees_of_freedom(m_agg)}")
    solver = SolverFactory("ipopt")
    results = solver.solve(m_agg)
    print(f"termination {results.solver.termination_condition}")
    # m_agg.fs.costing.display()
    # m_agg.fs.costing.LCOW.display()
    print(f"LCOW of ZLD System: {m_agg.fs.costing.LCOW()}")

    return m1, m_md, m_mec, m_en, m_agg


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
    b.base_currency = base_currency
    b.base_period = base_period
    b.capital_recovery_factor = Param(
        initialize=(value(models[0].fs.costing.capital_recovery_factor)),
        units=b.base_period**-1,
    )
    b.total_investment_factor = Param(
        initialize=(value(models[0].fs.costing.total_investment_factor)),
        units=b.base_period**-1,
    )
    b.utilization_factor = Param(
        initialize=value(models[0].fs.costing.utilization_factor)
    )
    b.electricity_cost = Param(initialize=value(models[0].fs.costing.electricity_cost))
    b.heat_cost = Param(initialize=value(models[0].fs.costing.heat_cost))
    b.maintenance_labor_chemical_factor = Param(
        initialize=value(models[0].fs.costing.maintenance_labor_chemical_factor)
    )

    for m in models:
        if m.name == "energy":
            continue
        assert value(m.fs.costing.capital_recovery_factor) == value(
            b.capital_recovery_factor
        )
        assert value(m.fs.costing.electricity_cost) == value(b.electricity_cost)
        assert value(m.fs.costing.heat_cost) == value(b.heat_cost)
        assert value(m.fs.costing.total_investment_factor) == value(
            b.total_investment_factor
        )
        assert value(m.fs.costing.maintenance_labor_chemical_factor) == value(
            b.maintenance_labor_chemical_factor
        )

    b.registered_unit_costing = Set()
    b.flow_types = Set()
    b.used_flows = Set()
    b.registered_flows = defaultdict(list)
    b.registered_flow_costs = defaultdict(list)
    b.used_flow_costs = dict()
    b.all_used_flow_costs = dict()

    for m in models:
        agg_flow_cost = getattr(m.fs.costing, "aggregate_flow_costs")

        for f in m.fs.costing.flow_types:
            if f not in b.flow_types:
                b.flow_types.add(f)

        b.all_used_flow_costs[m.name] = defaultdict(list)

        for f in m.fs.costing.used_flows:

            if f not in b.used_flows:
                b.used_flows.add(f)
                used_flow_cost = getattr(m.fs.costing, f"{f}_cost")
                b.used_flow_costs[f] = used_flow_cost

            if f in m.fs.costing._registered_flows.keys():
                b.registered_flows[f].extend(m.fs.costing._registered_flows[f])
                b.registered_flow_costs[f].append(value(agg_flow_cost[f]))

            used_flow_cost = getattr(m.fs.costing, f"{f}_cost")
            b.all_used_flow_costs[m.name][f].append(used_flow_cost)

        for ruc in m.fs.costing._registered_unit_costing:
            b.registered_unit_costing.add(ruc)

    b.total_capital_cost = Var(
        initialize=1e6, units=b.base_currency, bounds=(None, None)
    )

    @b.Constraint()
    def eq_total_capital_cost(blk):
        return blk.total_capital_cost == sum(
            value(
                pyunits.convert(
                    m.fs.costing.total_capital_cost, to_units=b.base_currency
                )
            )
            for m in models
        )

    b.total_operating_cost = Var(
        initialize=1e4, units=b.base_currency / b.base_period, bounds=(None, None)
    )

    @b.Constraint()
    def eq_total_operating_cost(blk):
        return blk.total_operating_cost == sum(
            value(
                pyunits.convert(
                    m.fs.costing.total_operating_cost,
                    to_units=b.base_currency / b.base_period,
                )
            )
            for m in models
        )

    b.aggregate_flow_costs = Var(b.used_flows, bounds=(None, None))

    @b.Constraint(b.used_flows)
    def eq_aggregate_flow_cost(blk, f):
        return blk.aggregate_flow_costs[f] == value(sum(blk.registered_flow_costs[f]))

    def agg_flow_rule(blk, f, funits):
        e = 0
        for flow in blk.registered_flows[f]:
            e += value(pyo.units.convert(flow, to_units=funits))
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
        pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / base_period)
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
            direct_capital_cost = pyo.units.convert(
                u.direct_capital_cost, to_units=base_currency
            )

            capital_cost = pyo.units.convert(u.capital_cost, to_units=base_currency)

            direct_capex_numerator += b.capital_recovery_factor * direct_capital_cost
            capex_numerator = b.capital_recovery_factor * (
                b.total_investment_factor * capital_cost
            )
            indirect_capex_numerator += capex_numerator - direct_capex_numerator

            total_capex_numerator += capex_numerator
            # total_capex_numerator += indirect_capex_numerator

            fixed_opex_numerator += b.maintenance_labor_chemical_factor * capital_cost
        if hasattr(u, "fixed_operating_cost"):
            fixed_opex_numerator += pyo.units.convert(
                u.fixed_operating_cost, to_units=base_currency / base_period
            )

        total_opex_numerator += fixed_opex_numerator

        if hasattr(u, "variable_operating_cost"):
            variable_opex_numerator += pyo.units.convert(
                u.variable_operating_cost, to_units=base_currency / base_period
            )
            total_opex_numerator += variable_opex_numerator

        if ".pump" in u.unit_model.name:
            unit_model_name = "pump"
        elif ".EC." in u.unit_model.name:
            unit_model_name = "EC"
        elif ".UF." in u.unit_model.name:
            unit_model_name = "UF"
        elif ".RO." in u.unit_model.name:
            unit_model_name = "RO"
        elif ".md." in u.unit_model.name:
            unit_model_name = "MD"
        elif ".MEC." in u.unit_model.name:
            unit_model_name = "MEC"
        elif ".pv" in u.unit_model.name:
            unit_model_name = "PV"
        elif ".trough" in u.unit_model.name:
            unit_model_name = "CST"
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
            flow_cost = pyo.units.convert(
                flow_expr * cost_var, to_units=base_currency / base_period
            )
            i = (flow_cost * b.utilization_factor) / denominator
            flow_lcows_check[ftype] += i


def build_agg_model(m1, m_md, m_mec, m_en):

    m1.name = "primary"
    m_md.name = "MD"
    m_mec.name = "MEC"
    m_en.name = "energy"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.feed_flows = [m1.fs.feed.properties[0].flow_vol_phase["Liq"]]
    m.fs.product_flows = [
        m1.fs.product.properties[0].flow_vol_phase["Liq"],
        m_md.fs.product.properties[0].flow_vol_phase["Liq"],
        m_mec.fs.product.properties[0].flow_vol_phase["Liq"],
    ]

    m.fs.costing = Block()

    m.fs.flow_in = Expression(expr=sum(value(f) for f in m.fs.feed_flows))
    m.fs.flow_treated = Expression(expr=sum(value(f) for f in m.fs.product_flows))
    m.fs.flow_waste = Expression(expr=m.fs.flow_in - m.fs.flow_treated)
    m.fs.system_recovery = Expression(expr=m.fs.flow_treated / m.fs.flow_in)

    build_agg_costing_blk(
        m.fs.costing,
        models=[m1, m_md, m_mec, m_en],
        # inlet_flows=m.fs.inlet_flows,
        product_flow=m.fs.flow_treated,
    )

    return m


def add_merge_cols(rd):
    for mc in merge_cols:
        rd[mc] = list()
    return rd


def append_merge_cols(rd, m1, m_md, m_mec, m_en, m_agg):
    rd["ro_recovery"].append(value(m1.fs.water_recovery))
    rd["md_recovery"].append(value(m_md.fs.md_recovery))
    rd["crystallization_yield"].append(value(m_mec.fs.crystallization_yield))
    rd["nacl_recovered_cost"].append(value(m_mec.fs.nacl_recovered_cost))
    rd["frac_elec_from_grid"].append(value(m_en.fs.frac_elec_from_grid))
    rd["frac_heat_from_grid"].append(value(m_en.fs.frac_heat_from_grid))
    rd["electricity_cost"].append(value(m_agg.fs.costing.electricity_cost))
    rd["heat_cost"].append(value(m_agg.fs.costing.heat_cost))
    return rd


def run_kbhdp_base_case():
    # base case conditions
    global rd_1, rd_md, rd_mec, rd_en, rd_agg

    primary_kwargs = dict(
        water_recovery=0.8,
        electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
        heat_cost=heat_cost_base,
    )
    md_kwargs = dict(
        water_recovery=0.7,
        electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
        heat_cost=heat_cost_base,
    )
    mec_op_kwargs = dict(
        crystallization_yield=0.9,
        nacl_recovered_cost=0,
        electricity_cost=electricity_cost_base,
        heat_cost=heat_cost_base,
    )
    en_kwargs = dict(
        frac_heat_from_grid=1,
        frac_elec_from_grid=1,
        pv_cost_per_watt_installed=1.6,
        cst_cost_per_capacity_capital=560,
        cst_cost_per_storage_capital=62,
    )

    m1, m_md, m_mec, m_en, m_agg = build_and_run_kbhdp_zld(
        primary_kwargs=primary_kwargs,
        md_kwargs=md_kwargs,
        mec_op_kwargs=mec_op_kwargs,
        en_kwargs=en_kwargs,
    )

    # build results dict and append baseline results

    rd_1 = build_results_dict(m1, skips=skips)
    rd_1 = add_merge_cols(rd_1)
    # rd_1 = results_dict_append(m1, rd_1)
    # rd_1 = append_merge_cols(rd_1, m1, m_md, m_mec, m_en, m_agg)

    rd_md = build_results_dict(m_md, skips=skips)
    rd_md = add_merge_cols(rd_md)
    # rd_md = results_dict_append(m_md, rd_md)
    # rd_md = append_merge_cols(rd_md, m1, m_md, m_mec, m_en, m_agg)

    rd_mec = build_results_dict(m_mec, skips=skips)
    rd_mec = add_merge_cols(rd_mec)
    # rd_mec = results_dict_append(m_mec, rd_mec)
    # rd_mec = append_merge_cols(rd_mec, m1, m_md, m_mec, m_en, m_agg)

    rd_en = build_results_dict(m_en, skips=skips)
    rd_en = add_merge_cols(rd_en)
    # rd_en = results_dict_append(m_en, rd_en)
    # rd_en = append_merge_cols(rd_en, m1, m_md, m_mec, m_en, m_agg)

    rd_agg = build_results_dict(m_agg, skips=skips)
    rd_agg = add_merge_cols(rd_agg)
    # rd_agg = results_dict_append(m_agg, rd_agg)
    # rd_agg = append_merge_cols(rd_agg, m1, m_md, m_mec, m_en, m_agg)

    return rd_1, rd_md, rd_mec, rd_en, rd_agg


def run_kbhdp_zld_sweep(
    ro_recovery=0.8,
    md_recovery=0.7,
    cryst_yields=0.9,
    nacl_recov_cost=-0,
    frac_elec_from_grid=0,
    frac_heat_from_grid=0,
    aluminum_cost=2.23,
    electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
    heat_cost=heat_cost_base,
    pv_cost=1.6,
    trough_cost=560,
    storage_cost=62,
    file_append="",
):

    # save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/kurby_reflo/case_studies/KBHDP/ZLD/zld_sweep_results"

    if not isinstance(ro_recovery, list):
        ro_recovery = [ro_recovery]
    if not isinstance(md_recovery, list):
        md_recovery = [md_recovery]
    if not isinstance(cryst_yields, list):
        cryst_yields = [cryst_yields]
    if not isinstance(nacl_recov_cost, list):
        nacl_recov_cost = [nacl_recov_cost]
    if not isinstance(frac_elec_from_grid, list):
        frac_elec_from_grid = [frac_elec_from_grid]
    if not isinstance(frac_heat_from_grid, list):
        frac_heat_from_grid = [frac_heat_from_grid]
    if not isinstance(electricity_cost, list):
        electricity_cost = [electricity_cost]
    if not isinstance(heat_cost, list):
        heat_cost = [heat_cost]
    if not isinstance(pv_cost, list):
        pv_cost = [pv_cost]
    if not isinstance(trough_cost, list):
        trough_cost = [trough_cost]
    if not isinstance(storage_cost, list):
        storage_cost = [storage_cost]
    if not isinstance(aluminum_cost, list):
        aluminum_cost = [aluminum_cost]

    # base case conditions

    primary_kwargs = dict(
        water_recovery=0.8,
        electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
        heat_cost=heat_cost_base,
    )
    md_kwargs = dict(
        water_recovery=0.7,
        electricity_cost=electricity_cost_base,  # USD_2018 / kWh, equivalent to 0.066 USD_2023/kWh
        heat_cost=heat_cost_base,
    )
    mec_op_kwargs = dict(
        crystallization_yield=0.9,
        nacl_recovered_cost=0,
        electricity_cost=electricity_cost_base,
        heat_cost=heat_cost_base,
    )
    en_kwargs = dict(
        frac_heat_from_grid=0,
        frac_elec_from_grid=0,
        pv_cost_per_watt_installed=1.6,
        cst_cost_per_capacity_capital=560,
        cst_cost_per_storage_capital=62,
    )

    m1, m_md, m_mec, m_en, m_agg = build_and_run_kbhdp_zld(
        primary_kwargs=primary_kwargs,
        md_kwargs=md_kwargs,
        mec_op_kwargs=mec_op_kwargs,
        en_kwargs=en_kwargs,
    )
    # rd_mec = build_results_dict(m_mec, skips=skips)
    # rd_mec = results_dict_append(m_mec, rd_mec)
    # for k, v in rd_mec.items():
    #     if "_recovered" in k:
    #         print(v)

    # build results dict and append baseline results

    rd_1 = build_results_dict(m1, skips=skips)
    rd_1 = add_merge_cols(rd_1)
    rd_1 = results_dict_append(m1, rd_1)
    rd_1 = append_merge_cols(rd_1, m1, m_md, m_mec, m_en, m_agg)

    rd_md = build_results_dict(m_md, skips=skips)
    rd_md = add_merge_cols(rd_md)
    rd_md = results_dict_append(m_md, rd_md)
    rd_md = append_merge_cols(rd_md, m1, m_md, m_mec, m_en, m_agg)

    rd_mec = build_results_dict(m_mec, skips=skips)
    rd_mec = add_merge_cols(rd_mec)
    rd_mec = results_dict_append(m_mec, rd_mec)
    rd_mec = append_merge_cols(rd_mec, m1, m_md, m_mec, m_en, m_agg)

    rd_en = build_results_dict(m_en, skips=skips)
    rd_en = add_merge_cols(rd_en)
    rd_en = results_dict_append(m_en, rd_en)
    rd_en = append_merge_cols(rd_en, m1, m_md, m_mec, m_en, m_agg)

    rd_agg = build_results_dict(m_agg, skips=skips)
    rd_agg = add_merge_cols(rd_agg)
    rd_agg = results_dict_append(m_agg, rd_agg)
    rd_agg = append_merge_cols(rd_agg, m1, m_md, m_mec, m_en, m_agg)

    for ro in ro_recovery:
        for md in md_recovery:
            for cy in cryst_yields:
                for nrc in nacl_recov_cost:
                    for fe in frac_elec_from_grid:
                        for fh in frac_heat_from_grid:
                            for ec in electricity_cost:
                                for hc in heat_cost:
                                    for pc in pv_cost:
                                        for tc in trough_cost:
                                            for sc in storage_cost:
                                                for al in aluminum_cost:

                                                    primary_kwargs = dict(
                                                        water_recovery=ro,
                                                        electricity_cost=ec,
                                                        heat_cost=hc,
                                                        aluminum_cost=al,
                                                    )
                                                    md_kwargs = dict(
                                                        water_recovery=md,
                                                        electricity_cost=ec,
                                                        heat_cost=hc,
                                                    )
                                                    mec_op_kwargs = dict(
                                                        crystallization_yield=cy,
                                                        nacl_recovered_cost=nrc,
                                                        electricity_cost=ec,
                                                        heat_cost=hc,
                                                    )
                                                    en_kwargs = dict(
                                                        frac_heat_from_grid=fh,
                                                        frac_elec_from_grid=fe,
                                                        pv_cost_per_watt_installed=pc,
                                                        cst_cost_per_capacity_capital=tc,
                                                        cst_cost_per_storage_capital=sc,
                                                    )

                                                    m1, m_md, m_mec, m_en, m_agg = (
                                                        build_and_run_kbhdp_zld(
                                                            primary_kwargs=primary_kwargs,
                                                            md_kwargs=md_kwargs,
                                                            mec_op_kwargs=mec_op_kwargs,
                                                            en_kwargs=en_kwargs,
                                                        )
                                                    )

                                                    rd_1 = results_dict_append(m1, rd_1)
                                                    rd_1 = append_merge_cols(
                                                        rd_1,
                                                        m1,
                                                        m_md,
                                                        m_mec,
                                                        m_en,
                                                        m_agg,
                                                    )
                                                    rd_md = results_dict_append(
                                                        m_md, rd_md
                                                    )
                                                    rd_md = append_merge_cols(
                                                        rd_md,
                                                        m1,
                                                        m_md,
                                                        m_mec,
                                                        m_en,
                                                        m_agg,
                                                    )
                                                    rd_mec = results_dict_append(
                                                        m_mec, rd_mec
                                                    )
                                                    rd_mec = append_merge_cols(
                                                        rd_mec,
                                                        m1,
                                                        m_md,
                                                        m_mec,
                                                        m_en,
                                                        m_agg,
                                                    )
                                                    rd_en = results_dict_append(
                                                        m_en, rd_en
                                                    )
                                                    rd_en = append_merge_cols(
                                                        rd_en,
                                                        m1,
                                                        m_md,
                                                        m_mec,
                                                        m_en,
                                                        m_agg,
                                                    )
                                                    rd_agg = results_dict_append(
                                                        m_agg, rd_agg
                                                    )
                                                    rd_agg = append_merge_cols(
                                                        rd_agg,
                                                        m1,
                                                        m_md,
                                                        m_mec,
                                                        m_en,
                                                        m_agg,
                                                    )

    # return rd_1, rd_md, rd_mec, rd_en, rd_agg, m1, m_md, m_mec, m_en, m_agg

    df_1 = pd.DataFrame.from_dict(rd_1)
    df_1.rename(
        columns={c: c.replace("fs", "m1.fs") for c in df_1.columns}, inplace=True
    )
    df_md = pd.DataFrame.from_dict(rd_md)
    df_md.rename(
        columns={c: c.replace("fs", "md.fs") for c in df_md.columns}, inplace=True
    )
    df_mec = pd.DataFrame.from_dict(rd_mec)
    df_mec.rename(
        columns={c: c.replace("fs", "m_mec.fs") for c in df_mec.columns}, inplace=True
    )
    df_en = pd.DataFrame.from_dict(rd_en)
    df_en.rename(
        columns={c: c.replace("fs", "m_en.fs") for c in df_en.columns}, inplace=True
    )
    df_agg = pd.DataFrame.from_dict(rd_agg)
    df_agg.rename(
        columns={c: c.replace("fs", "m_agg.fs") for c in df_agg.columns}, inplace=True
    )

    timestr = time.strftime("%Y%m%d-%H%M%S")

    df_1.to_csv(f"{save_path}/kbhdp_zld_1_{timestr}.csv", index=False)
    df_md.to_csv(f"{save_path}/kbhdp_zld_md_{timestr}.csv", index=False)
    df_mec.to_csv(f"{save_path}/kbhdp_zld_mec_{timestr}.csv", index=False)
    df_en.to_csv(f"{save_path}/kbhdp_zld_en_{timestr}.csv", index=False)
    df_agg.to_csv(f"{save_path}/kbhdp_zld_agg_{timestr}.csv", index=False)

    df_combo = pd.merge(df_1, df_md, on=merge_cols)
    df_combo = pd.merge(df_combo, df_mec, on=merge_cols)
    df_combo = pd.merge(df_combo, df_en, on=merge_cols)
    df_combo = pd.merge(df_combo, df_agg, on=merge_cols)

    df_combo.to_csv(
        f"{save_path}/kbhdp_zld_combo_{timestr}_{file_append}.csv", index=False
    )

    print("COMPLETED")
    return df_combo


def main():
    run_kbhdp_zld_sweep()


def run_kbhdp_zld_mp(
    primary_kwargs_input,
    md_kwargs_input,
    mec_kwargs_input,
    en_kwargs_input,
    rd_1,
    rd_md,
    rd_mec,
    rd_en,
    rd_agg,
    electricity_cost,
    heat_cost,
):

    primary_kwargs = dict(
        water_recovery=0.8,
        electricity_cost=electricity_cost,
        aluminum_cost=2.23,
        heat_cost=heat_cost,
    )
    md_kwargs = dict(
        water_recovery=0.7,
        electricity_cost=electricity_cost,
        heat_cost=heat_cost,
    )
    mec_op_kwargs = dict(
        crystallization_yield=0.9,
        nacl_recovered_cost=0,
        electricity_cost=electricity_cost,
        heat_cost=heat_cost,
    )
    en_kwargs = dict(
        frac_heat_from_grid=1,
        frac_elec_from_grid=1,
        pv_cost_per_watt_installed=1.6,
        cst_cost_per_capacity_capital=560,
        cst_cost_per_storage_capital=62,
    )

    for k, v in primary_kwargs.items():
        if k in primary_kwargs_input.keys():
            primary_kwargs[k] = primary_kwargs_input[k]

    for k, v in md_kwargs.items():
        if k in md_kwargs_input.keys():
            md_kwargs[k] = md_kwargs_input[k]

    for k, v in mec_op_kwargs.items():
        if k in mec_kwargs_input.keys():
            mec_op_kwargs[k] = mec_kwargs_input[k]

    for k, v in en_kwargs.items():
        if k in en_kwargs_input.keys():
            en_kwargs[k] = en_kwargs_input[k]

    m1, m_md, m_mec, m_en, m_agg = build_and_run_kbhdp_zld(
        primary_kwargs=primary_kwargs,
        md_kwargs=md_kwargs,
        mec_op_kwargs=mec_op_kwargs,
        en_kwargs=en_kwargs,
    )

    col_replace = ["m1.fs", "m_md.fs", "m_mec.fs", "m_en.fs", "m_agg.fs"]

    for i, (rd, m, cr) in enumerate(
        zip(
            [rd_1, rd_md, rd_mec, rd_en, rd_agg],
            [m1, m_md, m_mec, m_en, m_agg],
            col_replace,
        )
    ):

        rd = results_dict_append(m, rd)
        rd = append_merge_cols(rd, m1, m_md, m_mec, m_en, m_agg)

        if i == 0:
            df = pd.DataFrame.from_dict(rd)
            df.rename(
                columns={c: c.replace("fs", cr) for c in df.columns}, inplace=True
            )
        else:
            df1 = pd.DataFrame.from_dict(rd)
            df1.rename(
                columns={c: c.replace("fs", cr) for c in df1.columns}, inplace=True
            )

            df = pd.merge(df, df1, on=merge_cols)

    return df


def run_elec_heat_cost_sweep(num_pts=5):

    # base case treatment energy component sweep
    min_rel = 0.5
    max_rel = 1.25

    ### PV cost per watt
    al_cost_base = 2.23

    aluminum_cost = np.linspace(
        al_cost_base * min_rel, al_cost_base * max_rel, num_pts
    ).tolist()

    ### NaCl recovery cost
    nacl_recov_cost = [-0.05, -0.025, 0]

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

    ro_recovery = [0.6, 0.7, 0.8]
    md_recovery = [0.6, 0.7, 0.75]
    cryst_yield = [0.7, 0.8, 0.9]

    rd_1, rd_md, rd_mec, rd_en, rd_agg = run_kbhdp_base_case()

    primary_inputs = [{}]
    md_inputs = [{}]
    mec_inputs = [{}]
    en_inputs = [{}]

    # for ac in aluminum_cost:
    #     for nrc in nacl_recov_cost:
    # for ec in elec_cost:
    #     for hc in heat_cost:
    # for ro in ro_recovery:
    #     for ac in aluminum_cost:

    #         primary_inputs.append(
    #             dict(
    #                 water_recovery=ro,
    #                 aluminum_cost=ac,
    #                 # electricity_cost=ec,
    #                 # heat_cost=hc,
    #             )
    #         )

    # for md in md_recovery:

    #     md_inputs.append(
    #         dict(
    #             water_recovery=md,
    #             # electricity_cost=ec,
    #             # heat_cost=hc,
    #         )
    #     )

    # for cy in cryst_yield:
    #     for nrc in nacl_recov_cost:
    #         mec_inputs.append(
    #             dict(
    #                 crystallization_yield=cy,
    #                 nacl_recovered_cost=nrc,
    #                 # electricity_cost=ec,
    #                 # heat_cost=hc,
    #             )
    #         )

    # en_inputs.append({})
    d1 = {"a": 1}
    d2 = {"b": 1}
    d3 = {"c": 1}
    d4 = {"d": 1}
    d5 = {"e": 1}

    args = list(
        itertools.product(
            primary_inputs,
            md_inputs,
            mec_inputs,
            en_inputs,
            *[[rd] for rd in [rd_1, rd_md, rd_mec, rd_en, rd_agg]],
            # *[[rd] for rd in [d1, d2, d3, d4, d5]],
            *[elec_cost],
            *[heat_cost],
        )
    )
    # import pprint
    # import sys
    # pprint.pprint(primary_inputs)
    # pprint.pprint(args[:10])
    # print(len(args))
    # print(sys.getsizeof(args) * 1e-9)
    # print(multiprocessing.cpu_count())
    # assert False

    with multiprocessing.Pool(processes=8) as pool:
        results = pool.starmap(run_kbhdp_zld_mp, args)
    df_results = pd.concat(results)

    timestr = time.strftime("%Y%m%d-%H%M%S")
    # save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/ZLD_results"
    df_results.to_csv(
        f"{save_path}/kbhdp_zld_results_electricity_heat_cost_sweep_{timestr}.csv",
        index=False,
    )


def run_recov_chem_mat_cost_sweep(num_pts=5):

    min_rel = 0.5
    max_rel = 1.25

    ### electricity cost
    elec_cost = []
    elec_cost.append(electricity_cost_base)

    ### heat cost
    heat_cost = []
    heat_cost.append(heat_cost_base)

    ### PV cost per watt
    al_cost_base = 2.23

    aluminum_cost = np.linspace(
        al_cost_base * min_rel, al_cost_base * max_rel, num_pts
    ).tolist()

    ### NaCl recovery cost
    nacl_recov_cost = [-0.05, -0.025, 0]

    ro_recovery = [0.6, 0.7, 0.8]
    md_recovery = [0.5, 0.6, 0.7]
    cryst_yield = [0.7, 0.8, 0.9]

    rd_1, rd_md, rd_mec, rd_en, rd_agg = run_kbhdp_base_case()

    primary_inputs = [{}]
    md_inputs = [{}]
    mec_inputs = [{}]
    en_inputs = [{}]

    for ro in ro_recovery:
        for ac in aluminum_cost:

            primary_inputs.append(
                dict(
                    aluminum_cost=ac,
                    water_recovery=ro,
                )
            )

    for md in md_recovery:

        md_inputs.append(
            dict(
                water_recovery=md,
            )
        )

    for cy in cryst_yield:
        for nrc in nacl_recov_cost:
            mec_inputs.append(
                dict(
                    crystallization_yield=cy,
                    nacl_recovered_cost=nrc,
                )
            )

    # d1 = {"a": 1}
    # d2 = {"b": 1}
    # d3 = {"c": 1}
    # d4 = {"d": 1}
    # d5 = {"e": 1}

    args = list(
        itertools.product(
            primary_inputs,
            md_inputs,
            mec_inputs,
            en_inputs,
            *[[rd] for rd in [rd_1, rd_md, rd_mec, rd_en, rd_agg]],
            # *[[rd] for rd in [d1, d2, d3, d4, d5]],
            *[elec_cost],
            *[heat_cost],
        )
    )
    # import pprint
    # import sys
    # pprint.pprint(primary_inputs)
    # pprint.pprint(args)
    # print(len(args))
    # print(sys.getsizeof(args) * 1e-9)
    # print(multiprocessing.cpu_count())
    # assert False

    with multiprocessing.Pool(processes=9) as pool:
        results = pool.starmap(run_kbhdp_zld_mp, args)
    df_results = pd.concat(results)

    timestr = time.strftime("%Y%m%d-%H%M%S")
    # save_path = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/ZLD_results"
    df_results.to_csv(
        f"{save_path}/kbhdp_zld_results_recovery_alum_cost_nacl_recov_cost_sweep_{timestr}.csv",
        index=False,
    )


def run_energy_param_sweep(num_pts=4):

    # base case treatment energy component sweep
    min_rel = 0.5
    max_rel = 1.25

    ### PV cost per watt
    pv_cost_base = 1.6

    pv_cost = np.linspace(
        pv_cost_base * min_rel, pv_cost_base * max_rel, num_pts
    ).tolist()
    pv_cost.append(pv_cost_base)

    ### Trough cost per watt
    trough_cost_base = 560

    trough_cost = np.linspace(
        trough_cost_base * min_rel, trough_cost_base * max_rel, num_pts
    ).tolist()
    trough_cost.append(trough_cost_base)

    ### Trough storage cost per watt
    storage_cost_base = 62

    storage_cost = np.linspace(
        storage_cost_base * min_rel, storage_cost_base * max_rel, num_pts
    ).tolist()
    storage_cost.append(storage_cost_base)

    ### frac elec from grid
    frac_elec_from_grid = [0.5, 0.25, 0]
    ### frac heat from grid
    frac_heat_from_grid = [0.5, 0.25, 0]

    ### electricity cost
    # elec_cost = np.linspace(
    #     electricity_cost_base * min_rel, electricity_cost_base * max_rel, num_pts
    # ).tolist()
    # elec_cost.append(electricity_cost_base)
    # # elec_cost = list()

    # ### heat cost

    # heat_cost = np.linspace(
    #     heat_cost_base * min_rel, heat_cost_base * max_rel, num_pts
    # ).tolist()
    # heat_cost.append(heat_cost_base)
    # heat_cost = list()

    rd_1, rd_md, rd_mec, rd_en, rd_agg = run_kbhdp_base_case()

    primary_inputs = [{}]
    md_inputs = [{}]
    mec_inputs = [{}]
    en_inputs = [{}]

    for pc in pv_cost:
        for tc in trough_cost:
            for sc in storage_cost:
                for fe in frac_elec_from_grid:
                    for fh in frac_heat_from_grid:
                        en_inputs.append(
                            dict(
                                pv_cost_per_watt_installed=pc,
                                cst_cost_per_capacity_capital=tc,
                                cst_cost_per_storage_capital=sc,
                                frac_elec_from_grid=fe,
                                frac_heat_from_grid=fh,
                            )
                        )

    ### electricity cost
    elec_cost = []
    elec_cost.append(electricity_cost_base)

    ### heat cost
    heat_cost = []
    heat_cost.append(heat_cost_base)

    args = list(
        itertools.product(
            primary_inputs,
            md_inputs,
            mec_inputs,
            en_inputs,
            *[[rd] for rd in [rd_1, rd_md, rd_mec, rd_en, rd_agg]],
            *[elec_cost],
            *[heat_cost],
        )
    )

    # import pprint
    # import sys
    # # pprint.pprint(primary_inputs)
    # # pprint.pprint(args[:10])
    # print( len(args))
    # print(len(args) / 8 / 60)  # hr
    # print(sys.getsizeof(args) * 1e-9)
    # assert False

    with multiprocessing.Pool(processes=8) as pool:
        results = pool.starmap(run_kbhdp_zld_mp, args)
    df_results = pd.concat(results)

    timestr = time.strftime("%Y%m%d-%H%M%S")

    df_results.to_csv(
        f"{save_path}/kbhdp_zld_results_pv_trough_storage_cost_frac_from_grid_sweep-{timestr}.csv", index=False
    )


if __name__ == "__main__":
    # main()

    # run_elec_heat_cost_sweep()
    run_recov_chem_mat_cost_sweep()
    # run_energy_param_sweep()
