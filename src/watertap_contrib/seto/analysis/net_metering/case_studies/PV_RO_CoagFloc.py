import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker, cm
from matplotlib import ticker
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
import watertap_contrib.seto.analysis.net_metering.PV_RO_surrogate as PV_RO
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Param,
    Constraint,
    Expression,
    Block,
    Var,
    TransformationFactory,
    assert_optimal_termination,
    check_optimal_termination,
    value,
    units as pyunits,
)

from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models import Product, Feed, Translator
from idaes.core.util.model_statistics import *
from watertap.core.wt_database import Database

from idaes.core import UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import *
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)
from watertap.core.util.infeasible import *

from watertap_contrib.seto.analysis.net_metering.util import *
from watertap_contrib.seto.analysis.net_metering.case_studies.utils import *

from watertap_contrib.seto.costing import (
    TreatmentCosting,
    EnergyCosting,
    SETOSystemCosting,
)
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
from watertap_contrib.seto.solar_models.surrogate.pv import PVSurrogate
from watertap.unit_models.coag_floc_model import CoagulationFlocculation
from watertap.unit_models.zero_order.coag_and_floc_zo import CoagulationFlocculationZO
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting



def coag_floc_setup(m):
    """Builds the structure of the PV-RO system

    Returns:
        object: A Pyomo concrete optimization model and flowsheet
    """
    print(f'\n{"=======> BUILDING COAGULATION-FLOCCULATION <=======":^60}\n')

    # m = ConcreteModel()
    # m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = WaterParameterBlock(solute_list=["tds", "tss"])
    m.db = Database()
    
    m.fs.coag_and_floc = CoagulationFlocculationZO(
        property_package=m.fs.props, database=m.db
    )

    return m

def coag_floc_set_op_conds(m, Q, conc):
    m.fs.coag_and_floc.properties[0.0].flow_vol = Q * pyunits.m ** 3 / pyunits.s
    m.fs.coag_and_floc.inlet.flow_mass_comp[0, "H2O"].fix((1-conc/1000))
    m.fs.coag_and_floc.inlet.flow_mass_comp[0, "tss"].fix(0)
    m.fs.coag_and_floc.inlet.flow_mass_comp[0, "tds"].fix((conc/1000))
    # m.fs.coag_and_floc.inlet.flow_mass_comp[0, "H2O"].fix(1000*Q*(1-conc/1000))
    # m.fs.coag_and_floc.inlet.flow_mass_comp[0, "tss"].fix(0)
    # m.fs.coag_and_floc.inlet.flow_mass_comp[0, "tds"].fix(1000*Q*(conc/1000))

def coag_floc_load_params(m):
    data = m.db.get_unit_operation_parameters("coag_and_floc")
    m.fs.coag_and_floc.load_parameters_from_database()

def coag_floc_costing(m):
    m.fs.coag_floc_costing = ZeroOrderCosting()    
    m.fs.coag_and_floc.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.coag_floc_costing)

    m.fs.coag_and_floc.costing.LCOW = pyo.Var(
            initialize=0,
            units= pyunits.USD_2021 / pyunits.year,
            doc="Levelized cost of water",
        )
    
    m.fs.coag_and_floc.costing.LCOW_constraint = Constraint(
        expr=m.fs.coag_and_floc.costing.LCOW == (
        (pyunits.convert(m.fs.coag_and_floc.costing.capital_cost, to_units=pyunits.USD_2021)/
         m.fs.coag_floc_costing.plant_lifetime) /
        pyunits.convert(m.fs.coag_and_floc.properties[0.0].flow_vol, to_units=pyunits.m**3 / pyunits.year)
        )
    )

def coag_floc_init(m):
    print(degrees_of_freedom(m.fs.coag_and_floc))
    m.fs.coag_and_floc.initialize()
    
def coag_floc_solve(m):
    solver = get_solver()
    results = solver.solve(m)

    return results

def build_coag_floc(m, Q, conc):
    coag_floc_setup(m)
    coag_floc_set_op_conds(m, Q, conc)
    coag_floc_load_params(m)
    coag_floc_costing(m)
    coag_floc_init(m)

    return m

def build_system(Q_basis, conc_val, rec_val):
    m = PV_RO.model_setup(Q_basis, conc_val, rec_val)
    build_coag_floc(m, Q_basis, conc_val)

    return m

def solve_system(m):
    m, results = PV_RO.run(m)
    print_coag_floc_results(m)
    return m, results

def main():
    m = build_system(mgd_to_m3s(10), 35, 0.5)
    m, results = solve_system(m)

def param_sweep(Q_basis: int, x_bounds: list, y_bounds: list, x_res: int, y_res: int):
    conc_vals = np.linspace(*y_bounds, y_res)
    rec_vals = np.linspace(*x_bounds, x_res)
    temp_dict = {}
    count = 0
    solve_report = []
    for idx1, conc_val in enumerate(conc_vals):
        for idx2, rec_val in enumerate(rec_vals):
            try:
                m = build_system(mgd_to_m3s(Q_basis), conc_val, rec_val)

                try:
                    m, results = solve_system(m)
                    case = 'Optimal'
                    encoder = 1
                    lcow = m.fs.sys_costing.LCOW() + m.fs.coag_and_floc.costing.LCOW()
                    recovery = m.fs.treatment.ro.recovery_vol_phase[0, "Liq"]()
                    solve_report.append([f'Optimal solve for {Q_basis} MGD, {conc_val} g/L, {rec_val} recovery: LCOW={lcow}'])
                except:
                    print('Something went wrong during solving')
                    m, results = m, np.nan
                    lcow = np.nan
                    recovery = np.nan
                    case = 'Infeasible'
                    encoder = 0
                    solve_report.append([f'Infeasible solve for {Q_basis} MGD, {conc_val} g/L, {rec_val} recovery'])
                    # continue

            except:
                print('Something went wrong during initialization')
                m, results = m, np.nan
                encoder = 0
                lcow = np.nan
                recovery = np.nan
                solve_report.append([f'Failed Initialization for {Q_basis} MGD, {conc_val} g/L, {rec_val} recovery'])
                
            temp_dict['Scenario {}'.format(count)]={'idx1':idx1,'idx2':idx2,'var1':conc_val, 'var2':rec_val, 'model':m, 'results':results, 'Condition':case, 'Conditional':encoder, 'LCOW':lcow, 'Rec':recovery}
            count += 1

    return m, results, temp_dict, solve_report

if __name__ == "__main__":
    m, results, temp_dict, solve_report = param_sweep(Q_basis=10, x_bounds=[0.3, 0.5], y_bounds=[30, 60], x_res=3, y_res=3)
    x, y, z, grid = get_data(temp_dict, 'LCOW', gridshape = (3, 3))
    fig, ax, cbar = contour_LCOW(x, y, grid, x_label='Recovery (%)', y_label='TDS (g/L)', z_label='LCOW '+'USD/m$^3$', x_scale='linear', low=min(z), mid=(min(z)+max(z))/2, high=max(z), xlimits = [min(x), max(x)], ylimits = [min(y), max(y)], cmap_pallete="RdBu_r", contour_x_pos=0.4, contour_label_space=20, auto_ticks=np.linspace(round(min(z),1), round(max(z),2), 6).tolist())
    plt.show()