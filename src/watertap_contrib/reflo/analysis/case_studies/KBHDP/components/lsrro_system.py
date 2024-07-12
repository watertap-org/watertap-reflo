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
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
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
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)

# from analysisWaterTAP.utils.flowsheet_utils import *
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)

# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')


_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

def build_ro(m, blk, number_of_stages=1, prop_package=None) -> None:
    pass


def init_system(m, verbose=True, solver=None):
    pass


def init_ro_system(m, blk, verbose=True, solver=None):
    pass


def init_lsrro_stage(m, stage, solver=None):
    pass


def set_operating_conditions(
    m, Qin=None, Qout=None, Cin=None, water_recovery=None
):
    pass


def calc_scale(value):
    return math.floor(math.log(value, 10))


def set_lsrro_system_operating_conditions(
    m, blk, mem_area=100, RO_pump_pressure=15e5
):
    pass


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results


def get_sub_blocks(block, decend=False, report=False):
    blocks = []
    for v in block.component_data_objects(
        ctype=Block, active=True, descend_into=decend
    ):
        print(v)
        if report:
            try:
                table = v._get_stream_table_contents()
                for item in table:
                    print(table[item])
                # print(v._get_stream_table_contents())
            except:
                pass


def display_lsrro_system_build(m):
    get_sub_blocks(m.fs)
    get_sub_blocks(m.fs.lsrro)
    for stage in m.fs.lsrro.stage.values():
        get_sub_blocks(stage)
    print("\n")


def display_dof_breakdown(blk, decend=False, report=False):
    print(
        "\n\n-------------------- DEGREE OF FREEDOM BREAKDOWN --------------------\n\n"
    )
    print(f'{"BLOCK":<40s}{"DEGREES OF FREEDOM":<30s}')
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(f"{v.name:<40s}{degrees_of_freedom(v)}")


def display_inlet_conditions(blk):
    print("\n\n")
    # print(blk.feed.display())

    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    # assert False


def display_flow_table(blk):
    print("\n\n")
    print(
        f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
    )
    print(
        f'{"Feed":<34s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Product":<34s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )
    print(
        f'{"Disposal":<34s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<30.1f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
    )

    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )
    for idx, stage in blk.stage.items():
        print(
            f'{"RO Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.3f}'
        )


def report_LSRRO(m, blk):
    print(f"\n\n-------------------- RO Report --------------------\n")
    print(f'{"Recovery":<30s}{value(100*m.fs.water_recovery):<10.1f}{"%"}')
    print(f'{"RO Operating Pressure":<30s}{value(pyunits.convert(blk.pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}')


def build_system():
    m = ConcreteModel()

    return m


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    display_ro_system_build(m)
    set_operating_conditions(m, Qin=1, Cin=2.5)
    set_ro_system_operating_conditions(m, m.fs.ro, mem_area=10)
    init_system(m)
    solve(m)
    
    display_flow_table(m.fs.ro)
    report_RO(m, m.fs.ro)
    # print(m.fs.ro.stage[1].module.report())
    # print(m.fs.costing.display())

#FIX this flowsheet needs to get converted to a lsrro system