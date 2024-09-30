#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.testing import initialization_tester
from idaes.core.util.exceptions import InitializationError, ConfigurationError
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)

from watertap_contrib.reflo.unit_models import SolarStill
from watertap_contrib.reflo.costing import TreatmentCosting

# Get default solver for testing
solver = get_solver()


def build_ss():

    rho = 1000 * pyunits.kg / pyunits.m**3
    inlet_dict = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 31.4038218e-3},
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    water_yield_calc_dict = dict(
        input_weather_file_path="/Users/ksitterl/Documents/SETO/models/solar_still/SS_model-Sept2024/SS_Model/Data/TMY2 SAM CSV/data.csv",  # path to input weather file
        interval_day=150,  # interval at which to calculate water yield (e.g., every 15 days)
        salinity=20,  # salinity of influent water; g/L
        water_depth_basin=0.02,  # depth of water in solar still basin; m
        length_basin=0.6,  # length of each side of basin (length=width); m
    )

    tds_conc = 20 * pyunits.g / pyunits.liter
    daily_water_production = 100 * pyunits.m**3 / pyunits.day

    flow_mass_in = pyunits.convert(
        daily_water_production * rho, to_units=pyunits.kg / pyunits.s
    )
    flow_mass_tds = pyunits.convert(
        daily_water_production * tds_conc, to_units=pyunits.kg / pyunits.s
    )
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**inlet_dict)

    m.fs.unit = SolarStill(
        property_package=m.fs.properties,
        water_yield_calculation_args=water_yield_calc_dict,
    )

    m.fs.unit.properties_in[0].flow_vol_phase[...]
    m.fs.unit.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(value(flow_mass_in))
    m.fs.unit.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix(
        value(flow_mass_tds)
    )
    m.fs.unit.properties_in[0].pressure.fix(101325)
    m.fs.unit.properties_in[0].temperature.fix(293.15)

    return m


class TestSolarStill:
    @pytest.fixture(scope="class")
    def ss_frame(self):
        m = build_ss()
        return m
