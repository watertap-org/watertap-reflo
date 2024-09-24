#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from copy import deepcopy
from pyomo.environ import (
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    units as pyunits,
    Expr_if,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core.util.math import smooth_min, smooth_max, smooth_bound

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing.units.deep_well_injection import (
    cost_deep_well_injection,
)

__author__ = "Kurban Sitterley"

_log = idaeslog.getLogger(__name__)


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


@declare_process_block_class("DeepWellInjection")
class DeepWellInjectionData(InitializationMixin, UnitModelBlockData):
    """
    Zero deep well injection model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
        ),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    CONFIG.declare(
        "well_depth",
        ConfigValue(
            default=5000,
            domain=In([2500, 5000, 7500, 10000]),
            description="Depth of injection well. Costing is available for 2500, 5000, 7500, and 10000 ft depths.",
            doc="""Depth of injection well for costing purposes.""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        self.add_port(name="inlet", block=self.properties)

        self.pipe_diameter_coeff = Param(
            initialize=5.329,
            units=pyunits.dimensionless,
            doc="Base parameter for pipe diameter equation",
        )

        self.pipe_diameter_exponent = Param(
            initialize=0.4998,
            units=pyunits.dimensionless,
            doc="Exponent parameter for pipe diameter equation",
        )

        self.well_depth = Param(
            initialize=self.config.well_depth,
            units=pyunits.ft,
            doc="Depth of injection well",
        )

        self.flow_mgd = pyunits.convert(
            self.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.Mgallons / pyunits.day,
        )

        @self.Expression(doc="Injection pipe diameter in inches")
        def pipe_diameter(b):
            flow_mgd_dimensionless = pyunits.convert(
                pyunits.convert(
                    b.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.Mgallons / pyunits.day,
                )
                * pyunits.day
                * pyunits.Mgallons**-1,
                to_units=pyunits.dimensionless,
            )
            pipe_diameter = (
                b.pipe_diameter_coeff * flow_mgd_dimensionless**b.pipe_diameter_exponent
            )
            return smooth_bound(pipe_diameter, 2, 24) * pyunits.inches
            # return pipe_diameter
    
    @property
    def default_costing_method(self):
        return cost_deep_well_injection

if __name__ == "__main__":

    from pyomo.environ import (
        ConcreteModel,
        Expression,
        TerminationCondition,
        SolverStatus,
        value,
        Var,
        assert_optimal_termination,
        units as pyunits,
    )
    from pyomo.network import Port
    from idaes.core import FlowsheetBlock
    from pyomo.util.check_units import assert_units_consistent

    from watertap.core.solvers import get_solver

    from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
    from idaes.core.util.model_statistics import (
        degrees_of_freedom,
        number_variables,
        number_total_constraints,
        number_unused_variables,
    )
    from idaes.core.util.testing import initialization_tester
    from idaes.core.util.scaling import (
        calculate_scaling_factors,
        unscaled_variables_generator,
        badly_scaled_var_generator,
    )
    from idaes.core import UnitModelCostingBlock

    from watertap.property_models.multicomp_aq_sol_prop_pack import (
        MCASParameterBlock,
    )
    from watertap_contrib.reflo.costing import TreatmentCosting

    from pprint import pprint
    from watertap_contrib.reflo.costing import TreatmentCosting
    import watertap.property_models.unit_specific.cryst_prop_pack as props
    from watertap.property_models.water_prop_pack import WaterParameterBlock
    from watertap.core.util.model_diagnostics.infeasible import *
    from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
    from idaes.core.util.scaling import *
    from idaes.core.util.testing import initialization_tester
    from idaes.core.util.scaling import (
        calculate_scaling_factors,
        unscaled_variables_generator,
        badly_scaled_var_generator,
    )
    from watertap_contrib.reflo.kurby import print_unit_solutions, make_test_dict

    solver = get_solver()

    inlet_conc = {
        "Ca_2+": 1.43,
        "Mg_2+": 0.1814,
        "SiO2": 0.054,
        "Alkalinity_2-": 0.421,
    }

    flow_mgd = 5.08 * pyunits.Mgallons / pyunits.day
    flow_mgd = 21 * pyunits.Mgallons / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3

    flow_mass_water_in = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(solute_list=inlet_conc.keys())
    m.fs.dwi = dwi = DeepWellInjection(property_package=m.fs.properties)
    m.fs.costing = TreatmentCosting()
    m.fs.dwi.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(dwi.properties[0].flow_vol_phase["Liq"])

    prop_in = dwi.properties[0]
    prop_in.temperature.fix()
    prop_in.pressure.fix()
    flow_mass_phase_water = pyunits.convert(
        flow_mgd * rho, to_units=pyunits.kg / pyunits.s
    )
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
    for solute, conc in inlet_conc.items():
        mass_flow_solute = pyunits.convert(
            flow_mgd * conc * pyunits.kg / pyunits.m**3,
            to_units=pyunits.kg / pyunits.s,
        )
        prop_in.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / mass_flow_solute),
            index=("Liq", solute),
        )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / flow_mass_phase_water),
        index=("Liq", "H2O"),
    )

    results = solver.solve(m)

    print(f"pipe_diameter = {dwi.pipe_diameter()} {pyunits.get_units(dwi.pipe_diameter)}")
    print()
    print(f"DOF = {degrees_of_freedom(m)}")
    print(f"LCOW = {m.fs.costing.LCOW()}")
    # print(
    #     dwi.flow_mgd(),
    #     flow_mass_water_in(),
    #     dwi.properties[0].flow_mass_phase_comp["Liq", "H2O"](),
    #     dwi.properties[0].flow_vol_phase["Liq"](),
    # )
    # dwi.well_depth.display()

    # dwi.costing.display()
