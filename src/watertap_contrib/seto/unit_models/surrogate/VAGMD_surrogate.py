###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    value,
    Suffix,
    ConcreteModel,
    PositiveIntegers,
    Reference,
    Constraint,
    units as pyunits,
    exp,
    log,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

from pyomo.dae import ContinuousSet

# Import IDAES cores
from idaes.core import (
    FlowsheetBlock,
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

# Import Watertap packages
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import WaterParameterBlock

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)


@declare_process_block_class("VAGMDsurrogate")
class VAGMDData(UnitModelBlockData):
    """
    Vacuumed Membrane distillation (air-gapped) - batch operation model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be True",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
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
        "high_brine_salinity",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="""A ConfigBlock with arguments to be passed to indicate if the brine
            has a salinity higher than 175.3 g/L""",
        ),
    )
    CONFIG.declare(
        "module_type",
        ConfigValue(
            default="AS7C1.5L",
            domain=In(["AS7C1.5L", "AS26C7.2L"]),
            doc="""Selection of module type (7 for AS7C1.5L and 26 for AS26C7.2L)""",
        ),
    )
    CONFIG.declare(
        "cooling_system_type",
        ConfigValue(
            default="closed",
            domain=In(["open", "closed"]),
            doc="""Selection of cooling system type (open or closed)""",
        ),
    )
    CONFIG.declare(
        "number_cycles",
        ConfigValue(
            default=10,
            domain=int,
            doc="""Number of cycles required""",
        ),
    )

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        """
        MD Module type and corresponding area
        """
        if self.config.module_type == "AS7C1.5L":
            module_area = 7.2
            self.module_area = Param(
                initialize=7.2, units=pyunits.m**2, doc="Area of module AS7C1.5L"
            )
        else:  # module_type = 'AS26C7.2L'
            module_area = 25.92
            self.module_area = Param(
                initialize=25.92, units=pyunits.m**2, doc="Area of module AS26C7.2L"
            )

        """
        Add block for the batched feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of feed water",
            **tmp_dict
        )        

    """
    Equation to calculate pressure drop
    """
    def _get_pressure_drop(self, flow_rate, salinity):
        if self.config.module_type == 'AS7C1.5L':
            coefficients = [-158.2007422,
                            0.39402609,
                            0,
                            0.000585345,
                            8.93618e-5,
                            -0.000287828]
        else: # self.config.module_type == 'AS26C7.2L'
            coefficients = [-72.53793298,
                            0.110437201,
                            0,
                            0.000643495,
                            0.000189924,
                            -0.001111447]

        return   (coefficients[0] 
                + coefficients[1] * flow_rate
                + coefficients[2] * salinity
                + coefficients[3] * flow_rate * salinity
                + coefficients[4] * flow_rate**2
                + coefficients[5] * salinity**2)


    """
    Surrogate equations for membrane performance
    """
    def _get_membrane_performance(self, TEI, FFR, TCI, SgL, Ttank):
        # Model parameters
        PFluxAS26 = [
            0.798993148477908,
            0.314627216640160,
            0.559805181621833,
            -0.146236734128216,
            -0.659197144919924,
            0,
            0.185658514024503,
            -0.107221706014227,
            0,
            0,
            -0.187626469717738,
            0,
            0,
            0,
            0.128420664686447,
        ]
        PFluxAS7_high = [
            9.41014468300000,
            0,
            0,
            0,
            -0.0188989390000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        PFluxAS7_low = [
            4.82699491400000,
            1.37479276300000,
            1.91988177200000,
            -0.574212905000000,
            -0.641257664000000,
            0.399259954000000,
            0,
            0,
            0,
            0,
            0,
            0,
            -0.588924321000000,
            0,
            0,
        ]
        TCOAS26 = [
            65.1084685465240,
            9.15474718837607,
            -0.917460918908258,
            0.480070517276181,
            -1.06168979823129,
            0,
            0,
            0,
            0.142552983811052,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_high = [
            67.0068599900000,
            0,
            0,
            0,
            -0.0145469190000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_low = [
            58.8189910000000,
            8.79691366400000,
            -2.06741401400000,
            1.63187967600000,
            -0.914624645000000,
            -0.536574144000000,
            -0.249657477000000,
            0.398973861000000,
            -0.153760262000000,
            0.102355281000000,
            0.696768080000000,
            -0.300582958000000,
            -0.557410173000000,
        ]
        TEOAS26 = [
            29.2261439896435,
            0.569016152083381,
            0.824636694807529,
            4.62669502530487,
            1.37222105534565,
            0,
            0,
            0.220665657590258,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_high = [
            36.2497021800000,
            0,
            0,
            0,
            0.0126951860000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_low = [
            34.4309511000000,
            1.55140768500000,
            1.85928314600000,
            4.52887180500000,
            1.10791196800000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0.808447211000000,
        ]
        Coder = [
            [-5.68487382500000, 0.0705622560000000, 0.000152146000000000],
            [-1.58460599600000, 0.00102338700000000, 1.20000000000000e-06],
            [-4.27697973100000, 0.175533630000000, -0.000178178000000000],
            [-1.49331349500000, 0.0146627780000000, 5.62000000000000e-06],
            [-7, 0.100000000000000, 0],
            [-2.14285714285714, 0.00285714285714286, 0],
            [-5, 0.200000000000000, 0],
            [-1.33333333333333, 0.00950000000000000, 0],
        ]
        a1 = 0.983930048493388
        a2 = -4.8359231959954e-04
        S_c = a1 * SgL + a2 * SgL**2  # [g/kg]

        TEI -= 273.15
        TCI -= 273.15
        CoderVars = [
            [1, TEI, TEI**2],
            [1, FFR, FFR**2],
            [1, TCI, TCI**2],
            [1, S_c, S_c**2],
        ]

        # Model calculations
        if self.config.module_type == "AS7C1.5L":
            if self.config.high_brine_salinity:
                TEI = 0
                FFR = 0
                TCI = 0
                S_r = S_c

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high
            else:
                TEI = sum(CoderVars[0][j] * Coder[0][j] for j in range(len(Coder[0])))
                FFR = sum(CoderVars[1][j] * Coder[1][j] for j in range(len(Coder[1])))
                TCI = sum(CoderVars[2][j] * Coder[2][j] for j in range(len(Coder[2])))
                S_r = sum(CoderVars[3][j] * Coder[3][j] for j in range(len(Coder[3])))

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_low, TCOAS7_low, TEOAS7_low

            VarsAS7 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                TCI * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]
            VarsAS7_TCO = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                FFR**2,
                S_r**2,
                FFR**3,
            ]

            PFlux = sum(VarsAS7[j] * PFluxAS7[j] for j in range(len(VarsAS7)))
            TCO = sum(VarsAS7_TCO[j] * TCOAS7[j] for j in range(len(VarsAS7_TCO)))
            TEO = sum(VarsAS7[j] * TEOAS7[j] for j in range(len(VarsAS7)))

        else:
            TEI = sum(CoderVars[0][j] * Coder[4][j] for j in range(len(Coder[0])))
            FFR = sum(CoderVars[1][j] * Coder[5][j] for j in range(len(Coder[1])))
            TCI = sum(CoderVars[2][j] * Coder[6][j] for j in range(len(Coder[2])))
            S_r = sum(CoderVars[3][j] * Coder[7][j] for j in range(len(Coder[3])))

            VarsAS26 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                TCI * TEI,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                S_r * TCI,
                FFR * S_r,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]

            PFlux = sum(VarsAS26[j] * PFluxAS26[j] for j in range(len(VarsAS26)))
            TCO = sum(VarsAS26[j] * TCOAS26[j] for j in range(len(VarsAS26)))
            TEO = sum(VarsAS26[j] * TEOAS26[j] for j in range(len(VarsAS26)))

        return [PFlux, TCO, TEO, S_c]
