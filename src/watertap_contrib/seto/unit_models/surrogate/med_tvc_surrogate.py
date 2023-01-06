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
    NonNegativeReals,
    Reference,
    Constraint,
    units as pyunits,
    exp,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

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


@declare_process_block_class("MEDTVCSurrogate")
class MEDTVCData(UnitModelBlockData):
    """
    Multi-effect distillation with thermal vapor compressor model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
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
        "property_package2",
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

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        """
        Add system configurations
        """
        self.Nef = Param(
            initialize=12, units=pyunits.dimensionless, doc="Number of effects"
        )

        self.RR = Var(
            initialize=0.50,
            bounds=(0.30, 0.50),
            units=pyunits.dimensionless,
            doc="Recovery ratio",
        )

        self.Q_loss = Param(
            initialize=0.054, units=pyunits.dimensionless, doc="System thermal loss"
        )

        """
        Add block for feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict
        )

        # Set alias for feed water salinity and convert to mg/L for surrogate model
        self.Xf = pyunits.convert(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            to_units=pyunits.mg / pyunits.L,
        )

        # Set alias for feed water temperature and convert to degree C for surrogate model
        self.Tin = self.feed_props[0].temperature - 273.15 * pyunits.K

        # Set alias for feed water enthalpy and convert to kJ/kg
        h_sw = pyunits.convert(
            self.feed_props[0].enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
        )

        # Set alias for feed water mass flow rate (kg/s)
        m_f = sum(
            self.feed_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.feed_props.component_list
        )

        """
        Add block for distillate
        """
        self.distillate_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of distillate",
            **tmp_dict
        )

        # distillate salinity is 0
        @self.Constraint(doc="distillate salinity")
        def distillate_s(b):
            return b.distillate_props[0].conc_mass_phase_comp["Liq", "TDS"] == 0

        # Connect the plant capacity with the distillate mass flow rate
        self.Capacity = Var(
            initialize=2000,
            bounds=(2000, 100000),
            units=pyunits.m**3 / pyunits.d,
            doc="Capacity of the plant (m3/day)",
        )

        @self.Constraint(doc="distillate mass flow rate")
        def distillate_mfr(b):
            distillate_mfr = b.distillate_props[0].flow_vol_phase["Liq"]
            return b.Capacity == pyunits.convert(
                distillate_mfr, to_units=pyunits.m**3 / pyunits.d
            )

        # distillate temperature is the same as last effect vapor temperature, which is 10 deg higher than condenser inlet seawater temperature
        @self.Constraint(doc="distillate temperature")
        def distillate_temp(b):
            return (
                b.distillate_props[0].temperature
                == self.Tin + (273.15 + 10) * pyunits.K
            )

        # Set alias for distillate enthalpy and convert to kJ/kg
        h_d = pyunits.convert(
            self.distillate_props[0].enth_mass_phase["Liq"],
            to_units=pyunits.kJ / pyunits.kg,
        )

        # Set alias for distillate mass (kg/s) flow rate
        m_d = sum(
            self.distillate_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.distillate_props.component_list
        )

        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

        self.brine_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        # Relationship between brine salinity and feed salinity
        @self.Constraint(doc="Brine salinity")
        def s_b_cal(b):
            return b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"] == b.feed_props[
                0
            ].conc_mass_phase_comp["Liq", "TDS"] / (1 - b.RR)

        # Brine temperature
        @self.Constraint(doc="Brine temperature")
        def T_b_cal(b):
            return (
                b.brine_props[0].temperature
                == b.distillate_props[0].temperature
                + b.brine_props[0].boiling_point_elevation_phase["Liq"]
            )

        # Set alias for brine enthalpy and convert to kJ/kg
        h_b = pyunits.convert(
            self.brine_props[0].enth_mass_phase["Liq"], to_units=pyunits.kJ / pyunits.kg
        )
        # Set alias for brine volume (m3/hr) and mass (kg/s) flow rate
        m_b = sum(
            self.brine_props[0].flow_mass_phase_comp["Liq", j]
            for j in self.brine_props.component_list
        )
        q_b = pyunits.convert(
            self.brine_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        """
        Add block for reject cooling water
        """
        self.cooling_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of cooling reject",
            **tmp_dict
        )

        # Use the source water for cooling (Same salinity)
        @self.Constraint(doc="Cooling reject salinity")
        def s_cooling_cal(b):
            return (
                b.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
            )

        # Set alias for enthalpy and convert to kJ/kg
        h_cool = pyunits.convert(
            self.cooling_out_props[0].enth_mass_phase["Liq"],
            to_units=pyunits.kJ / pyunits.kg,
        )
        # Assumption: the temperature of cooling reject is 3 degC lower than in the condenser
        @self.Constraint(doc="Cooling reject temperature")
        def T_cool_cal(b):
            return (
                b.cooling_out_props[0].temperature
                == b.distillate_props[0].temperature - 3 * pyunits.K
            )

        # Set alias for cooling reject volume flow rate (m3/hr)
        self.q_cooling = pyunits.convert(
            self.cooling_out_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        """
        Add block for heating steam
        """
        tmp_dict["parameters"] = self.config.property_package2
        tmp_dict["defined_state"] = True

        self.steam_props = self.config.property_package2.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of heating steam",
            **tmp_dict
        )

        # Heating steam temperature is fixed at 70 C in this configuration
        @self.Constraint(doc="Heating steam temperature")
        def steam_temp(b):
            return b.steam_props[0].temperature == 70 + 273.15

        # All in steam
        @self.Constraint(doc="Inlet steam status")
        def steam_phase(b):
            return b.steam_props[0].flow_mass_phase_comp["Liq", "H2O"] == 0

        """
        Add block for motive steam
        """
        tmp_dict["parameters"] = self.config.property_package2
        tmp_dict["defined_state"] = True

        self.motive_props = self.config.property_package2.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of motive steam",
            **tmp_dict
        )

        # Set alias for motive steam pressure and convert to bar for surrogate model
        self.Pm = pyunits.convert(self.motive_props[0].pressure, to_units=pyunits.bar)

        # All in steam
        @self.Constraint(doc="Motive steam status")
        def motive_phase(b):
            return b.motive_props[0].flow_mass_phase_comp["Liq", "H2O"] == 0

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="distillate", block=self.distillate_props)
        self.add_port(name="brine", block=self.brine_props)

        """
        Add Vars for model outputs
        """
        self.P_req = Var(
            initialize=5000,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement (kW)",
        )

        self.sA = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**2 / (pyunits.m**3 / pyunits.d),
            doc="Specific area (m2/m3/day))",
        )

        self.STEC = Var(
            initialize=65,
            bounds=(0, None),
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific thermal power consumption (kWh/m3)",
        )

        self.GOR = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg / pyunits.kg,
            doc="Gained output ratio (kg of distillate water per kg of heating steam",
        )

        """
        Mass balances
        """
        # Feed flow rate calculation
        @self.Constraint(doc="Feed water volume flow rate")
        def qF_cal(b):
            return b.feed_props[0].flow_vol_phase["Liq"] == pyunits.convert(
                b.Capacity / b.RR, to_units=pyunits.m**3 / pyunits.s
            )

        # Set alias for feed volume flow rate and covert to m3/hr
        self.qF = pyunits.convert(
            self.feed_props[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.hour,
        )

        # Brine flow rate calculation
        @self.Constraint(doc="Brine volume flow rate")
        def q_b_cal(b):
            return (
                b.brine_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                - b.distillate_props[0].flow_vol_phase["Liq"]
            )

        """
        Add Vars for intermediate model variables
        """
        self.m_sw = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.kg / pyunits.s,
            doc="Feed and cooling water mass flow rate (kg/s)",
        )

        self.q_sw = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.m**3 / pyunits.h,
            doc="Feed and cooling water volume flow rate (m3/h)",
        )

        # Surrogate coefficients for calculating GOR
        Nef_vals = [8, 10, 12, 14, 16]  # Number of effects
        coeffs_set = range(21)  # Number of coefficients

        GOR_coeffs = {
            8: [
                -1.42e-06,
                7.764150858,
                -1.15e-05,
                0.582960027,
                1.24e-07,
                0.056555556,
                -2.57e-19,
                3.63e-25,
                3.23e-19,
                3.04e-21,
                0.104803303,
                -2.08e-09,
                0.012249322,
                -0.000299729,
                6.28e-22,
                -5.098563697,
                -0.001263974,
                -3.56e-25,
                -0.005858519,
                -8.896296296,
                -1.98e-11,
            ],
            10: [
                -2.28e-06,
                10.98911472,
                -1.31e-05,
                0.689378726,
                1.47e-07,
                0.080111111,
                -2.98e-19,
                -6.15e-25,
                4.06e-19,
                4.44e-21,
                0.124314705,
                -1.81e-10,
                0.018184282,
                -0.000334959,
                6.98e-22,
                -6.475922544,
                -0.001524489,
                -5.72e-25,
                -0.006804444,
                -12.28888889,
                -1.68e-11,
            ],
            12: [
                -3.05e-06,
                14.55682317,
                -1.30e-05,
                0.425476523,
                1.61e-07,
                0.097638889,
                6.39e-06,
                2.07e-13,
                -2.09e-06,
                -3.91e-08,
                0.143261656,
                6.78e-10,
                0.027642276,
                -0.00025,
                4.65e-08,
                -2.620877082,
                -0.001921861,
                -5.22e-11,
                -0.001659259,
                -15.9537037,
                -1.83e-11,
            ],
            14: [
                -5.77e-06,
                18.58044758,
                -9.55e-06,
                0.758813856,
                1.78e-07,
                0.148784722,
                2.27e-06,
                -4.93e-13,
                -3.10e-06,
                -3.26e-08,
                0.157319103,
                4.40e-09,
                0.033536585,
                -0.000335027,
                1.13e-08,
                -7.631229505,
                -0.002018089,
                -1.01e-11,
                -0.006560185,
                -21.38310185,
                -7.92e-12,
            ],
            16: [
                -8.38e-06,
                20.76389007,
                -3.04e-06,
                0.967432521,
                1.85e-07,
                0.216269841,
                1.58e-05,
                -6.41e-13,
                -2.27e-06,
                -4.83e-07,
                0.185976144,
                2.71e-09,
                0.046360821,
                -0.000913957,
                -4.66e-08,
                -10.82641732,
                -0.002230429,
                -1.06e-11,
                -0.009318519,
                -25.09448224,
                -4.22e-12,
            ],
        }

        @self.Constraint(doc="GOR surrogate equation")
        def GOR_cal(b):
            return (
                b.GOR
                == b.Xf * GOR_coeffs[b.Nef.value][0]
                + b.RR * GOR_coeffs[b.Nef.value][1]
                + b.Xf * b.RR * GOR_coeffs[b.Nef.value][2]
                + b.Tin * GOR_coeffs[b.Nef.value][3]
                + b.Tin * b.Xf * GOR_coeffs[b.Nef.value][4]
                + b.Tin * b.RR * GOR_coeffs[b.Nef.value][5]
                + b.Capacity * GOR_coeffs[b.Nef.value][6]
                + b.Capacity * b.Xf * GOR_coeffs[b.Nef.value][7]
                + b.Capacity * b.RR * GOR_coeffs[b.Nef.value][8]
                + b.Capacity * b.Tin * GOR_coeffs[b.Nef.value][9]
                + b.Pm * GOR_coeffs[b.Nef.value][10]
                + b.Pm * b.Xf * GOR_coeffs[b.Nef.value][11]
                + b.Pm * b.RR * GOR_coeffs[b.Nef.value][12]
                + b.Pm * b.Tin * GOR_coeffs[b.Nef.value][13]
                + b.Pm * b.Capacity * GOR_coeffs[b.Nef.value][14]
                + 1 * GOR_coeffs[b.Nef.value][15]
                + b.Pm**2 * GOR_coeffs[b.Nef.value][16]
                + b.Capacity**2 * GOR_coeffs[b.Nef.value][17]
                + b.Tin**2 * GOR_coeffs[b.Nef.value][18]
                + b.RR**2 * GOR_coeffs[b.Nef.value][19]
                + b.Xf**2 * GOR_coeffs[b.Nef.value][20]
            )

        # Surrogate coefficients for calculating sA
        sA_coeffs = {
            8: [
                -0.000113631,
                -13.11517461,
                0.000106181,
                -0.494480795,
                2.52e-06,
                0.214458333,
                9.34e-06,
                2.00e-11,
                1.72e-06,
                1.04e-07,
                -0.045572653,
                -6.57e-09,
                -0.000291328,
                6.29e-05,
                -7.80e-13,
                12.94569696,
                0.000878721,
                -1.10e-10,
                0.007633796,
                6.371296296,
                3.41e-10,
            ],
            10: [
                -0.000275164,
                -30.34222168,
                0.000245662,
                -0.964861655,
                5.70e-06,
                0.49675,
                -7.92e-07,
                1.92e-12,
                5.90e-07,
                5.40e-09,
                0.001209732,
                -2.56e-08,
                -0.005934959,
                -1.36e-06,
                1.14e-09,
                25.22049669,
                1.96e-05,
                3.20e-12,
                0.013141389,
                14.14444444,
                8.80e-10,
            ],
            12: [
                -0.00059378,
                -70.68831403,
                0.000537912,
                -2.002884776,
                1.10e-05,
                1.071916667,
                2.61e-05,
                -7.91e-11,
                -8.33e-06,
                -4.78e-07,
                -0.020075596,
                9.86e-08,
                0.014342818,
                0.000726897,
                -3.18e-08,
                52.69512155,
                -0.000184651,
                -6.33e-11,
                0.024886019,
                34.90462963,
                2.02e-09,
            ],
            14: [
                -0.00121164,
                -152.5537433,
                0.001101146,
                -3.461298325,
                2.14e-05,
                2.249461806,
                4.67e-05,
                -1.99e-10,
                -2.11e-05,
                -9.28e-07,
                -0.054534982,
                3.70e-07,
                0.068449356,
                0.001176152,
                -3.36e-08,
                97.79031075,
                -0.000342235,
                -6.33e-11,
                0.038331481,
                78.17853009,
                4.22e-09,
            ],
            16: [
                -0.003702813,
                -566.9002008,
                0.003490278,
                -10.24802298,
                5.87e-05,
                8.029940476,
                0.000180849,
                -1.67e-09,
                -0.000223884,
                -5.19e-06,
                -0.699017368,
                3.49e-06,
                0.438254936,
                0.010243428,
                -2.60e-06,
                312.1661769,
                0.00564047,
                1.39e-09,
                0.098896296,
                317.9402872,
                1.27e-08,
            ],
        }

        @self.Constraint(doc="sA surrogate equation")
        def sA_cal(b):
            return (
                b.sA
                == b.Xf * sA_coeffs[b.Nef.value][0]
                + b.RR * sA_coeffs[b.Nef.value][1]
                + b.Xf * b.RR * sA_coeffs[b.Nef.value][2]
                + b.Tin * sA_coeffs[b.Nef.value][3]
                + b.Tin * b.Xf * sA_coeffs[b.Nef.value][4]
                + b.Tin * b.RR * sA_coeffs[b.Nef.value][5]
                + b.Capacity * sA_coeffs[b.Nef.value][6]
                + b.Capacity * b.Xf * sA_coeffs[b.Nef.value][7]
                + b.Capacity * b.RR * sA_coeffs[b.Nef.value][8]
                + b.Capacity * b.Tin * sA_coeffs[b.Nef.value][9]
                + b.Pm * sA_coeffs[b.Nef.value][10]
                + b.Pm * b.Xf * sA_coeffs[b.Nef.value][11]
                + b.Pm * b.RR * sA_coeffs[b.Nef.value][12]
                + b.Pm * b.Tin * sA_coeffs[b.Nef.value][13]
                + b.Pm * b.Capacity * sA_coeffs[b.Nef.value][14]
                + 1 * sA_coeffs[b.Nef.value][15]
                + b.Pm**2 * sA_coeffs[b.Nef.value][16]
                + b.Capacity**2 * sA_coeffs[b.Nef.value][17]
                + b.Tin**2 * sA_coeffs[b.Nef.value][18]
                + b.RR**2 * sA_coeffs[b.Nef.value][19]
                + b.Xf**2 * sA_coeffs[b.Nef.value][20]
            )

        # Surrogate coefficients for calculating qs
        qs_coeffs = {
            8: [
                -7.59e-06,
                -34.3908768,
                9.69e-06,
                -0.131363139,
                1.48e-07,
                0.200569444,
                0.002087075,
                -1.31e-10,
                -0.000541286,
                -5.70e-06,
                -0.109758232,
                -3.05e-09,
                0.000630081,
                0.000228591,
                -9.38e-08,
                8.029964774,
                0.002098453,
                -2.59e-10,
                0.000841111,
                40.08888889,
                -1.42e-12,
            ],
            10: [
                1.59e-05,
                -29.09136792,
                -2.04e-05,
                -0.066787419,
                -7.21e-08,
                0.145083333,
                0.001764972,
                -1.33e-10,
                -0.000559573,
                -5.46e-06,
                0.005482395,
                -6.00e-08,
                -0.017926829,
                8.64e-06,
                -7.40e-08,
                5.722781223,
                5.07e-05,
                8.14e-12,
                0.00026963,
                36.86296296,
                -6.77e-11,
            ],
            12: [
                -4.48e-06,
                -38.28037765,
                1.02e-07,
                -0.413517098,
                1.86e-07,
                0.304472222,
                0.001604017,
                -1.50e-10,
                -0.00055003,
                -6.73e-06,
                -0.008556022,
                -3.24e-09,
                -0.002256098,
                0.000298408,
                -9.02e-08,
                12.82816936,
                3.81e-05,
                -3.04e-11,
                0.005071204,
                41.34814815,
                -1.16e-11,
            ],
            14: [
                -7.83e-07,
                -40.03225543,
                -7.58e-06,
                -0.359841547,
                1.82e-07,
                0.30859375,
                0.00148361,
                -1.61e-10,
                -0.000593416,
                -6.66e-06,
                -0.001763606,
                -3.31e-09,
                0.001147527,
                0.00014878,
                -6.89e-08,
                11.96567201,
                -5.86e-05,
                -7.44e-11,
                0.004328148,
                45.25318287,
                -2.16e-11,
            ],
            16: [
                3.35e-06,
                -38.68477421,
                -1.90e-05,
                -0.86243177,
                2.00e-07,
                0.233531746,
                0.00126754,
                -1.68e-10,
                -0.000627918,
                -4.45e-06,
                -0.077815222,
                -4.48e-09,
                -0.046564073,
                0.001702812,
                -8.14e-07,
                20.5725045,
                0.002389423,
                5.92e-10,
                0.012013148,
                49.04383976,
                -3.29e-11,
            ],
        }

        @self.Constraint(doc="qs surrogate equation")
        def qs_cal(b):
            return (
                b.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
                == b.Xf * qs_coeffs[b.Nef.value][0]
                + b.RR * qs_coeffs[b.Nef.value][1]
                + b.Xf * b.RR * qs_coeffs[b.Nef.value][2]
                + b.Tin * qs_coeffs[b.Nef.value][3]
                + b.Tin * b.Xf * qs_coeffs[b.Nef.value][4]
                + b.Tin * b.RR * qs_coeffs[b.Nef.value][5]
                + b.Capacity * qs_coeffs[b.Nef.value][6]
                + b.Capacity * b.Xf * qs_coeffs[b.Nef.value][7]
                + b.Capacity * b.RR * qs_coeffs[b.Nef.value][8]
                + b.Capacity * b.Tin * qs_coeffs[b.Nef.value][9]
                + b.Pm * qs_coeffs[b.Nef.value][10]
                + b.Pm * b.Xf * qs_coeffs[b.Nef.value][11]
                + b.Pm * b.RR * qs_coeffs[b.Nef.value][12]
                + b.Pm * b.Tin * qs_coeffs[b.Nef.value][13]
                + b.Pm * b.Capacity * qs_coeffs[b.Nef.value][14]
                + 1 * qs_coeffs[b.Nef.value][15]
                + b.Pm**2 * qs_coeffs[b.Nef.value][16]
                + b.Capacity**2 * qs_coeffs[b.Nef.value][17]
                + b.Tin**2 * qs_coeffs[b.Nef.value][18]
                + b.RR**2 * qs_coeffs[b.Nef.value][19]
                + b.Xf**2 * qs_coeffs[b.Nef.value][20]
            )

        # Surrogate coefficients for calculating qm
        qm_coeffs = {
            8: [
                1.49e-05,
                -29.83707466,
                2.86e-05,
                -2.345398834,
                -9.10e-07,
                0.267902778,
                0.002171651,
                4.00e-10,
                -0.000324519,
                -2.78e-05,
                -0.489303064,
                -7.25e-08,
                0.033783875,
                0.006929505,
                -4.33e-06,
                43.26022001,
                0.005568346,
                -1.65e-10,
                0.035492222,
                28.25555556,
                4.56e-11,
            ],
            10: [
                2.43e-05,
                -27.0360671,
                8.10e-06,
                -1.966013646,
                -8.70e-07,
                0.257819444,
                0.001822688,
                3.13e-10,
                -0.000341778,
                -2.35e-05,
                -0.349886659,
                -9.33e-08,
                0.025067751,
                0.00570271,
                -3.59e-06,
                36.5052648,
                0.003546979,
                4.13e-12,
                0.029556389,
                25.92222222,
                7.53e-12,
            ],
            12: [
                1.12e-05,
                -32.60635053,
                1.54e-05,
                -1.366729008,
                -6.01e-07,
                0.362472222,
                0.001598178,
                2.60e-10,
                -0.000344721,
                -2.04e-05,
                -0.312617943,
                -5.56e-08,
                0.035501355,
                0.004905556,
                -3.20e-06,
                28.5739492,
                0.003185042,
                2.08e-11,
                0.019064444,
                28.63333333,
                3.01e-11,
            ],
            14: [
                1.20e-05,
                -33.97061713,
                7.27e-06,
                -1.345914579,
                -4.99e-07,
                0.39046875,
                0.001457515,
                2.05e-10,
                -0.000376928,
                -1.86e-05,
                -0.283422499,
                -4.47e-08,
                0.04076897,
                0.004456843,
                -2.85e-06,
                28.05931057,
                0.002831752,
                2.13e-11,
                0.018738704,
                30.55700231,
                1.65e-11,
            ],
            16: [
                1.33e-05,
                -33.90192114,
                -1.17e-06,
                -1.544393645,
                -4.29e-07,
                0.364702381,
                0.001274023,
                1.65e-10,
                -0.000401789,
                -1.57e-05,
                -0.282901411,
                -2.79e-08,
                0.014909988,
                0.004355014,
                -3.09e-06,
                31.51108079,
                0.004025387,
                4.05e-10,
                0.021729259,
                33.11413454,
                8.59e-12,
            ],
        }

        @self.Constraint(doc="qm surrogate equation")
        def qm_cal(b):
            return (
                b.motive_props[0].flow_mass_phase_comp["Vap", "H2O"]
                == b.Xf * qm_coeffs[b.Nef.value][0]
                + b.RR * qm_coeffs[b.Nef.value][1]
                + b.Xf * b.RR * qm_coeffs[b.Nef.value][2]
                + b.Tin * qm_coeffs[b.Nef.value][3]
                + b.Tin * b.Xf * qm_coeffs[b.Nef.value][4]
                + b.Tin * b.RR * qm_coeffs[b.Nef.value][5]
                + b.Capacity * qm_coeffs[b.Nef.value][6]
                + b.Capacity * b.Xf * qm_coeffs[b.Nef.value][7]
                + b.Capacity * b.RR * qm_coeffs[b.Nef.value][8]
                + b.Capacity * b.Tin * qm_coeffs[b.Nef.value][9]
                + b.Pm * qm_coeffs[b.Nef.value][10]
                + b.Pm * b.Xf * qm_coeffs[b.Nef.value][11]
                + b.Pm * b.RR * qm_coeffs[b.Nef.value][12]
                + b.Pm * b.Tin * qm_coeffs[b.Nef.value][13]
                + b.Pm * b.Capacity * qm_coeffs[b.Nef.value][14]
                + 1 * qm_coeffs[b.Nef.value][15]
                + b.Pm**2 * qm_coeffs[b.Nef.value][16]
                + b.Capacity**2 * qm_coeffs[b.Nef.value][17]
                + b.Tin**2 * qm_coeffs[b.Nef.value][18]
                + b.RR**2 * qm_coeffs[b.Nef.value][19]
                + b.Xf**2 * qm_coeffs[b.Nef.value][20]
            )

        # Energy consumption
        @self.Constraint(doc="STEC calculation")
        def STEC_cal(b):
            return b.STEC == pyunits.convert(
                1
                / b.GOR
                * (
                    b.motive_props[0].enth_mass_phase["Vap"]
                    - b.steam_props[0].enth_mass_phase["Liq"]
                )
                * b.distillate_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(doc="Thermal power requirement calculation")
        def P_req_cal(b):
            return b.P_req == pyunits.convert(b.STEC * b.Capacity, to_units=pyunits.kW)

        # Mass flow rate
        @self.Constraint(doc="Feed and cooling water mass flow rate (kg/s)")
        def m_sw_cal(b):
            return b.m_sw * (h_cool - h_sw) == (
                (1 - b.Q_loss) * b.P_req - m_b * h_b - m_d * h_d + h_cool * m_f
            )

        # Volume flow rates
        @self.Constraint(doc="Feed and cooling water mass flow rate (m3/h)")
        def q_sw_cal(b):
            return b.q_sw == pyunits.convert(
                b.m_sw / b.feed_props[0].dens_mass_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hour,
            )

        @self.Constraint(doc="Cooling water mass flow rate (m3/h)")
        def q_cooling_cal(b):
            return b.q_cooling == b.q_sw - b.qF

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.RR) is None:
            iscale.set_scaling_factor(self.RR, 1e1)

        if iscale.get_scaling_factor(self.Capacity) is None:
            iscale.set_scaling_factor(self.Capacity, 1e-3)

        if iscale.get_scaling_factor(self.STEC) is None:
            iscale.set_scaling_factor(self.STEC, 1e-3)

        if iscale.get_scaling_factor(self.P_req) is None:
            iscale.set_scaling_factor(
                self.P_req,
                iscale.get_scaling_factor(self.Capacity)
                * iscale.get_scaling_factor(self.STEC),
            )

        if iscale.get_scaling_factor(self.sA) is None:
            iscale.set_scaling_factor(self.sA, 1e-1)

        if iscale.get_scaling_factor(self.GOR) is None:
            iscale.set_scaling_factor(self.GOR, 1e-1)

        if (
            iscale.get_scaling_factor(
                self.feed_props[0].flow_mass_phase_comp["Liq", "H2O"]
            )
            is None
        ):
            iscale.set_scaling_factor(
                self.feed_props[0].flow_mass_phase_comp["Liq", "H2O"], 1e1
            )

        if iscale.get_scaling_factor(self.m_sw) is None:
            iscale.set_scaling_factor(self.m_sw, 1e-3)

        if iscale.get_scaling_factor(self.q_sw) is None:
            iscale.set_scaling_factor(self.q_sw, 1e-3)

        # Transforming constraints
        sf = iscale.get_scaling_factor(
            self.distillate_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.distillate_s, sf)

        sf = iscale.get_scaling_factor(self.Capacity)
        iscale.constraint_scaling_transform(self.distillate_mfr, sf)

        sf = iscale.get_scaling_factor(self.distillate_props[0].temperature)
        iscale.constraint_scaling_transform(self.distillate_temp, sf)

        sf = iscale.get_scaling_factor(self.cooling_out_props[0].temperature)
        iscale.constraint_scaling_transform(self.T_cool_cal, sf)

        sf = iscale.get_scaling_factor(
            self.cooling_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.s_cooling_cal, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].temperature)
        iscale.constraint_scaling_transform(self.T_b_cal, sf)

        sf = iscale.get_scaling_factor(
            self.brine_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.s_b_cal, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].temperature)
        iscale.constraint_scaling_transform(self.steam_temp, sf)

        sf = iscale.get_scaling_factor(
            self.steam_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.steam_phase, sf)

        sf = iscale.get_scaling_factor(
            self.motive_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.motive_phase, sf)

        sf = iscale.get_scaling_factor(self.feed_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.qF_cal, sf)

        sf = iscale.get_scaling_factor(self.brine_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.q_b_cal, sf)

        sf = iscale.get_scaling_factor(self.GOR)
        iscale.constraint_scaling_transform(self.GOR_cal, sf)

        sf = iscale.get_scaling_factor(self.sA)
        iscale.constraint_scaling_transform(self.sA_cal, sf)

        sf = iscale.get_scaling_factor(
            self.steam_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        iscale.constraint_scaling_transform(self.qs_cal, sf)

        sf = iscale.get_scaling_factor(
            self.motive_props[0].flow_mass_phase_comp["Vap", "H2O"]
        )
        iscale.constraint_scaling_transform(self.qm_cal, sf)

        sf = iscale.get_scaling_factor(self.STEC)
        iscale.constraint_scaling_transform(self.STEC_cal, sf)

        sf = iscale.get_scaling_factor(self.P_req)
        iscale.constraint_scaling_transform(self.P_req_cal, sf)

        sf = (
            iscale.get_scaling_factor(self.m_sw)
            * iscale.get_scaling_factor(
                self.cooling_out_props[0].enth_mass_phase["Liq"]
            )
            * 1e3
        )
        iscale.constraint_scaling_transform(self.m_sw_cal, sf)

        sf = iscale.get_scaling_factor(self.q_sw)
        iscale.constraint_scaling_transform(self.q_sw_cal, sf)

        sf = (
            iscale.get_scaling_factor(self.cooling_out_props[0].flow_vol_phase["Liq"])
            / 3600
        )
        iscale.constraint_scaling_transform(self.q_cooling_cal, sf)
