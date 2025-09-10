#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.math import smooth_min, smooth_max
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.units.chemical_softening import (
    cost_chemical_softening,
)

__author__ = "Mukta Hardikar, Abdiel Lugo, Kurban Sitterley, Zachary Binger"


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


@declare_process_block_class("ChemicalSoftening")
class ChemicalSofteningData(InitializationMixin, UnitModelBlockData):
    """
    Chemical softening model
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
        "softening_procedure_type",
        ConfigValue(
            default=SofteningProcedureType.single_stage_lime,
            domain=In(SofteningProcedureType),
            description="Procedure of softening system",
        ),
    )

    CONFIG.declare(
        "silica_removal",
        ConfigValue(
            default=False,
            domain=bool,
            description="Include silica removal in softening process -- assumes conditions are proper",
        ),
    )

    def build(self):
        super().build()

        """
        1.  Crittenden, J. C., & Montgomery Watson Harza (Firm). (2012). Water treatment principles and design. Hoboken, N.J: J.Wiley.
        2.  Davis, M. L. (2010). Water and wastewater engineering: Design principles and practice.
        3.  Baruth. (2005). Water treatment plant design / American Water Works Association, American Society of Civil Engineers; Edward E. Baruth, technical editor. (Fourth edition.). McGraw-Hill.
        4.  Edzwald, J. K., & American Water Works Association. (2011). Water quality & treatment: A handbook on drinking water. New York: McGraw-Hill.
        5.  R.O. Mines Environmental Engineering: Principles and Practice, 1st Ed, John Wiley & Sons
        6.  Lee, C. C., & Lin, S. D. (2007). Handbook of environmental engineering calculations. New York: McGraw Hill.
        """

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        required_comps = ["Ca_2+", "Mg_2+", "Alkalinity_2-"]

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        comps = self.config.property_package.solute_set
        if not all(rc in comps for rc in required_comps):
            raise ConfigurationError(
                "ChemicalSoftening requires Ca_2+, Mg_2+, and Alkalinity_2- as solutes in inlet stream but not all were provided."
            )
        non_hardness_comps = [
            j
            for j in self.config.property_package.solute_set
            if j not in ["Ca_2+", "Mg_2+", "Alkalinity_2-"]
        ]

        # MW of components
        self.CaO_mw = Param(
            initialize=56,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of CaO (dry lime)",
        )

        self.Ca_mw = Param(
            initialize=40,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of Ca",
        )

        self.CaCO3_mw = Param(
            initialize=100,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of CaCO3",
        )

        self.Na2CO3_mw = Param(
            initialize=106,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of soda",
        )

        self.CO2_mw = Param(
            initialize=44,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of carbon dioxide",
        )

        self.Mg_mw = Param(
            initialize=24.3,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of Mg2+",
        )

        self.MgCl2_mw = Param(
            initialize=95.2,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of MgCl2",
        )

        self.MgOH2_mw = Param(
            initialize=58.3,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of Mg(OH)2",
        )

        self.SiO2_check_conv = Param(
            initialize=2.35,
            units=pyunits.dimensionless,
            doc="Silica conversion factor to compare with Mg2+ concentration",
        )

        self.MgCl2_SiO2_ratio = Param(
            initialize=5.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ratio of Mg2+ to SiO2",
        )

        self.Ca_hardness_CaCO3_sludge_factor = Param(
            initialize=2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Sludge produced per kg of Ca in CaCO3 hardness",
        )

        self.Mg_carbonate_hardness_sludge_factor = Param(
            initialize=2.6,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Sludge produced per kg of Mg in CaCO3 hardness",
        )

        self.Mg_noncarbonate_hardness_sludge_factor = Param(
            initialize=1.6,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Sludge produced per kg of Mg in non-CaCO3 hardness",
        )

        self.excess_CaO_coeff = Param(
            initialize=0.05,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Multiplication factor to calculate excess CaO dose",
        )

        self.Ca_CaCO3_conv = Param(
            initialize=2.5,
            units=pyunits.dimensionless,
            doc="Conversion factor for Ca to equivalent CaCO3",
        )

        self.Mg_CaCO3_conv = Param(
            initialize=4.12,
            units=pyunits.dimensionless,
            doc="Conversion factor for Mg to equivalent CaCO3",
        )

        self.number_mixers = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of mixers",
        )

        self.number_floc = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of flocculators",
        )

        self.eps = Param(
            initialize=1e-18,
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Smoothing factor",
        )

        self.ca_eff_target = Var(
            initialize=0.20,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Target Ca concentration in the effluent",
        )

        self.mg_eff_target = Var(
            initialize=0.10,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Target Mg concentration in the effluent",
        )

        # Creating default effluent composition of all solutes other than Ca and Mg
        removal_eff_dict = dict(
            zip(
                non_hardness_comps,
                [
                    0.85 if j == "TSS" else 0.7 if j != "TDS" else (1 - 1e-3)
                    for j in non_hardness_comps
                ],
            )
        )

        self.removal_efficiency = Var(
            non_hardness_comps,
            initialize=removal_eff_dict,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Removal efficiency of all solutes other than Ca and Mg",
        )

        # System design variables

        self.retention_time_mixer = Var(
            initialize=1,
            bounds=(0.1, 5),
            units=pyunits.minutes,
            doc="Retention time for mixer",
        )

        self.retention_time_floc = Var(
            initialize=12,
            bounds=(10, 45),
            units=pyunits.minutes,
            doc="Retention time for flocculator",
        )

        self.retention_time_sed = Var(
            initialize=150,
            bounds=(120, 240),
            units=pyunits.minutes,
            doc="Retention time for sedimentation basin",
        )

        self.retention_time_recarb = Var(
            initialize=18,
            bounds=(15, 30),
            units=pyunits.minutes,
            doc="Retention time for recarbonation",
        )

        self.sedimentation_overflow = Var(
            initialize=90,
            bounds=(30, 100),
            units=pyunits.m / pyunits.day,
            doc="Sedimentation basin overflow rate",
        )

        self.volume_mixer = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of single mixer",
        )

        self.volume_floc = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of single flocculator",
        )

        self.volume_sed = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of sedimentation basin",
        )

        self.volume_recarb = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Volume of recarbonator",
        )

        self.vel_gradient_mix = Var(
            initialize=500,
            units=pyunits.s**-1,
            bounds=(300, 1000),
            doc="Velocity gradient of rapid mixer",
        )

        self.vel_gradient_floc = Var(
            initialize=30,
            units=pyunits.s**-1,
            bounds=(20, 80),
            doc="Velocity gradient of flocculator",
        )

        self.frac_mass_water_recovery = Var(
            initialize=0.99,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Fractional recovery of water on mass basis",
        )

        self.CaO_dosing = Var(
            initialize=1e5,
            bounds=(0, None),
            units=pyunits.kg / pyunits.day,
            doc="Lime requirements",
        )

        self.Na2CO3_dosing = Var(
            initialize=1e5,
            bounds=(0, None),
            units=pyunits.kg / pyunits.day,
            doc="Soda ash requirements",
        )

        self.CO2_first_basin = Var(
            initialize=1e2,
            bounds=(0, None),
            units=pyunits.kg / pyunits.day,
            doc="CO2 flow required for recarbonation with single basin",
        )

        self.CO2_second_basin = Var(
            initialize=1e2,
            bounds=(0, None),
            units=pyunits.kg / pyunits.day,
            doc="CO2 flow required for recarbonation in second basin only in excess lime scenario",
        )

        self.excess_CaO = Var(
            initialize=0,
            bounds=(0, None),  # typically 30-70 mg/L, MWH
            units=pyunits.kg / pyunits.m**3,
            doc="Excess lime requiremenent",
        )

        self.pH = Var(
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 14),
            doc="Incoming water pH",
        )

        self.CO2_CaCO3 = Var(
            initialize=0.1,
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Carbonic acid in mg/L as CaCO3",
        )

        self.sludge_prod = Var(
            initialize=1e6,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Sludge production rate in kg/day",
        )

        self.total_hardness = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Total hardness",
        )

        self.carbonate_hardness = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Total carbonate hardness",
        )

        self.noncarbonate_hardness = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Total noncarbonate hardness",
        )

        # Calculate the CO2 in the inlet

        @self.Constraint(doc="Equation to calculate the inlet CO2 in equivalent CaCO3")
        def eq_CO2_CaCO3(b):
            dimensionless_temp = pyunits.convert(
                b.properties_in[0].temperature * pyunits.degK**-1,
                to_units=pyunits.dimensionless,
            )
            K1 = 10 ** (
                14.8435 - 3404.71 / dimensionless_temp - 0.032786 * dimensionless_temp
            )
            K2 = 10 ** (
                6.498 - 2909.39 / dimensionless_temp - 0.02379 * dimensionless_temp
            )

            alpha = (
                1 / ((10 ** (-b.pH)) / K1 + 1 + K2 / (10 ** (-b.pH)))
            ) * pyunits.dimensionless
            CT = (
                b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"] / alpha
            )  # mg CaCO3 /L
            CO2 = (
                CT - b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
            )  # mg CaCO3 /L

            return b.CO2_CaCO3 == CO2

        # Add Ca,Mg carbonate hardness

        @self.Expression(doc="Calcium in influent converted to equivalent CaCO3")
        def Ca_CaCO3(b):
            return pyunits.convert(
                b.properties_in[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                * b.Ca_CaCO3_conv,
                to_units=pyunits.kg / pyunits.m**3,
            )

        @self.Expression(doc="Magnesium in influent converted to equivalent CaCO3")
        def Mg_CaCO3(b):
            return pyunits.convert(
                b.properties_in[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                * b.Mg_CaCO3_conv,
                to_units=pyunits.kg / pyunits.m**3,
            )

        @self.Constraint(doc="Total Hardness in equivalent CaCO3")
        def eq_total_hardness(b):
            return b.total_hardness == pyunits.convert(
                b.Ca_CaCO3 + b.Mg_CaCO3, to_units=pyunits.kg / pyunits.m**3
            )

        @self.Constraint(doc="Carbonate hardness in CaCO3")
        def eq_carbonate_hardness(b):
            return b.carbonate_hardness == pyunits.convert(
                smooth_min(
                    b.total_hardness,
                    b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"],
                    b.eps,
                ),
                to_units=pyunits.kg / pyunits.m**3,
            )

        @self.Constraint(doc="Noncarbonate hardness in equivalent CaCO3")
        def eq_noncarbonate_hardness(b):
            return b.noncarbonate_hardness == pyunits.convert(
                b.total_hardness - b.carbonate_hardness,
                to_units=pyunits.kg / pyunits.m**3,
            )

        @self.Expression(doc="Calculate Calcium carbonate hardness in equivalent CaCO3")
        def Ca_carbonate_hardness_CaCO3(b):
            return smooth_min(b.Ca_CaCO3, b.carbonate_hardness, b.eps)

        @self.Expression(doc="Calculate Calcium non carbonate hardness")
        def Ca_noncarbonate_hardness_CaCO3(b):
            return smooth_max(
                b.Ca_CaCO3 - b.carbonate_hardness,
                0 * pyunits.kg / pyunits.m**3,
                b.eps,
            )

        @self.Expression(
            doc="Calculate Magnesium carbonate hardness in equivalent CaCO3"
        )
        def Mg_carbonate_hardness_CaCO3(b):
            return smooth_max(
                b.carbonate_hardness - b.Ca_CaCO3,
                0 * pyunits.kg / pyunits.m**3,
                b.eps,
            )

        @self.Expression(doc="Calculate Magnesium non carbonate hardness")
        def Mg_noncarbonate_hardness_CaCO3(b):
            return smooth_min(
                b.Mg_CaCO3,
                b.total_hardness - b.carbonate_hardness,
                b.eps,
            )

        @self.Expression(doc="Mg removal efficiency")
        def mg_removal_eff(b):
            mass_flow_mg_eff = pyunits.convert(
                b.mg_eff_target * b.properties_out[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            )
            return (
                1
                - mass_flow_mg_eff
                / b.properties_in[0].flow_mass_phase_comp["Liq", "Mg_2+"]
            )

        @self.Expression(doc="Ca removal efficiency")
        def ca_removal_eff(b):
            mass_flow_ca_eff = pyunits.convert(
                b.ca_eff_target * b.properties_out[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            )
            return (
                1
                - mass_flow_ca_eff
                / b.properties_in[0].flow_mass_phase_comp["Liq", "Ca_2+"]
            )

        @self.Expression(
            doc="Lime concentration for inlet stream to achieve dosing target."
        )
        def lime_concentration(b):
            return pyunits.convert(
                b.CaO_dosing / b.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.gram / pyunits.liter,
            )

        @self.Expression(
            doc="Soda ash concentration for inlet stream to achieve dosing target."
        )
        def soda_concentration(b):
            return pyunits.convert(
                b.Na2CO3_dosing / b.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.gram / pyunits.liter,
            )

        @self.Expression(
            doc="CO2 concentration for inlet stream in first basin to achieve dosing target."
        )
        def co2_first_basin_concentration(b):
            return pyunits.convert(
                b.CO2_first_basin / b.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.gram / pyunits.liter,
            )

        @self.Expression(
            doc="CO2 concentration for inlet stream in second basin to achieve dosing target."
        )
        def co2_second_basin_concentration(b):
            return pyunits.convert(
                b.CO2_second_basin / b.properties_in[0].flow_vol_phase["Liq"],
                to_units=pyunits.gram / pyunits.liter,
            )

        @self.Constraint(doc="Isothermal outlet")
        def eq_isothermal_outlet(b):
            return b.properties_in[0].temperature == b.properties_out[0].temperature

        @self.Constraint(doc="Isothermal waste")
        def eq_isothermal_waste(b):
            return b.properties_in[0].temperature == b.properties_waste[0].temperature

        @self.Constraint(doc="Isobaric outlet")
        def eq_isobaric_outlet(b):
            return b.properties_in[0].pressure == b.properties_out[0].pressure

        @self.Constraint(doc="Isobaric waste")
        def eq_isobaric_waste(b):
            return b.properties_in[0].pressure == b.properties_waste[0].pressure

        # Calculating chemical dosing

        if (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime
        ):
            # These variables aren't otherwise constrained for this SofteningProcedureType
            self.excess_CaO.fix(0)
            self.CO2_second_basin.fix(0)
            self.Na2CO3_dosing.fix(0)

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                    )
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                co2_required_expr = (
                    pyunits.convert(
                        (
                            b.properties_in[0].conc_mass_phase_comp[
                                "Liq", "Alkalinity_2-"
                            ]
                            - b.Ca_CaCO3
                            + b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                            * b.Ca_CaCO3_conv
                        )
                        * b.properties_in[0].flow_vol_phase["Liq"],
                        to_units=pyunits.kg / pyunits.d,
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                )
                return b.CO2_first_basin == Expr_if(
                    co2_required_expr <= 0 * pyunits.kg / pyunits.day,
                    0 * pyunits.kg / pyunits.day,
                    co2_required_expr,
                )

        elif self.config.softening_procedure_type is SofteningProcedureType.excess_lime:

            # These variables aren't otherwise constrained for this SofteningProcedureType
            self.CO2_second_basin.fix(0)
            self.Na2CO3_dosing.fix(0)

            @self.Constraint(doc="Excess lime addition in CaCO3 basis")
            def eq_excess_CaO(b):
                return (
                    b.excess_CaO
                    == (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                    )
                    * b.excess_CaO_coeff
                )

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                        + b.excess_CaO
                    )
                    * b.properties_in[0].flow_vol_phase["Liq"]
                    * b.CaO_mw
                    / b.CaCO3_mw,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                co2_required_expr = pyunits.convert(
                    (
                        b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
                        - b.total_hardness
                        + b.excess_CaO
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                        * b.Ca_CaCO3_conv
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                        * b.Mg_CaCO3_conv
                    )
                    * b.properties_in[0].flow_vol_phase["Liq"]
                    * b.CO2_mw
                    / b.CaCO3_mw,
                    to_units=pyunits.kg / pyunits.d,
                )
                return b.CO2_first_basin == Expr_if(
                    co2_required_expr <= 0 * pyunits.kg / pyunits.day,
                    0 * pyunits.kg / pyunits.day,
                    co2_required_expr,
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime_soda
        ):

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                    )
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="Soda dosing")
            def eq_Na2CO3_dosing(b):
                return b.Na2CO3_dosing == pyunits.convert(
                    b.noncarbonate_hardness
                    * b.Na2CO3_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                co2_required_expr = pyunits.convert(
                    (
                        b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
                        + b.noncarbonate_hardness
                        - b.Ca_CaCO3
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                        * b.Ca_CaCO3_conv
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )
                return b.CO2_first_basin == Expr_if(
                    co2_required_expr <= 0 * pyunits.kg / pyunits.day,
                    0 * pyunits.kg / pyunits.day,
                    co2_required_expr,
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.excess_lime_soda
        ):

            @self.Constraint(doc="Excess lime addition in CaCO3 basis")
            def eq_excess_CaO(b):
                return (
                    b.excess_CaO
                    == (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                    )
                    * b.excess_CaO_coeff
                )

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (
                        b.CO2_CaCO3
                        + b.properties_in[0].conc_mass_phase_comp[
                            "Liq", "Alkalinity_2-"
                        ]
                        + b.Mg_CaCO3
                        + b.excess_CaO
                    )
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="Soda dosing")
            def eq_Na2CO3_dosing(b):
                return b.Na2CO3_dosing == pyunits.convert(
                    b.noncarbonate_hardness
                    * b.Na2CO3_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                co2_required_expr = pyunits.convert(
                    (
                        b.excess_CaO
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                        * b.Mg_CaCO3_conv
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )
                return b.CO2_first_basin == Expr_if(
                    co2_required_expr <= 0 * pyunits.kg / pyunits.day,
                    0 * pyunits.kg / pyunits.day,
                    co2_required_expr,
                )

            @self.Constraint(doc="CO2 for second basin")
            def eq_CO2_second_basin(b):
                return b.CO2_second_basin == pyunits.convert(
                    (
                        b.properties_in[0].conc_mass_phase_comp["Liq", "Alkalinity_2-"]
                        + b.noncarbonate_hardness
                        - b.total_hardness
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                        * b.Ca_CaCO3_conv
                        + b.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                        * b.Mg_CaCO3_conv
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.d,
                )

        if self.config.silica_removal and ["SiO2"] in comps:

            @self.Expression(doc="MgCl2 dosing constraint returns Mg2+ concentration")
            def MgCl2_dosing(b):
                return Expr_if(
                    (
                        (
                            b.properties_in[0].conc_mass_phase_comp["Liq", "SiO2"]
                            * b.SiO2_check_conv
                            > b.properties_in[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                        )
                    ),
                    b.MgCl2_SiO2_ratio
                    * b.properties_in[0].conc_mass_phase_comp["Liq", "SiO2"],
                    1e-15 * pyunits.kg / pyunits.m**3,
                )

        else:
            self.MgCl2_dosing = Var(
                initialize=0,
                units=pyunits.kg / pyunits.m**3,
                bounds=(0, None),
                doc="Magnesium Chloride requirements dosing in Mg2+ concentration",
            )

            if not self.config.silica_removal:
                # No MgCl2 is used if not targeting silica removal
                self.MgCl2_dosing.fix(0)

        @self.Constraint(doc="Water recovery")
        def eq_water_recovery(b):
            return (
                b.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                * b.frac_mass_water_recovery
            )

        @self.Constraint(doc="Water mass balance")
        def eq_water_mass_balance(b):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
                + b.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(non_hardness_comps, doc="Non-hardness component mass balance")
        def eq_non_hardness_comp_mass_balance(b, j):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", j]
                == b.properties_out[0].flow_mass_phase_comp["Liq", j]
                + b.properties_waste[0].flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(non_hardness_comps, doc="Non-hardness component Removal")
        def eq_non_hardness_comp_removal(b, j):
            return b.properties_out[0].flow_mass_phase_comp[
                "Liq", j
            ] == b.properties_in[0].flow_mass_phase_comp["Liq", j] * (
                1 - b.removal_efficiency[j]
            )

        @self.Constraint(doc="Ca in effluent")
        def eq_effluent_ca(b):
            return b.properties_out[0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ] == pyunits.convert(
                b.ca_eff_target * b.properties_out[0].flow_vol_phase["Liq"],
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(doc="Ca mass balance")
        def eq_mass_balance_ca(b):
            return (
                b.properties_waste[0].flow_mass_phase_comp["Liq", "Ca_2+"]
                == b.properties_in[0].flow_mass_phase_comp["Liq", "Ca_2+"]
                + pyunits.convert(
                    b.excess_CaO
                    * (b.Ca_mw / b.CaCO3_mw)
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )
                - b.properties_out[0].flow_mass_phase_comp["Liq", "Ca_2+"]
            )

        @self.Constraint(doc="Mg in effluent")
        def eq_effluent_mg(b):
            return b.properties_out[0].flow_mass_phase_comp[
                "Liq", "Mg_2+"
            ] == pyunits.convert(
                (b.mg_eff_target * b.properties_out[0].flow_vol_phase["Liq"]),
                to_units=pyunits.kg / pyunits.s,
            )

        @self.Constraint(doc="Mg mass balance")
        def eq_mass_balance_mg(b):
            return b.properties_waste[0].flow_mass_phase_comp[
                "Liq", "Mg_2+"
            ] == b.properties_in[0].flow_mass_phase_comp[
                "Liq", "Mg_2+"
            ] - b.properties_out[
                0
            ].flow_mass_phase_comp[
                "Liq", "Mg_2+"
            ] + pyunits.convert(
                (b.MgCl2_dosing * b.properties_out[0].flow_vol_phase["Liq"]),
                to_units=pyunits.kg / pyunits.s,
            )

        # Alkalinity removal for each softening type is:
        # single stage lime (example in book Crittenden): source water alkalinity - Ca hardness + residual Ca hardness
        # excess lime (example in book Crittenden): source water alkalinity - Ca hardness - excess lime dose + residual Ca hardness
        # single stage lime soda (only hydroxide alkalinity): source water alkalinity + soda ash - total hardness + residual hardness
        # excess lime soda (only hydroxide alkalinity): source water alkalinity + soda ash - total hardness + residual hardness

        if (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime
        ):

            @self.Constraint(doc="Alkalinity mass balance")
            def eq_effluent_alk(b):
                return b.properties_out[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] == b.properties_in[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] - pyunits.convert(
                    b.Ca_CaCO3 * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) + pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                    * b.Ca_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        elif self.config.softening_procedure_type is SofteningProcedureType.excess_lime:

            @self.Constraint(doc="Alkalinity mass balance")
            def eq_effluent_alk(b):
                return b.properties_out[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] == b.properties_in[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] - pyunits.convert(
                    b.total_hardness * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) - pyunits.convert(
                    b.excess_CaO * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) + pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                    * b.Ca_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) + pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                    * b.Ca_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) + pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                    * b.Mg_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime_soda
        ):

            @self.Constraint(doc="Alkalinity mass balance")
            def eq_effluent_alk(b):
                return b.properties_out[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] == pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                    * b.Ca_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.excess_lime_soda
        ):

            @self.Constraint(doc="Alkalinity mass balance")
            def eq_effluent_alk(b):
                return b.properties_out[0].flow_mass_phase_comp[
                    "Liq", "Alkalinity_2-"
                ] == pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Ca_2+"]
                    * b.Ca_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                ) + pyunits.convert(
                    b.properties_out[0].conc_mass_phase_comp["Liq", "Mg_2+"]
                    * b.Mg_CaCO3_conv
                    * b.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        @self.Constraint(doc="Alkalinity mass balance")
        def eq_mass_balance_alk(b):
            return (
                b.properties_waste[0].flow_mass_phase_comp["Liq", "Alkalinity_2-"]
                == b.properties_in[0].flow_mass_phase_comp["Liq", "Alkalinity_2-"]
                - b.properties_out[0].flow_mass_phase_comp["Liq", "Alkalinity_2-"]
            )

        if ["TSS"] in comps:

            @self.Constraint(doc="Sludge production")
            def eq_sludge_prod(b):
                return b.sludge_prod == pyunits.convert(
                    (
                        b.Ca_hardness_CaCO3_sludge_factor
                        * b.Ca_carbonate_hardness_CaCO3
                        + b.Mg_carbonate_hardness_sludge_factor
                        * b.Mg_carbonate_hardness_CaCO3
                        + b.Ca_noncarbonate_hardness_CaCO3
                        + b.Mg_noncarbonate_hardness_sludge_factor
                        * b.Mg_noncarbonate_hardness_CaCO3
                        + b.excess_CaO
                        + b.properties_in[0].conc_mass_phase_comp["Liq", "TSS"]
                        + b.MgCl2_dosing
                        * b.MgOH2_mw
                        / b.Mg_mw  # to convert to Mg(OH)2 solid
                    )
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        else:

            @self.Constraint(doc="Sludge production")
            def eq_sludge_prod(b):
                return b.sludge_prod == pyunits.convert(
                    (
                        b.Ca_hardness_CaCO3_sludge_factor
                        * b.Ca_carbonate_hardness_CaCO3
                        + b.Mg_carbonate_hardness_sludge_factor
                        * b.Mg_carbonate_hardness_CaCO3
                        + b.Ca_noncarbonate_hardness_CaCO3
                        + b.Mg_noncarbonate_hardness_sludge_factor
                        * b.Mg_noncarbonate_hardness_CaCO3
                        + b.excess_CaO
                        + b.MgCl2_dosing
                        * b.MgOH2_mw
                        / b.Mg_mw  # to convert to Mg(OH)2 solid
                    )
                    * b.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.kg / pyunits.s,
                )

        @self.Constraint(doc="Volume of mixer")
        def eq_volume_mixer(b):
            return b.volume_mixer == pyunits.convert(
                b.properties_in[0].flow_vol_phase["Liq"] * b.retention_time_mixer,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Volume of flocculator")
        def eq_volume_floc(b):
            return b.volume_floc == pyunits.convert(
                b.properties_in[0].flow_vol_phase["Liq"] * b.retention_time_floc,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Volume of sedimentation basin")
        def eq_volume_sed(b):
            return b.volume_sed == pyunits.convert(
                b.properties_in[0].flow_vol_phase["Liq"] * b.retention_time_sed,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Volume of recarbonation basin")
        def eq_volume_recarb(b):
            return b.volume_recarb == pyunits.convert(
                b.properties_in[0].flow_vol_phase["Liq"] * b.retention_time_recarb,
                to_units=pyunits.m**3,
            )

    def initialize_build(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = blk.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            blk.state_args = state_args = {}
            state_dict = blk.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args_out = deepcopy(state_args)

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_waste = deepcopy(state_args)

        blk.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        blk.state_args_out = state_args_out
        blk.state_args_waste = state_args_waste

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.ca_eff_target) is None:
            sf = iscale.get_scaling_factor(
                self.properties_in[0].conc_mass_phase_comp["Liq", "Ca_2+"]
            )
            iscale.set_scaling_factor(self.ca_eff_target, sf)

        if iscale.get_scaling_factor(self.mg_eff_target) is None:
            sf = iscale.get_scaling_factor(
                self.properties_in[0].conc_mass_phase_comp["Liq", "Mg_2+"]
            )
            iscale.set_scaling_factor(self.mg_eff_target, 1)

        if iscale.get_scaling_factor(self.removal_efficiency) is None:
            iscale.set_scaling_factor(self.removal_efficiency, 1)

        if iscale.get_scaling_factor(self.retention_time_mixer) is None:
            iscale.set_scaling_factor(self.retention_time_mixer, 1)

        if iscale.get_scaling_factor(self.retention_time_floc) is None:
            iscale.set_scaling_factor(self.retention_time_floc, 0.1)

        if iscale.get_scaling_factor(self.retention_time_sed) is None:
            iscale.set_scaling_factor(self.retention_time_sed, 1e-2)

        if iscale.get_scaling_factor(self.retention_time_recarb) is None:
            iscale.set_scaling_factor(self.retention_time_recarb, 0.1)

        if iscale.get_scaling_factor(self.sedimentation_overflow) is None:
            iscale.set_scaling_factor(self.sedimentation_overflow, 0.1)

        if iscale.get_scaling_factor(self.number_mixers) is None:
            iscale.set_scaling_factor(self.number_mixers, 1)

        if iscale.get_scaling_factor(self.number_floc) is None:
            iscale.set_scaling_factor(self.number_floc, 1)

        if iscale.get_scaling_factor(self.volume_mixer) is None:
            iscale.set_scaling_factor(self.volume_mixer, 1e-1)

        if iscale.get_scaling_factor(self.volume_floc) is None:
            iscale.set_scaling_factor(self.volume_floc, 1e-2)

        if iscale.get_scaling_factor(self.volume_sed) is None:
            iscale.set_scaling_factor(self.volume_sed, 1e-2)

        if iscale.get_scaling_factor(self.volume_recarb) is None:
            iscale.set_scaling_factor(self.volume_recarb, 1e-2)

        if iscale.get_scaling_factor(self.vel_gradient_mix) is None:
            iscale.set_scaling_factor(self.vel_gradient_mix, 1e-2)

        if iscale.get_scaling_factor(self.vel_gradient_floc) is None:
            iscale.set_scaling_factor(self.vel_gradient_floc, 1e-2)

        if iscale.get_scaling_factor(self.frac_mass_water_recovery) is None:
            iscale.set_scaling_factor(self.frac_mass_water_recovery, 1)

        if iscale.get_scaling_factor(self.CaO_dosing) is None:
            iscale.set_scaling_factor(self.CaO_dosing, 1e-4)

        if iscale.get_scaling_factor(self.excess_CaO) is None:
            iscale.set_scaling_factor(self.excess_CaO, 10)

        if iscale.get_scaling_factor(self.CO2_first_basin) is None:
            iscale.set_scaling_factor(self.CO2_first_basin, 1e-2)

        if iscale.get_scaling_factor(self.CO2_second_basin) is None:
            iscale.set_scaling_factor(self.CO2_second_basin, 1e-2)

        if iscale.get_scaling_factor(self.Na2CO3_dosing) is None:
            iscale.set_scaling_factor(self.Na2CO3_dosing, 1e-5)

        if iscale.get_scaling_factor(self.CO2_CaCO3) is None:
            iscale.set_scaling_factor(self.CO2_CaCO3, 1)

        if iscale.get_scaling_factor(self.sludge_prod) is None:
            iscale.set_scaling_factor(self.sludge_prod, 1)

        if iscale.get_scaling_factor(self.total_hardness) is None:
            iscale.set_scaling_factor(self.total_hardness, 1e-1)

        if iscale.get_scaling_factor(self.carbonate_hardness) is None:
            iscale.set_scaling_factor(self.carbonate_hardness, 1e-1)

        if iscale.get_scaling_factor(self.noncarbonate_hardness) is None:
            iscale.set_scaling_factor(self.noncarbonate_hardness, 1e-1)

        if iscale.get_scaling_factor(self.pH) is None:
            iscale.set_scaling_factor(self.pH, 1)

        if isinstance(self.MgCl2_dosing, Var):
            if iscale.get_scaling_factor(self.MgCl2_dosing) is None:
                iscale.set_scaling_factor(self.MgCl2_dosing, 10)

        if (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime_soda
        ):
            iscale.constraint_scaling_transform(self.Na2CO3_dosing, 1e-3)

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.excess_lime_soda
        ):
            iscale.constraint_scaling_transform(self.eq_Na2CO3_dosing, 1e-4)
            iscale.constraint_scaling_transform(self.eq_CO2_second_basin, 1e-2)

        iscale.constraint_scaling_transform(self.eq_CaO_dosing, 1e-3)

        iscale.constraint_scaling_transform(self.eq_CO2_first_basin, 1e-3)

        iscale.constraint_scaling_transform(self.eq_effluent_ca, 1)

        iscale.constraint_scaling_transform(self.eq_effluent_mg, 10)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Waste Outlet": self.waste,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        var_dict["Lime dose"] = self.CaO_dosing
        var_dict["Soda dose"] = self.Na2CO3_dosing
        var_dict["CO2 in first basin"] = self.CO2_first_basin
        var_dict["CO2 in second basin"] = self.CO2_second_basin
        if (
            self.config.silica_removal
            and ["SiO2"] in self.config.property_package.solute_set
        ):
            expr_dict["MgCl2 dose"] = self.MgCl2_dosing
        else:
            var_dict["MgCl2 dose"] = self.MgCl2_dosing
        expr_dict["Sludge produced"] = self.sludge_prod

        return {"vars": var_dict, "exprs": expr_dict}

    @property
    def default_costing_method(self):
        return cost_chemical_softening
