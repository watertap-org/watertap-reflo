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

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum

from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin

__author__ = "Mukta Hardikar, Abdiel Lugo"

_log = idaeslog.getLogger(__name__)


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


@declare_process_block_class("ChemicalSofteningZO")
class ChemicalSofteningZOData(InitializationMixin, UnitModelBlockData):
    """
    Zero order chemical softening model
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
        1.	Crittenden, J. C., & Montgomery Watson Harza (Firm). (2012). Water treatment principles and design. Hoboken, N.J: J.Wiley.
        2.	Davis, M. L. (2010). Water and wastewater engineering: Design principles and practice.
        3.	Baruth. (2005). Water treatment plant design / American Water Works Association, American Society of Civil Engineers; Edward E. Baruth, technical editor. (Fourth edition.). McGraw-Hill.
        4.	Edzwald, J. K., & American Water Works Association. (2011). Water quality & treatment: A handbook on drinking water. New York: McGraw-Hill.
        5.	R.O. Mines Environmental Engineering: Principles and Practice, 1st Ed, John Wiley & Sons
        6.	Lee, C. C., & Lin, S. D. (2007). Handbook of environmental engineering calculations. New York: McGraw Hill.
        """

        # This creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        # Add ports
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_waste = self.properties_waste[0]
        comps = self.config.property_package.solute_set
        non_hardness_comps = [
            j
            for j in self.config.property_package.solute_set
            if j not in ["Ca_2+", "Mg_2+"]
        ]

        # MW of components

        self.CaO_mw = Param(
            initialize=56,
            units=pyunits.g / pyunits.mol,
            doc="Molecular weight of CaO (dry lime)",
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
            initialize=44, units=pyunits.g / pyunits.mol, doc="Molecular weight of soda"
        )

        self.MgCl2_SiO2_ratio = Param(initialize=9.2, doc="Ratio of MgCl2 to SiO2")

        self.Ca_hardness_CaCO3_sludge_factor = Param(
            initialize=2, doc="Sludge produced per kg Ca in CaCO3 hardness"
        )

        self.Mg_hardness_CaCO3_sludge_factor = Param(
            initialize=2.6, doc="Sludge produced per kg Mg in CaCO3 hardness"
        )

        self.Mg_hardness_nonCaCO3_sludge_prod_factor = Param(
            initialize=1.6, doc="Sludge produced per kg Mg in non-CaCO3 hardness"
        )

        self.excess_CaO_coeff = Param(
            initialize=0.15, doc="Multiplication factor to calculate excess CaO"
        )

        self.ca_eff_target = Var(
            initialize=0.020,
            units=pyunits.kg / pyunits.m**3,
            doc="Target Ca in CaCO3 equivalents in the effluent",
        )

        self.mg_eff_target = Var(
            initialize=0.010,
            units=pyunits.kg / pyunits.m**3,
            doc="Target Mg in CaCO3 equivalents in the effluent",
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
            doc="Removal efficiency of all solutes other than Ca and Mg",
        )

        # System design variables

        self.retention_time_mixer = Var(
            initialize=1,
            units=pyunits.minutes,
            bounds=(0.1, 5),
            doc="Retention time for mixer",
        )

        self.retention_time_floc = Var(
            initialize=12,
            units=pyunits.minutes,
            bounds=(10, 45),
            doc="Retention time for flocculator",
        )

        self.retention_time_sed = Var(
            initialize=150,
            units=pyunits.minutes,
            bounds=(120, 240),
            doc="Retention time for sedimentation basin",
        )

        self.retention_time_recarb = Var(
            initialize=18,
            units=pyunits.minutes,
            bounds=(15, 30),
            doc="Retention time for recarbonation",
        )

        self.sedimentation_overflow = Var(
            initialize=90,
            units=pyunits.m / pyunits.day,
            bounds=(30, 100),
            doc="Sedimentation basin overflow rate",
        )

        self.no_of_mixer = Var(
            initialize=1, units=pyunits.dimensionless, doc="Volume of mixer"
        )

        self.no_of_floc = Var(
            initialize=1, units=pyunits.dimensionless, doc="Volume of mixer"
        )

        self.volume_mixer = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of mixer"
        )

        self.volume_floc = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of flocculator"
        )

        self.volume_sed = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of sedimentation basin"
        )

        self.volume_recarb = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of recarbonator"
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

        self.Ca_CaCO3_conv = Param(
            initialize=2.5,
            units=pyunits.dimensionless,
            doc="Conversion factor for Ca to equivalent CaCO3",
        )

        self.Mg_CaCO3_conv = Param(
            initialize=4.12,
            units=pyunits.dimensionless,
            doc="Conversion factor for Ca to equivalent CaCO3",
        )

        self.frac_vol_recovery = Var(
            initialize=0.99,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Fractional volumetric recovery of water",
        )

        self.CaO_dosing = Var(
            initialize=1e5,
            units=pyunits.kg / pyunits.day,
            bounds=(0, None),
            doc="Lime requirements",
        )

        self.Na2CO3_dosing = Var(
            initialize=1e5,
            units=pyunits.kg / pyunits.day,
            bounds=(0, None),
            doc="Soda ash requirements",
        )

        self.CO2_first_basin = Var(
            initialize=1e2,
            units=pyunits.kg / pyunits.day,
            bounds=(0, None),
            doc="CO2 concentration required for recarbonation with single basin",
        )

        self.CO2_second_basin = Var(
            initialize=1e2,
            units=pyunits.kg / pyunits.day,
            bounds=(0, None),
            doc="CO2 concentration required for recarbonation in second basin only in excess lime scenario",
        )

        self.excess_CaO = Var(
            initialize=0,
            units=pyunits.kg / pyunits.m**3,
            doc="Excess lime requiremenent",
        )

        self.CO2_CaCO3 = Var(
            initialize=115,
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

        # Add Ca,Mg carbonate hardness --> Created as temporary variables so that basic property package can be used
        if "Ca_2+" in comps:

            @self.Expression(doc="Calcium in influent converted to equivalent CaCO3")
            def Ca_CaCO3(b):
                return pyunits.convert(
                    prop_in.conc_mass_comp["Ca_2+"] * b.Ca_CaCO3_conv,
                    to_units=pyunits.kg / pyunits.m**3,
                )

        else:
            self.Ca_CaCO3 = Param(
                initialize=0,
                units=pyunits.kg / pyunits.m**3,
                doc="Calcium in influent converted to kg/m3 as CaCO3",
            )

        if "Mg_2+" in comps:

            @self.Expression(doc="Magnesium in influent converted to kg/m3 as CaCO3")
            def Mg_CaCO3(b):
                return pyunits.convert(
                    prop_in.conc_mass_comp["Mg_2+"] * b.Mg_CaCO3_conv,
                    to_units=pyunits.kg / pyunits.m**3,
                )

        else:
            self.Mg_CaCO3 = Param(
                initialize=0,
                units=pyunits.kg / pyunits.m**3,
                doc="Magnesium in influent converted to kg/m3 as CaCO3",
            )

        if "Alkalinity_2-" in comps:

            @self.Expression(doc="Calculate Calcium carbonate hardness")
            def Ca_hardness_CaCO3(b):
                return Expr_if(
                    b.Ca_CaCO3 >= prop_in.conc_mass_comp["Alkalinity_2-"],
                    prop_in.conc_mass_comp["Alkalinity_2-"],
                    b.Ca_CaCO3,
                )

            @self.Expression(doc="Calculate Calcium non carbonate hardness")
            def Ca_hardness_nonCaCO3(b):
                return Expr_if(
                    b.Ca_CaCO3 >= prop_in.conc_mass_comp["Alkalinity_2-"],
                    b.Ca_CaCO3 - prop_in.conc_mass_comp["Alkalinity_2-"],
                    1e-15 * pyunits.kg / pyunits.m**3,
                )

            @self.Expression(doc="Calculate Magnesium carbonate hardness")
            def Mg_hardness_CaCO3(b):
                return Expr_if(
                    b.Ca_CaCO3 >= prop_in.conc_mass_comp["Alkalinity_2-"],
                    1e-15 * pyunits.kg / pyunits.m**3,
                    prop_in.conc_mass_comp["Alkalinity_2-"] - b.Ca_CaCO3,
                )

            @self.Expression(doc="Calculate Magnesium non carbonate hardness")
            def Mg_hardness_nonCaCO3(b):
                return Expr_if(
                    b.Ca_CaCO3 >= prop_in.conc_mass_comp["Alkalinity_2-"],
                    b.Mg_CaCO3,
                    b.Mg_CaCO3 + b.Ca_CaCO3 - prop_in.conc_mass_comp["Alkalinity_2-"],
                )

        # Calculating chemical dosing

        if (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime
        ):

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (b.CO2_CaCO3 + b.Ca_hardness_CaCO3)
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                return (
                    b.CO2_first_basin
                    == pyunits.convert(
                        (
                            prop_in.conc_mass_comp["Alkalinity_2-"]
                            - b.Ca_CaCO3
                            + prop_out.conc_mass_comp["Ca_2+"] * b.Ca_CaCO3_conv
                        )
                        * prop_in.flow_vol,
                        to_units=pyunits.kg / pyunits.d,
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                )

        elif self.config.softening_procedure_type is SofteningProcedureType.excess_lime:

            @self.Constraint(doc="Excess lime addition")
            def eq_excess_CaO(b):
                return (
                    b.excess_CaO
                    == (
                        b.CO2_CaCO3
                        + prop_in.conc_mass_comp["Alkalinity_2-"]
                        + b.Mg_CaCO3
                    )
                    * b.excess_CaO_coeff
                )

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return (
                    b.CaO_dosing
                    == (
                        b.CO2_CaCO3
                        + prop_in.conc_mass_comp["Alkalinity_2-"]
                        + b.Mg_CaCO3
                        + b.excess_CaO
                    )
                    * prop_in.flow_vol
                    * b.CaO_mw
                    / b.CaCO3_mw
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                return (
                    b.CO2_first_basin
                    == (
                        prop_in.conc_mass_comp["Alkalinity_2-"]
                        - (b.Ca_CaCO3 + b.Ca_CaCO3)
                        + b.excess_CaO
                        + prop_out.conc_mass_comp["Ca_2+"] * b.Ca_CaCO3_conv
                        + prop_out.conc_mass_comp["Mg_2+"] * b.Mg_CaCO3_conv
                    )
                    * prop_in.flow_vol
                    * b.CO2_mw
                    / b.CaCO3_mw
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime_soda
        ):

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (b.CO2_CaCO3 + b.Ca_hardness_CaCO3)
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="Soda dosing")
            def eq_Na2CO3_dosing(b):
                return b.Na2CO3_dosing == pyunits.convert(
                    (b.Ca_hardness_nonCaCO3 + b.Mg_hardness_nonCaCO3)
                    * b.Na2CO3_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                return b.CO2_first_basin == pyunits.convert(
                    (
                        prop_in.conc_mass_comp["Alkalinity_2-"]
                        + (b.Ca_hardness_nonCaCO3 + b.Mg_hardness_nonCaCO3)
                        - b.Ca_CaCO3
                        + prop_out.conc_mass_comp["Ca_2+"] * b.Ca_CaCO3_conv
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.excess_lime_soda
        ):

            @self.Constraint(doc="Excess lime addition")
            def eq_excess_CaO(b):
                return (
                    b.excess_CaO
                    == (
                        b.CO2_CaCO3
                        + prop_in.conc_mass_comp["Alkalinity_2-"]
                        + b.Mg_CaCO3
                    )
                    * b.excess_CaO_coeff
                )

            @self.Constraint(doc="Lime dosing")
            def eq_CaO_dosing(b):
                return b.CaO_dosing == pyunits.convert(
                    (
                        b.CO2_CaCO3
                        + b.Ca_hardness_CaCO3
                        + 2 * b.Mg_hardness_CaCO3
                        + b.Mg_hardness_nonCaCO3
                        + b.excess_CaO
                    )
                    * b.CaO_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="Soda dosing")
            def eq_Na2CO3_dosing(b):
                return b.Na2CO3_dosing == pyunits.convert(
                    (b.Ca_hardness_nonCaCO3 + b.Mg_hardness_nonCaCO3)
                    * b.Na2CO3_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for first basin")
            def eq_CO2_first_basin(b):
                return b.CO2_first_basin == pyunits.convert(
                    (b.excess_CaO + prop_out.conc_mass_comp["Mg_2+"] * b.Mg_CaCO3_conv)
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

            @self.Constraint(doc="CO2 for second basin")
            def eq_CO2_second_basin(b):
                return b.CO2_second_basin == pyunits.convert(
                    (
                        prop_in.conc_mass_comp["Alkalinity_2-"]
                        + (b.Ca_hardness_nonCaCO3 + b.Mg_hardness_nonCaCO3)
                        - (b.Ca_CaCO3 + b.Ca_CaCO3)
                        + prop_out.conc_mass_comp["Ca_2+"] * b.Ca_CaCO3_conv
                        + prop_out.conc_mass_comp["Mg_2+"] * b.Mg_CaCO3_conv
                    )
                    * b.CO2_mw
                    / b.CaCO3_mw
                    * prop_in.flow_vol,
                    to_units=pyunits.kg / pyunits.d,
                )

        if self.config.silica_removal and ["SiO2"] in comps:

            @self.Expression(doc="MgCl2 dosing constraint")
            def MgCl2_dosing(b):
                return Expr_if(
                    (
                        (
                            prop_in.conc_mass_comp["SiO2"] * 2.35
                            > prop_in.conc_mass_comp["Mg_2+"]
                        )
                    ),
                    b.MgCl2_SiO2_ratio * prop_in.conc_mass_comp["SiO2"],
                    1e-15 * pyunits.kg / pyunits.m**3,
                )

        else:
            self.MgCl2_dosing = Var(
                initialize=0,
                units=pyunits.kg / pyunits.m**3,
                bounds=(0, None),
                doc="Magnesium Chloride requirments dosing",
            )

        # value(prop_in.conc_mass_comp["SiO2"]>0.05) and value

        @self.Constraint(doc="Recovery")
        def eq_recovery(b):
            return prop_in.flow_vol == prop_out.flow_vol + prop_waste.flow_vol

        @self.Constraint(doc="Waste flow")
        def eq_waste_flow(b):
            return prop_out.flow_vol == prop_in.flow_vol * b.frac_vol_recovery

        @self.Constraint(non_hardness_comps, doc="Mass balance")
        def eq_mass_balance(b, j):
            return (
                prop_in.flow_mass_comp[j]
                == prop_out.flow_mass_comp[j] + prop_waste.flow_mass_comp[j]
            )

        @self.Constraint(non_hardness_comps, doc="Component Removal")
        def eq_component_removal(b, j):
            return prop_waste.flow_mass_comp[j] == prop_in.flow_mass_comp[j] * (
                b.removal_efficiency[j]
            )

        @self.Constraint(doc="Volume of mixer")
        def eq_volume_mixer(b):
            return (
                b.volume_mixer
                == pyunits.convert(
                    b.properties_in[0].flow_vol * b.retention_time_mixer,
                    to_units=pyunits.m**3,
                )
                * b.no_of_mixer
            )

        @self.Constraint(doc="Volume of flocculator")
        def eq_volume_floc(b):
            return (
                b.volume_floc
                == pyunits.convert(
                    b.properties_in[0].flow_vol * b.retention_time_floc,
                    to_units=pyunits.m**3,
                )
                * b.no_of_floc
            )

        @self.Constraint(doc="Volume of sedimentation basin")
        def eq_volume_sed(b):
            return b.volume_sed == pyunits.convert(
                b.properties_in[0].flow_vol * b.retention_time_sed,
                to_units=pyunits.m**3,
            )

        @self.Constraint(doc="Volume of recarbonation basin")
        def eq_volume_recarb(b):
            return b.volume_recarb == pyunits.convert(
                b.properties_in[0].flow_vol * b.retention_time_recarb,
                to_units=pyunits.m**3,
            )

        if ["Ca_2+"] in comps:

            @self.Constraint(doc="Ca in effluent")
            def eq_effluent_ca(b):
                return (
                    prop_out.conc_mass_comp["Ca_2+"]
                    == b.ca_eff_target / b.Ca_CaCO3_conv
                )

            @self.Constraint(doc="Ca mass balance")
            def eq_mass_balance_ca(b):
                return prop_waste.flow_mass_comp["Ca_2+"] == prop_in.flow_mass_comp[
                    "Ca_2+"
                ] - prop_out.flow_mass_comp["Ca_2+"] + pyunits.convert(
                    b.excess_CaO * prop_in.flow_vol, to_units=pyunits.kg / pyunits.s
                )

        if ["Mg_2+"] in comps:

            @self.Constraint(doc="Mg in effluent")
            def eq_effluent_mg(b):
                return (
                    prop_out.conc_mass_comp["Mg_2+"]
                    == b.mg_eff_target / b.Mg_CaCO3_conv
                )

            @self.Constraint(doc="Mg mass balance")
            def eq_mass_balance_mg(b):
                return prop_waste.flow_mass_comp["Mg_2+"] == prop_in.flow_mass_comp[
                    "Mg_2+"
                ] - prop_out.flow_mass_comp["Mg_2+"] + pyunits.convert(
                    b.MgCl2_dosing * prop_in.flow_vol, to_units=pyunits.kg / pyunits.s
                )

        if ["TSS"] in comps:

            @self.Constraint(doc="Sludge production")
            def eq_sludge_prod(b):
                return (
                    b.sludge_prod
                    == (
                        b.Ca_hardness_CaCO3_sludge_factor * b.Ca_hardness_CaCO3
                        + b.Mg_hardness_CaCO3_sludge_factor * b.Mg_hardness_CaCO3
                        + b.Ca_hardness_nonCaCO3
                        + b.Mg_hardness_nonCaCO3_sludge_prod_factor
                        * b.Mg_hardness_nonCaCO3
                        + b.excess_CaO
                        + prop_in.conc_mass_comp["TSS"]
                        + b.MgCl2_dosing
                    )
                    * prop_in.flow_vol
                )

        else:

            @self.Constraint(doc="Sludge production")
            def eq_sludge_prod(b):
                return (
                    b.sludge_prod
                    == (
                        b.Ca_hardness_CaCO3_sludge_factor * b.Ca_hardness_CaCO3
                        + b.Mg_hardness_CaCO3_sludge_factor * b.Mg_hardness_CaCO3
                        + b.Ca_hardness_nonCaCO3
                        + b.Mg_hardness_nonCaCO3_sludge_prod_factor
                        * b.Mg_hardness_nonCaCO3
                        + b.excess_CaO
                        + b.MgCl2_dosing
                    )
                    * prop_in.flow_vol
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
        iscale.set_scaling_factor(self.ca_eff_target, 1)
        iscale.set_scaling_factor(self.mg_eff_target, 1)
        iscale.set_scaling_factor(self.removal_efficiency, 1)
        iscale.set_scaling_factor(self.retention_time_mixer, 1)
        iscale.set_scaling_factor(self.retention_time_floc, 0.1)
        iscale.set_scaling_factor(self.retention_time_sed, 1e-2)
        iscale.set_scaling_factor(self.retention_time_recarb, 0.1)
        iscale.set_scaling_factor(self.sedimentation_overflow, 0.1)
        iscale.set_scaling_factor(self.no_of_mixer, 1)
        iscale.set_scaling_factor(self.no_of_floc, 1)
        iscale.set_scaling_factor(self.volume_mixer, 1e-1)
        iscale.set_scaling_factor(self.volume_floc, 1e-2)
        iscale.set_scaling_factor(self.volume_sed, 1e-2)
        iscale.set_scaling_factor(self.volume_recarb, 1e-2)
        iscale.set_scaling_factor(self.vel_gradient_mix, 1e-2)
        iscale.set_scaling_factor(self.vel_gradient_floc, 1e-2)
        iscale.set_scaling_factor(self.frac_vol_recovery, 1)
        iscale.set_scaling_factor(self.CaO_dosing, 1e-4)
        iscale.set_scaling_factor(self.excess_CaO, 1)
        iscale.set_scaling_factor(self.CO2_first_basin, 1e-2)
        iscale.set_scaling_factor(self.CO2_second_basin, 1e-2)
        iscale.set_scaling_factor(self.Na2CO3_dosing, 1e-5)
        iscale.set_scaling_factor(self.CO2_CaCO3, 1)
        iscale.set_scaling_factor(self.sludge_prod, 1e-6)

        comps = self.config.property_package.solute_set
        non_hardness_comps = [
            j
            for j in self.config.property_package.solute_set
            if j not in ["Ca_2+", "Mg_2+"]
        ]

        if self.config.softening_procedure_type is SofteningProcedureType.excess_lime:

            sf = iscale.get_scaling_factor(
                self.properties_in[0].conc_mass_comp["Alkalinity_2-"]
            )
            iscale.constraint_scaling_transform(self.eq_excess_CaO, sf)

        if (
            self.config.softening_procedure_type
            is SofteningProcedureType.single_stage_lime_soda
        ):
            sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol) * 3600 * 24
            iscale.constraint_scaling_transform(self.Na2CO3_dosing, sf)

        elif (
            self.config.softening_procedure_type
            is SofteningProcedureType.excess_lime_soda
        ):

            sf = iscale.get_scaling_factor(
                self.properties_in[0].conc_mass_comp["Alkalinity_2-"]
            )
            iscale.constraint_scaling_transform(self.eq_excess_CaO, sf)

            sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol) * 3600 * 24
            iscale.constraint_scaling_transform(self.eq_Na2CO3_dosing, sf)

            sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol) * 3600 * 24
            iscale.constraint_scaling_transform(self.eq_CO2_second_basin, sf)

        # Scaling constraints
        sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol) * 3600 * 24
        iscale.constraint_scaling_transform(self.eq_CaO_dosing, sf)

        sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol) * 3600 * 24
        iscale.constraint_scaling_transform(self.eq_CO2_first_basin, sf)

        if self.config.silica_removal and ["SiO2"] in comps:
            iscale.set_scaling_factor(self.MgCl2_dosing, 1)

        sf = iscale.get_scaling_factor(self.properties_in[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_recovery, sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_waste_flow, sf)

        for c in non_hardness_comps:
            sf = iscale.get_scaling_factor(self.properties_waste[0].flow_mass_comp[c])
            iscale.constraint_scaling_transform(self.eq_mass_balance[c], sf)

        for c in non_hardness_comps:
            sf = iscale.get_scaling_factor(self.properties_waste[0].flow_mass_comp[c])
            iscale.constraint_scaling_transform(self.eq_component_removal[c], sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_volume_mixer, sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_volume_floc, sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_volume_sed, sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_volume_recarb, sf)

        iscale.constraint_scaling_transform(self.eq_effluent_ca, 10)

        iscale.constraint_scaling_transform(self.eq_effluent_mg, 10)

        sf = iscale.get_scaling_factor(self.properties_in[0].conc_mass_comp["Ca_2+"])
        iscale.constraint_scaling_transform(self.eq_mass_balance_ca, sf)

        sf = iscale.get_scaling_factor(self.properties_in[0].conc_mass_comp["Mg_2+"])
        iscale.constraint_scaling_transform(self.eq_mass_balance_mg, sf)

        sf = iscale.get_scaling_factor(self.properties_out[0].flow_vol)
        iscale.constraint_scaling_transform(self.eq_sludge_prod, sf)

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
        var_dict["Lime dose"] = self.CaO_dosing.value
        var_dict["Soda dose"] = self.Na2CO3_dosing.value
        var_dict["CO2 in first basin"] = self.CO2_first_basin.value
        var_dict["CO2 in second basin"] = self.CO2_second_basin.value
        if (
            self.config.silica_removal
            and ["SiO2"] in self.config.property_package.solute_set
        ):
            var_dict["MgCl2 dose"] = self.MgCl2_dosing()
        else:
            var_dict["MgCl2 dose"] = self.MgCl2_dosing.value
        var_dict["Sludge produced"] = self.sludge_prod.value

        return {"vars": var_dict, "exprs": expr_dict}
