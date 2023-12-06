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

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    Constraint,
    check_optimal_termination,
    Param,
    Suffix,
    Block,
    log,
    log10,
    exp,
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from enum import Enum, auto

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum, extract_data
from idaes.core.util.exceptions import InitializationError, ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin

__author__ = "Kurban Sitterley"


_log = idaeslog.getLogger(__name__)


class PackingMaterial(StrEnum):
    PVC = "PVC"
    Ceramic = "ceramic"
    StainlessSteel = "stainless_steel"


@declare_process_block_class("AirStripping0D")
class AirStripping0DData(InitializationMixin, UnitModelBlockData):
    """
    Zero dimensional air stripping model
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
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentPhase,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.none,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.none.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
        **default** - MomentumBalanceType.pressureTotal.
        **Valid values:** {
        **MomentumBalanceType.none** - exclude momentum balances,
        **MomentumBalanceType.pressureTotal** - single pressure balance for material,
        **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
        **MomentumBalanceType.momentumTotal** - single momentum balance for material,
        **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "target",
        ConfigValue(
            default="X",
            domain=str,
            description="Designates targeted species for removal",
        ),
    )

    CONFIG.declare(
        "packing_material",
        ConfigValue(
            default=PackingMaterial.PVC,
            domain=In(PackingMaterial),
            description="Material used for product. ",
        ),
    )

    def build(self):
        super().build()

        target = self.config.target
        self.target_set = Set(initialize=[target])
        comps = self.config.property_package.component_list
        solutes = self.config.property_package.solute_set
        phase_set = self.config.property_package.phase_list
        self.phase_target_set = phase_set * self.target_set

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.process_flow = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.process_flow.add_state_blocks(has_phase_equilibrium=True)
        self.process_flow.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.process_flow.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=False
        )
        # self.process_flow.add_phase_component_balance()
        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        prop_in = self.process_flow.properties_in[0]

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.air_water_ratio_param = Param(
            initialize=3.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor multiplied by minimum air-to-water flow ratio for operating air-to-water flow ratio",
        )

        @self.Expression(doc="Operational air-to-water ratio")
        def air_water_ratio_op(b):
            return prop_in.flow_vol_phase["Vap"] / prop_in.flow_vol_phase["Liq"]

        self.target_reduction_frac = Param(
            self.target_set,
            initialize=0.9,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Fractional reduction for target component",
        )
        self.overall_mass_transfer_coeff_sf = Param(
            initialize=0.7,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Safety factor for overall mass transfer coeff",
        )

        @self.Expression(self.target_set, doc="Remaining fraction of target component")
        def target_remaining_frac(b, j):
            return 1 - b.target_reduction_frac[j]

        # self.conc_mass_interface_comp = Var(
        #     self.target_set,
        #     initialize=1,
        #     bounds=(0, None),
        #     units=pyunits.kg / pyunits.m**3,
        #     doc="Concentration at air-water interface for target component",
        # )

        self.packing_surface_area_total = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.m**-1,
            doc="Total specific surface area of packing.",
        )

        self.packing_surface_area_wetted = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.m**-1,
            doc="Wetted specific surface area of packing.",
        )

        self.packing_diam_nominal = Var(
            initialize=0.1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Nominal diameter of packing material.",
        )

        self.packing_factor = Var(
            initialize=0.1,
            bounds=(0, None),
            units=pyunits.m**-1,
            doc="Packing factor.",
        )

        self.packing_surf_tension = Var(
            initialize=0.033,
            bounds=(0, None),
            units=pyunits.kg * pyunits.s**-2,
            doc="Surface tension of packing",
        )

        self.surf_tension_water = Var(
            initialize=0.05,
            bounds=(0, None),
            units=pyunits.kg * pyunits.s**-2,
            doc="Surface tension of water",
        )

        self.stripping_factor = Var(
            self.target_set,
            initialize=2,
            bounds=(1, 20),
            units=pyunits.dimensionless,
            doc="Stripping factor",
        )

        self.air_water_ratio_min = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Minumum air-to-water ratio",
        )

        self.tower_height = Var(
            initialize=1,
            bounds=(0, 14),
            units=pyunits.m,
            doc="Height of packed tower",
        )

        self.mass_transfer_coeff = Var(
            self.phase_target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Mass transfer coefficient in tower",
        )

        mass_load_bounds = {"Liq": (0.75, 44), "Vap": (0.013, 1.8)}  # Edzvald, 2011

        self.mass_loading_rate = Var(
            phase_set,
            initialize=1,
            bounds=extract_data(mass_load_bounds),
            units=pyunits.kg / (pyunits.m**2 * pyunits.s),
            doc="Mass loading rate in tower",
        )

        self.height_transfer_unit = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Height of one transfer unit",
        )

        self.number_transfer_unit = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.N_Re = Var(
            phase_set,
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Reynolds number",
        )

        self.N_Fr = Var(
            initialize=0.1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Froude number",
        )

        self.N_We = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Weber number",
        )

        self.N_Sc = Var(
            self.phase_target_set,
            initialize=700,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.wetability_parameter = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Wettability parameter",
        )

        @self.Expression()
        def packing_efficiency_number(b):
            return b.packing_surface_area_total * b.packing_diam_nominal

        # self.packing_efficiency_number = Var(
        #     initialize=30,
        #     bounds=(0, None),
        #     units=pyunits.dimensionless,
        #     doc="Packing efficiency number",
        # )

        self.overall_mass_transfer_coeff = Var(
            self.target_set,
            initialize=30,
            bounds=(0, None),
            units=pyunits.s**-1,
            doc="Overall mass transfer coeff (K_L*a)",
        )

        self.build_oto()

        self.pressure_drop_gradient = Var(
            initialize=100,
            bounds=(20, 1200),
            units=pyunits.Pa / pyunits.m,
            doc="Pressure drop per length of packed bed",
        )

        @self.Constraint()
        def eq_oto_F(b):
            # return b.pressure_drop_gradient == 10**b.oto_F
            return log10(b.pressure_drop_gradient) == b.oto_F

        @self.Expression()
        def pressure_drop(b):
            return b.pressure_drop_gradient * b.tower_height

        @self.Expression(doc="Cross sectional area of tower")
        def tower_area(b):
            return prop_in.flow_mass_phase["Liq"] / b.mass_loading_rate["Liq"]

        @self.Expression(doc="Diameter of tower")
        def tower_diam(b):
            return ((4 * b.tower_area) / Constants.pi) ** 0.5

        @self.Expression(doc="Diameter of tower")
        def tower_volume(b):
            return b.tower_area * b.tower_height

        @self.Constraint(phase_set, doc="Reynolds number by phase")
        def eq_Re(b, p):
            if p == "Liq":
                surf_area = b.packing_surface_area_wetted
            if p == "Vap":
                surf_area = b.packing_surface_area_total
            return b.N_Re[p] == b.mass_loading_rate[p] / (
                surf_area * prop_in.visc_d_phase[p]
            )

        @self.Constraint(self.phase_target_set, doc="Schmidt number by phase")
        def eq_Sc(b, p, j):
            return (
                b.N_Sc[p, j]
                * (prop_in.dens_mass_phase[p] * prop_in.diffus_phase_comp[p, j])
                == prop_in.visc_d_phase[p]
            )

        @self.Constraint(doc="Sherwood number")
        def eq_Sh(b):
            return (
                b.N_Sh * prop_in.visc_d_phase["Liq"] * Constants.acceleration_gravity
                == prop_in.dens_mass_phase["Liq"]
            )

        @self.Constraint(doc="Froude number")
        def eq_Fr(b):
            return (
                b.N_Fr
                * prop_in.dens_mass_phase["Liq"] ** 2
                * Constants.acceleration_gravity
                == b.mass_loading_rate["Liq"] ** 2 * b.packing_surface_area_total
            )

        @self.Constraint(doc="Weber number")
        def eq_We(b):
            return (
                b.N_We
                * prop_in.dens_mass_phase["Liq"]
                * b.packing_surface_area_total
                * b.surf_tension_water
                == b.mass_loading_rate["Liq"] ** 2
            )

        @self.Constraint(phase_set)
        def eq_mass_loading_rate(b, p):
            if p == "Vap":
                dens_vap = prop_in.dens_mass_phase[p]
                dens_liq = prop_in.dens_mass_phase["Liq"]
                visc_liq = prop_in.visc_d_phase["Liq"]
                M = b.oto_M
                return (
                    b.mass_loading_rate[p]
                    == (
                        (M * dens_vap * (dens_liq - dens_vap))
                        / (b.packing_factor * (visc_liq**0.1))
                    )
                    ** 0.5
                )
            if p == "Liq":
                return (
                    b.mass_loading_rate[p]
                    == (b.mass_loading_rate["Vap"] * prop_in.flow_mass_phase[p])
                    / prop_in.flow_mass_phase["Vap"]
                )

        @self.Constraint(self.target_set, doc="Stripping factor equation")
        def eq_stripping_factor(b, j):
            return (
                b.stripping_factor[j]
                == b.air_water_ratio_op * prop_in.henry_constant_comp[j]
            )

        @self.Constraint(self.target_set, doc="Minimum air-to-water ratio")
        def eq_air_water_ratio_min(b, j):
            c0 = prop_in.conc_mass_phase_comp["Liq", j]
            return b.air_water_ratio_min * (c0 * prop_in.henry_constant_comp[j]) == (
                c0 - c0 * b.target_remaining_frac[j]
            )

        @self.Constraint(doc="Operational air-to-water ratio")
        def eq_air_water_ratio_op(b):
            return (
                prop_in.flow_vol_phase["Vap"]
                == (b.air_water_ratio_min * b.air_water_ratio_param)
                * prop_in.flow_vol_phase["Liq"]
            )

        # @self.Constraint(self.target_set, doc="Concentration at air-water interface")
        # def eq_conc_mass_interface_comp(b, j):
        #     c0 = prop_in.conc_mass_phase_comp["Liq", j]
        #     return b.conc_mass_interface_comp[j] * prop_in.henry_constant_comp[
        #         j
        #     ] == b.air_water_ratio_op * (c0 - c0 * b.target_remaining_frac[j])

        @self.Constraint(phase_set, doc="Iosthermal constraint")
        def eq_isothermal(b, p):
            return (
                b.process_flow.properties_in[0].temperature[p]
                == b.process_flow.properties_out[0].temperature[p]
            )

        # @self.Constraint(self.phase_target_set, doc="Effluent concentration")
        # def eq_conc_out(b, p, j):
        #     return (
        #         b.process_flow.properties_in[0].conc_mass_phase_comp[p, j]
        #         * b.target_remaining_frac[j]
        #         == b.process_flow.properties_out[0].conc_mass_phase_comp[p, j]
        #     )

        # @self.Constraint(self.target_set, doc="Mass transfer term Liq >> Vap")
        # def eq_mass_transfer_cv(b, j):
        #     return (
        #         b.process_flow.mass_transfer_term[0, "Vap", j]
        #         == -b.process_flow.mass_transfer_term[0, "Liq", j]
        #     )

        self.process_flow.mass_transfer_term[0, "Liq", "H2O"].fix(0)
        self.process_flow.mass_transfer_term[0, "Liq", "Air"].fix(0)
        self.process_flow.mass_transfer_term[0, "Vap", "H2O"].fix(0)
        self.process_flow.mass_transfer_term[0, "Vap", "Air"].fix(0)

    def build_oto(self):
        """

        OTO = Onda, Takeuchi, Okumoto
        Method to build parameters, variables, and constraints for OTO Model.

        Onda, K., Takeuchi, H., & Okumoto, Y. (1968).
        Mass Transfer Coefficients between Gas and Liquid Phases in Packed Columns.
        Journal of Chemical Engineering of Japan, 1(1), 56-62. doi:10.1252/jcej.1.56

        """
        prop_in = self.process_flow.properties_in[0]
        phase_set = self.config.property_package.phase_list

        self.oto_E = E = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO E parameter",
        )

        self.oto_F = F = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO F parameter",
        )

        self.oto_M = M = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO M parameter",
        )

        # a0 = -6.6599 + 4.3077*F - 1.3503*F^2 + 0.15931*F^3
        self.oto_a0_param1 = a01 = Param(
            initialize=-6.6599,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a0 term, first parameter",
        )
        self.oto_a0_param2 = a02 = Param(
            initialize=4.3077,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a0 term, second parameter",
        )
        self.oto_a0_param3 = a03 = Param(
            initialize=-1.3503,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a0 term, third parameter",
        )
        self.oto_a0_param4 = a04 = Param(
            initialize=0.15931,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a0 term, fourth parameter",
        )

        self.oto_a0 = a0 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a0 term",
        )

        @self.Constraint(doc="OTO a0 equation")
        def eq_oto_a0(b):
            return a0 == a01 + a02 * F + a03 * F**2 + a04 * F**3

        # a1 = 3.0945 - 4.3512*F + 1.6240*F^2 - 0.20855*F^3
        self.oto_a1_param1 = a11 = Param(
            initialize=3.0945,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a1 term, first parameter",
        )
        self.oto_a1_param2 = a12 = Param(
            initialize=-4.3512,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a1 term, second parameter",
        )
        self.oto_a1_param3 = a13 = Param(
            initialize=1.6240,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a1 term, third parameter",
        )
        self.oto_a1_param4 = a14 = Param(
            initialize=-0.20855,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a1 term, fourth parameter",
        )
        self.oto_a1 = a1 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a1 term",
        )

        @self.Constraint(doc="OTO a1 equation")
        def eq_oto_a1(b):
            return a1 == a11 + a12 * F + a13 * F**2 + a14 * F**3

        # a2 = 1.7611 - 2.3394*F + 0.89914*F^2 - 0.11597*F^3
        self.oto_a2_param1 = a21 = Param(
            initialize=1.7611,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a2 term, first parameter",
        )
        self.oto_a2_param2 = a22 = Param(
            initialize=-2.3394,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a2 term, second parameter",
        )
        self.oto_a2_param3 = a23 = Param(
            initialize=0.89914,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a2 term, third parameter",
        )
        self.oto_a2_param4 = a24 = Param(
            initialize=-0.115971,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a2 term, fourth parameter",
        )

        self.oto_a2 = a2 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO correlation: Pressure drop a2 term",
        )

        @self.Constraint(doc="OTO a2 equation")
        def eq_oto_a2(b):
            return a2 == a21 + a22 * F + a23 * F**2 + a24 * F**3

        @self.Constraint(doc="OTO E Parameter")
        def eq_oto_E(b):
            dens_vap = b.process_flow.properties_in[0].dens_mass_phase["Vap"]
            dens_liq = b.process_flow.properties_in[0].dens_mass_phase["Liq"]
            return E == -1 * (
                log10(
                    (b.air_water_ratio_op)
                    * (dens_vap / dens_liq - (dens_vap / dens_liq) ** 2) ** 0.5
                )
            )

        @self.Constraint(doc="OTO M Parameter")
        def eq_oto_M(b):
            return M == 10 ** (a0 + a1 * E + a2 * E**2)

        self.oto_aw_param = aw_param = Param(
            initialize=-1.45,
            units=pyunits.dimensionless,
            doc="OTO wetted surface area of packing correlation parameter",
        )
        self.oto_aw_exp1 = aw_exp1 = Param(
            initialize=0.75,
            units=pyunits.dimensionless,
            doc="OTO wetted surface area of packing correlation - exponent 1",
        )
        self.oto_aw_exp2 = aw_exp2 = Param(
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="OTO wetted surface area of packing correlation - exponent 2",
        )
        self.oto_aw_exp3 = aw_exp3 = Param(
            initialize=-0.05,
            units=pyunits.dimensionless,
            doc="OTO wetted surface area of packing correlation - exponent 3",
        )
        self.oto_aw_exp4 = aw_exp4 = Param(
            initialize=0.2,
            units=pyunits.dimensionless,
            doc="OTO wetted surface area of packing correlation - exponent 4",
        )

        @self.Constraint(doc="OTO equation for wetted surface area of packing material")
        def eq_packing_surf_area_wetted(b):
            a_t = b.packing_surface_area_total
            a_w = b.packing_surface_area_wetted
            sigma_c = b.packing_surf_tension
            sigma_w = b.surf_tension_water
            Lm = b.mass_loading_rate["Liq"]
            visc_liq = prop_in.visc_d_phase["Liq"]
            dens_liq = prop_in.dens_mass_phase["Liq"]
            g = Constants.acceleration_gravity

            exp_term = (
                aw_param
                * (sigma_c / sigma_w) ** aw_exp1
                * (Lm / (a_t * visc_liq)) ** aw_exp2
                * ((Lm**2 * a_t) / (dens_liq**2 * g)) ** aw_exp3
                * (Lm**2 / (dens_liq * a_t * sigma_w)) ** aw_exp4
            )
            return a_w / a_t == 1 - exp(exp_term)

        self.oto_liq_mass_xfr_param = kl_param = Param(
            initialize=0.0051,
            units=pyunits.m / pyunits.s,
            doc="OTO liquid mass transfer correlation parameter",
        )
        self.oto_liq_mass_xfr_exp1 = kl_exp1 = Param(
            initialize=0.667,
            units=pyunits.dimensionless,
            doc="OTO liquid mass transfer correlation Re exponent",
        )
        self.oto_liq_mass_xfr_exp2 = kl_exp2 = Param(
            initialize=-0.5,
            units=pyunits.dimensionless,
            doc="OTO liquid mass transfer correlation Sc exponent",
        )
        self.oto_liq_mass_xfr_exp3 = kl_exp3 = Param(
            initialize=0.4,
            units=pyunits.dimensionless,
            doc="OTO liquid mass transfer correlation Er exponent",
        )
        self.oto_liq_mass_xfr_exp4 = kl_exp4 = Param(
            initialize=-0.3334,
            units=pyunits.dimensionless,
            doc="OTO liquid mass transfer correlation Sh exponent",
        )

        self.oto_gas_mass_xfr_param = kg_param = Param(
            initialize=5.23,
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation parameter",
        )
        self.oto_gas_mass_xfr_exp1 = kg_exp1 = Param(
            initialize=0.7,
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation Re exponent",
        )
        self.oto_gas_mass_xfr_exp2 = kg_exp2 = Param(
            initialize=0.3334,
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation Sc exponent",
        )
        self.oto_gas_mass_xfr_exp3 = kg_exp3 = Param(
            initialize=-2,
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation Er exponent",
        )

        @self.Constraint(
            self.phase_target_set, doc="OTO model mass transfer coefficient equation"
        )
        def eq_mass_transfer_coeff(b, p, j):
            if p == "Liq":
                return (
                    b.mass_transfer_coeff[p, j]
                    == kl_param
                    * b.N_Re[p] ** kl_exp1
                    * b.N_Sc[p, j] ** kl_exp2
                    * b.packing_efficiency_number**kl_exp3
                    * b.N_Sh**kl_exp4
                )
            if p == "Vap":
                return (
                    b.mass_transfer_coeff[p, j]
                    == kg_param
                    * (b.packing_surface_area_total * prop_in.diffus_phase_comp[p, j])
                    * b.N_Re[p] ** kg_exp1
                    * b.N_Sc[p, j] ** kg_exp2
                    * b.packing_efficiency_number**kg_exp3
                )

        @self.Constraint(self.target_set, doc="OTO model - overall mass transfer coeff")
        def eq_overall_mass_transfer_coeff(b, j):
            KLa = b.overall_mass_transfer_coeff[j]
            KLa_sf = b.overall_mass_transfer_coeff_sf
            kl_aw = b.mass_transfer_coeff["Liq", j] * b.packing_surface_area_wetted
            kg_aw_H = (
                b.mass_transfer_coeff["Vap", j]
                * b.packing_surface_area_wetted
                * prop_in.henry_constant_comp[j]
            )
            return KLa == (kl_aw * kg_aw_H) / (kl_aw + kg_aw_H) * KLa_sf

    def initialize_build(
        self,
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
        pass
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.process_flow.properties_in.initialize(
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
            self.state_args = state_args = {}
            state_dict = self.process_flow.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        self.state_args_out = state_args_out = deepcopy(state_args)

        for p, j in self.process_flow.properties_out.phase_component_set:
            if p == "Liq" and j in self.target_set:
                state_args_out["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, j)] * value(self.target_remaining_frac[j])
                )

        self.process_flow.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        # self.process_flow.properties_in[0].flow_mass_phase_comp["Vap", "Air"].unfix()

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # if iscale.get_scaling_factor(self.conc_mass_interface_comp) is None:
        #     iscale.set_scaling_factor(self.conc_mass_interface_comp, 1e2)

        if iscale.get_scaling_factor(self.packing_surface_area_total) is None:
            iscale.set_scaling_factor(self.packing_surface_area_total, 1e-3)

        if iscale.get_scaling_factor(self.packing_surface_area_wetted) is None:
            iscale.set_scaling_factor(self.packing_surface_area_wetted, 1e-3)

        if iscale.get_scaling_factor(self.packing_diam_nominal) is None:
            iscale.set_scaling_factor(self.packing_diam_nominal, 10)

        if iscale.get_scaling_factor(self.packing_factor) is None:
            iscale.set_scaling_factor(self.packing_factor, 1e-2)

        if iscale.get_scaling_factor(self.packing_surf_tension) is None:
            iscale.set_scaling_factor(self.packing_surf_tension, 10)

        if iscale.get_scaling_factor(self.surf_tension_water) is None:
            iscale.set_scaling_factor(self.surf_tension_water, 10)

        if iscale.get_scaling_factor(self.stripping_factor) is None:
            iscale.set_scaling_factor(self.stripping_factor, 0.1)

        if iscale.get_scaling_factor(self.air_water_ratio_min) is None:
            iscale.set_scaling_factor(self.air_water_ratio_min, 0.1)

        if iscale.get_scaling_factor(self.tower_height) is None:
            iscale.set_scaling_factor(self.tower_height, 0.1)

        if iscale.get_scaling_factor(self.mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.mass_transfer_coeff, 1e3)

        if iscale.get_scaling_factor(self.mass_loading_rate) is None:
            for p, v in self.mass_loading_rate.items():
                if p == "Liq":
                    sf = 1e-2
                if p == "Vap":
                    sf = 10
                iscale.set_scaling_factor(v, sf)

        if iscale.get_scaling_factor(self.height_transfer_unit) is None:
            iscale.set_scaling_factor(self.height_transfer_unit, 1)

        if iscale.get_scaling_factor(self.number_transfer_unit) is None:
            iscale.set_scaling_factor(self.number_transfer_unit, 0.1)

        if iscale.get_scaling_factor(self.N_Re) is None:
            iscale.set_scaling_factor(self.N_Re, 1e2)

        if iscale.get_scaling_factor(self.N_Fr) is None:
            iscale.set_scaling_factor(self.N_Fr, 10)

        if iscale.get_scaling_factor(self.N_We) is None:
            iscale.set_scaling_factor(self.N_We, 10)

        if iscale.get_scaling_factor(self.N_Sc) is None:
            for (p, j), v in self.N_Sc.items():
                if p == "Liq" and j in self.target_set:
                    sf = 1e-4
                if p == "Vap" and j in self.target_set:
                    sf = 0.1
                iscale.set_scaling_factor(v, sf)

        if iscale.get_scaling_factor(self.N_Sh) is None:
            iscale.set_scaling_factor(self.N_Sh, 1e-5)

        if iscale.get_scaling_factor(self.wetability_parameter) is None:
            iscale.set_scaling_factor(self.wetability_parameter, 1e-2)

        if iscale.get_scaling_factor(self.overall_mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.overall_mass_transfer_coeff, 1e2)

        if iscale.get_scaling_factor(self.oto_E) is None:
            iscale.set_scaling_factor(self.oto_E, 0.1)

        if iscale.get_scaling_factor(self.oto_F) is None:
            iscale.set_scaling_factor(self.oto_F, 0.1)

        if iscale.get_scaling_factor(self.oto_M) is None:
            iscale.set_scaling_factor(self.oto_M, 1e2)

        if iscale.get_scaling_factor(self.oto_a0) is None:
            iscale.set_scaling_factor(self.oto_a0, 0.1)

        if iscale.get_scaling_factor(self.oto_a1) is None:
            iscale.set_scaling_factor(self.oto_a1, 10)

        if iscale.get_scaling_factor(self.oto_a2) is None:
            iscale.set_scaling_factor(self.oto_a2, 10)

        if iscale.get_scaling_factor(self.pressure_drop_gradient) is None:
            iscale.set_scaling_factor(self.pressure_drop_gradient, 0.1)

        # if iscale.get_scaling_factor(self.conc_mass_interface_comp) is None:
        #     iscale.set_scaling_factor(self.conc_mass_interface_comp, 1e2)

        # if iscale.get_scaling_factor(self.conc_mass_interface_comp) is None:
        #     iscale.set_scaling_factor(self.conc_mass_interface_comp, 1e2)

        # if iscale.get_scaling_factor(self.conc_mass_interface_comp) is None:
        #     iscale.set_scaling_factor(self.conc_mass_interface_comp, 1e2)


    def _get_stream_table_contents(self, time_point=0):
        pass

    #     return create_stream_table_dataframe(
    #         {},
    #         time_point=time_point,
    #     )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        pass
