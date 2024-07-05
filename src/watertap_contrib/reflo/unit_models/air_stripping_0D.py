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
import itertools

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    log10,
    exp,
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum, extract_data
from idaes.core.util.exceptions import InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.core import ControlVolume0DBlock, InitializationMixin

from watertap_contrib.reflo.costing.units.air_stripping import cost_air_stripping

__author__ = "Kurban Sitterley"


"""
REFERENCES:

Onda, K., Takeuchi, H., & Okumoto, Y. (1968). 
Mass Transfer Coefficients between Gas and Liquid Phases in Packed Columns. 
Journal of Chemical Engineering of Japan, 1(1), 56-62. doi:10.1252/jcej.1.56

Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). 
Chap. 7, 14 in MWH's Water Treatment: Principles and Design (3rd ed.). doi:10.1002/9781118131473

Edzvald, J. (2011). Chapter 6: Gas-Liquid Processes: Principles and Applications. 
In Water Quality & Treatment: A Handbook on Drinking Water (6 ed.): American Water Works Association.
ISBN 9780071630115

"""


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
        **MomentumBalanceType.momentumTotal** - single momentum balance for material,
        **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "target",
        ConfigValue(
            default="VOC",
            domain=str,
            description="Designates targeted species for removal.",
        ),
    )

    CONFIG.declare(
        "packing_material",
        ConfigValue(
            default=PackingMaterial.PVC,
            domain=In(PackingMaterial),
            description="Packing material used in tower.",
        ),
    )

    def build(self):
        super().build()

        target = self.config.target
        self.target_set = Set(initialize=[target])
        comps = self.config.property_package.component_list
        solutes = self.config.property_package.solute_set
        phase_set = self.config.property_package.phase_list
        self.liq_target_set = Set(initialize=["Liq"]) * self.target_set
        phase_target_idx = list(itertools.product(["Liq", "Vap"], self.target_set))
        self.phase_target_set = Set(initialize=phase_target_idx)

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

        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True,
        )

        prop_in = self.process_flow.properties_in[0]
        prop_out = self.process_flow.properties_out[0]

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)

        self.process_flow.mass_transfer_term[0, "Liq", "H2O"].fix(0)
        self.process_flow.mass_transfer_term[0, "Vap", "Air"].fix(0)

        self.air_water_ratio_param = Param(
            initialize=3.5,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor multiplied by minimum air-to-water flow ratio for operating air-to-water flow ratio",
        )

        self.pressure_drop_tower_param = Param(
            initialize=275,
            mutable=True,
            units=(pyunits.newton * pyunits.s**2) / pyunits.m**4,
            doc="Factor to calculate pressure drop through tower",
        )

        self.tower_height_factor = Param(
            initialize=1.2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Factor to calculate tower height",
        )

        self.tower_port_diameter = Param(
            initialize=6,
            mutable=True,
            units=pyunits.inches,
            doc="Diameter of tower access ports: 2-24 inches",
        )

        self.tower_pipe_diameter = Param(
            initialize=6,
            mutable=True,
            units=pyunits.inches,
            doc="Diameter of tower inlet and outlet piping: 2-24 inches",
        )

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

        self.blower_efficiency = Param(
            initialize=0.4,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Blower efficiency",
        )

        self.pump_efficiency = Param(
            initialize=0.85,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Pump efficiency",
        )

        self.power_blower_denom_coeff = Param(
            initialize=0.283,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Blower power equation denominator coefficient",
        )

        self.power_blower_exponent = Param(
            initialize=0.283,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Blower power equation exponent",
        )

        self.blower_power = Var(
            initialize=50,
            bounds=(0, None),
            units=pyunits.kilowatt,
            doc="Air blower power requirement",
        )

        self.pump_power = Var(
            initialize=50,
            bounds=(0, None),
            units=pyunits.kilowatt,
            doc="Water pump power requirement",
        )

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

        self.packing_height = Var(
            initialize=1,
            bounds=(0, 25),
            units=pyunits.m,
            doc="Height of packed tower",
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
            self.target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Height of one transfer unit",
        )

        self.number_transfer_unit = Var(
            self.target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.pressure_drop_gradient = Var(
            initialize=100,
            bounds=(20, 1200),
            units=pyunits.Pa / pyunits.m,
            doc="Pressure drop per length of packed bed",
        )

        self.overall_mass_transfer_coeff = Var(
            self.target_set,
            initialize=30,
            bounds=(0, None),
            units=pyunits.s**-1,
            doc="Overall mass transfer coefficient (K_L*a)",
        )

        @self.Expression(doc="Operational air-to-water ratio")
        def air_water_ratio_op(b):
            return prop_in.flow_vol_phase["Vap"] / prop_in.flow_vol_phase["Liq"]

        @self.Expression(doc="Packing efficiency number")
        def packing_efficiency_number(b):
            return b.packing_surface_area_total * b.packing_diam_nominal

        @self.Expression(doc="Cross sectional area of tower")
        def tower_area(b):
            return prop_in.flow_mass_phase["Liq"] / b.mass_loading_rate["Liq"]

        @self.Expression(doc="Diameter of tower")
        def tower_diam(b):
            return ((4 * b.tower_area) / Constants.pi) ** 0.5

        @self.Expression(doc="Height of tower")
        def tower_height(b):
            return b.packing_height * b.tower_height_factor

        @self.Expression(doc="Volume of tower")
        def tower_volume(b):
            return b.tower_area * b.tower_height

        @self.Expression(doc="Volume of packing")
        def packing_volume(b):
            return b.tower_area * b.packing_height

        @self.Expression(self.target_set, doc="Remaining fraction of target component")
        def target_remaining_frac(b, j):
            return 1 - b.target_reduction_frac[j]

        @self.Expression(doc="Pressure drop")
        def pressure_drop(b):
            return b.pressure_drop_gradient * b.tower_height

        @self.Expression(
            doc="Pressure drop through demister, packing support plate, duct work, and inlet/outlet"
        )
        def pressure_drop_tower(b):
            return pyunits.convert(
                (prop_in.flow_vol_phase["Vap"] / b.tower_area) ** 2
                * b.pressure_drop_tower_param,
                to_units=pyunits.Pa,
            )

        @self.Expression(self.phase_target_set, doc="Schmidt Number")
        def N_Sc(b, p, j):
            return pyunits.convert(
                prop_in.visc_d_phase[p]
                / (prop_in.dens_mass_phase[p] * prop_in.diffus_phase_comp[p, j]),
                to_units=pyunits.dimensionless,
            )

        @self.Expression(doc="Reynolds number")
        def N_Re(b):
            return pyunits.convert(
                b.mass_loading_rate["Liq"]
                / (b.packing_surface_area_total * prop_in.visc_d_phase["Liq"]),
                to_units=pyunits.dimensionless,
            )

        @self.Expression(doc="Froude number")
        def N_Fr(b):
            return (b.mass_loading_rate["Liq"] ** 2 * b.packing_surface_area_total) / (
                prop_in.dens_mass_phase["Liq"] ** 2 * Constants.acceleration_gravity
            )

        @self.Expression(doc="Weber number")
        def N_We(b):
            return b.mass_loading_rate["Liq"] ** 2 / (
                prop_in.dens_mass_phase["Liq"]
                * b.packing_surface_area_total
                * b.surf_tension_water
            )

        self.build_oto()

        @self.Constraint(phase_set)
        def eq_mass_loading_rate(b, p, doc="Mass loading rate equation per phase"):
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
            ce = c0 * b.target_remaining_frac[j]
            return b.air_water_ratio_min * (c0 * prop_in.henry_constant_comp[j]) == (
                c0 - ce
            )

        @self.Constraint(self.target_set, doc="Overall mass transfer coeff")
        def eq_overall_mass_transfer_coeff(b, j):
            KLa = b.overall_mass_transfer_coeff[j]
            KLa_sf = b.overall_mass_transfer_coeff_sf
            kl_aw = b.oto_mass_transfer_coeff["Liq", j] * b.packing_surface_area_wetted
            kg_aw_H = (
                b.oto_mass_transfer_coeff["Vap", j]
                * b.packing_surface_area_wetted
                * prop_in.henry_constant_comp[j]
            )
            return KLa == (kl_aw * kg_aw_H) / (kl_aw + kg_aw_H) * KLa_sf

        @self.Constraint(phase_set, doc="Iosthermal constraint")
        def eq_isothermal(b, p):
            return prop_in.temperature[p] == prop_out.temperature[p]

        @self.Constraint(self.target_set)
        def eq_liq_to_vap(b, j):
            return (
                prop_out.flow_mass_phase_comp["Vap", j]
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        @self.Constraint(doc="Control volume deltaP")
        def eq_deltaP(b):
            return b.process_flow.deltaP[0] == (b.pressure_drop + b.pressure_drop_tower)

        @self.Constraint(["Liq"], self.target_set, doc="Effluent concentration")
        def eq_conc_out(b, p, j):
            return (
                b.process_flow.properties_in[0].conc_mass_phase_comp[p, j]
                * b.target_remaining_frac[j]
                == b.process_flow.properties_out[0].conc_mass_phase_comp[p, j]
            )

        @self.Constraint(self.target_set, doc="Height of transfer unit")
        def eq_htu(b, j):
            return b.height_transfer_unit[j] == prop_in.flow_vol_phase["Liq"] / (
                b.tower_area * b.overall_mass_transfer_coeff[j]
            )

        @self.Constraint(self.target_set, doc="Number of transfer unit")
        def eq_ntu(b, j):
            S = b.stripping_factor[j]
            log_term = (
                1
                + (
                    prop_in.conc_mass_phase_comp["Liq", j]
                    / prop_out.conc_mass_phase_comp["Liq", j]
                )
                * (S - 1)
            ) / S
            return b.number_transfer_unit[j] == (S / (S - 1)) * log(log_term)

        @self.Constraint(self.target_set, doc="Packing height calculation")
        def eq_packing_height(b, j):
            return (
                b.packing_height
                == b.height_transfer_unit[j] * b.number_transfer_unit[j]
            )

        @self.Constraint(doc="Pumping power required.")
        def eq_pump_power(b):
            return (
                b.pump_power
                == pyunits.convert(
                    prop_in.flow_mass_phase["Liq"]
                    * b.tower_height
                    * Constants.acceleration_gravity,
                    to_units=pyunits.kilowatt,
                )
                / b.pump_efficiency
            )

        @self.Constraint(doc="Blower power required.")
        def eq_blower_power(b):
            pressure_ambient = 101325 * pyunits.Pa
            return b.blower_power == pyunits.convert(
                (
                    prop_in.flow_mass_phase["Vap"]
                    * Constants.gas_constant
                    * prop_in.temperature["Vap"]
                )
                / (
                    prop_in.mw_comp["Air"]
                    * b.power_blower_denom_coeff
                    * b.blower_efficiency
                )
                * (
                    (prop_out.pressure / pressure_ambient) ** b.power_blower_exponent
                    - 1
                ),
                to_units=pyunits.kilowatt,
            )

    def build_oto(self):
        """

        OTO = Onda, Takeuchi, Okumoto
        Method to build parameters, variables, and constraints for OTO Model.

        Onda, K., Takeuchi, H., & Okumoto, Y. (1968).
        Mass Transfer Coefficients between Gas and Liquid Phases in Packed Columns.
        Journal of Chemical Engineering of Japan, 1(1), 56-62. doi:10.1252/jcej.1.56

        """
        prop_in = self.process_flow.properties_in[0]

        self.oto_E = E = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: E parameter",
        )

        @self.Constraint(doc="OTO model E calculation")
        def eq_oto_E(b):
            dens_vap = b.process_flow.properties_in[0].dens_mass_phase["Vap"]
            dens_liq = b.process_flow.properties_in[0].dens_mass_phase["Liq"]
            return E == -1 * (
                log10(
                    (b.air_water_ratio_op)
                    * (dens_vap / dens_liq - (dens_vap / dens_liq) ** 2) ** 0.5
                )
            )

        self.oto_F = F = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: F parameter",
        )

        @self.Constraint(doc="OTO model F calculation")
        def eq_oto_F(b):
            return log10(b.pressure_drop_gradient * (pyunits.m / pyunits.Pa)) == b.oto_F

        # a0 = -6.6599 + 4.3077*F - 1.3503*F^2 + 0.15931*F^3
        self.oto_a0_param1 = a01 = Param(
            initialize=-6.6599,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a0 term, first parameter",
        )
        self.oto_a0_param2 = a02 = Param(
            initialize=4.3077,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a0 term, second parameter",
        )
        self.oto_a0_param3 = a03 = Param(
            initialize=-1.3503,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a0 term, third parameter",
        )
        self.oto_a0_param4 = a04 = Param(
            initialize=0.15931,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a0 term, fourth parameter",
        )

        self.oto_a0 = a0 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a0 term",
        )

        @self.Constraint(doc="OTO a0 equation")
        def eq_oto_a0(b):
            return a0 == a01 + a02 * F + a03 * F**2 + a04 * F**3

        # a1 = 3.0945 - 4.3512*F + 1.6240*F^2 - 0.20855*F^3
        self.oto_a1_param1 = a11 = Param(
            initialize=3.0945,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a1 term, first parameter",
        )
        self.oto_a1_param2 = a12 = Param(
            initialize=-4.3512,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a1 term, second parameter",
        )
        self.oto_a1_param3 = a13 = Param(
            initialize=1.6240,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a1 term, third parameter",
        )
        self.oto_a1_param4 = a14 = Param(
            initialize=-0.20855,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a1 term, fourth parameter",
        )
        self.oto_a1 = a1 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a1 term",
        )

        @self.Constraint(doc="OTO a1 equation")
        def eq_oto_a1(b):
            return a1 == a11 + a12 * F + a13 * F**2 + a14 * F**3

        # a2 = 1.7611 - 2.3394*F + 0.89914*F^2 - 0.11597*F^3
        self.oto_a2_param1 = a21 = Param(
            initialize=1.7611,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a2 term, first parameter",
        )
        self.oto_a2_param2 = a22 = Param(
            initialize=-2.3394,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a2 term, second parameter",
        )
        self.oto_a2_param3 = a23 = Param(
            initialize=0.89914,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a2 term, third parameter",
        )
        self.oto_a2_param4 = a24 = Param(
            initialize=-0.115971,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a2 term, fourth parameter",
        )

        self.oto_a2 = a2 = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: Pressure drop a2 term",
        )

        @self.Constraint(doc="OTO a2 equation")
        def eq_oto_a2(b):
            return a2 == a21 + a22 * F + a23 * F**2 + a24 * F**3

        self.oto_M = M = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="OTO model: M parameter",
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

            exp_term = (
                aw_param
                * (sigma_c / sigma_w) ** aw_exp1
                * (Lm / (a_t * visc_liq)) ** aw_exp2
                * (b.N_Fr) ** aw_exp3
                * (b.N_We) ** aw_exp4
            )
            return a_w / a_t == 1 - exp(exp_term)

        self.oto_liq_mass_xfr_param = kl_param = Param(
            initialize=0.0051,
            units=pyunits.m / pyunits.s,
            doc="OTO liquid mass transfer correlation parameter",
        )
        self.oto_liq_mass_xfr_exp1 = kl_exp1 = Param(
            initialize=(2 / 3),
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
            doc="OTO liquid mass transfer correlation packing efficiency exponent",
        )
        self.oto_liq_mass_xfr_exp4 = kl_exp4 = Param(
            initialize=-(1 / 3),
            units=pyunits.dimensionless,
            doc="OTO liquid mass transfer correlation fourth exponent",
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
            initialize=(1 / 3),
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation Sc exponent",
        )
        self.oto_gas_mass_xfr_exp3 = kg_exp3 = Param(
            initialize=-2,
            units=pyunits.dimensionless,
            doc="OTO gas mass transfer correlation packing efficiency exponent",
        )

        self.oto_mass_transfer_coeff = Var(
            self.phase_target_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="OTO model: phase mass transfer coefficient in tower",
        )

        @self.Expression()
        def oto_kl_term(b):
            return prop_in.dens_mass_phase["Liq"] / (
                prop_in.visc_d_phase["Liq"] * Constants.acceleration_gravity
            )

        @self.Constraint(
            self.phase_target_set, doc="OTO model mass transfer coefficient equation"
        )
        def eq_oto_mass_transfer_coeff(b, p, j):
            if p == "Liq":
                term1 = b.mass_loading_rate[p] / (
                    b.packing_surface_area_wetted * prop_in.visc_d_phase[p]
                )
                return (
                    b.oto_mass_transfer_coeff[p, j]
                    == kl_param
                    * term1**kl_exp1
                    * b.N_Sc[p, j] ** kl_exp2
                    * b.packing_efficiency_number**kl_exp3
                    * b.oto_kl_term**kl_exp4
                )
            if p == "Vap":
                return (
                    b.oto_mass_transfer_coeff[p, j]
                    == kg_param
                    * (b.packing_surface_area_total * prop_in.diffus_phase_comp[p, j])
                    * b.N_Re**kg_exp1
                    * b.N_Sc[p, j] ** kg_exp2
                    * b.packing_efficiency_number**kg_exp3
                )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        calc_from_constr=None,
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
            if j == self.config.target:
                if p == "Vap":
                    state_args_out["flow_mass_phase_comp"][(p, j)] = state_args[
                        "flow_mass_phase_comp"
                    ][("Liq", j)] * value(self.target_remaining_frac[j])

        self.process_flow.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        self.calc_from_constr_dict = {
            "oto_F": "eq_oto_F",
            "oto_a0": "eq_oto_a0",
            "oto_a1": "eq_oto_a1",
            "oto_a2": "eq_oto_a2",
            "oto_E": "eq_oto_E",
            "oto_M": "eq_oto_M",
        }

        if calc_from_constr is not None:
            if not isinstance(calc_from_constr, dict):
                raise InitializationError(
                    "calc_from_constr must be a dict with var, constraint pairs"
                )
            for k, v in calc_from_constr:
                self.calc_from_constr_dict[k] = v

        for v, c in self.calc_from_constr_dict.items():
            axv = getattr(self, v)
            axc = getattr(self, c)
            for i, x in axc.items():
                if i:
                    calculate_variable_from_constraint(axv[i], x)
                else:
                    calculate_variable_from_constraint(axv, x)
        init_log.info("Initialization Step 1c Complete.")

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

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        pf = self.process_flow
        prop_in = pf.properties_in[0]

        if (
            iscale.get_scaling_factor(
                pf.mass_transfer_term[0, "Vap", self.config.target]
            )
            is None
        ):
            iscale.set_scaling_factor(
                pf.mass_transfer_term[0, "Vap", self.config.target], 1e5
            )

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

        if iscale.get_scaling_factor(self.packing_height) is None:
            iscale.set_scaling_factor(self.packing_height, 0.1)

        if iscale.get_scaling_factor(self.oto_mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.oto_mass_transfer_coeff, 1e3)

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

        if iscale.get_scaling_factor(self.blower_power) is None:
            iscale.set_scaling_factor(self.blower_power, 0.1)

        if iscale.get_scaling_factor(self.pump_power) is None:
            iscale.set_scaling_factor(self.pump_power, 0.1)

        iscale.constraint_scaling_transform(
            self.eq_deltaP, iscale.get_scaling_factor(prop_in.pressure)
        )

        for (p, j), c in self.eq_oto_mass_transfer_coeff.items():
            iscale.constraint_scaling_transform(c, 1e6)

        for j, c in self.eq_liq_to_vap.items():
            sf = iscale.get_scaling_factor(pf.mass_transfer_term[0, "Liq", j])
            iscale.constraint_scaling_transform(c, sf)

        for (p, j), c in self.eq_conc_out.items():
            sf = iscale.get_scaling_factor(prop_in.conc_mass_phase_comp[p, j])
            iscale.constraint_scaling_transform(c, sf)

    def _get_stream_table_contents(self, time_point=0):

        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        target = self.config.target
        var_dict = dict()

        var_dict["Wetted surface area of packing"] = self.packing_surface_area_wetted
        var_dict["Nominal diameter of packing"] = self.packing_diam_nominal
        var_dict["Packing factor"] = self.packing_factor
        var_dict[f"Stripping factor [{target}]"] = self.stripping_factor[target]
        var_dict["Air-to-water flow ratio, minimum"] = self.air_water_ratio_min
        var_dict["Packing height"] = self.packing_height
        var_dict["Mass loading rate, liquid"] = self.mass_loading_rate["Liq"]
        var_dict["Mass loading rate, vapor"] = self.mass_loading_rate["Vap"]
        var_dict[f"Height transfer unit [{target}]"] = self.height_transfer_unit[target]
        var_dict[f"Number transfer units [{target}]"] = self.number_transfer_unit[
            target
        ]
        var_dict[
            f"Overall mass transfer coeff (KL_a) [{target}]"
        ] = self.overall_mass_transfer_coeff[target]
        var_dict["Pressure drop gradient"] = self.pressure_drop_gradient
        var_dict["OTO Model - M parameter"] = self.oto_M
        var_dict["OTO Model - E parameter"] = self.oto_E
        var_dict["OTO Model - F parameter"] = self.oto_F
        var_dict[
            f"OTO Model - liquid mass transfer coeff [{target}]"
        ] = self.oto_mass_transfer_coeff["Liq", target]
        var_dict[
            f"OTO Model - vapor mass transfer coeff [{target}]"
        ] = self.oto_mass_transfer_coeff["Vap", target]
        var_dict[
            f"CV mass transfer term [{target}]"
        ] = self.process_flow.mass_transfer_term[time_point, "Liq", target]
        var_dict[f"CV delta P"] = self.process_flow.deltaP[time_point]
        var_dict[f"Blower power required"] = self.blower_power
        var_dict[f"Pump power required"] = self.pump_power

        expr_dict = dict()

        expr_dict["Air-to-water flow ratio, minimum"] = self.air_water_ratio_op
        expr_dict["Pressure drop through tower"] = self.pressure_drop_tower
        expr_dict["Tower height"] = self.tower_height
        expr_dict["Tower diameter"] = self.tower_diam
        expr_dict["Tower volume"] = self.tower_volume
        expr_dict["Tower area"] = self.tower_area
        expr_dict["Packing volume"] = self.packing_volume
        expr_dict[f"Schmidt number, liquid [{target}]"] = self.N_Sc["Liq", target]
        expr_dict[f"Schmidt number, gas [{target}]"] = self.N_Sc["Vap", target]
        expr_dict[f"Reynolds number"] = self.N_Re
        expr_dict[f"Froude number"] = self.N_Fr
        expr_dict[f"Weber number"] = self.N_We

        return {"vars": var_dict, "exprs": expr_dict}

    @property
    def default_costing_method(self):
        return cost_air_stripping
