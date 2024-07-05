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
    Var,
    Param,
    Suffix,
    units as pyunits,
    check_optimal_termination,
    exp,
    log,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import WaterTAP cores
from watertap.core import InitializationMixin

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.util.misc import StrEnum
import idaes.logger as idaeslog
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing.units.vagmd_surrogate import cost_vagmd_surrogate

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


class ModuleType(StrEnum):
    AS7C15L = "AS7C1.5L"
    AS26C72L = "AS7C7.2L"


class CoolingType(StrEnum):
    open = "open"
    closed = "closed"


@declare_process_block_class("VAGMDSurrogateBase")
class VAGMDBaseData(InitializationMixin, UnitModelBlockData):
    """
    Vacuum Air-Gap Membrane Distillation - simulation of one Aquastill module
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
        "property_package_seawater",
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
        "property_package_water",
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

    def build(self):
        super().build()
        self._tech_type = "vagmd_surrogate"

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = (
            self.config.property_package_seawater.get_metadata().get_derived_units
        )

        """
        Model parameters
        """
        self.pump_efficiency = Param(
            initialize=0.6, units=pyunits.dimensionless, doc="Pump efficiency"
        )

        self.heat_exchanger_area = Param(
            initialize=1.34,
            units=pyunits.m**2,
            doc="Effective heat transfer area",
        )

        self.cooling_flow_rate = Param(
            initialize=1265,
            units=pyunits.L / pyunits.h,
            doc="Cooling water volumetric flow rate, fixed to 1265L/h, to maintain vacuum pressure inside the MD module around 200 mbar",
        )

        self.thermal_heat_transfer_coeff = Param(
            initialize=3168,
            units=pyunits.W / pyunits.m**2 / pyunits.K,
            doc="Overall heat transfer coefficient",
        )

        self.cooling_flow_pressure_drop = Param(
            initialize=170, units=pyunits.mbar, doc="Cooling flow pressure drop"
        )

        """
        MD Module type and corresponding area
        """
        if self.config.module_type == ModuleType.AS7C15L:
            self.module_area = Param(
                initialize=7.2, units=pyunits.m**2, doc="Area of module AS7C1.5L"
            )
        else:  # module_type = ModuleType.AS26C72L
            self.module_area = Param(
                initialize=25.92, units=pyunits.m**2, doc="Area of module AS26C7.2L"
            )

        """
        Intermediate variables
        """
        self.permeate_flux = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.L / pyunits.h / pyunits.m**2,
            doc="Permeate flux",
        )

        self.feed_flow_pressure_drop = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.mbar,
            doc="Feed flow pressure drop",
        )

        self.log_mean_temp_dif = Var(
            initialize=0,
            bounds=(0, None),
            units=pyunits.K,
            doc="Log mean temperature difference in the MD module",
        )

        self.thermal_resistance_hot = Var(
            initialize=5e3,
            bounds=(0, None),
            units=pyunits.W / pyunits.K,
            doc="Thermal resistance on the hot side",
        )

        self.thermal_resistance_cold = Var(
            initialize=1e3,
            bounds=(0, None),
            units=pyunits.W / pyunits.K,
            doc="Thermal resistance on the cold side",
        )

        self.number_transfer_units = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="The number of transfer units",
        )

        self.effectiveness_heat_exchanger = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="The effectiveness of the heat exchanger",
        )

        """
        Output variables
        """
        self.thermal_power = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal power requirement",
        )

        self.feed_pump_power_elec = Var(
            initialize=0.005,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electric power for pumping feed water",
        )

        self.cooling_pump_power_elec = Var(
            initialize=0.005,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electric power for pumping cooling water",
        )

        self.cooling_power_thermal = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Cooling power requirement",
        )

        """
        Add block for the batched feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_seawater
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict,
        )

        """
        Add block for the permeated water
        """
        tmp_dict["defined_state"] = False

        self.permeate_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of permeated water",
            **tmp_dict,
        )
        # salinity in permeate is zero
        self.permeate_props[0].flow_mass_phase_comp["Liq", "TDS"].fix(0)

        # permeate temperature is assumed the same as feed water
        @self.Constraint(doc="distillate temperature")
        def eq_distillate_temp(b):
            return b.permeate_props[0].temperature == b.feed_props[0].temperature

        """
        Add block for water at the evaporator inlet
        """
        tmp_dict["defined_state"] = True

        self.evaporator_in_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of evaporator inlet water",
                **tmp_dict,
            )
        )

        """
        Add block for water at the evaporator outlet
        """
        tmp_dict["defined_state"] = False

        self.evaporator_out_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of evaporator outlet water",
                **tmp_dict,
            )
        )
        # calculate brine volumetric flow rate
        @self.Constraint(doc="brine volumetric flow rate")
        def eq_brine_volumetric_flow_rate(b):
            return b.evaporator_out_props[0].flow_vol_phase["Liq"] == b.feed_props[
                0
            ].flow_vol_phase["Liq"] - pyunits.convert(
                b.permeate_flux * b.module_area, to_units=pyunits.m**3 / pyunits.s
            )

        # calculate brine salinity
        @self.Constraint(doc="brine salinity")
        def eq_brine_salinity(b):
            return (
                b.evaporator_out_props[0].conc_mass_phase_comp["Liq", "TDS"]
                == b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
                * b.feed_props[0].flow_vol_phase["Liq"]
                / b.evaporator_out_props[0].flow_vol_phase["Liq"]
            )

        # make brine pressure at 1 bar
        self.evaporator_out_props[0].pressure.fix(101325)

        """
        Add block for water at the condenser inlet
        """
        tmp_dict["defined_state"] = True

        self.condenser_in_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of condenser inlet water",
                **tmp_dict,
            )
        )

        """
        Add block for water at the condenser outlet
        """
        tmp_dict["defined_state"] = False

        self.condenser_out_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of condenser outlet water",
                **tmp_dict,
            )
        )

        """
        Add block for the inlet cooling water 
        """
        self.cooling_in_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of the inlet cooling water",
            **tmp_dict,
        )

        # Specify the concentration and temperature with any value for initialization,
        # because only volumetric flow rate is known and the temperature will be released and calculated afterwards
        self.cooling_in_props.calculate_state(
            var_args={
                # Cooling water volumetric flow rate, fixed to 1265L/h, to maintain vacuum pressure inside the MD module around 200 mbar
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    1265 * pyunits.L / pyunits.h, to_units=pyunits.m**3 / pyunits.s
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): 50,
                ("temperature", None): 25 + 273.15,
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        self.cooling_in_props[0].temperature.unfix()

        """
        Add block for the outlet cooling water
        """
        self.cooling_out_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of the outlet cooling water",
                **tmp_dict,
            )
        )
        # Specify the concentration and temperature with any value for initialization,
        # because only volumetric flow rate is known and the temperature will be released and calculated afterwards
        self.cooling_out_props.calculate_state(
            var_args={
                # Cooling water volumetric flow rate, fixed to 1265L/h, to maintain vacuum pressure inside the MD module around 200 mbar
                ("flow_vol_phase", "Liq"): pyunits.convert(
                    1265 * pyunits.L / pyunits.h, to_units=pyunits.m**3 / pyunits.s
                ),
                ("conc_mass_phase_comp", ("Liq", "TDS")): 50,
                ("temperature", None): 25 + 273.15,
                ("pressure", None): 101325,
            },
            hold_state=True,
        )
        self.cooling_out_props[0].temperature.unfix()

        """
        Add block for the average status in the cooler
        """
        self.avg_cooling_props = (
            self.config.property_package_seawater.state_block_class(
                self.flowsheet().config.time,
                doc="Average properties in the cooler",
                **tmp_dict,
            )
        )

        """
        Add block for the average status in the heater
        """
        self.avg_feed_props = self.config.property_package_seawater.state_block_class(
            self.flowsheet().config.time,
            doc="Average properties of the feed water",
            **tmp_dict,
        )

        tmp_dict["parameters"] = self.config.property_package_water

        """
        Add block for the average status in the condenser
        """
        self.avg_condenser_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Average properties of water in the condenser",
            **tmp_dict,
        )
        self.avg_condenser_props[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
        self.avg_condenser_props[0].pressure.fix(101325)

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="brine", block=self.evaporator_out_props)

        """
        Constraint equations
        """
        # Cooling system specification
        if self.config.cooling_system_type == CoolingType.closed:

            @self.Constraint(doc="Calculate cooling power requirement")
            def eq_cooling_power_thermal(b):
                return b.cooling_power_thermal == pyunits.convert(
                    b.thermal_resistance_hot
                    * (
                        b.feed_props[0].temperature
                        - b.condenser_in_props[0].temperature
                    ),
                    to_units=pyunits.kW,
                )

            @self.Constraint(
                doc="Calculate inlet cooling water temperature in closed cooling mode"
            )
            def eq_cooling_in_temp(b):
                return (
                    b.cooling_in_props[0].temperature
                    == b.feed_props[0].temperature
                    - pyunits.convert(b.cooling_power_thermal, to_units=pyunits.W)
                    / b.effectiveness_heat_exchanger
                    / b.thermal_resistance_hot
                )

            @self.Constraint(
                doc="Calculate outlet cooling water temperature in closed cooling mode"
            )
            def eq_cooling_out_temp(b):
                return b.cooling_out_props[0].temperature == b.cooling_in_props[
                    0
                ].temperature + b.thermal_resistance_hot / b.thermal_resistance_cold * (
                    b.feed_props[0].temperature - b.condenser_in_props[0].temperature
                )

        else:  # self.config.cooling_system_type == CoolingType.open

            @self.Constraint(doc="Calculate cooling power requirment")
            def eq_cooling_power_thermal(b):
                return b.cooling_power_thermal == pyunits.convert(
                    b.effectiveness_heat_exchanger
                    * b.thermal_resistance_hot
                    * (b.feed_props[0].temperature - b.cooling_in_props[0].temperature),
                    to_units=pyunits.kW,
                )

            @self.Constraint(
                doc="Calculate condenser inlet water temperature in open cooling mode"
            )
            def eq_condenser_in_temp(b):
                return (
                    b.condenser_in_props[0].temperature
                    == b.feed_props[0].temperature
                    - pyunits.convert(b.cooling_power_thermal, to_units=pyunits.W)
                    / b.thermal_resistance_hot
                )

            @self.Constraint(
                doc="Calculate outlet cooling water temperature in open cooling mode"
            )
            def eq_cooling_out_temp(b):
                return (
                    b.cooling_out_props[0].temperature
                    == b.cooling_in_props[0].temperature
                    + pyunits.convert(b.cooling_power_thermal, to_units=pyunits.W)
                    / b.thermal_resistance_cold
                )

        # Set alias for state properties. The acronyms below are:
        #   FFR: Feed flow rate
        #   S  : Feed salinity
        #   TEI: Evaporator inlet temperature
        #   TCI: Condenser inlet temperature
        #   Ttank: Tank temperature
        #   TEO: Evaporator outlet temperature
        #   TCO: Condenser outlet temperature

        FFR = pyunits.convert(
            self.feed_props[0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.h
        )
        S = self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        TEI = self.evaporator_in_props[0].temperature
        TCI = self.condenser_in_props[0].temperature
        Ttank = self.feed_props[0].temperature
        TEO = self.evaporator_out_props[0].temperature
        TCO = self.condenser_out_props[0].temperature

        @self.Constraint(doc="Average temperature of the feed flow in the heat source")
        def eq_avg_temp_heat_tank(b):
            return b.avg_feed_props[0].temperature == (TEI + TCO) / 2

        @self.Constraint(doc="Average temperature of the flow in the condenser")
        def eq_avg_temp_condenser(b):
            return b.avg_condenser_props[0].temperature == (TCI + TCO) / 2

        @self.Constraint(doc="Average temperature in the cooloer")
        def eq_avg_temp_cooling(b):
            return b.avg_cooling_props[0].temperature == (TCI + Ttank) / 2

        @self.Constraint(doc="Average salinity of the feed flow")
        def eq_avg_salinity_feed_tank(b):
            return (
                b.avg_feed_props[0].mass_frac_phase_comp["Liq", "TDS"]
                == b._get_membrane_performance(TEI, FFR, TCI, S)[3] / 1000
            )

        @self.Constraint(doc="Average salinity in the cooler")
        def eq_avg_salinity_cooler(b):
            return b.avg_cooling_props[0].mass_frac_phase_comp["Liq", "TDS"] == 0

        @self.Constraint(doc="Flowrate of the average feed flow block")
        def eq_feed_volumetric_flow_rate(b):
            return (
                b.avg_feed_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="Flowrate in the cooler")
        def eq_cooling_volumetric_flow_rate(b):
            return (
                b.avg_cooling_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="Permeate flux")
        def eq_permeate_flux(b):
            return b.permeate_flux == b._get_membrane_performance(TEI, FFR, TCI, S)[0]

        @self.Constraint(doc="Feed flow pressure drop")
        def eq_feed_flow_pressure_drop(b):
            return b.feed_flow_pressure_drop == b._get_pressure_drop(
                pyunits.convert(
                    b.feed_props[0].flow_vol_phase["Liq"],
                    to_units=pyunits.L / pyunits.h,
                ),
                b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"],
            )

        @self.Constraint(doc="Evaporatore outlet temperature")
        def eq_evaporator_outlet_temp(b):
            return (
                b.evaporator_out_props[0].temperature
                == b._get_membrane_performance(TEI, FFR, TCI, S)[2] + 273.15
            )

        @self.Constraint(doc="Condenser outlet temperature")
        def eq_condenser_outlet_temp(b):
            return (
                b.condenser_out_props[0].temperature
                == b._get_membrane_performance(TEI, FFR, TCI, S)[1] + 273.15
            )

        @self.Constraint(doc="Permeate flow rate")
        def eq_permeate_volumetric_flow_rate(b):
            return b.permeate_props[0].flow_vol_phase["Liq"] == pyunits.convert(
                b.permeate_flux * self.module_area,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Constraint(doc="Condenser mass flow rate")
        def eq_condenser_mass_flow_rate(b):
            return (
                b.avg_condenser_props[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.permeate_props[0].flow_mass_phase_comp["Liq", "H2O"]
            )

        @self.Constraint(doc="initial log mean temperature difference")
        def eq_log_mean_temp_dif(b):
            return b.log_mean_temp_dif == ((TEI - TCO) - (TEO - TCI)) / log(
                (TEI - TCO) / (TEO - TCI + 1e-6)
            )

        @self.Constraint(doc="Thermal power")
        def eq_thermal_power(b):
            CpF = b.avg_feed_props[0].cp_mass_phase["Liq"]
            RhoF = b.avg_feed_props[0].dens_mass_phase["Liq"]
            return b.thermal_power == pyunits.convert(
                (FFR * CpF * (TEI - TCO)) * (RhoF), to_units=pyunits.kW
            )

        @self.Constraint(doc="Electric power for pumping feed water")
        def eq_feed_pump_power_elec(b):
            return b.feed_pump_power_elec == pyunits.convert(
                b.feed_flow_pressure_drop
                / b.pump_efficiency
                * b.feed_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.kW,
            )

        @self.Constraint(doc="Electric power for pumping cooling water")
        def eq_cooling_pump_power_elec(b):
            return b.cooling_pump_power_elec == pyunits.convert(
                b.cooling_flow_pressure_drop
                / b.pump_efficiency
                * b.cooling_in_props[0].flow_vol_phase["Liq"],
                to_units=pyunits.kW,
            )

        @self.Constraint(doc="Calculate hot side thermal resistance")
        def eq_thermal_resistance_hot(b):
            return b.thermal_resistance_hot == pyunits.convert(
                b.feed_props[0].flow_vol_phase["Liq"]
                * b.avg_feed_props[0].dens_mass_phase["Liq"]
                * b.avg_feed_props[0].cp_mass_phase["Liq"],
                to_units=pyunits.W / pyunits.K,
            )

        @self.Constraint(doc="Calculate cold side thermal resistance")
        def eq_thermal_resistance_cold(b):
            return b.thermal_resistance_cold == pyunits.convert(
                b.cooling_in_props[0].flow_vol_phase["Liq"]
                * b.avg_cooling_props[0].dens_mass_phase["Liq"]
                * b.avg_cooling_props[0].cp_mass_phase["Liq"],
                to_units=pyunits.W / pyunits.K,
            )

        @self.Constraint(doc="Calculate the number of transfer units")
        def eq_number_transfer_units(b):
            return (
                b.number_transfer_units
                == b.thermal_heat_transfer_coeff
                * b.heat_exchanger_area
                / b.thermal_resistance_hot
            )

        @self.Constraint(doc="Calculate the effectiveness of the heat exchanger")
        def eq_effectiveness_heat_exchanger(b):
            return b.effectiveness_heat_exchanger == (
                1
                - exp(
                    -(1 - b.thermal_resistance_hot / b.thermal_resistance_cold)
                    * b.number_transfer_units
                )
            ) / (
                1
                - b.thermal_resistance_hot
                / b.thermal_resistance_cold
                * exp(
                    -(1 - b.thermal_resistance_hot / b.thermal_resistance_cold)
                    * b.number_transfer_units
                )
                + 1e-8
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        self.config.property_package_seawater.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        self.config.property_package_seawater.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "TDS")
        )
        self.config.property_package_water.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )

        if iscale.get_scaling_factor(self.permeate_flux) is None:
            iscale.set_scaling_factor(self.permeate_flux, 1)

        if iscale.get_scaling_factor(self.log_mean_temp_dif) is None:
            iscale.set_scaling_factor(self.log_mean_temp_dif, 1e-1)

        if iscale.get_scaling_factor(self.thermal_power) is None:
            iscale.set_scaling_factor(self.thermal_power, 1e-1)

        if iscale.get_scaling_factor(self.cooling_power_thermal) is None:
            iscale.set_scaling_factor(self.cooling_power_thermal, 1e-1)

        if iscale.get_scaling_factor(self.feed_flow_pressure_drop) is None:
            iscale.set_scaling_factor(self.feed_flow_pressure_drop, 1e-1)

        if iscale.get_scaling_factor(self.feed_pump_power_elec) is None:
            iscale.set_scaling_factor(self.feed_pump_power_elec, 1e3)

        if iscale.get_scaling_factor(self.cooling_pump_power_elec) is None:
            iscale.set_scaling_factor(self.cooling_pump_power_elec, 1e3)

        if iscale.get_scaling_factor(self.thermal_resistance_hot) is None:
            iscale.set_scaling_factor(self.thermal_resistance_hot, 1e-3)

        if iscale.get_scaling_factor(self.thermal_resistance_cold) is None:
            iscale.set_scaling_factor(self.thermal_resistance_cold, 1e-3)

        if iscale.get_scaling_factor(self.number_transfer_units) is None:
            iscale.set_scaling_factor(self.number_transfer_units, 1e-1)

        if iscale.get_scaling_factor(self.effectiveness_heat_exchanger) is None:
            iscale.set_scaling_factor(self.effectiveness_heat_exchanger, 1e0)

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")
        init_log.info_low("Starting initialization...")
        # Iniitialize feed properties

        flags = self.feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )

        init_log.info_low("Initialization for feed Completed.")

        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.feed_props[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value
        init_log.info_low("Starting initialization of evaporator_in_props ")
        self.evaporator_in_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of evaporator_out_props ")
        self.evaporator_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_low("Starting initialization of condenser_in_props ")
        self.condenser_in_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of condenser_out_props ")
        self.condenser_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of avg_feed_props ")
        self.avg_feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of avg_cooling_props ")
        self.avg_cooling_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of cooling_in_props ")
        self.cooling_in_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_low("Starting initialization of cooling_out_props ")
        self.cooling_out_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        state_args_water = {}
        state_dict_water = self.avg_condenser_props[
            self.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_water.keys():
            if state_dict_water[k].is_indexed():
                state_args_water[k] = {}
                for m in state_dict_water[k].keys():
                    state_args_water[k][m] = state_dict_water[k][m].value
            else:
                state_args_water[k] = state_dict_water[k].value

        for p, j in self.avg_condenser_props.phase_component_set:
            if p == "Vap" and j == "H2O":
                state_args_water["flow_mass_phase_comp"][(p, j)] = 0

        init_log.info_low("Starting initialization of avg_condenser_props ")
        self.avg_condenser_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_water,
        )

        state_args_permeate = deepcopy(state_args)
        for p, j in self.permeate_props.phase_component_set:
            if j == "TDS":
                state_args_permeate["flow_mass_phase_comp"][(p, j)] = 0

        init_log.info_low("Starting initialization of permeate_props ")
        self.permeate_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,
        )

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"Degrees of freedom of {self.name} is {degrees_of_freedom(self)}"
                f"before the initialization of VAGMD unit model."
            )

        opt = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization status {}.".format(idaeslog.condition(res)))

    def _determine_salinity_mode(
        self,
        feed_flow_rate,
        evap_inlet_temp,
        cond_inlet_temp,
        module_type,
        high_brine_salinity,
        cooling_system_type,
    ):
        # This function rewrites the operation parameters and cooling type
        # when the brine salinity is high (>175.3 g/L) for module AS7C1.5L
        if module_type == ModuleType.AS7C15L and high_brine_salinity:
            feed_flow_rate = 1100  # L/h
            evap_inlet_temp = 80  # deg C
            cond_inlet_temp = 25  # deg C
            self.config.cooling_system_type = CoolingType.closed

        return (feed_flow_rate, evap_inlet_temp, cond_inlet_temp, cooling_system_type)

    """
    Equation to calculate pressure drop
    """

    def _get_pressure_drop(self, flow_rate, salinity):
        if self.config.module_type == ModuleType.AS7C15L:
            coefficients = [
                -158.2007422,
                0.39402609,
                0,
                0.000585345,
                8.93618e-5,
                -0.000287828,
            ]
        else:  # self.config.module_type == ModuleType.AS26C72L
            coefficients = [
                -72.53793298,
                0.110437201,
                0,
                0.000643495,
                0.000189924,
                -0.001111447,
            ]

        return (
            coefficients[0]
            + coefficients[1] * flow_rate
            + coefficients[2] * salinity
            + coefficients[3] * flow_rate * salinity
            + coefficients[4] * flow_rate**2
            + coefficients[5] * salinity**2
        )

    """
    Surrogate equations for membrane performance
    """

    def _get_membrane_performance(self, TEI, FFR, TCI, SgL):
        # Function arguments:
        #   FFR: Feed flow rate
        #   TEI: Evaporator inlet temperature
        #   TCI: Condenser inlet temperature
        #   SgL: Feed salinity

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
        coefficients = [
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

        surrogate_vars = [
            [1, TEI, TEI**2],
            [1, FFR, FFR**2],
            [1, TCI, TCI**2],
            [1, S_c, S_c**2],
        ]

        # Model calculations
        if self.config.module_type == ModuleType.AS7C15L:
            if self.config.high_brine_salinity:
                S_r = S_c

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high
            else:
                TEI = sum(
                    surrogate_vars[0][j] * coefficients[0][j]
                    for j in range(len(coefficients[0]))
                )
                FFR = sum(
                    surrogate_vars[1][j] * coefficients[1][j]
                    for j in range(len(coefficients[1]))
                )
                TCI = sum(
                    surrogate_vars[2][j] * coefficients[2][j]
                    for j in range(len(coefficients[2]))
                )
                S_r = sum(
                    surrogate_vars[3][j] * coefficients[3][j]
                    for j in range(len(coefficients[3]))
                )

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

        else:  # self.config.module_type == ModuleType.AS26C72L
            TEI = sum(
                surrogate_vars[0][j] * coefficients[4][j]
                for j in range(len(coefficients[0]))
            )
            FFR = sum(
                surrogate_vars[1][j] * coefficients[5][j]
                for j in range(len(coefficients[1]))
            )
            TCI = sum(
                surrogate_vars[2][j] * coefficients[6][j]
                for j in range(len(coefficients[2]))
            )
            S_r = sum(
                surrogate_vars[3][j] * coefficients[7][j]
                for j in range(len(coefficients[3]))
            )

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

    @property
    def default_costing_method(self):
        return cost_vagmd_surrogate
