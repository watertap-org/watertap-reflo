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


@declare_process_block_class("VAGMDbatch_surrogate")
class VAGMDbatchData(UnitModelBlockData):
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
        Model parameters
        """
        self.discharge_volume = Param(
            initialize=3.2175, units=pyunits.L, doc="Discharge volume (L)"
        )

        self.min_batch_volume = Param(
            initialize=50, units=pyunits.L, doc="Void volume of an Ags module"
        )

        self.min_salinity = Param(
            initialize=35,
            units=pyunits.g / pyunits.L,
            doc="Minimum feed salinity (g/L)",
        )

        # self.max_salinity = Param(
        #     self.module_type,
        #     initialize={7: 292.2, 26: 245.5},
        #     units=pyunits.g / pyunits.L,
        #     doc="Maximum feed salinity for the selected module (g/L)",
        # )

        self.heat_transfer_coeff = Param(
            initialize=3168,
            units=pyunits.W / pyunits.m**2 / pyunits.C,
            doc="Overall heat transfer coefficient (W/m2/C)",
        )

        self.heat_transfer_area = Param(
            initialize=1.34,
            units=pyunits.m**2,
            doc="Effective heat transfer surface area (m2)",
        )

        self.cooling_flow_rate = Param(
            initialize=1265, units=pyunits.L / pyunits.h, doc="Cooling flow rate (L/h)"
        )
        # Cooling flow rate is fixed to 1265 to maintain vacuum pressure inside the MD module around 200 mbar.

        self.pres_drop_cooling = Param(
            initialize=170, units=pyunits.mbar, doc="Pressure drop of cooling water"
        )

        self.pump_eff = Param(initialize=0.6, doc="Cooling water pump efficiency")

        # TCoolIn is only applicable for open cooling system
        if self.config.cooling_system_type == "open":
            self.TCoolIn = Var(
                initialize=25,
                bounds=(20, 30),
                units=pyunits.K,
                doc="Initial temperature of the cooling water (C)",
            )

        """
        Input variables
        """
        self.initial_batch_volume = Var(
            initialize=50,
            bounds=(self.min_batch_volume, None),
            units=pyunits.L,
            doc="Initial batch volume (L))",
        )

        self.dt = Var(
            initialize=30,
            bounds=(0, None),
            units=pyunits.s,
            doc="Time step of the simulation (s)",
        )

        self.cycles = Set(
            initialize=[i for i in range(self.config.number_cycles)],
            doc="Number of cycles",
        )

        """
        Intermediate variables
        """
        # self.max_recovery = Var(
        #     initialize=1,
        #     bounds=(0, 1),
        #     doc="Maximum recovery rate (%)"
        # )

        self.permeate_flux = Var(
            self.cycles,
            initialize=10,
            units=pyunits.L / pyunits.h / pyunits.m**2,
            doc="Permeate flux (L/h/m2)",
        )

        self.acc_distillate_volume = Var(
            self.cycles,
            initialize=0,
            units=pyunits.L,
            doc="Accumulated volume of distillate(L)",
        )

        self.acc_recovery_ratio = Var(
            self.cycles,
            initialize=0,
            units=pyunits.dimensionless,
            doc="Accumulated recovery rate after each recirculation",
        )

        self.log_mean_temp_dif = Var(
            self.cycles,
            initialize=0,
            units=pyunits.C,
            doc="Log mean temperature difference in the MD module",
        )

        """
        Add variables for output
        """
        self.thermal_power = Var(
            self.cycles,
            initialize=10,
            units=pyunits.kW,
            doc="Thermal power during each timestep (kW)",
        )

        self.thermal_energy = Var(
            self.cycles,
            initialize=0,
            units=pyunits.kW,
            doc="Thermal energy consumption during each timestep (kWh)",
        )

        self.specific_energy_consumption_thermal = Var(
            self.cycles,
            initialize=100,
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific thermal power consumption (kWh/m3)",
        )

        self.gain_output_ratio = Var(
            self.cycles,
            initialize=5,
            units=pyunits.kW / pyunits.kW,
            doc="Gained output ratio",
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
        Add block for the permeated water
        """
        tmp_dict["defined_state"] = False

        self.permeate_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of permeated water",
            **tmp_dict
        )
        # salinity in permeate is zero
        @self.Constraint(self.cycles, doc="permeate salinity")
        def eq_permeate_salinity(b, t):
            return b.permeate_props[0, t].flow_mass_phase_comp["Liq", "TDS"] == 0

        # permeate temperature is assumed the same as feed water, for model initialization purpose
        @self.Constraint(self.cycles, doc="distillate temperature")
        def eq_distillate_temp(b, t):
            return b.permeate_props[0, t].temperature == b.feed_props[0, t].temperature

        """
        Add block for water at the condenser inlet
        """
        tmp_dict["defined_state"] = True

        self.condenser_in_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of condenser inlet water",
            **tmp_dict
        )

        """
        Add block for water at the condenser outlet
        """
        tmp_dict["defined_state"] = False

        self.condenser_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of condenser outlet water",
            **tmp_dict
        )

        """
        Add block for water at the evaporator inlet
        """
        tmp_dict["defined_state"] = True

        self.evaporator_in_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of evaporator inlet water",
            **tmp_dict
        )

        """
        Add block for water at the evaporator outlet
        """
        tmp_dict["defined_state"] = False

        self.evaporator_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        """
        Add block for inlet cooling water
        """
        self.cooling_in_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of inlet cooling water",
            **tmp_dict
        )

        """
        Add block for outlet cooling water
        """
        self.cooling_out_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of outlet cooling water",
            **tmp_dict
        )

        """
        Add block for the average status of the feed flow
        """
        self.avg_feed_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        """
        Add block for the average status of the cooling flow
        """
        self.avg_cool_props = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            self.cycles,
            doc="Material properties of evaporator outlet water",
            **tmp_dict
        )

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="permeate", block=self.permeate_props)

        """
        System configurations
        """
        # TODO: raise error if exceeding maximum recovery rate
        # Maximum recovery rate (calculated based on the max allowed salinity)
        # @self.Constraint(doc="Maximum recovery rate")
        # def max_recovery_cal(b):
        #     return b.max_recovery == (
        #         1
        #         - (
        #             b.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        #             / b.max_salinity[b.config.module_type]
        #         )
        #     )

        """
        Set up for the beginning of the first circulation
        """
        S0 = self.feed_props[0, 0].conc_mass_phase_comp["Liq", "TDS"]
        FFR0 = pyunits.convert(
            self.feed_props[0, 0].flow_vol_phase["Liq"], to_units=pyunits.L / pyunits.h
        )
        Ttank0 = self.feed_props[0, 0].temperature - 273.15 * pyunits.K
        TEI = self.evaporator_in_props[0, 0].temperature
        TCI = self.condenser_in_props[0, 0].temperature

        @self.Constraint(doc="initial permeate flux")
        def eq_permeate_flux_0(b):
            return (
                b.permeate_flux[0]
                == b._get_membrane_performance(TEI, FFR0, TCI, S0, Ttank0)[0]
            )

        @self.Constraint(doc="initial permeate flow rate")
        def eq_permeate_volumetric_flow_rate_0(b):
            return b.permeate_props[0, 0].flow_vol_phase["Liq"] == pyunits.convert(
                b.permeate_flux[0] * self.module_area,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Constraint(doc="initial condenser outlet temperature")
        def eq_condenser_outlet_temp_0(b):
            return (
                b.condenser_out_props[0, 0].temperature
                == b._get_membrane_performance(TEI, FFR0, TCI, S0, Ttank0)[1] + 273.15
            )

        @self.Constraint(doc="initial evaporatore outlet temperature")
        def eq_evaporator_outlet_temp_0(b):
            return (
                b.evaporator_out_props[0, 0].temperature
                == b._get_membrane_performance(TEI, FFR0, TCI, S0, Ttank0)[2] + 273.15
            )

        @self.Constraint(doc="initial accumulated volume of distillate")
        def eq_acc_distillate_volume_0(b):
            return b.acc_distillate_volume[0] == 0

        @self.Constraint(doc="initial recovery rate")
        def eq_recovery_rate_0(b):
            return b.acc_recovery_ratio[0] == 0

        @self.Constraint(doc="initial log mean temperature difference")
        def eq_log_mean_temp_dif_0(b):
            TCO = b.condenser_out_props[0, 0].temperature - 273.15 * pyunits.K
            TEO = b.evaporator_out_props[0, 0].temperature - 273.15 * pyunits.K
            TEI = b.evaporator_in_props[0, 0].temperature - 273.15 * pyunits.K
            TCI = b.condenser_in_props[0, 0].temperature - 273.15 * pyunits.K
            return b.log_mean_temp_dif[0] == ((TEI - TCO) - (TEO - TCI)) / log(
                (TEI - TCO) / (TEO - TCI + 1e-6)
            )

        @self.Constraint(
            doc="initial average temperature of the feed flow in the heat source"
        )
        def eq_avg_temp_heat_tank_0(b):
            return (
                b.avg_feed_props[0, 0].temperature
                == (
                    b.evaporator_in_props[0, 0].temperature
                    + b.condenser_out_props[0, 0].temperature
                )
                / 2
            )

        @self.Constraint(doc="initial average salinity of the feed flow")
        def eq_avg_salinity_feed_tank_0(b):
            return (
                b.avg_feed_props[0, 0].mass_frac_phase_comp["Liq", "TDS"]
                == b._get_membrane_performance(TEI, FFR0, TCI, S0, Ttank0)[3] / 1000
            )

        @self.Constraint(doc="initial flowrate of the feed flow")
        def eq_feed_volumetric_flow_rate_0(b):
            return (
                b.avg_feed_props[0, 0].flow_vol_phase["Liq"]
                == b.feed_props[0, 0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="initial thermal power")
        def eq_thermal_power_0(b):
            CpF = b.avg_feed_props[0, 0].cp_mass_phase["Liq"]
            RhoF = b.avg_feed_props[0, 0].dens_mass_phase["Liq"]
            TCO = b.condenser_out_props[0, 0].temperature - 273.15 * pyunits.K
            TEI = b.evaporator_in_props[0, 0].temperature - 273.15 * pyunits.K
            return b.thermal_power[0] == pyunits.convert(
                (FFR0 * CpF * (TEI - TCO)) * (RhoF), to_units=pyunits.kW
            )

        @self.Constraint(doc="initial thermal energy consumption")
        def eq_thermal_energy_0(b):
            return b.thermal_energy[0] == pyunits.convert(
                b.thermal_power[0] * b.dt, to_units=pyunits.kWh
            )

        @self.Constraint(doc="initial inlet cooling water temperature")
        def eq_cooling_in_temp_0(b):
            return b.cooling_in_props[0, 0].temperature == Ttank0 + 273.15 * pyunits.K

        @self.Constraint(doc="initial outlet cooling water temperature")
        def eq_cooling_out_temp_0(b):
            return b.cooling_out_props[0, 0].temperature == Ttank0 + 273.15 * pyunits.K

        @self.Constraint(doc="initial specific thermal energy consumption")
        def eq_specific_thermal_energy_consumption_0(b):
            return b.specific_energy_consumption_thermal[0] == 0

        @self.Constraint(doc="initial gain output ratio")
        def eq_gain_output_ratio_0(b):
            return b.gain_output_ratio[0] == 0

        """
        Calculation for each circulation (time step)
        """

        @self.Constraint(
            self.cycles,
            doc="Accumulated volume of distillate of each time step",
        )
        def eq_acc_distillate_volume_t(b, t):
            if t > 0:
                return b.acc_distillate_volume[t] == b.acc_distillate_volume[
                    t - 1
                ] + pyunits.convert(
                    b.permeate_props[0, t - 1].flow_vol_phase["Liq"] * b.dt,
                    to_units=pyunits.L,
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Salinity of the batched feed water at each time step",
        )
        def eq_feed_salinity_t(b, t):
            if t > 0:
                return b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"] * (
                    b.initial_batch_volume - b.acc_distillate_volume[t]
                ) == b.feed_props[0, t - 1].conc_mass_phase_comp["Liq", "TDS"] * (
                    b.initial_batch_volume - b.acc_distillate_volume[t - 1]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Feed flow rate remains the same during the recirculations",
        )
        def eq_feed_flow_rate_t(b, t):
            if t > 0:
                return (
                    b.feed_props[0, t].flow_vol_phase["Liq"]
                    == b.feed_props[0, t - 1].flow_vol_phase["Liq"]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Temperature of the batched feed water at each time step",
        )
        def eq_feed_temp_t(b, t):
            if t > 0:
                return b.feed_props[0, t].temperature == (
                    pyunits.convert(
                        b.feed_props[0, t - 1].flow_vol_phase["Liq"] * b.dt,
                        to_units=pyunits.L,
                    )
                    * b.evaporator_out_props[0, t - 1].temperature
                    + (b.initial_batch_volume - b.acc_distillate_volume[t - 1])
                    * b.feed_props[0, t - 1].temperature
                ) / (
                    pyunits.convert(
                        b.feed_props[0, t - 1].flow_vol_phase["Liq"] * b.dt,
                        to_units=pyunits.L,
                    )
                    + b.initial_batch_volume
                    - b.acc_distillate_volume[t - 1]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Accumulated recovery ratio after each circulation",
        )
        def eq_acc_recovery_ratio_t(b, t):
            if t > 0:
                return (
                    b.acc_recovery_ratio[t]
                    == 1
                    - b.feed_props[0, 0].conc_mass_phase_comp["Liq", "TDS"]
                    / b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Permeate flux at each time step",
        )
        def eq_permeate_flux_t(b, t):
            if t > 0:
                return (
                    b.permeate_flux[t]
                    == b._get_membrane_performance(
                        b.evaporator_in_props[0, 0].temperature,
                        pyunits.convert(
                            self.feed_props[0, 0].flow_vol_phase["Liq"],
                            to_units=pyunits.L / pyunits.h,
                        ),
                        b.condenser_in_props[0, 0].temperature,
                        b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"],
                        b.feed_props[0, t].temperature,
                    )[0]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Permeate flow rate at each time step",
        )
        def eq_permeate_volumetric_flow_rate_t(b, t):
            if t > 0:
                return b.permeate_props[0, t].flow_vol_phase["Liq"] == pyunits.convert(
                    b.permeate_flux[t] * self.module_area,
                    to_units=pyunits.m**3 / pyunits.s,
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Condenser outlet temperature at each time step",
        )
        def eq_condenser_outlet_temp_t(b, t):
            if t > 0:
                return (
                    b.condenser_out_props[0, t].temperature
                    == b._get_membrane_performance(
                        b.evaporator_in_props[0, 0].temperature,
                        pyunits.convert(
                            self.feed_props[0, 0].flow_vol_phase["Liq"],
                            to_units=pyunits.L / pyunits.h,
                        ),
                        b.condenser_in_props[0, 0].temperature,
                        b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"],
                        b.feed_props[0, t].temperature,
                    )[1]
                    + 273.15
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Evaporator outlet temperature at each time step",
        )
        def eq_evaporator_outlet_temp_t(b, t):
            if t > 0:
                return (
                    b.evaporator_out_props[0, t].temperature
                    == b._get_membrane_performance(
                        b.evaporator_in_props[0, 0].temperature,
                        pyunits.convert(
                            self.feed_props[0, 0].flow_vol_phase["Liq"],
                            to_units=pyunits.L / pyunits.h,
                        ),
                        b.condenser_in_props[0, 0].temperature,
                        b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"],
                        b.feed_props[0, t].temperature,
                    )[2]
                    + 273.15
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Log mean temperature difference at each time step",
        )
        def eq_log_mean_temp_dif_t(b, t):
            if t > 0:
                TCO = b.condenser_out_props[0, t].temperature - 273.15 * pyunits.K
                TEO = b.evaporator_out_props[0, t].temperature - 273.15 * pyunits.K
                TEI = b.evaporator_in_props[0, 0].temperature - 273.15 * pyunits.K
                TCI = b.condenser_in_props[0, 0].temperature - 273.15 * pyunits.K
                return b.log_mean_temp_dif[t] == ((TEI - TCO) - (TEO - TCI)) / log(
                    (TEI - TCO) / (TEO - TCI + 1e-6)
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Average temperature of the feed flow in the heat source at each time step",
        )
        def eq_avg_temp_heat_tank_t(b, t):
            if t > 0:
                return (
                    b.avg_feed_props[0, t].temperature
                    == (
                        b.evaporator_in_props[0, 0].temperature
                        + b.condenser_out_props[0, t].temperature
                    )
                    / 2
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Salinity of the feed flow in the heat source at each time step",
        )
        def eq_avg_salinity_heat_tank_t(b, t):
            if t > 0:
                return (
                    b.avg_feed_props[0, t].mass_frac_phase_comp["Liq", "TDS"]
                    == b._get_membrane_performance(
                        b.evaporator_in_props[0, 0].temperature,
                        pyunits.convert(
                            self.feed_props[0, 0].flow_vol_phase["Liq"],
                            to_units=pyunits.L / pyunits.h,
                        ),
                        b.condenser_in_props[0, 0].temperature,
                        b.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"],
                        b.feed_props[0, t].temperature,
                    )[3]
                    / 1000
                )

            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Feed flow rate in the heat source at each time step",
        )
        def eq_avg_flow_rate_heat_tank_t(b, t):
            if t > 0:
                return (
                    b.avg_feed_props[0, t].flow_vol_phase["Liq"]
                    == b.feed_props[0, 0].flow_vol_phase["Liq"]
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Thermal power at each time step",
        )
        def eq_thermal_power_t(b, t):
            if t > 0:
                FFR = b.feed_props[0, t].flow_vol_phase["Liq"]
                CpF = b.avg_feed_props[0, t].cp_mass_phase["Liq"]
                RhoF = b.avg_feed_props[0, t].dens_mass_phase["Liq"]
                TCO = b.condenser_out_props[0, t].temperature - 273.15 * pyunits.K
                TEI = b.evaporator_in_props[0, 0].temperature - 273.15 * pyunits.K
                return b.thermal_power[t] == pyunits.convert(
                    (FFR * CpF * (TEI - TCO)) * (RhoF), to_units=pyunits.kW
                )
            else:
                return Constraint.Skip

        @self.Constraint(
            self.cycles,
            doc="Thermal energy at each time step",
        )
        def eq_thermal_energy_t(b, t):
            if t > 0:
                return b.thermal_energy[t] == pyunits.convert(
                    b.thermal_power[t] * b.dt, to_units=pyunits.kWh
                )
            else:
                return Constraint.Skip

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if self.config.cooling_system_type == "open":
            if iscale.get_scaling_factor(self.TCoolIn) is None:
                iscale.set_scaling_factor(self.TCoolIn, 1e-1)

        if iscale.get_scaling_factor(self.initial_batch_volume) is None:
            iscale.set_scaling_factor(self.initial_batch_volume, 1)

        if iscale.get_scaling_factor(self.dt) is None:
            iscale.set_scaling_factor(self.dt, 1e-1)

        for t in self.cycles:
            if iscale.get_scaling_factor(self.permeate_flux) is None:
                iscale.set_scaling_factor(self.permeate_flux, 1)

            if iscale.get_scaling_factor(self.acc_distillate_volume) is None:
                iscale.set_scaling_factor(self.acc_distillate_volume, 1e-1)

            if iscale.get_scaling_factor(self.acc_recovery_ratio) is None:
                iscale.set_scaling_factor(self.acc_recovery_ratio, 10)

            if iscale.get_scaling_factor(self.log_mean_temp_dif) is None:
                iscale.set_scaling_factor(self.log_mean_temp_dif, 1e-1)

            if iscale.get_scaling_factor(self.thermal_power) is None:
                iscale.set_scaling_factor(self.thermal_power, 1)

            if iscale.get_scaling_factor(self.thermal_energy) is None:
                iscale.set_scaling_factor(self.thermal_energy, 1e2)

            if (
                iscale.get_scaling_factor(self.specific_energy_consumption_thermal)
                is None
            ):
                iscale.set_scaling_factor(
                    self.specific_energy_consumption_thermal, 1e-2
                )

            if iscale.get_scaling_factor(self.gain_output_ratio) is None:
                iscale.set_scaling_factor(self.gain_output_ratio, 1e-1)

        # Transforming constraints
        sf = iscale.get_scaling_factor(self.permeate_flux[0])
        iscale.constraint_scaling_transform(self.eq_permeate_flux_0, sf)

        sf = iscale.get_scaling_factor(self.permeate_props[0, 0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_permeate_volumetric_flow_rate_0, sf)

        sf = iscale.get_scaling_factor(self.condenser_out_props[0, 0].temperature)
        iscale.constraint_scaling_transform(self.eq_condenser_outlet_temp_0, sf)

        sf = iscale.get_scaling_factor(self.evaporator_out_props[0, 0].temperature)
        iscale.constraint_scaling_transform(self.eq_evaporator_outlet_temp_0, sf)

        sf = iscale.get_scaling_factor(self.acc_distillate_volume[0])
        iscale.constraint_scaling_transform(self.eq_acc_distillate_volume_0, sf)

        sf = iscale.get_scaling_factor(self.acc_recovery_ratio[0])
        iscale.constraint_scaling_transform(self.eq_recovery_rate_0, sf)

        sf = iscale.get_scaling_factor(self.log_mean_temp_dif[0])
        iscale.constraint_scaling_transform(self.eq_log_mean_temp_dif_0, sf)

        sf = iscale.get_scaling_factor(self.avg_feed_props[0, 0].temperature)
        iscale.constraint_scaling_transform(self.eq_avg_temp_heat_tank_0, sf)

        sf = iscale.get_scaling_factor(
            self.avg_feed_props[0, 0].mass_frac_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.eq_avg_salinity_feed_tank_0, sf)

        sf = iscale.get_scaling_factor(self.avg_feed_props[0, 0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_feed_volumetric_flow_rate_0, sf)

        sf = iscale.get_scaling_factor(self.thermal_power[0])
        iscale.constraint_scaling_transform(self.eq_thermal_power_0, sf)

        sf = iscale.get_scaling_factor(self.thermal_energy[0])
        iscale.constraint_scaling_transform(self.eq_thermal_energy_0, sf)

        sf = iscale.get_scaling_factor(self.cooling_in_props[0, 0].temperature)
        iscale.constraint_scaling_transform(self.eq_cooling_in_temp_0, sf)

        sf = iscale.get_scaling_factor(self.cooling_out_props[0, 0].temperature)
        iscale.constraint_scaling_transform(self.eq_cooling_out_temp_0, sf)

        sf = iscale.get_scaling_factor(self.specific_energy_consumption_thermal[0])
        iscale.constraint_scaling_transform(
            self.eq_specific_thermal_energy_consumption_0, sf
        )

        sf = iscale.get_scaling_factor(self.gain_output_ratio[0])
        iscale.constraint_scaling_transform(self.eq_gain_output_ratio_0, sf)

        for t in self.cycles:
            sf = iscale.get_scaling_factor(
                self.permeate_props[0, t].flow_mass_phase_comp["Liq", "TDS"]
            )
            iscale.constraint_scaling_transform(self.eq_permeate_salinity[t], sf)

            sf = iscale.get_scaling_factor(self.permeate_props[0, t].temperature)
            iscale.constraint_scaling_transform(self.eq_distillate_temp[t], sf)

            if t > 0:
                sf = iscale.get_scaling_factor(self.acc_distillate_volume[t])
                iscale.constraint_scaling_transform(
                    self.eq_acc_distillate_volume_t[t], sf
                )

                sf1 = iscale.get_scaling_factor(
                    self.feed_props[0, t].conc_mass_phase_comp["Liq", "TDS"]
                )
                sf2 = iscale.get_scaling_factor(self.initial_batch_volume)
                sf3 = iscale.get_scaling_factor(self.acc_distillate_volume[t])
                iscale.constraint_scaling_transform(
                    self.eq_feed_salinity_t[t], sf1 * (sf2 - sf3)
                )

                sf = iscale.get_scaling_factor(
                    self.feed_props[0, t].flow_vol_phase["Liq"]
                )
                iscale.constraint_scaling_transform(self.eq_feed_flow_rate_t[t], sf)

                sf = iscale.get_scaling_factor(self.feed_props[0, t].temperature)
                iscale.constraint_scaling_transform(self.eq_feed_temp_t[t], sf)

                sf = iscale.get_scaling_factor(self.acc_recovery_ratio[t])
                iscale.constraint_scaling_transform(self.eq_acc_recovery_ratio_t[t], sf)

                sf = iscale.get_scaling_factor(self.permeate_flux[t])
                iscale.constraint_scaling_transform(self.eq_permeate_flux_t[t], sf)

                sf = iscale.get_scaling_factor(
                    self.permeate_props[0, t].flow_vol_phase["Liq"]
                )
                iscale.constraint_scaling_transform(
                    self.eq_permeate_volumetric_flow_rate_t[t], sf
                )

                sf = iscale.get_scaling_factor(
                    self.condenser_out_props[0, t].temperature
                )
                iscale.constraint_scaling_transform(
                    self.eq_condenser_outlet_temp_t[t], sf
                )

                sf = iscale.get_scaling_factor(
                    self.evaporator_out_props[0, t].temperature
                )
                iscale.constraint_scaling_transform(
                    self.eq_evaporator_outlet_temp_t[t], sf
                )

                sf = iscale.get_scaling_factor(self.log_mean_temp_dif[t])
                iscale.constraint_scaling_transform(self.eq_log_mean_temp_dif_t[t], sf)

                sf = iscale.get_scaling_factor(self.avg_feed_props[0, t].temperature)
                iscale.constraint_scaling_transform(self.eq_avg_temp_heat_tank_t[t], sf)

                sf = iscale.get_scaling_factor(
                    self.avg_feed_props[0, t].mass_frac_phase_comp["Liq", "TDS"]
                )
                iscale.constraint_scaling_transform(
                    self.eq_avg_salinity_heat_tank_t[t], sf
                )

                sf = iscale.get_scaling_factor(
                    self.avg_feed_props[0, t].flow_vol_phase["Liq"]
                )
                iscale.constraint_scaling_transform(
                    self.eq_avg_flow_rate_heat_tank_t[t], sf
                )

                sf = iscale.get_scaling_factor(self.thermal_power[t])
                iscale.constraint_scaling_transform(self.eq_thermal_power_t[t], sf)

                sf = iscale.get_scaling_factor(self.thermal_energy[t])
                iscale.constraint_scaling_transform(self.eq_thermal_energy_t[t], sf)

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
