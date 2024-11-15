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

from pyomo.common.config import ConfigValue
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.util.scaling import get_scaling_factor, set_scaling_factor
from idaes.core.util.misc import add_object_reference

from watertap.costing.watertap_costing_package import (
    WaterTAPCostingData,
    WaterTAPCostingBlockData,
)
from watertap.costing.zero_order_costing import _load_case_study_definition

from watertap_contrib.reflo.core import PySAMWaterTAP
from watertap_contrib.reflo.solar_models.surrogate.pv.pv_surrogate import (
    PVSurrogateData,
)
from watertap_contrib.reflo.costing.tests.dummy_costing_units import (
    DummyElectricityUnit,
)


@declare_process_block_class("REFLOCosting")
class REFLOCostingData(WaterTAPCostingData):

    CONFIG = WaterTAPCostingData.CONFIG()
    CONFIG.declare(
        "case_study_definition",
        ConfigValue(
            default=None,
            doc="Path to YAML file defining global parameters for case study. If "
            "not provided, WaterTAP-REFLO values are used.",
        ),
    )

    def build_global_params(self):

        super().build_global_params()

        # Override WaterTAP default value of USD_2018
        self.base_currency = pyo.units.USD_2021

        self.sales_tax_frac = pyo.Param(
            initialize=0.05,
            mutable=True,
            doc="Sales tax as fraction of capital costs",
            units=pyo.units.dimensionless,
        )

        self.heat_cost = pyo.Var(
            initialize=0.0,
            doc="Heat cost",
            units=self.base_currency / pyo.units.kWh,
        )

        self.defined_flows["heat"] = self.heat_cost

        self.heat_cost.fix(0.0)
        self.electricity_cost.fix(0.0)
        self.plant_lifetime.fix(20)
        self.utilization_factor.fix(1)

        # This should override default values
        if self.config.case_study_definition is not None:
            self.case_study_def = _load_case_study_definition(self)
            # Register currency and conversion rates
            if "currency_definitions" in self.case_study_def:
                pyo.units.load_definitions_from_strings(
                    self.case_study_def["currency_definitions"]
                )
            # If currency definition is defined in case study yaml,
            # we should be able to set it here.
            if "base_currency" in self.case_study_def:
                self.base_currency = getattr(
                    pyo.units, self.case_study_def["base_currency"]
                )
            if "base_period" in self.case_study_def:
                self.base_period = getattr(
                    pyo.units, self.case_study_def["base_period"]
                )
            # Define expected flows
            for f, v in self.case_study_def["defined_flows"].items():
                value = v["value"]
                units = getattr(pyo.units, v["units"])
                if self.component(f + "_cost") is not None:
                    self.component(f + "_cost").fix(value * units)
                else:
                    self.defined_flows[f] = value * units


@declare_process_block_class("TreatmentCosting")
class TreatmentCostingData(REFLOCostingData):
    def build_global_params(self):
        super().build_global_params()

    def build_process_costs(self):
        super().build_process_costs()

    def add_specific_electric_energy_consumption(
        self, flow_rate, name="specific_electric_energy_consumption"
    ):
        """
        Add specific electric energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific electric energy consumption
        """

        specific_electric_energy_consumption = pyo.Var(
            initialize=1,
            units=pyo.units.kilowatt * pyo.units.hr * pyo.units.m**-3,
            doc=f"Specific electric energy consumption based on flow {flow_rate.name}",
        )

        self.add_component(name, specific_electric_energy_consumption)

        specific_electric_energy_consumption_constraint = pyo.Constraint(
            expr=specific_electric_energy_consumption
            == pyo.units.convert(
                self.aggregate_flow_electricity / flow_rate,
                to_units=pyo.units.kilowatt * pyo.units.hr * pyo.units.m**-3,
            )
        )

        self.add_component(
            "specific_electric_energy_consumption_constraint",
            specific_electric_energy_consumption_constraint,
        )

    def add_specific_thermal_energy_consumption(
        self, flow_rate, name="specific_thermal_energy_consumption"
    ):
        """
        Add specific thermal energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific thermal energy consumption
        """

        specific_thermal_energy_consumption = pyo.Var(
            initialize=1,
            units=pyo.units.kilowatt * pyo.units.hr * pyo.units.m**-3,
            doc=f"Specific thermal energy consumption based on flow {flow_rate.name}",
        )

        self.add_component(name, specific_thermal_energy_consumption)

        specific_thermal_energy_consumption_constraint = pyo.Constraint(
            expr=specific_thermal_energy_consumption
            == pyo.units.convert(
                self.aggregate_flow_heat / flow_rate,
                to_units=pyo.units.kilowatt * pyo.units.hr * pyo.units.m**-3,
            )
        )

        self.add_component(
            "specific_thermal_energy_consumption_constraint",
            specific_thermal_energy_consumption_constraint,
        )


@declare_process_block_class("EnergyCosting")
class EnergyCostingData(REFLOCostingData):
    def build_global_params(self):
        super().build_global_params()

        # If creating an energy unit that generates electricity,
        # set this flag to True in costing package.
        # See PV costing package for example.
        self.has_electricity_generation = False

        self.base_energy_units = pyo.units.kilowatt * pyo.units.hour

        self.plant_lifetime_set = pyo.Set(
            initialize=range(pyo.value(self.plant_lifetime) + 1)
        )

        self.annual_electrical_system_degradation = pyo.Param(
            initialize=0.005,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Yearly performance degradation of electric energy system",
        )

        self.annual_heat_system_degradation = pyo.Param(
            initialize=0.005,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Yearly performance degradation of electric energy system",
        )

        self.yearly_electricity_production = pyo.Var(
            self.plant_lifetime_set,
            initialize=1e4,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kilowatt * pyo.units.hour,
            doc="Yearly electricity production over facility lifetime",
        )

        self.lifetime_electricity_production = pyo.Var(
            initialize=1e6,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kilowatt * pyo.units.hour,
            doc="Total electricity production over facility lifetime",
        )

        self.yearly_heat_production = pyo.Var(
            self.plant_lifetime_set,
            initialize=1e4,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kilowatt * pyo.units.hour,
            doc="Yearly heat production over facility lifetime",
        )

        self.lifetime_heat_production = pyo.Var(
            initialize=1e6,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kilowatt * pyo.units.hour,
            doc="Total heat production over facility lifetime",
        )

    def build_process_costs(self):
        super().build_process_costs()

    def build_LCOE_params(self):

        def rule_yearly_electricity_production(b, y):
            if y == 0:
                return b.yearly_electricity_production[y] == pyo.units.convert(
                    self.aggregate_flow_electricity * -1 * pyo.units.year,
                    to_units=pyo.units.kilowatt * pyo.units.hour,
                )
            else:
                return b.yearly_electricity_production[
                    y
                ] == b.yearly_electricity_production[y - 1] * (
                    1 - b.annual_electrical_system_degradation
                )

        self.yearly_electricity_production_constraint = pyo.Constraint(
            self.plant_lifetime_set, rule=rule_yearly_electricity_production
        )

        def rule_lifetime_electricity_production(b):
            return (
                b.lifetime_electricity_production
                == sum(b.yearly_electricity_production[y] for y in b.plant_lifetime_set)
                * b.utilization_factor
            )

        self.lifetime_electricity_production_constraint = pyo.Constraint(
            rule=rule_lifetime_electricity_production
        )

        if get_scaling_factor(self.yearly_electricity_production) is None:
            set_scaling_factor(self.yearly_electricity_production, 1e-4)

        if get_scaling_factor(self.lifetime_electricity_production) is None:
            set_scaling_factor(self.lifetime_electricity_production, 1e-4)

    def build_LCOH_params(self):

        def rule_yearly_heat_production(b, y):
            if y == 0:
                return b.yearly_heat_production[y] == pyo.units.convert(
                    self.aggregate_flow_heat * -1 * pyo.units.year,
                    to_units=pyo.units.kilowatt * pyo.units.hour,
                )
            else:
                return b.yearly_heat_production[y] == b.yearly_heat_production[
                    y - 1
                ] * (1 - b.annual_heat_system_degradation)

        self.yearly_heat_production_constraint = pyo.Constraint(
            self.plant_lifetime_set, rule=rule_yearly_heat_production
        )

        def rule_lifetime_heat_production(b):
            return (
                b.lifetime_heat_production
                == sum(b.yearly_heat_production[y] for y in b.plant_lifetime_set)
                * b.utilization_factor
            )

        self.lifetime_heat_production_constraint = pyo.Constraint(
            rule=rule_lifetime_heat_production
        )

        if get_scaling_factor(self.yearly_heat_production) is None:
            set_scaling_factor(self.yearly_heat_production, 1e-4)

        if get_scaling_factor(self.lifetime_heat_production) is None:
            set_scaling_factor(self.lifetime_heat_production, 1e-4)

    def add_LCOE(self):
        """
        Add Levelized Cost of Energy (LCOE) to costing block.
        """

        # https://www.nrel.gov/analysis/tech-lcoe-documentation.html

        self.build_LCOE_params()

        numerator = pyo.units.convert(
            (
                self.total_capital_cost * self.capital_recovery_factor
                + self.aggregate_fixed_operating_cost
            )
            * self.plant_lifetime,
            to_units=self.base_currency,
        )

        LCOE_expr = pyo.Expression(
            expr=pyo.units.convert(
                numerator / self.lifetime_electricity_production,
                to_units=self.base_currency / self.base_energy_units,
            )
        )

        self.add_component("LCOE", LCOE_expr)

    def add_LCOH(self):
        """
        Add Levelized Cost of Heat (LCOH) to costing block.
        """

        # https://www.nrel.gov/analysis/tech-lcoe-documentation.html

        self.build_LCOH_params()

        numerator = pyo.units.convert(
            (
                self.total_capital_cost * self.capital_recovery_factor
                + self.aggregate_fixed_operating_cost
            )
            * self.plant_lifetime,
            to_units=self.base_currency,
        )

        LCOH_expr = pyo.Expression(
            expr=pyo.units.convert(
                numerator / self.lifetime_heat_production,
                to_units=self.base_currency / self.base_energy_units,
            )
        )

        self.add_component("LCOH", LCOH_expr)

    def add_LCOW(self, *args, **kwargs):

        raise ValueError("Can't add LCOW to EnergyCosting package.")


@declare_process_block_class("REFLOSystemCosting")
class REFLOSystemCostingData(WaterTAPCostingBlockData):

    def build_global_params(self):
        super().build_global_params()

        self.base_currency = pyo.units.USD_2021

        # Fix the parameters
        self.electricity_cost.fix(0.0)
        self.plant_lifetime.fix(20)
        self.utilization_factor.fix(1)

        self.electricity_cost_buy = pyo.Param(
            mutable=True,
            initialize=0.07,
            doc="Electricity cost to buy",
            units=self.base_currency / pyo.units.kWh,
        )

        self.electricity_cost_sell = pyo.Param(
            mutable=True,
            initialize=0.05,
            doc="Electricity cost to sell",
            units=self.base_currency / pyo.units.kWh,
        )

        self.heat_cost_buy = pyo.Param(
            mutable=True,
            initialize=0.01,
            doc="Heat cost to buy",
            units=self.base_currency / pyo.units.kWh,
        )

        self.heat_cost_sell = pyo.Param(
            mutable=True,
            initialize=0.01,
            doc="Heat cost to sell",
            units=self.base_currency / pyo.units.kWh,
        )

        # Build the integrated system costs
        self.build_integrated_costs()

    def build_integrated_costs(self):

        treat_cost = self._get_treatment_cost_block()
        energy_cost = self._get_energy_cost_block()

        # Check if all parameters are equivalent
        self._check_common_param_equivalence(treat_cost, energy_cost)

        # Add all treatment and energy units to _registered_unit_costing
        # so aggregated costs can be calculated at system level.
        for b in [treat_cost, energy_cost]:
            for u in b._registered_unit_costing:
                self._registered_unit_costing.append(u)

        self.total_capital_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Total capital cost for integrated system",
            units=self.base_currency,
        )

        self.total_operating_cost = pyo.Var(
            initialize=1e3,
            doc="Total operating cost for integrated system",
            units=self.base_currency / self.base_period,
        )

        self.aggregate_flow_electricity = pyo.Var(
            initialize=1e3,
            doc="Aggregated system electricity flow",
            units=pyo.units.kW,
        )

        self.aggregate_flow_heat = pyo.Var(
            initialize=1e3,
            doc="Aggregated system heat flow",
            units=pyo.units.kW,
        )

        self.total_electric_operating_cost = pyo.Var(
            initialize=1e3,
            doc="Total electricity related operating cost",
            units=self.base_currency / self.base_period,
        )

        self.total_heat_operating_cost = pyo.Var(
            initialize=1e3,
            doc="Total heat related operating cost",
            units=self.base_currency / self.base_period,
        )

        self.frac_elec_from_grid = pyo.Var(
            initialize=0.1,
            domain=pyo.NonNegativeReals,
            bounds=(0, 1.00001),
            doc="Fraction of electricity from grid",
            units=pyo.units.dimensionless,
        )

        self.aggregate_flow_electricity_purchased = pyo.Var(
            initialize=100,
            domain=pyo.NonNegativeReals,
            doc="Aggregated electricity consumed",
            units=pyo.units.kW,
        )

        self.aggregate_flow_electricity_sold = pyo.Var(
            initialize=0,
            domain=pyo.NonNegativeReals,
            doc="Aggregated electricity produced",
            units=pyo.units.kW,
        )

        self.aggregate_flow_heat_purchased = pyo.Var(
            initialize=100,
            domain=pyo.NonNegativeReals,
            doc="Aggregated heat consumed",
            units=pyo.units.kW,
        )

        self.aggregate_flow_heat_sold = pyo.Var(
            initialize=0,
            domain=pyo.NonNegativeReals,
            doc="Aggregated heat produced",
            units=pyo.units.kW,
        )

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == pyo.units.convert(
                treat_cost.total_capital_cost + energy_cost.total_capital_cost,
                to_units=self.base_currency,
            )
        )

        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == pyo.units.convert(
                treat_cost.total_operating_cost
                + energy_cost.total_operating_cost
                + self.total_electric_operating_cost
                + self.total_heat_operating_cost,
                to_units=self.base_currency / self.base_period,
            )
        )
        # energy producer's electricity flow is negative
        self.aggregate_electricity_balance = pyo.Constraint(
            expr=(
                self.aggregate_flow_electricity_purchased
                + -1 * energy_cost.aggregate_flow_electricity
                == treat_cost.aggregate_flow_electricity
                + self.aggregate_flow_electricity_sold
            )
        )

        # Calculate fraction of electricity from grid when an electricity generating unit is present
        if energy_cost.has_electricity_generation:
            elec_gen_unit = self._get_electricity_generation_unit()
            self.frac_elec_from_grid_constraint = pyo.Constraint(
                expr=(
                    self.frac_elec_from_grid
                    == 1
                    - (
                        elec_gen_unit.electricity
                        / (
                            elec_gen_unit.electricity
                            + self.aggregate_flow_electricity_purchased
                        )
                    )
                )
            )

        else:
            self.frac_elec_from_grid.fix(1)

        self.aggregate_electricity_complement = pyo.Constraint(
            expr=self.aggregate_flow_electricity_purchased
            * self.aggregate_flow_electricity_sold
            == 0
        )

        if all(hasattr(b, "aggregate_flow_heat") for b in [treat_cost, energy_cost]):

            # treatment block is consuming heat and energy block is generating it
            self.has_heat_flows = True
            self.frac_heat_from_grid = pyo.Var(
                initialize=0,
                domain=pyo.NonNegativeReals,
                bounds=(0, 1.00001),
                doc="Fraction of heat from grid",
                units=pyo.units.dimensionless,
            )

            # energy producer's heat flow is negative
            self.aggregate_heat_balance = pyo.Constraint(
                expr=(
                    self.aggregate_flow_heat_purchased
                    + -1 * energy_cost.aggregate_flow_heat
                    == treat_cost.aggregate_flow_heat + self.aggregate_flow_heat_sold
                )
            )

            self.frac_heat_from_grid_constraint = pyo.Constraint(
                expr=(
                    self.frac_heat_from_grid
                    == 1
                    - (
                        -1
                        * energy_cost.aggregate_flow_heat
                        / treat_cost.aggregate_flow_heat
                    )
                )
            )

            self.aggregate_heat_complement = pyo.Constraint(
                expr=self.aggregate_flow_heat_purchased * self.aggregate_flow_heat_sold
                == 0
            )

        elif hasattr(treat_cost, "aggregate_flow_heat"):

            # treatment block is consuming heat but energy block isn't generating
            # we still want to cost the heat consumption

            self.has_heat_flows = True
            self.aggregate_flow_heat_sold.fix(0)

            self.aggregate_heat_balance = pyo.Constraint(
                expr=(
                    self.aggregate_flow_heat_purchased == treat_cost.aggregate_flow_heat
                )
            )

        else:
            # treatment block isn't consuming heat and energy block isn't generating heat
            self.has_heat_flows = False
            self.aggregate_flow_heat_purchased.fix(0)
            self.aggregate_flow_heat_sold.fix(0)

        # positive is for cost and negative for revenue
        self.total_electric_operating_cost_constraint = pyo.Constraint(
            expr=self.total_electric_operating_cost
            == pyo.units.convert(
                (
                    pyo.units.convert(
                        self.aggregate_flow_electricity_purchased,
                        to_units=pyo.units.kWh / pyo.units.year,
                    )
                    * self.electricity_cost_buy
                    - pyo.units.convert(
                        self.aggregate_flow_electricity_sold,
                        to_units=pyo.units.kWh / pyo.units.year,
                    )
                    * self.electricity_cost_sell
                ),
                to_units=self.base_currency / self.base_period,
            )
        )

        # positive is for cost and negative for revenue
        self.total_heat_operating_cost_constraint = pyo.Constraint(
            expr=self.total_heat_operating_cost
            == pyo.units.convert(
                (
                    pyo.units.convert(
                        self.aggregate_flow_heat_purchased,
                        to_units=pyo.units.kWh / pyo.units.year,
                    )
                    * self.heat_cost_buy
                    - pyo.units.convert(
                        self.aggregate_flow_heat_sold,
                        to_units=pyo.units.kWh / pyo.units.year,
                    )
                    * self.heat_cost_sell
                ),
                to_units=self.base_currency / self.base_period,
            )
        )

        # positive is for consumption
        self.aggregate_flow_electricity_constraint = pyo.Constraint(
            expr=self.aggregate_flow_electricity
            == self.aggregate_flow_electricity_purchased
            - self.aggregate_flow_electricity_sold
        )

        self.aggregate_flow_heat_constraint = pyo.Constraint(
            expr=self.aggregate_flow_heat
            == self.aggregate_flow_heat_purchased - self.aggregate_flow_heat_sold
        )

    def initialize_build(self):

        energy_cost = self._get_energy_cost_block()

        self.aggregate_flow_electricity_sold.fix(0)
        self.aggregate_electricity_complement.deactivate()

        calculate_variable_from_constraint(
            self.aggregate_flow_electricity_purchased,
            self.aggregate_electricity_balance,
        )

        # Commented code remains as a PSA:
        # If you send a fixed variable to calculate_variable_from_constraint,
        # the variable will come out a different value but still be fixed!

        # if hasattr(self, "frac_elec_from_grid_constraint"):
        #     calculate_variable_from_constraint(
        #         self.frac_elec_from_grid, self.frac_elec_from_grid_constraint
        #     )

        calculate_variable_from_constraint(
            self.total_electric_operating_cost,
            self.total_electric_operating_cost_constraint,
        )

        calculate_variable_from_constraint(
            self.aggregate_flow_electricity,
            self.aggregate_flow_electricity_constraint,
        )

        self.aggregate_flow_electricity_sold.unfix()
        self.aggregate_electricity_complement.activate()

        if not self.has_heat_flows:
            self.total_heat_operating_cost.fix(0)
            self.total_heat_operating_cost_constraint.deactivate()
            self.aggregate_flow_heat.fix(0)
            self.aggregate_flow_heat_constraint.deactivate()

        else:
            if hasattr(self, "aggregate_heat_complement"):

                self.aggregate_flow_heat_sold.fix(0)
                self.aggregate_heat_complement.deactivate()

                if not self.aggregate_flow_heat_purchased.is_fixed():
                    calculate_variable_from_constraint(
                        self.aggregate_flow_heat_purchased,
                        self.aggregate_heat_balance,
                    )

            calculate_variable_from_constraint(
                self.total_heat_operating_cost,
                self.total_heat_operating_cost_constraint,
            )
            calculate_variable_from_constraint(
                self.aggregate_flow_heat,
                self.aggregate_flow_heat_constraint,
            )

            if hasattr(self, "aggregate_heat_complement"):

                self.aggregate_flow_heat_sold.unfix()
                self.aggregate_heat_complement.activate()

        super().initialize_build()

        if hasattr(self, "LCOT"):
            calculate_variable_from_constraint(
                self.LCOT,
                self.LCOT_constraint,
            )

        if hasattr(self, "LCOE"):
            for y, c in energy_cost.yearly_electricity_production_constraint.items():

                calculate_variable_from_constraint(
                    energy_cost.yearly_electricity_production[y],
                    c,
                )

            calculate_variable_from_constraint(
                energy_cost.lifetime_electricity_production,
                energy_cost.lifetime_electricity_production_constraint,
            )
            calculate_variable_from_constraint(
                energy_cost.LCOE,
                energy_cost.LCOE_constraint,
            )

        if hasattr(self, "LCOH"):

            for y, c in energy_cost.yearly_heat_production_constraint.items():

                calculate_variable_from_constraint(
                    energy_cost.yearly_heat_production[y],
                    c,
                )

            calculate_variable_from_constraint(
                energy_cost.lifetime_heat_production,
                energy_cost.lifetime_heat_production_constraint,
            )
            calculate_variable_from_constraint(
                energy_cost.LCOH,
                energy_cost.LCOH_constraint,
            )

    def calculate_scaling_factors(self):

        if get_scaling_factor(self.total_capital_cost) is None:
            set_scaling_factor(self.total_capital_cost, 1e-3)

        if get_scaling_factor(self.total_operating_cost) is None:
            set_scaling_factor(self.total_operating_cost, 1e-3)

        if get_scaling_factor(self.total_electric_operating_cost) is None:
            set_scaling_factor(self.total_electric_operating_cost, 1e-2)

        if get_scaling_factor(self.total_heat_operating_cost) is None:
            set_scaling_factor(self.total_heat_operating_cost, 1)

        if get_scaling_factor(self.aggregate_flow_electricity) is None:
            set_scaling_factor(self.aggregate_flow_electricity, 0.1)

        if get_scaling_factor(self.aggregate_flow_heat) is None:
            set_scaling_factor(self.aggregate_flow_heat, 0.1)

        if get_scaling_factor(self.aggregate_flow_electricity_purchased) is None:
            sf = get_scaling_factor(self.aggregate_flow_electricity)
            set_scaling_factor(self.aggregate_flow_electricity_purchased, sf)

        if get_scaling_factor(self.aggregate_flow_electricity_sold) is None:
            set_scaling_factor(self.aggregate_flow_electricity_sold, 1)

        if get_scaling_factor(self.aggregate_flow_heat_purchased) is None:
            sf = get_scaling_factor(self.aggregate_flow_heat)
            set_scaling_factor(self.aggregate_flow_heat_purchased, sf)

        if get_scaling_factor(self.aggregate_flow_electricity_sold) is None:
            set_scaling_factor(self.aggregate_flow_electricity_sold, 1)

        if get_scaling_factor(self.frac_elec_from_grid) is None:
            set_scaling_factor(self.frac_elec_from_grid, 1)

        if hasattr(self, "frac_heat_from_grid"):
            if get_scaling_factor(self.frac_heat_from_grid) is None:
                set_scaling_factor(self.frac_heat_from_grid, 1)

    def build_process_costs(self):
        """
        Not used in place of build_integrated_costs
        """
        pass

    def add_LCOT(self, flow_rate):
        """
        Add Levelized Cost of Treatment (LCOT) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOT
        """

        LCOT = pyo.Var(
            doc=f"Levelized Cost of Treatment based on flow {flow_rate.name}",
            units=self.base_currency / pyo.units.m**3,
        )
        self.add_component("LCOT", LCOT)

        LCOT_constraint = pyo.Constraint(
            expr=LCOT
            == (
                self.total_capital_cost * self.capital_recovery_factor
                + self.total_operating_cost
            )
            / (
                pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / self.base_period)
                * self.utilization_factor
            ),
            doc=f"Constraint for Levelized Cost of Treatment based on flow {flow_rate.name}",
        )
        self.add_component("LCOT_constraint", LCOT_constraint)

    def add_LCOE(self):
        """
        Add Levelized Cost of Energy (LCOE) to costing block.
        """

        energy_cost = self._get_energy_cost_block()
        if not hasattr(energy_cost, "LCOE"):
            energy_cost.add_LCOE()

        add_object_reference(self, "LCOE", energy_cost.LCOE)

    def add_LCOW(self, flow_rate, name="LCOW"):
        """
        Add Levelized Cost of Water (LCOW) to costing block.
        """

        treat_cost = self._get_treatment_cost_block()

        if not hasattr(treat_cost, "LCOW"):
            treat_cost.add_LCOW(flow_rate, name="LCOW")

        add_object_reference(self, name, getattr(treat_cost, name))

    def add_LCOE(self):
        """
        Add Levelized Cost of Energy (LCOE) to costing block.
        """

        energy_cost = self._get_energy_cost_block()
        if not hasattr(energy_cost, "LCOE"):
            energy_cost.add_LCOE()

        add_object_reference(self, "LCOE", energy_cost.LCOE)

    def add_LCOH(self):
        """
        Add Levelized Cost of Heat (LCOH) to costing block.
        """

        energy_cost = self._get_energy_cost_block()
        if not hasattr(energy_cost, "LCOH"):
            energy_cost.add_LCOH()

        add_object_reference(self, "LCOH", energy_cost.LCOH)

    def add_specific_electric_energy_consumption(self, flow_rate, name="SEEC"):
        """
        Add specific electric energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific electric energy consumption
        """
        treat_cost = self._get_treatment_cost_block()

        if not hasattr(treat_cost, "specific_electric_energy_consumption_constraint"):
            treat_cost.add_specific_electric_energy_consumption(flow_rate, name=name)

        add_object_reference(self, name, getattr(treat_cost, name))

    def add_specific_thermal_energy_consumption(self, flow_rate, name="STEC"):
        """
        Add specific thermal energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific thermal energy consumption
        """
        treat_cost = self._get_treatment_cost_block()

        if not hasattr(treat_cost, "specific_thermal_energy_consumption_constraint"):
            treat_cost.add_specific_thermal_energy_consumption(flow_rate, name=name)

        add_object_reference(self, name, getattr(treat_cost, name))

    def _check_common_param_equivalence(self, treat_cost, energy_cost):
        """
        Check if the common costing parameters across all three costing packages
        (treatment, energy, and system) have the same value.
        """

        common_params = [
            "electricity_cost",
            "heat_cost",
            "electrical_carbon_intensity",
            "maintenance_labor_chemical_factor",
            "plant_lifetime",
            "utilization_factor",
            "base_currency",
            "base_period",
            "sales_tax_frac",
            "TIC",
            "TPEC",
            "wacc",
        ]

        for cp in common_params:
            tp = getattr(treat_cost, cp)
            ep = getattr(energy_cost, cp)
            if (isinstance(tp, pyo.Var)) or isinstance(tp, pyo.Param):
                param_is_equivalent = pyo.value(tp) == pyo.value(ep)
            else:
                param_is_equivalent = tp == ep
            if not param_is_equivalent:
                err_msg = f"The common costing parameter {cp} was found to have a different value "
                err_msg += f"on the energy and treatment costing blocks. "
                err_msg += "Common costing parameters must be equivalent across all costing blocks "
                err_msg += "to use REFLOSystemCosting."
                raise ValueError(err_msg)

            if hasattr(self, cp):
                # if REFLOSystemCosting has this parameter,
                # we fix it to the treatment costing block value
                p = getattr(self, cp)
                if isinstance(p, pyo.Var):
                    p.fix(pyo.value(tp))
                elif isinstance(p, pyo.Param):
                    p.set_value(pyo.value(tp))

            if cp == "base_currency":
                self.base_currency = treat_cost.base_currency
            if cp == "base_period":
                self.base_period = treat_cost.base_period

    def _get_treatment_cost_block(self):
        """
        Get the TreatmentCosting block, if present.
        """
        tb = None
        for b in self.model().component_objects(pyo.Block):
            if isinstance(b, TreatmentCostingData):
                tb = b
        if tb is None:
            err_msg = "REFLOSystemCosting package requires a TreatmentCosting block"
            err_msg += " but one was not found."
            raise ValueError(err_msg)
        else:
            return tb

    def _get_energy_cost_block(self):
        """
        Get the EnergyCosting block, if present.
        """
        eb = None
        for b in self.model().component_objects(pyo.Block):
            if isinstance(b, EnergyCostingData):
                eb = b
        if eb is None:
            err_msg = "REFLOSystemCosting package requires a EnergyCosting block"
            err_msg += " but one was not found."
            raise ValueError(err_msg)
        else:
            return eb

    def _get_electricity_generation_unit(self):
        """
        Get the electricity generating unit on the flowsheet, if present.
        """
        elec_gen_unit = None
        for b in self.model().component_objects(pyo.Block):
            if isinstance(
                b, PVSurrogateData
            ):  # PV is only electricity generation model currently
                elec_gen_unit = b
            if isinstance(b, DummyElectricityUnit):  # only used for testing
                elec_gen_unit = b
        if elec_gen_unit is None:
            err_msg = (
                f"{self.name} indicated an electricity generation model was present "
            )
            err_msg += "on the flowsheet, but none was found."
            raise ValueError(err_msg)
        else:
            return elec_gen_unit

    def _get_pysam(self):
        """
        Get the PySAMWaterTAP block on flowsheet.
        """
        pysam_block_test_lst = []
        for k, v in vars(self.model()).items():
            if isinstance(v, PySAMWaterTAP):
                pysam_block_test_lst.append(k)

        if len(pysam_block_test_lst) != 1:
            raise ValueError("There is no instance of PySAMWaterTAP on this model.")

        else:
            pysam = getattr(self.model(), pysam_block_test_lst[0])
            return pysam
