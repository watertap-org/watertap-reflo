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

import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.models.unit_models import Mixer

from watertap.costing.watertap_costing_package import (
    WaterTAPCostingData,
)

from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

from watertap_contrib.reflo.core import PySAMWaterTAP


@declare_process_block_class("REFLOCosting")
class REFLOCostingData(WaterTAPCostingData):
    def build_global_params(self):
        super().build_global_params()

        self.base_currency = pyo.units.USD_2021
        self.plant_lifetime = pyo.Var(
            initialize=20, units=self.base_period, doc="Plant lifetime"
        )

        self.sales_tax_frac = pyo.Param(
            initialize=0.05,
            mutable=True,
            doc="Sales tax as fraction of capital costs",
            units=pyo.units.dimensionless,
        )

        self.heat_cost = pyo.Param(
            mutable=True,
            initialize=0.01,
            doc="Heat cost",
            units=pyo.units.USD_2018 / pyo.units.kWh,
        )
        self.register_flow_type("heat", self.heat_cost)

        self.plant_lifetime.fix()
        self.utilization_factor.fix(1)

    def build_process_costs(self):
        # super().build_process_costs()

        self.total_capital_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Total capital cost",
            units=self.base_currency,
        )
        self.maintenance_labor_chemical_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Maintenance-labor-chemical operating cost",
            units=self.base_currency / self.base_period,
        )
        self.total_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.Reals,
            doc="Total operating cost",
            units=self.base_currency / self.base_period,
        )

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.factor_total_investment * self.aggregate_capital_cost
        )
        self.maintenance_labor_chemical_operating_cost_constraint = pyo.Constraint(
            expr=self.maintenance_labor_chemical_operating_cost
            == self.factor_maintenance_labor_chemical * self.total_capital_cost
        )

        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == self.maintenance_labor_chemical_operating_cost
            + self.aggregate_fixed_operating_cost
            + self.aggregate_variable_operating_cost
            + sum(self.aggregate_flow_costs.values()) * self.utilization_factor
        )


@declare_process_block_class("TreatmentCosting")
class TreatmentCostingData(REFLOCostingData):
    def build_global_params(self):
        super().build_global_params()

    def build_process_costs(self):
        super().build_process_costs()


@declare_process_block_class("EnergyCosting")
class EnergyCostingData(REFLOCostingData):
    def build_global_params(self):
        super().build_global_params()

    def build_process_costs(self):
        super().build_process_costs()


@declare_process_block_class("REFLOSystemCosting")
class REFLOSystemCostingData(FlowsheetCostingBlockData):
    def build(self):
        super().build()

        self._registered_LCOWs = {}

    def build_global_params(self):
        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        self.base_currency = pyo.units.USD_2021

        self.base_period = pyo.units.year

        self.utilization_factor = pyo.Var(
            initialize=1,
            doc="Plant capacity utilization [fraction of uptime]",
            units=pyo.units.dimensionless,
        )

        self.plant_lifetime = pyo.Var(
            initialize=20, units=self.base_period, doc="Plant lifetime"
        )

        self.factor_total_investment = pyo.Var(
            initialize=1,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.factor_maintenance_labor_chemical = pyo.Var(
            initialize=0.03,
            doc="Maintenance-labor-chemical factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )

        self.wacc = pyo.Param(
            initialize=0.05,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Weighted Average Cost of Capital [WACC]",
        )

        self.electricity_cost = pyo.Param(
            mutable=True,
            initialize=0.0718,  # From EIA for 2021
            doc="Electricity cost",
            units=self.base_currency / pyo.units.kWh,
        )

        self.register_flow_type("electricity", self.electricity_cost)

        self.electrical_carbon_intensity = pyo.Param(
            mutable=True,
            initialize=0.475,
            doc="Grid carbon intensity [kgCO2_eq/kWh]",
            units=pyo.units.kg / pyo.units.kWh,
        )

        self.factor_capital_annualization = pyo.Expression(
            expr=(
                (
                    self.wacc
                    * (1 + self.wacc) ** (self.plant_lifetime / self.base_period)
                )
                / (((1 + self.wacc) ** (self.plant_lifetime / self.base_period)) - 1)
                / self.base_period
            )
        )
        # fix the parameters
        self.fix_all_vars()
        # Build the integrated system costs
        self.build_integrated_costs()

    def build_process_costs(self):
        pass

    def build_integrated_costs(self):
        treat_cost = self._get_treatment_cost_block()
        en_cost = self._get_energy_cost_block()

        self.total_capital_cost = pyo.Var(
            initialize=1e3,
            # domain=pyo.NonNegativeReals,
            doc="Total capital cost for integrated system",
            units=self.base_currency,
        )
        self.total_operating_cost = pyo.Var(
            initialize=1e3,
            # domain=pyo.NonNegativeReals,
            doc="Total operating cost for integrated system",
            units=self.base_currency / self.base_period,
        )
        self.aggregate_flow_electricity = pyo.Var(
            initialize=1e3,
            # domain=pyo.NonNegativeReals,
            doc="Aggregated electricity flow",
            units=pyo.units.kW,
        )

        # if all("heat" in b.defined_flows for b in [treat_cost, en_cost]):
        if all(hasattr(b, "aggregate_flow_heat") for b in [treat_cost, en_cost]):
            self.aggregate_flow_heat = pyo.Var(
                initialize=1e3,
                # domain=pyo.NonNegativeReals,
                doc="Aggregated heat flow",
                units=pyo.units.kW,
            )

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == pyo.units.convert(
                treat_cost.total_capital_cost + en_cost.total_capital_cost,
                to_units=self.base_currency,
            )
        )

        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == pyo.units.convert(
                treat_cost.total_operating_cost + en_cost.total_operating_cost,
                to_units=self.base_currency / self.base_period,
            )
        )

        self.aggregate_flow_electricity_constraint = pyo.Constraint(
            expr=self.aggregate_flow_electricity
            == treat_cost.aggregate_flow_electricity
            + en_cost.aggregate_flow_electricity
        )

        # if all("heat" in b.defined_flows for b in [treat_cost, en_cost]):
        if all(hasattr(b, "aggregate_flow_heat") for b in [treat_cost, en_cost]):
            self.aggregate_flow_heat_constraint = pyo.Constraint(
                expr=self.aggregate_flow_heat
                == treat_cost.aggregate_flow_heat + en_cost.aggregate_flow_heat
            )

    def add_LCOW(self, flow_rate, name="LCOW"):
        """
        Add Levelized Cost of Water (LCOW) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOW
            name (optional) - name for the LCOW variable (default: LCOW)
        """

        LCOW = pyo.Var(
            doc=f"Levelized Cost of Water based on flow {flow_rate.name}",
            units=self.base_currency / pyo.units.m**3,
        )
        self.add_component(name, LCOW)

        LCOW_constraint = pyo.Constraint(
            expr=LCOW
            == (
                self.total_capital_cost * self.factor_capital_annualization
                + self.total_operating_cost
            )
            / (
                pyo.units.convert(
                    flow_rate, to_units=pyo.units.m**3 / self.base_period
                )
                * self.utilization_factor
            ),
            doc=f"Constraint for Levelized Cost of Water based on flow {flow_rate.name}",
        )
        self.add_component(name + "_constraint", LCOW_constraint)

        self._registered_LCOWs[name] = (LCOW, LCOW_constraint)

    def add_LCOE(self, e_model="pysam"):
        """
        Add Levelized Cost of Energy (LCOE) to costing block.
        Args:
            e_model - energy modeling approach used (PySAM or surrogate)
        """

        if e_model == "pysam":
            pysam = self._get_pysam()

            if not pysam._has_been_run:
                raise Exception(
                    f"PySAM model {pysam._pysam_model_name} has not yet been run, so there is no annual_energy data available."
                    "You must run the PySAM model before adding LCOE metric."
                )

            en_cost = self._get_energy_cost_block()

            self.annual_energy_generated = pyo.Param(
                initialize=pysam.annual_energy,
                units=pyo.units.kWh / pyo.units.year,
                doc=f"Annual energy generated by {pysam._pysam_model_name}",
            )
            LCOE_expr = pyo.Expression(
                expr=(
                    en_cost.total_capital_cost * self.factor_capital_annualization
                    + (
                        en_cost.aggregate_fixed_operating_cost
                        + en_cost.aggregate_variable_operating_cost
                    )
                )
                / self.annual_energy_generated
                * self.utilization_factor
            )
            self.add_component("LCOE", LCOE_expr)

        if e_model == "surrogate":
            raise NotImplementedError(
                "add_LCOE for surrogate models not available yet."
            )

    def add_specific_electric_energy_consumption(self, flow_rate):
        """
        Add specific electric energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific energy consumption
        """

        specific_electric_energy_consumption = pyo.Var(
            initialize=100,
            doc=f"Specific electric energy consumption based on flow {flow_rate.name}",
        )

        self.add_component(
            "specific_electric_energy_consumption", specific_electric_energy_consumption
        )

        specific_electric_energy_consumption_constraint = pyo.Constraint(
            expr=specific_electric_energy_consumption
            == self.aggregate_flow_electricity
            / pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / pyo.units.hr)
        )

        self.add_component(
            "specific_electric_energy_consumption_constraint",
            specific_electric_energy_consumption_constraint,
        )

    def add_specific_thermal_energy_consumption(self, flow_rate):
        """
        Add specific thermal energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific energy consumption
        """

        specific_thermal_energy_consumption = pyo.Var(
            initialize=100,
            doc=f"Specific thermal energy consumption based on flow {flow_rate.name}",
        )

        self.add_component(
            "specific_thermal_energy_consumption", specific_thermal_energy_consumption
        )

        specific_thermal_energy_consumption_constraint = pyo.Constraint(
            expr=specific_thermal_energy_consumption
            == self.aggregate_flow_heat
            / pyo.units.convert(flow_rate, to_units=pyo.units.m**3 / pyo.units.hr)
        )

        self.add_component(
            "specific_thermal_energy_consumption_constraint",
            specific_thermal_energy_consumption_constraint,
        )

    def add_defined_flow(self, flow_name, flow_cost):
        """
        This method adds a defined flow to the costing block.

        NOTE: Use this method to add `defined_flows` to the costing block
              to ensure updates to `flow_cost` get propagated in the model.
              See https://github.com/IDAES/idaes-pse/pull/1014 for details.

        Args:
            flow_name: string containing the name of the flow to register
            flow_cost: Pyomo expression that represents the flow unit cost

        Returns:
            None
        """
        flow_cost_name = flow_name + "_cost"
        current_flow_cost = self.component(flow_cost_name)
        if current_flow_cost is None:
            self.add_component(flow_cost_name, pyo.Expression(expr=flow_cost))
            self.defined_flows._setitem(flow_name, self.component(flow_cost_name))
        elif current_flow_cost is flow_cost:
            self.defined_flows._setitem(flow_name, current_flow_cost)
        else:
            # if we get here then there's an attribute named
            # flow_cost_name on the block, which is an error
            raise RuntimeError(
                f"Attribute {flow_cost_name} already exists "
                f"on the costing block, but is not {flow_cost}"
            )

    def _get_treatment_cost_block(self):
        for b in self.model().component_objects(pyo.Block):
            if isinstance(b, TreatmentCostingData):
                return b

    def _get_energy_cost_block(self):
        for b in self.model().component_objects(pyo.Block):
            if isinstance(b, EnergyCostingData):
                return b

    def _get_pysam(self):
        pysam_block_test_lst = []
        for k, v in vars(self.model()).items():
            if isinstance(v, PySAMWaterTAP):
                pysam_block_test_lst.append(k)

        if len(pysam_block_test_lst) != 1:
            raise Exception("There is no instance of PySAMWaterTAP on this model.")

        else:
            pysam = getattr(self.model(), pysam_block_test_lst[0])
            return pysam