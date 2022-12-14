from watertap.core.zero_order_costing import ZeroOrderCosting, ZeroOrderCostingData

from idaes.core import declare_process_block_class

from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

import pyomo.environ as pyo

from watertap_contrib.seto.solar_models.zero_order import PhotovoltaicZO


@declare_process_block_class("SETOZeroOrderCosting")
class SETOZeroOrderCostingData(ZeroOrderCostingData):
    """
    General costing package for zero-order processes modified for SETO project.
    """

    CONFIG = ZeroOrderCostingData.CONFIG()
    # CONFIG.declare()

    # def build_global_params(self):
    #     super().build_global_params()

    def build_process_costs(self):
        super().build_process_costs()

        self.total_variable_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Total variable operating cost",
            units=self.base_currency / self.base_period,
        )

    # def initialize_build(self):
    #     super().initialize_build()

    def cost_pv(blk):

        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of solar energy generation system",
        )

        blk.direct_cost = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct costs of PV system",
        )

        blk.indirect_cost = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Indirect costs of PV system",
        )

        blk.sales_tax = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Sales tax for PV system",
        )

        blk.system_capacity = pyo.Var(
            initialize=0,
            units=pyo.units.watt,
            bounds=(0, None),
            doc="DC system capacity for PV system",
        )

        blk.cost_per_watt_module = pyo.Var(
            initialize=0.41,
            units=blk.config.flowsheet_costing_block.base_currency / pyo.units.watt,
            bounds=(0, None),
            doc="Cost per watt for solar module",
        )

        blk.cost_per_watt_inverter = pyo.Var(
            initialize=0.05,
            units=blk.config.flowsheet_costing_block.base_currency / pyo.units.watt,
            bounds=(0, None),
            doc="Cost per watt for inverter",
        )

        blk.cost_per_watt_other = pyo.Var(
            initialize=0.1,
            units=blk.config.flowsheet_costing_block.base_currency / pyo.units.watt,
            bounds=(0, None),
            doc="Cost per watt for other equipment, installation, and margin/overhead",
        )

        blk.cost_per_watt_indirect = pyo.Var(
            initialize=0.13,
            units=blk.config.flowsheet_costing_block.base_currency / pyo.units.watt,
            bounds=(0, None),
            doc="Cost per watt for permitting, environmental studies, engineering, land prep, and grid interconnection",
        )

        blk.land_cost_per_acre = pyo.Var(
            initialize=4000,
            units=blk.config.flowsheet_costing_block.base_currency / pyo.units.acre,
            bounds=(0, None),
            doc="Land cost per acre required",
        )

        blk.land_area = pyo.Var(
            initialize=0,
            units=pyo.units.acre,
            bounds=(0, None),
            doc="Land area required for PV system",
        )

        blk.contingency_frac_direct_cost = pyo.Var(
            initialize=0.03,
            units=pyo.units.dimensionless,
            bounds=(0, 1),
            doc="Fraction of direct costs for contingency",
        )

        blk.tax_frac_direct_cost = pyo.Var(
            initialize=0.05,
            units=pyo.units.dimensionless,
            bounds=(0, 1),
            doc="Fraction of direct costs for sales tax",
        )

        blk.direct_cost_constraint = pyo.Constraint(
            expr=blk.direct_cost
            == blk.system_capacity
            * (
                blk.cost_per_watt_module
                + blk.cost_per_watt_inverter
                + blk.cost_per_watt_other
            )
            + (
                blk.system_capacity
                * (
                    blk.cost_per_watt_module
                    + blk.cost_per_watt_inverter
                    + blk.cost_per_watt_other
                )
            )
            * blk.contingency_frac_direct_cost
        )

        blk.indirect_cost_constraint = pyo.Constraint(
            expr=blk.indirect_cost
            == (blk.system_capacity * blk.cost_per_watt_indirect)
            + (blk.land_area * blk.land_cost_per_acre)
        )

        blk.sales_tax_constraint = pyo.Constraint(
            expr=blk.sales_tax == blk.direct_cost * blk.tax_frac_direct_cost
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.direct_cost + blk.indirect_cost + blk.sales_tax
        )

        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of solar system",
        )

        blk.fixed_operating_by_capacity = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency
            / (pyo.units.kW * blk.config.flowsheet_costing_block.base_period),
            bounds=(0, None),
            doc="Fixed operating cost of solar system per kW generated",
        )

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == blk.fixed_operating_by_capacity
            * pyo.units.convert(blk.system_capacity, to_units=pyo.units.kW)
        )

        blk.variable_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Variable operating cost of solar",
        )

        blk.annual_generation = pyo.Var(
            initialize=0,
            units=pyo.units.MWh,
            bounds=(0, None),
            doc="Annual electricity generation of solar system",
        )

        blk.variable_operating_by_generation = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency
            / (pyo.units.MWh * blk.config.flowsheet_costing_block.base_period),
            bounds=(0, None),
            doc="Annual operating cost of solar system per MWh generated",
        )

        blk.variable_operating_cost_constraint = pyo.Constraint(
            expr=blk.variable_operating_cost
            == blk.variable_operating_by_generation * blk.annual_generation
        )
        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity, "electricity"
        )

    unit_mapping = {PhotovoltaicZO: cost_pv}
