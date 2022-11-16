from watertap.core.zero_order_costing import ZeroOrderCosting, ZeroOrderCostingData

from idaes.core import declare_process_block_class

from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

import pyomo.environ as pyo

from watertap_contrib.seto.solar_models.zero_order import SolarEnergyZO


@declare_process_block_class("SETOZeroOrderCosting")
class SETOZeroOrderCostingData(ZeroOrderCostingData):
    """
    General costing package for zero-order processes modified for SETO project.
    """

    CONFIG = ZeroOrderCostingData.CONFIG()
    # CONFIG.declare()

    def build_global_params(self):
        super().build_global_params()

    def build_process_costs(self):
        super().build_process_costs()

    def initialize_build(self):
        super().initialize_build()

    def cost_solar_energy(blk):

        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of solar energy generation system",
        )

        blk.installed_cost = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of solar energy generation system",
        )

        blk.installed_cost.fix()

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.installed_cost
        )

        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of solar system",
        )

        blk.system_capacity = pyo.Var(
            initialize=0,
            units=pyo.units.kW,
            bounds=(0, None),
            doc="Capacity of solar system",
        )

        blk.system_capacity.fix()

        blk.fixed_operating_by_capacity = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency
            / (pyo.units.kW * blk.config.flowsheet_costing_block.base_period),
            bounds=(0, None),
            doc="Fixed operating cost of solar system per kW generated",
        )

        blk.fixed_operating_by_capacity.fix()

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == blk.fixed_operating_by_capacity * blk.system_capacity
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

        blk.annual_generation.fix()

        blk.variable_operating_by_generation = pyo.Var(
            initialize=0,
            units=blk.config.flowsheet_costing_block.base_currency
            / (pyo.units.MWh * blk.config.flowsheet_costing_block.base_period),
            bounds=(0, None),
            doc="Annual operating cost of solar system per MWh generated",
        )

        blk.variable_operating_by_generation.fix()

        blk.variable_operating_cost_constraint = pyo.Constraint(
            expr=blk.variable_operating_cost
            == blk.variable_operating_by_generation * blk.annual_generation
        )
        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity, "electricity"
        )

    unit_mapping = {SolarEnergyZO: cost_solar_energy}
