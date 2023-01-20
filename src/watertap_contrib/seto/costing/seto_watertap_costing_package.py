import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.models.unit_models import Mixer

from watertap.costing.watertap_costing_package import WaterTAPCostingData
from watertap.unit_models import (
    ReverseOsmosis0D,
    ReverseOsmosis1D,
    NanoFiltration0D,
    NanofiltrationZO,
    PressureExchanger,
    Crystallization,
    Ultraviolet0D,
    Pump,
    EnergyRecoveryDevice,
    Electrodialysis0D,
    Electrodialysis1D,
    IonExchange0D,
    GAC,
)

from watertap.costing.units.crystallizer import cost_crystallizer
from watertap.costing.units.electrodialysis import cost_electrodialysis
from watertap.costing.units.energy_recovery_device import cost_energy_recovery_device
from watertap.costing.units.gac import cost_gac
from watertap.costing.units.ion_exchange import cost_ion_exchange
from watertap.costing.units.nanofiltration import cost_nanofiltration
from watertap.costing.units.mixer import cost_mixer
from watertap.costing.units.pressure_exchanger import cost_pressure_exchanger
from watertap.costing.units.pump import cost_pump
from watertap.costing.units.reverse_osmosis import cost_reverse_osmosis
from watertap.costing.units.uv_aop import cost_uv_aop

from watertap_contrib.seto.solar_models.zero_order import SolarEnergyZO, PhotovoltaicZO
from watertap_contrib.seto.costing.solar.photovoltaic import cost_pv
from watertap_contrib.seto.unit_models.surrogate import LTMEDSurrogate
from watertap_contrib.seto.costing.units.lt_med_surrogate import cost_lt_med_surrogate


@declare_process_block_class("SETOWaterTAPCosting")
class SETOWaterTAPCostingData(WaterTAPCostingData):
    unit_mapping = {
        SolarEnergyZO: cost_pv,  # Keeping this for now
        LTMEDSurrogate: cost_lt_med_surrogate,
        PhotovoltaicZO: cost_pv,
        Mixer: cost_mixer,
        Pump: cost_pump,
        EnergyRecoveryDevice: cost_energy_recovery_device,
        PressureExchanger: cost_pressure_exchanger,
        ReverseOsmosis0D: cost_reverse_osmosis,
        ReverseOsmosis1D: cost_reverse_osmosis,
        NanoFiltration0D: cost_nanofiltration,
        NanofiltrationZO: cost_nanofiltration,
        Crystallization: cost_crystallizer,
        Ultraviolet0D: cost_uv_aop,
        Electrodialysis0D: cost_electrodialysis,
        Electrodialysis1D: cost_electrodialysis,
        IonExchange0D: cost_ion_exchange,
        GAC: cost_gac,
    }

    def build_global_params(self):
        super().build_global_params()

        if "USD_2021" not in pyo.units._pint_registry:
            pyo.units.load_definitions_from_strings(["USD_2021 = 500/708.0 * USD_CE500"])

        self.base_currency = pyo.units.USD_2021
        self.plant_lifetime = pyo.Var(
            initialize=20, units=self.base_period, doc="Plant lifetime"
        )

        self.heat_cost = pyo.Param(
            mutable=True,
            initialize=0.01,
            doc="Heat cost",
            units=pyo.units.USD_2018 / pyo.units.kWh,
        )
        self.add_defined_flow("heat", self.heat_cost)

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
