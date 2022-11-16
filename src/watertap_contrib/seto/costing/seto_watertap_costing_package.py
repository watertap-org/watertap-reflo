# from watertap.costing import ROType, PumpType, EnergyRecoveryDeviceType, MixerType, CrystallizerCostType
from watertap.costing.watertap_costing_package import WaterTAPCostingData

# from enum import Enum

import pyomo.environ as pyo

# from pyomo.util.calc_var_value import calculate_variable_from_constraint

from watertap_contrib.seto.solar_models.zero_order import SolarEnergyZO
from idaes.core import declare_process_block_class

from idaes.models.unit_models import Mixer
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

@declare_process_block_class("SETOWaterTAPCosting")
class SETOWaterTAPCostingData(WaterTAPCostingData):

    def build(self):
        super().build()
    
    def build_global_params(self):
        super().build_global_params()

        self.base_currency = pyo.units.USD_2020
    
    @staticmethod
    def cost_solar_energy(blk):
        pass 
     
SETOWaterTAPCostingData.unit_mapping = {SolarEnergyZO: SETOWaterTAPCostingData.cost_solar_energy,
                        Mixer: WaterTAPCostingData.cost_mixer,
                        Pump: WaterTAPCostingData.cost_pump,
                        EnergyRecoveryDevice: WaterTAPCostingData.cost_energy_recovery_device,
                        PressureExchanger: WaterTAPCostingData.cost_pressure_exchanger,
                        ReverseOsmosis0D: WaterTAPCostingData.cost_reverse_osmosis,
                        ReverseOsmosis1D: WaterTAPCostingData.cost_reverse_osmosis,
                        NanoFiltration0D: WaterTAPCostingData.cost_nanofiltration,
                        NanofiltrationZO: WaterTAPCostingData.cost_nanofiltration,
                        Crystallization: WaterTAPCostingData.cost_crystallizer,
                        Ultraviolet0D: WaterTAPCostingData.cost_uv_aop,
                        Electrodialysis0D: WaterTAPCostingData.cost_electrodialysis,
                        Electrodialysis1D: WaterTAPCostingData.cost_electrodialysis,
                        IonExchange0D: WaterTAPCostingData.cost_ion_exchange,
                        GAC: WaterTAPCostingData.cost_gac,}

def make_capital_cost_var(blk):
    blk.capital_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit capital cost",
    )