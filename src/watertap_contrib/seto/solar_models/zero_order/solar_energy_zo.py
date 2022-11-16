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
"""
This module contains a zero-order representation of a solar energy
operation.
"""


from idaes.core import declare_process_block_class
from idaes.core.util.misc import StrEnum

from pyomo.environ import Constraint, Var, units as pyunits
from pyomo.common.config import ConfigValue, In

from watertap.core import build_pt, ZeroOrderBaseData
from watertap_contrib.seto.energy import solar_energy

import PySAM.Pvsamv1 as pv
import PySAM.Grid as grid
import PySAM.Utilityrate5 as utilityrate
import PySAM.Singleowner as singleowner

import json
import numpy as np

__author__ = "Kurban Sitterley"


class SolarEnergyZOType(StrEnum):
    PV = "PV"


@declare_process_block_class("SolarEnergyZO")
class SolarEnergyZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a solar energy unit model.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()
    CONFIG.declare(
        "solar_energy_type",
        ConfigValue(
            default="PV",
            domain=In(SolarEnergyZOType),
            description="Indicates type of solar energy ZO model",
        ),
    )

    def build(self):
        super().build()

        self._tech_type = "solar_energy"

        build_pt(self)
        solar_energy(self)

    def pysam_pv(self):
        """
        Something something PySAM
        """

        def spe_power(spe_eff_level, spe_rad_level, spe_area):
            return spe_eff_level / 100 * spe_rad_level * spe_area

        def size_pv_array(tech_model, desired_size, desired_dcac_ratio):
            # Sizing rules
            # 1. Voc < Vdcmax
            # 2. Vmp > Vmin
            # 3. Vmp < Vmax
            # 4. num series * num_parallel is about desired array size (num_parallel = desired / (num series * mod_power)
            # 5. num inverters is about desired array size (num_inv = num_series * num_parallel * mod_power) / inv_power

            m = flatten_dict(tech_model.export())

            # module parameters
            module_model = int(m["module_model"])
            m["spe_imp"] = (
                spe_power(m["spe_eff4"], m["spe_rad4"], m["spe_area"]) / m["spe_vmp"]
            )  # 4 = reference conditions
            mod_imp = [
                m["spe_imp"],
                m["cec_i_mp_ref"],
                m["sixpar_imp"],
                m["snl_impo"],
                m["sd11par_Imp0"],
            ][module_model]
            mod_vmp = [
                m["spe_vmp"],
                m["cec_v_mp_ref"],
                m["sixpar_vmp"],
                m["snl_vmpo"],
                m["sd11par_Vmp0"],
            ][module_model]
            mod_voc = [
                m["spe_voc"],
                m["cec_v_oc_ref"],
                m["sixpar_voc"],
                m["snl_voco"],
                m["sd11par_Voc0"],
            ][module_model]
            mod_power = mod_vmp * mod_imp

            # inverter parameters
            inv_vmin = m["mppt_low_inverter"]
            inv_vmax = m["mppt_hi_inverter"]
            inverter_model = int(m["inverter_model"])
            m["inv_ds_pdco"] = m["inv_ds_paco"] / (m["inv_ds_eff"] / 100)
            inv_vdcmax = [
                m["inv_snl_vdcmax"],
                m["inv_ds_vdcmax"],
                m["inv_pd_vdcmax"],
                m["inv_cec_cg_vdcmax"],
            ][
                inverter_model
            ]  # Vdcmax
            inv_power = [
                m["inv_snl_paco"],
                m["inv_ds_paco"],
                m["inv_pd_paco"],
                m["inv_cec_cg_paco"],
            ][
                inverter_model
            ]  # Paco
            inv_dc_power = [
                m["inv_snl_pdco"],
                m["inv_ds_pdco"],
                m["inv_pd_pdco"],
                m["inv_cec_cg_pdco"],
            ][
                inverter_model
            ]  # Pdco

            # DC-connected battery parameters (assumed to use common inverter)
            batt_max_power_dc = 0
            if m["en_batt"] and m["batt_ac_or_dc"] == 0:
                batt_max_power_dc = m["batt_max_power"]

            if mod_vmp > 0:
                num_series = 0.5 * (inv_vmin + inv_vmax) / mod_vmp
                if inv_vdcmax > 0:
                    while num_series > 0 and (num_series * mod_voc) > inv_vdcmax:
                        num_series -= 1

            num_series = max(1, round(num_series))
            num_parallel = desired_size * 1000 / (num_series * mod_power)
            num_parallel = max(1, round(num_parallel))
            if desired_dcac_ratio > 0:
                inverters = ((num_series * num_parallel * mod_power)) / (
                    desired_dcac_ratio * inv_power
                )
                # round inverters for best DC-AC ratio
                if inverters - np.floor(inverters) < 0.5:
                    num_inverters = np.floor(inverters)
                else:
                    num_inverters = np.ceil(inverters)
            else:
                num_inverters = np.ceil(
                    ((num_series * num_parallel * mod_power)) / inv_power
                )
            num_inverters = max(1, num_inverters)
            total_modules = num_series * num_parallel
            total_ac_capacity = inv_power * num_inverters / 1000
            total_dc_inverter_capacity = inv_dc_power * num_inverters / 1000

            # check that the sizing was close to the desired sizes, otherwise error out
            nameplate_dc = total_modules * mod_power / 1000
            proposed_ratio = nameplate_dc / (num_inverters * inv_power / 1000)
            if abs(nameplate_dc - desired_size) / desired_size > 0.2:
                num_inverters = None
                num_series = None
                num_parallel = None
                total_modules = None
                nameplate_dc = None
                total_ac_capacity = None
                total_dc_inverter_capacity = None

            return {
                "inverter_count": num_inverters,
                "subarray1_modules_per_string": num_series,
                "subarray1_nstrings": num_parallel,
                "total_modules": total_modules,
                "system_capacity": nameplate_dc,
                "total_inverter_capacity": total_ac_capacity,
                "total_dc_inverter_capacity": total_dc_inverter_capacity,
            }

        def setup_pv_singleowner(model_name, config_files, weather_file):
            tech_model = pv.new()  # .default()
            grid_model = grid.from_existing(tech_model, model_name)
            rate_model = utilityrate.from_existing(tech_model, model_name)
            cash_model = singleowner.from_existing(tech_model, model_name)
            modules = [tech_model, grid_model, rate_model, cash_model]

            load_config_files(config_files, modules)
            tech_model.SolarResource.solar_resource_file = weather_file

            return {
                "tech_model": tech_model,
                "grid_model": grid_model,
                "rate_model": rate_model,
                "cash_model": cash_model,
            }  # ORDER IS IMPORTANT - must run models

        def get_pv_energy_gen(desired_size, model_name=None, dcac_ratio=1.2):
            if model_name is None:
                model_name = "FlatPlatePVSingleOwner"

            # These config files are gotten from the SAM UI -> Generate code... -> PySAM JSON :
            config_files = [
                "pysam/untitled_pvsamv1.json",
                "pysam/untitled_grid.json",
                "pysam/untitled_utilityrate5.json",
                "pysam/untitled_singleowner.json",
            ]

            weather_file = "pysam/phoenix_az_33.450495_-111.983688_psmv3_60_tmy.csv"

            modules = setup_pv_singleowner(model_name, config_files, weather_file)

            tech_model = modules["tech_model"]
            cash_model = modules["cash_model"]
            cash_model.value("om_fixed", [1e4])
            cash_model.value("om_production", [20])
            results = size_pv_array(tech_model, desired_size, dcac_ratio)
            tech_model.value("inverter_count", results["inverter_count"])
            tech_model.value("subarray1_nstrings", results["subarray1_nstrings"])
            nameplate_dc_capacity = results["system_capacity"]
            total_ac_capacity = results["total_inverter_capacity"]
            for module in modules.values():
                module.execute()
            return cash_model, tech_model, results


def flatten_dict(d):
    def get_key_values(d):
        for key, value in d.items():
            if isinstance(value, dict):
                yield from get_key_values(value)
            else:
                yield key, value

    return {key: value for (key, value) in get_key_values(d)}


def load_config_files(file_names, modules):
    for file_name, module in zip(file_names, modules):
        with open(file_name, "r") as file:
            data = json.load(file)
            missing_values = []  # for debugging
            for k, v in data.items():
                if k != "number_inputs":
                    try:
                        module.value(k, v)
                    except:
                        missing_values.append(k)
            pass
