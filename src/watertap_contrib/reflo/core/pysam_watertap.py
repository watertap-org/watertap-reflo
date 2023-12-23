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

import json
from os.path import join, dirname
from math import floor, ceil
import PySAM.Pvsamv1 as pv
import PySAM.Grid as grid
import PySAM.Utilityrate5 as utilityrate
import PySAM.Singleowner as singleowner

__author__ = "Kurban Sitterley, Matthew Boyd"


class PySAMWaterTAP:
    def __init__(
        self,
        pysam_model=None,
        tech_config_file=None,
        grid_config_file=None,
        rate_config_file=None,
        cash_config_file=None,
        weather_file=None,
        dcac_ratio=1.2,
    ):

        if pysam_model is None:
            raise Exception("PySAM module is required.")
        if not all(
            [tech_config_file, grid_config_file, rate_config_file, cash_config_file]
        ):
            raise Exception("One of the PySAM model configuration files is missing.")
        if weather_file is None:
            raise Exception("Weather file is required.")

        self._pysam_model = pysam_model
        self._tech_config_file = tech_config_file
        self._grid_config_file = grid_config_file
        self._rate_config_file = rate_config_file
        self._cash_config_file = cash_config_file
        self._weather_file = weather_file
        self._dcac_ratio = dcac_ratio
        self._config_files = [
            self._tech_config_file,
            self._grid_config_file,
            self._rate_config_file,
            self._cash_config_file,
        ]

        if self._pysam_model == "pv":
            self._pysam_model_name = "FlatPlatePVSingleOwner"
            self.setup_pv_single_owner()

    def setup_pv_single_owner(self):
        print(f"\nBuilding PySAM model {self._pysam_model_name}...\n")
        self.tech_model = pv.new()
        self.grid_model = grid.from_existing(self.tech_model, self._pysam_model_name)
        self.rate_model = utilityrate.from_existing(
            self.tech_model, self._pysam_model_name
        )
        self.cash_model = singleowner.from_existing(
            self.tech_model, self._pysam_model_name
        )
        self._modules = [
            self.tech_model,
            self.grid_model,
            self.rate_model,
            self.cash_model,
        ]
        self._load_config_files()
        self.tech_model.SolarResource.solar_resource_file = self._weather_file

        self._has_been_run = False

    def run_pv_single_owner(
        self,
        desired_size=50,
        desired_dcac_ratio=1.2,
        tech_model_kwargs={},
        cash_model_kwargs={},
    ):
        self._size_pv_array(
            desired_size=desired_size, desired_dcac_ratio=desired_dcac_ratio
        )
        print(
            f"\nRunning PySAM model {self._pysam_model_name} for desired size = {self.desired_size} kW...\n"
        )
        self.tech_model.value("inverter_count", self.num_inverters)
        self.tech_model.value("subarray1_nstrings", self.num_parallel)
        for param, val in tech_model_kwargs.items():
            # if not isinstance(val, list):
            #     val = [val]
            self.tech_model.value(param, val)
        for param, val in cash_model_kwargs.items():
            if not isinstance(val, list):
                val = [val]
            self.cash_model.value(param, val)
        for mod in self._modules:
            mod.execute()
        print("PySAM run finished.\n")
        self.tech_model = self._modules[0]
        self.cash_model = self._modules[3]
        self.annual_energy = self.tech_model.Outputs.annual_energy
        self.hourly_energy = self.tech_model.Outputs.gen
        self.hourly_energy = [x if x > 0 else 0 for x in self.hourly_energy]
        self.dc_capacity_factor = self.tech_model.Outputs.capacity_factor
        self.ac_capacity_factor = self.tech_model.Outputs.capacity_factor_ac
        self.lcoe_real = self.cash_model.Outputs.lcoe_real
        self._has_been_run = True

    def _size_pv_array(self, desired_size=50, desired_dcac_ratio=1.2):
        # Sizing rules
        # 1. Voc < Vdcmax
        # 2. Vmp > Vmin
        # 3. Vmp < Vmax
        # 4. num series * num_parallel is about desired array size (num_parallel = desired / (num series * mod_power)
        # 5. num inverters is about desired array size (num_inv = num_series * num_parallel * mod_power) / inv_power

        if not hasattr(self, "desired_size") or desired_size is not None:
            self.desired_size = desired_size
        if not hasattr(self, "desired_dcac_ratio") or desired_dcac_ratio is not None:
            self.desired_dcac_ratio = desired_dcac_ratio

        self.m = m = self._flatten_dict(self.tech_model.export())

        # module parameters
        module_model = int(m["module_model"])
        m["spe_imp"] = (
            self._spe_power(m["spe_eff4"], m["spe_rad4"], m["spe_area"]) / m["spe_vmp"]
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
        num_parallel = self.desired_size * 1000 / (num_series * mod_power)
        num_parallel = max(1, round(num_parallel))
        if self.desired_dcac_ratio > 0:
            inverters = ((num_series * num_parallel * mod_power)) / (
                self.desired_dcac_ratio * inv_power
            )
            # round inverters for best DC-AC ratio
            if inverters - floor(inverters) < 0.5:
                num_inverters = floor(inverters)
            else:
                num_inverters = ceil(inverters)
        else:
            num_inverters = ceil(((num_series * num_parallel * mod_power)) / inv_power)
        num_inverters = max(1, num_inverters)
        total_modules = num_series * num_parallel
        total_ac_capacity = inv_power * num_inverters / 1000
        total_dc_inverter_capacity = inv_dc_power * num_inverters / 1000

        # check that the sizing was close to the desired sizes, otherwise error out
        nameplate_dc = total_modules * mod_power / 1000
        proposed_ratio = nameplate_dc / (num_inverters * inv_power / 1000)
        if abs(nameplate_dc - self.desired_size) / self.desired_size > 0.2:
            num_inverters = None
            num_series = None
            num_parallel = None
            total_modules = None
            nameplate_dc = None
            total_ac_capacity = None
            total_dc_inverter_capacity = None
        else:
            self.num_inverters = num_inverters
            self.num_series = num_series
            self.num_parallel = num_parallel
            self.total_modules = total_modules
            self.nameplate_dc = nameplate_dc
            self.total_ac_capacity = total_ac_capacity
            self.total_dc_inverter_capacity = total_dc_inverter_capacity

        self.size_pv_array = {
            "inverter_count": num_inverters,
            "subarray1_modules_per_string": num_series,
            "subarray1_nstrings": num_parallel,
            "total_modules": total_modules,
            "system_capacity": nameplate_dc,
            "total_inverter_capacity": total_ac_capacity,
            "total_dc_inverter_capacity": total_dc_inverter_capacity,
        }

    def _flatten_dict(self, d):
        def get_key_values(d):
            for key, value in d.items():
                if isinstance(value, dict):
                    yield from get_key_values(value)
                else:
                    yield key, value

        return {key: value for (key, value) in get_key_values(d)}

    def _load_config_files(self):
        for file_name, module in zip(self._config_files, self._modules):
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

    def _spe_power(self, spe_eff_level, spe_rad_level, spe_area):
        return spe_eff_level / 100 * spe_rad_level * spe_area
