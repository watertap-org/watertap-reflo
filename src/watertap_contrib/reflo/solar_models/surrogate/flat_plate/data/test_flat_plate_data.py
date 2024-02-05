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

import pytest
from pathlib import Path
import PySAM.Swh as swh
from .pysam_run_flat_plate import *


class TestFlatPlateData:
    @pytest.fixture(scope="class")
    def pvdata_frame(self):
        temperatures = {
            "T_cold": 20,
            "T_hot": 70,  # this will be overwritten by temperature_hot value
            "T_amb": 18,
        }
        param_file = Path(__file__).parent / "swh-reflo.json"
        weather_file = (
            Path(__file__).parent / "tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
        )
        config_data = read_module_datafile(param_file)
        if "solar_resource_file" in config_data:
            del config_data["solar_resource_file"]
        tech_model = setup_model(
            temperatures=temperatures,
            weather_file=str(weather_file),
            config_data=config_data,
        )
        return tech_model

    @pytest.mark.unit
    def test_run_model(self, pvdata_frame):
        tech_model = pvdata_frame
        result = run_model(
            tech_model, heat_load_mwt=1000, hours_storage=1, temperature_hot=70
        )
        assert result["heat_annual"] == pytest.approx(2163996973, rel=1e-2)
        assert result["electricity_annual"] == pytest.approx(72998038, rel=1e-2)

        result = run_model(
            tech_model, heat_load_mwt=100, hours_storage=26, temperature_hot=100
        )
        assert result["heat_annual"] == pytest.approx(232954184, rel=1e-2)
        assert result["electricity_annual"] == pytest.approx(4944004, rel=1e-2)
