#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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

from pyomo.environ import value

from watertap_contrib.reflo.flowsheets.KBHDP.components import (
    CST,
    DWI,
    EC,
    FPC,
    MEC,
    MD,
    LTMED,
    PV,
    ro_system,
    softener,
    UF,
)


class TestKBHDPComponents:

    @pytest.mark.component
    def test_CST_component(self):
        m = CST.main()
        assert pytest.approx(value(m.fs.cst.unit.heat), rel=1e-3) == 31237
        assert pytest.approx(value(m.fs.costing.LCOH), rel=1e-3) == 0.074473
        assert (
            pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 176170761
        )

    @pytest.mark.component
    def test_DWI_component(self):
        m = DWI.main()
        assert (
            pytest.approx(
                value(m.fs.DWI.unit.costing.variable_operating_cost), rel=1e-3
            )
            == 1671474
        )
        assert value(m.fs.costing.deep_well_injection.dwi_lcow) == 0.58

    @pytest.mark.component
    def test_EC_component(self):
        m = EC.main()
        assert pytest.approx(value(m.fs.EC.unit.conductivity), rel=1e-3) == 2.41627
        assert (
            pytest.approx(
                value(m.fs.costing.aggregate_flow_costs["aluminum"]), rel=1e-3
            )
            == 1248380
        )
        assert (
            pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 3281147
        )

    @pytest.mark.component
    def test_FPC_component(self):
        m = FPC.main()
        assert pytest.approx(value(m.fs.fpc.unit.heat), rel=1e-3) == 7039
        assert pytest.approx(value(m.fs.costing.LCOH), rel=1e-3) == 0.17035
        assert (
            pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 86666006
        )

    @pytest.mark.component
    def test_LTMED_component(self):
        m = LTMED.main()
        assert (
            pytest.approx(value(m.fs.LTMED.unit.thermal_power_requirement), rel=1e-3)
            == 18138
        )
        assert (
            pytest.approx(value(m.fs.LTMED.unit.gain_output_ratio), rel=1e-3) == 9.81539
        )
        assert (
            pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 19302723
        )

    @pytest.mark.component
    def test_MD_component(self):
        m = MD.main()
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.0945

    @pytest.mark.component
    def test_MEC_component(self):
        m = MEC.main()
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.66278
        assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 1.49032
        assert pytest.approx(value(m.fs.costing.SEC_th), rel=1e-3) == 180.7096

    @pytest.mark.component
    def test_PV_component(self):
        m = PV.main()
        assert pytest.approx(value(m.fs.pv.unit.electricity), rel=1e-3) == 258.3
        assert pytest.approx(value(m.fs.pv.unit.land_req), rel=1e-3) == 6.496
        assert pytest.approx(value(m.fs.costing.LCOE), rel=1e-3) == 0.0928

    @pytest.mark.component
    def test_ro_system_component(self):
        m = ro_system.main()
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 0.10502
        assert pytest.approx(value(m.fs.water_recovery), rel=1e-3) == 0.33857

    @pytest.mark.component
    def test_softener_component(self):
        m = softener.main()
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.00431
        assert pytest.approx(value(m.fs.softener.unit.CaO_dosing), rel=1e-3) == 7160.9

    @pytest.mark.component
    def test_UF_component(self):
        m = UF.main()
        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 0.0718
        assert pytest.approx(value(m.fs.UF.unit.electricity[0]), rel=1e-3) == 31.199
