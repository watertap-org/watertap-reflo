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

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import FlowsheetBlock
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.unit_models.zero_order.forward_osmosis_zo import (
    ForwardOsmosisZO,
)
from watertap_contrib.reflo.property_models.fo_draw_solution_properties import (
    FODrawSolutionParameterBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestFO:
    @pytest.fixture(scope="class")
    def FO_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.water_prop = SeawaterParameterBlock()
        m.fs.draw_solution_prop = FODrawSolutionParameterBlock()
        m.fs.fo = ForwardOsmosisZO(
            property_package_water=m.fs.water_prop,
            property_package_draw_solution=m.fs.draw_solution_prop,
        )

        fo = m.fs.fo
        strong_draw = fo.strong_draw_props[0]
        product = fo.product_props[0]

        # System specifications
        recovery_ratio = 0.3  # Assumed FO recovery ratio
        nanofiltration_recovery_ratio = 0.8  # Nanofiltration recovery ratio
        dp_brine = 0  # Required pressure over brine osmotic pressure (Pa)
        heat_mixing = 75.6  # Heat of mixing in the membrane (MJ/m3 product)
        reneration_temp = 90  # Separation temperature of the draw solution (C)
        separator_temp_loss = 1  # Temperature loss in the separator (K)
        feed_temperature = 13  # Feed water temperature (C)
        feed_vol_flow = 3.704  # Feed water volumetric flow rate (m3/s)
        feed_TDS_mass = 0.035  # TDS mass fraction of feed
        strong_draw_temp = 20  # Strong draw solution inlet temperature (C)
        strong_draw_mass_frac = 0.8  # Strong draw solution mass fraction
        product_draw_mas_frac = 0.01  # Mass fraction of draw in the product water

        fo.recovery_ratio.fix(recovery_ratio)
        fo.nanofiltration_recovery_ratio.fix(nanofiltration_recovery_ratio)
        fo.dp_brine.fix(dp_brine)
        fo.heat_mixing.fix(heat_mixing)
        fo.regeneration_temp.fix(reneration_temp + 273.15)
        fo.separator_temp_loss.fix(separator_temp_loss)

        # Specify strong draw solution properties
        fo.strong_draw_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_vol_flow,
                (
                    "mass_frac_phase_comp",
                    ("Liq", "DrawSolution"),
                ): strong_draw_mass_frac,
                ("temperature", None): strong_draw_temp + 273.15,
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        strong_draw.flow_mass_phase_comp["Liq", "DrawSolution"].unfix()

        # Specify product water properties
        fo.product_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_vol_flow,
                (
                    "mass_frac_phase_comp",
                    ("Liq", "DrawSolution"),
                ): product_draw_mas_frac,
                ("temperature", None): reneration_temp - separator_temp_loss + 273.15,
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        product.flow_mass_phase_comp["Liq", "H2O"].unfix()
        product.temperature.unfix()

        # Specify feed properties
        fo.feed_props.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): feed_vol_flow,
                ("mass_frac_phase_comp", ("Liq", "TDS")): feed_TDS_mass,
                ("temperature", None): feed_temperature + 273.15,
                ("pressure", None): 101325,
            },
            hold_state=True,
        )

        # Set scaling factors for mass flow rates
        m.fs.water_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
        )
        m.fs.water_prop.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "TDS")
        )
        m.fs.draw_solution_prop.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.draw_solution_prop.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "DrawSolution")
        )

        return m

    @pytest.mark.unit
    def test_build(self, FO_frame):
        m = FO_frame

        # test ports
        port_lst = ["feed", "brine", "strong_draw", "weak_draw", "reg_draw", "product"]
        for port_str in port_lst:
            port = getattr(m.fs.fo, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 170
        assert number_total_constraints(m) == 56
        assert number_unused_variables(m) == 70  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, FO_frame):
        m = FO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, FO_frame):
        m = FO_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, FO_frame):
        m = FO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_initialize(self, FO_frame):
        m = FO_frame
        initialization_tester(m, unit=m.fs.fo, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_solve(self, FO_frame):
        m = FO_frame

        # Unfix the state variables and fix mass fraction of draw solution state blocks
        strong_draw_mass = 0.8  # Strong draw solution mass fraction
        product_draw_mass = 0.01  # Mass fraction of draw in the product water
        m.fs.fo.unfix_and_fix_freedom(strong_draw_mass, product_draw_mass)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, FO_frame):
        m = FO_frame
        strong_draw = m.fs.fo.strong_draw_props[0]
        weak_draw = m.fs.fo.weak_draw_props[0]
        brine = m.fs.fo.brine_props[0]
        product = m.fs.fo.product_props[0]
        reg_draw = m.fs.fo.reg_draw_props[0]

        assert value(
            strong_draw.flow_mass_phase_comp["Liq", "DrawSolution"]
        ) == pytest.approx(1114.80, rel=1e-3)
        assert value(strong_draw.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            278.70, rel=1e-3
        )
        assert value(strong_draw.flow_vol_phase["Liq"]) == pytest.approx(
            1.284, rel=1e-3
        )
        assert value(
            weak_draw.flow_mass_phase_comp["Liq", "DrawSolution"]
        ) == pytest.approx(1114.80, rel=1e-3)
        assert value(weak_draw.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            1379.35, rel=1e-3
        )
        assert value(weak_draw.flow_vol_phase["Liq"]) == pytest.approx(2.349, rel=1e-3)
        assert value(weak_draw.temperature) == pytest.approx(293.25, rel=1e-3)
        assert value(
            weak_draw.mass_frac_phase_comp["Liq", "DrawSolution"]
        ) == pytest.approx(0.447, rel=1e-3)
        assert value(brine.mass_frac_phase_comp["Liq", "TDS"]) == pytest.approx(
            0.04956, rel=1e-3
        )
        assert value(brine.conc_mass_phase_comp["Liq", "TDS"]) == pytest.approx(
            51.32, rel=1e-3
        )
        assert value(brine.pressure_osm_phase["Liq"]) == pytest.approx(
            3689472, rel=1e-3
        )
        assert value(brine.temperature) == pytest.approx(293.25, rel=1e-3)
        assert value(product.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            1377.37, rel=1e-3
        )
        assert value(
            product.flow_mass_phase_comp["Liq", "DrawSolution"]
        ) == pytest.approx(13.913, rel=1e-3)
        assert value(product.flow_vol_phase["Liq"]) == pytest.approx(1.389, rel=1e-3)
        assert value(reg_draw.flow_mass_phase_comp["Liq", "H2O"]) == pytest.approx(
            278.80, rel=1e-3
        )
        assert value(
            reg_draw.flow_mass_phase_comp["Liq", "DrawSolution"]
        ) == pytest.approx(1114.80, rel=1e-3)
        assert value(m.fs.fo.delta_temp_membrane) == pytest.approx(5.90, rel=1e-3)
        assert value(m.fs.fo.membrane_temp) == pytest.approx(287.351, rel=1e-3)
        assert value(m.fs.fo.heat_transfer_to_weak_draw) == pytest.approx(
            30.78, rel=1e-3
        )
        assert value(m.fs.fo.heat_transfer_to_brine) == pytest.approx(44.82, rel=1e-3)
