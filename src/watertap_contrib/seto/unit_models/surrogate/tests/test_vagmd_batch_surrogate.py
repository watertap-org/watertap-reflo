import pytest
from pyomo.environ import (
    ConcreteModel,
    Set,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.seto.unit_models.surrogate import VAGMDbatch_surrogate

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting
from idaes.core.util.testing import initialization_tester
from watertap.core.util.initialization import assert_no_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

import idaes.logger as idaeslog
import numpy as np

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestVAGMD:
    @pytest.fixture(scope="class")
    def VAGMD_frame(self):

        module_type = "AS7C1.5L"  # Module type
        recovery_ratio = 0.5  # Target recovery rate
        feed_salinity = 35  # Feed salinity (g/L)
        temp_evap_in = 80  # Evaporator inelt temperature (degC)
        feed_flow_rate = 600  # Feed flow rate (L/h)
        temp_cond_in = 25  # Condensor inlet temperature (degC)
        temp_tank = 25  # Initial temperature of the feed (degC)
        batch_volume = 50  # Initial batch volume

        # Calculate final salinity temporarily
        final_salinity = feed_salinity / (1 - recovery_ratio)

        # Estimate the initial membrane performance
        (
            permeate_flux_init,
            temp_cond_out_init,
            temp_evap_out_init,
            Area,
        ) = self._get_membrane_performance(
            temp_evap_in,
            feed_flow_rate,
            temp_cond_in,
            feed_salinity,
            temp_tank,
            module_type,
            final_salinity,
        )
        permeate_flow_rate_init = (
            permeate_flux_init * Area
        )  # Initial permeate flow rate (L/h)

        dt = 20352.55 / feed_flow_rate  # Time step (s) # TODO: move it to unit model
        V_init = permeate_flow_rate_init * dt / 3600  # Initial permeate volume (L)
        N = (
            int(batch_volume * recovery_ratio / V_init) + 2
        )  # TODO: Update model to calculate N

        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            dynamic=False,
        )
        m.fs.properties = SeawaterParameterBlock()
        m.fs.vagmd = VAGMDbatch_surrogate(
            property_package=m.fs.properties,
            module_type=module_type,
            high_brine_salinity=final_salinity > 175.3,
            cooling_system_type="closed",
            number_cycles=N + 1,
        )

        m.fs.vagmd.feed_props[0, 0].conc_mass_phase_comp["Liq", "TDS"].fix(
            feed_salinity
        )
        m.fs.vagmd.feed_props[0, 0].temperature.fix(temp_tank + 273.15)
        m.fs.vagmd.feed_props[0, 0].flow_vol_phase["Liq"].fix(
            feed_flow_rate / 3600 / 1000
        )  # TODO: convert unit here
        m.fs.vagmd.evaporator_in_props[0, 0].temperature.fix(temp_evap_in + 273.15)
        m.fs.vagmd.condenser_in_props[0, 0].temperature.fix(temp_cond_in + 273.15)

        if m.fs.vagmd.config.cooling_system_type == "open":
            m.fs.vagmd.TCoolIn.fix(TCoolIn)
        m.fs.vagmd.dt.fix(dt)
        m.fs.vagmd.initial_batch_volume.fix(batch_volume)

        return m

    @pytest.mark.unit
    def test_config(self, VAGMD_frame):
        m = VAGMD_frame
        # check unit config arguments
        assert len(m.fs.vagmd.config) == 8

        assert not m.fs.vagmd.config.dynamic
        assert not m.fs.vagmd.config.has_holdup
        assert m.fs.vagmd.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, VAGMD_frame):
        m = VAGMD_frame

        # test ports
        port_lst = ["feed", "permeate"]
        for port_str in port_lst:
            port = getattr(m.fs.vagmd, port_str)
            assert isinstance(port, Port)
            assert len(port.vars) == 3

        # test statistics
        assert number_variables(m) == 4850
        assert number_total_constraints(m) == 2521
        assert (
            number_unused_variables(m) == 2300
        )  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, VAGMD_frame):
        m = VAGMD_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_var_scaling(self, VAGMD_frame):
        m = VAGMD_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, VAGMD_frame):
        m = VAGMD_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_frame):
        m = VAGMD_frame
        N = m.fs.vagmd.config.number_cycles - 1

        # Print time series to check solution
        # print(
        #     "t  ",
        #     "Timestamp",
        #     "  PFlux",
        #     "   AccVd",
        #     "     TCO",
        #     "        TEO",
        #     "         ATml",
        #     "    ThEnergy",
        #     "     S",
        #     "      Ttank",
        #     "         R",
        #     "     ThPower",
        # )
        # for i in range(len(m.fs.vagmd.permeate_flux)):
        #     print(
        #         "{:<3} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} ".format(
        #             round(i),
        #             round(value(m.fs.vagmd.dt * i / 60, 4), 2),
        #             round(value(m.fs.vagmd.permeate_flux[i]), 4),
        #             round(value(m.fs.vagmd.acc_distillate_volume[i]), 4),
        #             round(
        #                 value(m.fs.vagmd.condenser_out_props[0,i].temperature - 273.15), 4
        #             ),
        #             round(
        #                 value(m.fs.vagmd.evaporator_out_props[0,i].temperature - 273.15),
        #                 4,
        #             ),
        #             round(value(m.fs.vagmd.log_mean_temp_dif[i]), 4),
        #             round(value(m.fs.vagmd.thermal_energy[i]), 4),
        #             round(
        #                 value(
        #                     m.fs.vagmd.feed_props[0,i].conc_mass_phase_comp["Liq", "TDS"]
        #                 ),
        #                 4,
        #             ),
        #             round(value(m.fs.vagmd.feed_props[0,i].temperature - 273.15), 4),
        #             round(value(m.fs.vagmd.acc_recovery_ratio[i]), 6),
        #             round(value(m.fs.vagmd.thermal_power[i]), 4),
        #         )
        #     )

        assert pytest.approx(0.503, rel=1e-3) == value(m.fs.vagmd.acc_recovery_ratio[N])

    """
    Calculation models (copied here for initializing the model and run test cases at the end)
    """

    def _get_membrane_performance(self, TEI, FFR, TCI, SgL, Ttank, k, Sf):
        # Model parameters
        PFluxAS26 = [
            0.798993148477908,
            0.314627216640160,
            0.559805181621833,
            -0.146236734128216,
            -0.659197144919924,
            0,
            0.185658514024503,
            -0.107221706014227,
            0,
            0,
            -0.187626469717738,
            0,
            0,
            0,
            0.128420664686447,
        ]
        PFluxAS7_high = [
            9.41014468300000,
            0,
            0,
            0,
            -0.0188989390000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        PFluxAS7_low = [
            4.82699491400000,
            1.37479276300000,
            1.91988177200000,
            -0.574212905000000,
            -0.641257664000000,
            0.399259954000000,
            0,
            0,
            0,
            0,
            0,
            0,
            -0.588924321000000,
            0,
            0,
        ]
        TCOAS26 = [
            65.1084685465240,
            9.15474718837607,
            -0.917460918908258,
            0.480070517276181,
            -1.06168979823129,
            0,
            0,
            0,
            0.142552983811052,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_high = [
            67.0068599900000,
            0,
            0,
            0,
            -0.0145469190000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TCOAS7_low = [
            58.8189910000000,
            8.79691366400000,
            -2.06741401400000,
            1.63187967600000,
            -0.914624645000000,
            -0.536574144000000,
            -0.249657477000000,
            0.398973861000000,
            -0.153760262000000,
            0.102355281000000,
            0.696768080000000,
            -0.300582958000000,
            -0.557410173000000,
        ]
        TEOAS26 = [
            29.2261439896435,
            0.569016152083381,
            0.824636694807529,
            4.62669502530487,
            1.37222105534565,
            0,
            0,
            0.220665657590258,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_high = [
            36.2497021800000,
            0,
            0,
            0,
            0.0126951860000000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        TEOAS7_low = [
            34.4309511000000,
            1.55140768500000,
            1.85928314600000,
            4.52887180500000,
            1.10791196800000,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0.808447211000000,
        ]
        Coder = [
            [-5.68487382500000, 0.0705622560000000, 0.000152146000000000],
            [-1.58460599600000, 0.00102338700000000, 1.20000000000000e-06],
            [-4.27697973100000, 0.175533630000000, -0.000178178000000000],
            [-1.49331349500000, 0.0146627780000000, 5.62000000000000e-06],
            [-7, 0.100000000000000, 0],
            [-2.14285714285714, 0.00285714285714286, 0],
            [-5, 0.200000000000000, 0],
            [-1.33333333333333, 0.00950000000000000, 0],
        ]
        a1 = 0.983930048493388
        a2 = -4.8359231959954e-04
        S_c = a1 * SgL + a2 * SgL**2  # [g/kg]

        CoderVars = [
            [1, TEI, TEI**2],
            [1, FFR, FFR**2],
            [1, TCI, TCI**2],
            [1, S_c, S_c**2],
        ]

        # Model calculations
        if k == "AS7C1.5L":
            A = 7.2  # Membrane Area [m2]
            if Sf > 175.3:
                TEI = 0
                FFR = 0
                TCI = 0
                S_r = S_c

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_high, TCOAS7_high, TEOAS7_high
            else:
                TEI = np.dot(CoderVars[0], Coder[0])
                FFR = np.dot(CoderVars[1], Coder[1])
                TCI = np.dot(CoderVars[2], Coder[2])
                S_r = np.dot(CoderVars[3], Coder[3])

                PFluxAS7, TCOAS7, TEOAS7 = PFluxAS7_low, TCOAS7_low, TEOAS7_low

            VarsAS7 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                TCI * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]
            VarsAS7_TCO = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                FFR * S_r,
                S_r * TCI,
                FFR**2,
                S_r**2,
                FFR**3,
            ]

            PFlux = np.dot(VarsAS7, PFluxAS7)
            TCO = np.dot(VarsAS7_TCO, TCOAS7)
            TEO = np.dot(VarsAS7, TEOAS7)

        else:
            A = 25.92
            TEI = np.dot(CoderVars[0], Coder[4])
            FFR = np.dot(CoderVars[1], Coder[5])
            TCI = np.dot(CoderVars[2], Coder[6])
            S_r = np.dot(CoderVars[3], Coder[7])

            VarsAS26 = [
                1,
                TEI,
                FFR,
                TCI,
                S_r,
                TCI * TEI,
                FFR * TEI,
                S_r * TEI,
                FFR * TCI,
                S_r * TCI,
                FFR * S_r,
                TEI**2,
                FFR**2,
                TCI**2,
                S_r**2,
            ]

            PFlux = np.dot(VarsAS26, PFluxAS26)
            TCO = np.dot(VarsAS26, TCOAS26)
            TEO = np.dot(VarsAS26, TEOAS26)

        return [PFlux, TCO, TEO, A]
