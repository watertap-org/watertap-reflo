import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap_contrib.seto.multiperiod.VAGMD_batch_flowsheet_multiperiod import (
    create_multiperiod_vagmd_batch_model,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog

solver = get_solver()

from watertap.core.util.model_diagnostics.infeasible import *


class TestVAGMDbatch:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame(self):
        mp = create_multiperiod_vagmd_batch_model(
            # n_time_points=34,
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity= 50,
            recovery_ratio=0.72,
            initial_batch_volume=50,
            module_type="AS7C1.5L",
            cooling_system_type="open",
            cooling_inlet_temp = 25, # not required if cooling system type is "closed"
        )

        return mp

    @pytest.mark.component
    def test_solve(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        print("multi dof", degrees_of_freedom(mp))
        assert degrees_of_freedom(mp) == 0
        results = solver.solve(mp)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        blks = mp.get_active_process_blocks()
        print_close_to_bounds(mp)
        print_infeasible_constraints(mp)
        print(len(blks)) 
        print(
            "{:<3}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}".format(
                "t  ",
                "t_minute",
                "S",
                "Pflux",
                "AccVd",
                "Ttank",
                "TEO",
                "TCO",
                "RR",
                "ThPower",
                "AccThE",
                "STEC",
                "GOR",
                "Cpower",
                "CEnergy",
                "TCI",
                "TCoolIn",
                "TCoolOut",
            )
        )
        for i in range(len(blks)):
            print(
                "{:<3}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}".format(
                    round(i),
                    round(value(blks[i].fs.dt) * i / 60, 4),
                    round(
                        value(
                            blks[i]
                            .fs.vagmd.feed_props[0]
                            .conc_mass_phase_comp["Liq", "TDS"]
                        ),
                        4,
                    ),
                    round(value(blks[i].fs.vagmd.permeate_flux), 4),
                    round(value(blks[i].fs.acc_distillate_volume), 4),
                    round(
                        value(blks[i].fs.vagmd.feed_props[0].temperature - 273.15), 4
                    ),
                    round(
                        value(
                            blks[i].fs.vagmd.evaporator_out_props[0].temperature
                            - 273.15
                        ),
                        4,
                    ),
                    round(
                        value(
                            blks[i].fs.vagmd.condenser_out_props[0].temperature - 273.15
                        ),
                        4,
                    ),
                    round(value(blks[i].fs.acc_recovery_ratio) * 100, 4),
                    round(value(blks[i].fs.vagmd.thermal_power), 4),
                    round(value(blks[i].fs.acc_thermal_energy), 4),
                    round(value(blks[i].fs.specific_energy_consumption_thermal), 4),
                    round(value(blks[i].fs.gain_output_ratio), 4),
                    round(value(blks[i].fs.vagmd.cooling_power_thermal), 4),
                    round(value(blks[i].fs.acc_cooling_energy), 4),
                    round(
                        value(blks[i].fs.vagmd.condenser_in_props[0].temperature)
                        - 273.15,
                        4,
                    ),
                    round(
                        value(blks[i].fs.vagmd.cooling_in_props[0].temperature)
                        - 273.15,
                        4,
                    ),
                    round(
                        value(blks[i].fs.vagmd.cooling_out_props[0].temperature)
                        - 273.15,
                        4,
                    ),
                )
            )

        # Check final recovery rate
        assert pytest.approx(0.502548, rel=1e-3) == value(
            blks[-1].fs.acc_recovery_ratio
        )

        assert False
