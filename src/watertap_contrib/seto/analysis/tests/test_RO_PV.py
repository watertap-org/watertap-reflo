import pytest
import os
from pyomo.environ import Var, assert_optimal_termination, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
import watertap_contrib.seto.analysis.net_metering.PV_RO as PV_RO
from watertap_contrib.seto.core import SETODatabase, PySAMWaterTAP
from idaes.core.util.testing import initialization_tester
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

import idaes.core.util.model_statistics as stats

from idaes.core.util.scaling import (
    unscaled_variables_generator,
    unscaled_constraints_generator,
    calculate_scaling_factors,
)
import idaes.logger as idaeslog

@pytest.mark.unit
def test_config():
    m = PV_RO.build_ro_pv()

    assert hasattr(m.fs, "energy")
    assert hasattr(m.fs, "treatment")
    
    for component in ['feed','product','disposal','p1', 'ro', 'erd']:
        assert hasattr(m.fs.treatment, component)
        unit = getattr(m.fs.treatment, component)
        assert not unit.config.dynamic
        assert not unit.config.has_holdup
        if component == 'feed':
            assert unit.properties[0].flow_vol_phase
        if component == 'ro':
            assert unit.config.property_package is m.fs.properties
            assert unit.config.material_balance_type == MaterialBalanceType.useDefault
            assert unit.config.energy_balance_type == EnergyBalanceType.useDefault
            assert unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
            assert unit.config.pressure_change_type == PressureChangeType.calculated
            assert unit.config.concentration_polarization_type == ConcentrationPolarizationType.calculated
            assert unit.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
            assert unit.config.has_pressure_change
            assert isinstance(unit.feed_side.deltaP, Var)
            assert isinstance(unit.deltaP, Var)

@pytest.mark.unit
def test_set_operating_conditions():
    m = PV_RO.build_ro_pv()
    PV_RO.set_operating_conditions(m)

    assert stats.number_unused_variables(m) <= 3
    assert stats.number_fixed_variables(m.fs.treatment.p1.control_volume.properties_out[0]) == 1
    assert stats.number_total_blocks(m.fs.treatment) == 29
    assert stats.number_total_blocks(m.fs.energy) == 3

@pytest.mark.unit
def test_initialization():
    m = PV_RO.build_ro_pv()
    PV_RO.set_operating_conditions(m)
    solver = get_solver()
    optarg = solver.options
    m.fs.treatment.feed.initialize(optarg=optarg)
    PV_RO.initialize_treatment(m)
    PV_RO.initialize_energy(m)

    assert stats.degrees_of_freedom(m) == 2
    assert stats.number_variables(m) == 212
    assert stats.number_total_constraints(m) == 174
    assert stats.number_unused_variables(m) == 3

    for component in ['feed','product','disposal','p1', 'ro', 'erd']:
        assert hasattr(m.fs.treatment, component)
        unit = getattr(m.fs.treatment, component)
        assert_units_consistent(unit)

    initialization_tester(m, unit=m.fs.treatment.ro, dof=2, outlvl=idaeslog.DEBUG)

@pytest.mark.component
def test_solve():
    m = PV_RO.build_ro_pv()
    PV_RO.set_operating_conditions(m)
    PV_RO.initialize_sys(m)
    PV_RO.add_costing(m)
    PV_RO.fix_pv_costing(m)
    PV_RO.fix_treatment_global_params(m)
    PV_RO.optimize_setup(m, m.fs.sys_costing.LCOW)
    results = PV_RO.solve(m)

    # Check for optimal solution
    assert_optimal_termination(results)