import pytest

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port

from watertap_contrib.seto.solar_models.physical.flat_plate import IrradianceModelIsoSky

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)


class TestIrradianceModelIsoSky:
    @pytest.fixture(scope="class")
    def isosky_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.isosky = IrradianceModelIsoSky()

        return m

    @pytest.mark.unit
    def test_build(self, isosky_frame):
        m = isosky_frame

        assert len(m.fs.isosky.config) == 3
        assert not m.fs.isosky.config.dynamic
        assert not m.fs.isosky.config.has_holdup
        assert m.fs.isosky._tech_type == "irradiance_model_iso_sky"

        no_ports = list()
        for c in m.fs.isosky.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)
        assert len(no_ports) == 0
        assert number_variables(m.fs.isosky) == 13
        assert number_unused_variables(m.fs.isosky) == 2
        assert number_total_constraints(m.fs.isosky) == 11

    @pytest.mark.unit
    def test_tilted_radiation_model(self, isosky_frame):
        """Duffie and Beckman 4th Edition Example 2.15.1, page 90"""
        m = isosky_frame
        m.fs.isosky.phi.set_value(40)  # Latitude of collector
        m.fs.isosky.lon.set_value(89.4)  # Longitude of collector
        m.fs.isosky.std_meridian.set_value(
            90
        )  # Standard meridian corresponding to longitude
        m.fs.isosky.standard_time.set_value(9.695)  # Standard time, 9:41 AM
        m.fs.isosky.beta.set_value(60)  # Tilt angle of collector
        m.fs.isosky.rho_g.set_value(0.6)  # Ground reflectance
        m.fs.isosky.day_of_year.set_value(51)  # Day of year (Feb 20th)
        m.fs.isosky.G_bn.set_value(305.40)  # Beam normal radiation
        m.fs.isosky.G_d.set_value(796)  # Diffuse radiation on horizontal surface

        assert degrees_of_freedom(m) == 0

        solver = SolverFactory("ipopt")
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(-14.04, rel=1e-3) == value(m.fs.isosky.eqn_of_time)
        assert pytest.approx(9.501, rel=1e-3) == value(m.fs.isosky.solar_time)
        assert pytest.approx(-37.5, rel=1e-3) == value(m.fs.isosky.omega)
        assert pytest.approx(62.20, rel=1e-3) == value(m.fs.isosky.theta_z)
        assert pytest.approx(-11.58, rel=1e-3) == value(m.fs.isosky.delta)
        assert pytest.approx(1.7144, rel=1e-3) == value(m.fs.isosky.R_b)
        assert pytest.approx(36.97, rel=1e-3) == value(m.fs.isosky.theta)
        assert pytest.approx(243.9, rel=1e-3) == value(m.fs.isosky.G_b)
        assert pytest.approx(1040, rel=1e-3) == value(m.fs.isosky.G)
        assert pytest.approx(1171.24, rel=1e-3) == value(m.fs.isosky.G_T)
