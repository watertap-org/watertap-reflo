import pytest

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port

from watertap_contrib.seto.solar_models.physical.flat_plate import FlatPlatePhysical

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

class TestFlatPlatePhysical:
    @pytest.fixture(scope="class")
    def flatplate_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.flatplate = FlatPlatePhysical()

        return m

    @pytest.mark.unit
    def test_build(self, flatplate_frame):
        m = flatplate_frame

        assert len(m.fs.flatplate.config) == 3
        assert not m.fs.flatplate.config.dynamic
        assert not m.fs.flatplate.config.has_holdup
        assert m.fs.flatplate._tech_type == "flat_plate"

        no_ports = list()
        for c in m.fs.flatplate.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)
        assert len(no_ports) == 0
        assert number_variables(m.fs.flatplate) == 13
        assert number_unused_variables(m.fs.flatplate) == 2
        assert number_total_constraints(m.fs.flatplate) == 11

    # @pytest.mark.unit
    # def test_physical_variable_bounds(self, flatplate_frame):
    #     m = flatplate_frame
    #     assert m.fs.flatplate.heat_load.bounds == tuple([100, 1000])

    @pytest.mark.unit
    def test_tilted_radiation_model(self, flatplate_frame):
        """Duffie and Beckman 4th Edition Example 2.15.1, page 90"""
        m = flatplate_frame
        m.fs.flatplate.phi.set_value(40)                # Latitude of collector
        m.fs.flatplate.lon.set_value(89.4)              # Longitude of collector
        m.fs.flatplate.std_meridian.set_value(90)       # Standard meridian corresponding to longitude
        m.fs.flatplate.standard_time.set_value(9.695)   # Standard time, 9:41 AM
        m.fs.flatplate.beta.set_value(60)               # Tilt angle of collector
        m.fs.flatplate.rho_g.set_value(0.6)             # Ground reflectance
        m.fs.flatplate.day_of_year.set_value(51)        # Day of year (Feb 20th)
        m.fs.flatplate.G_bn.set_value(305.40)           # Beam normal radiation
        m.fs.flatplate.G_d.set_value(796)               # Diffuse radiation on horizontal surface

        assert degrees_of_freedom(m) == 0
        
        solver = SolverFactory("ipopt")
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(-14.04, rel=1e-3) == value(m.fs.flatplate.eqn_of_time)
        assert pytest.approx(9.501, rel=1e-3) == value(m.fs.flatplate.solar_time)
        assert pytest.approx(-37.5, rel=1e-3) == value(m.fs.flatplate.omega)
        assert pytest.approx(62.20, rel=1e-3) == value(m.fs.flatplate.theta_z)
        assert pytest.approx(-11.58, rel=1e-3) == value(m.fs.flatplate.delta)
        assert pytest.approx(1.7144, rel=1e-3) == value(m.fs.flatplate.R_b)
        assert pytest.approx(36.97, rel=1e-3) == value(m.fs.flatplate.theta)
        assert pytest.approx(243.9, rel=1e-3) == value(m.fs.flatplate.G_b)
        assert pytest.approx(1040, rel=1e-3) == value(m.fs.flatplate.G)
        assert pytest.approx(1171.24, rel=1e-3) == value(m.fs.flatplate.G_T)
