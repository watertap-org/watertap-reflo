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
        assert number_variables(m.fs.flatplate) == 6
        assert number_unused_variables(m.fs.flatplate) == 2
        assert number_total_constraints(m.fs.flatplate) == 4

    @pytest.mark.unit
    def test_flat_plate_physical_model(self, flatplate_frame):
        """docstring"""
        m = flatplate_frame
        m.fs.flatplate.area_coll.set_value(2.98)
        m.fs.flatplate.FRta.set_value(0.689)
        m.fs.flatplate.FRUL.set_value(3.85)
        m.fs.flatplate.iam.set_value(0.2)
        m.fs.flatplate.mdot_test.set_value(0.045528)
        m.fs.flatplate.cp_test.set_value(3400)  # specific heat of glycol [J/kg-K]
        m.fs.flatplate.cp_use.set_value(3400)   # specific heat of glycol [J/kg-K]
        m.fs.flatplate.ncoll.set_value(2)
        m.fs.flatplate.pump_watts.set_value(45)
        m.fs.flatplate.pump_eff.set_value(0.85)
        m.fs.flatplate.T_amb.set_value(12)      # default SAM model at noon on Jan. 1
        m.fs.flatplate.T_in.set_value(38.2)     # default SAM model at noon on Jan. 1
        m.fs.flatplate.G_trans.set_value(540)   # default SAM model at noon on Jan. 1

        assert degrees_of_freedom(m) == 0

        solver = SolverFactory("ipopt")
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(4.00, rel=1e-3) == value(m.fs.flatplate.FprimeUL)
        assert pytest.approx(1.00, rel=1e-3) == value(m.fs.flatplate.r)
        assert pytest.approx(1616.29, rel=1e-3) == value(m.fs.flatplate.Q_useful)   # [W], with no pipe or heat exchanger losses
        assert pytest.approx(52.94, rel=1e-3) == value(m.fs.flatplate.P_pump)       # [W]
