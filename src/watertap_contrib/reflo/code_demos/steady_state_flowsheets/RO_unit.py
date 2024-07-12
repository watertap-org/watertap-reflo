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
from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
)
from pyomo.network import Arc

# Imports from IDAES
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state

# Import function to check degrees of freedom
from idaes.core.util.model_statistics import degrees_of_freedom

# Import utility function for calculating scaling factors
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.models.unit_models import Product, Feed

# Imports from WaterTAP
from watertap.core.solvers import get_solver
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
)


def create_RO_unit():
    # Create a Pyomo concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.pump = Pump(property_package=m.fs.properties)

    # Add an RO unit to the flowsheet.
    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_pressure_change=False,
    )

    m.fs.erd = EnergyRecoveryDevice(property_package=m.fs.properties)

    m.fs.a1 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.a2 = Arc(source=m.fs.pump.outlet, destination=m.fs.RO.inlet)
    m.fs.a3 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
    m.fs.a4 = Arc(source=m.fs.RO.retentate, destination=m.fs.erd.inlet)
    m.fs.a5 = Arc(source=m.fs.erd.outlet, destination=m.fs.disposal.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_RO_operating_conditions(m, flow_in=1e-2, conc_in=30, water_recovery=0.5):
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): flow_in,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "NaCl")): conc_in * 1e-3,
            ("pressure", None): 101325,
            ("temperature", None): 298.15,
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    m.fs.pump.control_volume.properties_out[0].pressure.fix(50e5)
    m.fs.pump.efficiency_pump.fix(0.8)

    m.fs.RO.area.fix(50)  # membrane area (m^2)
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability (m/Pa/s)
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability (m/s)
    m.fs.RO.permeate.pressure[0].fix(101325)  # permeate pressure (Pa)


def set_RO_scaling(m):
    # Set scaling factors for component mass flowrates.
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    # Set scaling factor for membrane area.
    set_scaling_factor(m.fs.RO.area, 1e-2)

    # Calculate scaling factors for all other variables.
    calculate_scaling_factors(m)


def initialize_RO(m):
    print(degrees_of_freedom(m))
    m.fs.feed.initialize()
    propagate_state(m.fs.a1)
    m.fs.pump.initialize()
    propagate_state(m.fs.a2)
    m.fs.RO.initialize()
    propagate_state(m.fs.a4)
    m.fs.erd.efficiency_pump.fix(0.95)
    m.fs.erd.control_volume.properties_out[0].pressure.fix(101325)
    m.fs.erd.initialize()
    propagate_state(m.fs.a3)
    propagate_state(m.fs.a5)


def solve_RO(m):
    # Setup solver
    solver = get_solver()
    simulation_results = solver.solve(m)


def print_RO_results(m):
    m.fs.RO.report()
    print(m.fs.RO)


def main():
    m = create_RO_unit()
    set_RO_operating_conditions(m)
    set_RO_scaling(m)
    initialize_RO(m)
    solve_RO(m)
    print_RO_results(m)


if __name__ == "__main__":
    main()
