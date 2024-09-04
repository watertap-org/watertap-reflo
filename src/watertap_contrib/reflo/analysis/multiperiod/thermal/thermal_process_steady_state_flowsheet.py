from pyomo.environ import (
    Var,
    ConcreteModel,
    Objective,
    Param,
    Expression,
    Constraint,
    Block,
    log10,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
    NonNegativeReals,
)

from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.util.model_statistics import *
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.solvers.get_solver import get_solver
from idaes.models.unit_models import Product, Feed
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import (
    set_scaling_factor,
    calculate_scaling_factors,
    constraint_scaling_transform,
)
from idaes.core import UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog

from watertap_contrib.reflo.unit_models.thermal_energy_storage import ThermalEnergyStorage
from watertap_contrib.reflo.solar_models.zero_order.flat_plate_physical import FlatPlatePhysical
from watertap_contrib.reflo.core import SolarModelType

from idaes.models.unit_models.heat_exchanger import(  
    HeatExchanger,
    delta_temperature_lmtd_callback,
    delta_temperature_amtd_callback,
    HeatExchangerFlowPattern,

)
# Using a slightly modified version tp use the same property package across all the other models
from watertap_contrib.reflo.unit_models.heat_exchanger_ntu import (
    HeatExchangerNTU, 
    HXNTUInitializer,
)
from watertap.core.util.model_diagnostics.infeasible import *

from idaes.models.unit_models import Heater

from watertap.property_models.water_prop_pack import WaterParameterBlock

from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
    list_badly_scaled_variables,
    constraint_autoscale_large_jac
)


def propagate_state(arc):
    _prop_state(arc)


def build_thermal_flowsheet(m = None):

    '''
    This function adds the relevant components to the flowsheet
    '''

    m.fs = FlowsheetBlock(dynamic = False)
    m.fs.properties = WaterParameterBlock()

    # Include FPC model
    m.fs.fpc =  FlatPlatePhysical(
        property_package = m.fs.properties, solar_model_type=SolarModelType.physical
    )

    # Include TES model
    m.fs.tes = ThermalEnergyStorage(
        property_package = m.fs.properties
    )

    # HX between FPC and TES
    m.fs.hx_solar = HeatExchangerNTU(
                hot_side_name = 'fpc',
                cold_side_name = 'tes',
                fpc = {"property_package": m.fs.properties},
                tes = {"property_package": m.fs.properties}
        )

    # Grid heater to meet process inlet temperature set point
    m.fs.grid_heater = Heater(
        property_package = m.fs.properties)
    
    return m
    
    
def create_coupling_variables(blk):
    # Create coupling variables

    blk.fs.previous_hx_solar_hot_outlet_temperature = Var(
        domain = NonNegativeReals,
        initialize = 35+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc = 'Outlet temperature from HX solar hot side from the previous time step'
    )

    blk.fs.previous_fpc_outlet_temperature = Var(
        domain = NonNegativeReals,
        initialize = 40+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc = 'Outlet temperature from FPC from the previous time step'
    )

    blk.fs.previous_tes_tank_temp = Var(
        domain = NonNegativeReals,
        initialize = 30+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc= 'Temperature of the thermal storage tank from the previous time step'
    )

    blk.fs.previous_hx_solar_cold_outlet_temperature = Var(
        domain = NonNegativeReals,
        initialize = 40+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc = 'Outlet temperature from HX solar cold side from the previous time step'
    )


    blk.fs.previous_process_outlet_temperature = Var(
        domain = NonNegativeReals,
        initialize = 35+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc = 'Outlet temperature from process from the previous time step'
    )

    blk.fs.previous_grid_duty = Var(
        initialize = 0,
        units = pyunits.W,
        doc = 'Grid heat duty from previous step'
    )

    blk.fs.acc_grid_duty = Var(
        initialize = 0,
        units = pyunits.W,
        doc = 'Accumulated heat duty at each step'
    )

    blk.fs.previous_acc_grid_duty = Var(
        initialize = 0,
        units = pyunits.W,
        doc = 'Accumulated heat duty from previous step'
    )


def create_feed_streams(m,
                        mass_fr_fpc = 0.05,
                        mass_fr_tes_hx_solar=0.1,
                        mass_fr_tes_process=0.05,
                        ):

    # Creating state blocks for the outlet stream from TES with the TES tank temperature from the previous step

    # Outlet stream from the FPC going back to the solar HX hot inlet
    m.fs.fpc_outlet = Feed(property_package = m.fs.properties)
    m.fs.fpc_outlet.properties[0].temperature.fix(m.fs.previous_fpc_outlet_temperature())
    m.fs.fpc_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_fpc)
    m.fs.fpc_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.fpc_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from the TES going back to the solar HX cold inlet
    m.fs.tes_hx_outlet = Feed(property_package = m.fs.properties)
    m.fs.tes_hx_outlet.properties[0].temperature.fix(m.fs.previous_tes_tank_temp)
    m.fs.tes_hx_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_hx_solar)
    m.fs.tes_hx_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.tes_hx_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from the solar HX hot outlet to the FPC inlet
    m.fs.hx_solar_hot_outlet = Feed(property_package = m.fs.properties)
    m.fs.hx_solar_hot_outlet.properties[0].temperature.fix(m.fs.previous_hx_solar_hot_outlet_temperature())
    m.fs.hx_solar_hot_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_fpc)
    m.fs.hx_solar_hot_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.hx_solar_hot_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from the solar HX cold outlet to the TES
    m.fs.hx_solar_cold_outlet = Feed(property_package = m.fs.properties)
    m.fs.hx_solar_cold_outlet.properties[0].temperature.fix(m.fs.previous_hx_solar_cold_outlet_temperature())
    m.fs.hx_solar_cold_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_hx_solar)
    m.fs.hx_solar_cold_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.hx_solar_cold_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from process outlet to the TES
    m.fs.process_outlet = Feed(property_package = m.fs.properties)
    m.fs.process_outlet.properties[0].temperature.fix(m.fs.previous_process_outlet_temperature())
    m.fs.process_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_hx_solar)
    m.fs.process_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.process_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from the TES going back to the grid heater inlet
    m.fs.tes_process_outlet = Feed(property_package = m.fs.properties)
    m.fs.tes_process_outlet.properties[0].temperature.fix(m.fs.previous_tes_tank_temp())
    m.fs.tes_process_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_process)
    m.fs.tes_process_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.tes_process_outlet.properties[0].pressure.fix(101325)

    return m


def create_arcs(m):

    m.fs.fpc_hx_solar = Arc(source = m.fs.fpc_outlet.outlet, 
                            destination = m.fs.hx_solar.hot_side_inlet,
                            doc = 'Connect FPC outlet to the solar HX hotside inlet')

    m.fs.hx_solar_fpc = Arc(source = m.fs.hx_solar_hot_outlet.outlet,
                            destination = m.fs.fpc.inlet,
                            doc = 'Connect solar hx hot side outlet back to fpc')

    m.fs.tes_hx_solar = Arc(source = m.fs.tes_hx_outlet.outlet,
                            destination = m.fs.hx_solar.cold_side_inlet,
                            doc = 'Connect TES from previous time step to the solar HX cold side inlet')

    m.fs.hx_solar_tes = Arc(source = m.fs.hx_solar_cold_outlet.outlet,
                            destination = m.fs.tes.tes_hx_inlet,
                            doc = 'Connect solar HX cold side outlet back to TES inlet')

    m.fs.tes_gridHtr = Arc(source = m.fs.tes_process_outlet.outlet,
                            destination = m.fs.grid_heater.inlet,
                            doc = 'Connect TES to the grid heater inlet')


    m.fs.process_tes = Arc(source = m.fs.process_outlet.outlet,
                            destination = m.fs.tes.tes_process_inlet,
                            doc = 'Connect process outlet to the TES')

    TransformationFactory("network.expand_arcs").apply_to(m)


def fix_dof_and_initialize(
    blk,
    dt,
    GHI,
    mass_fr_tes_hx_solar, 
    mass_fr_tes_process,
    process_inlet_temp,
    outlvl=idaeslog.WARNING,
):
    """Fix degrees of freedom and initialize the flowsheet

    This function fixes the degrees of freedom of each unit and initializes the entire flowsheet.

    Args:
        m: Pyomo `Block` or `ConcreteModel` containing the flowsheet
        outlvl: Logger (default: idaeslog.WARNING)
    """

    solver = get_solver()
    optarg = solver.options

    # blk.fs.previous_hx_solar_hot_outlet_temperature.fix(41+273.15)
    # blk.fs.previous_fpc_outlet_temperature.fix(41+273.15)
    # blk.fs.previous_tes_tank_temp.fix(39+273.15)
    # blk.fs.previous_hx_solar_cold_outlet_temperature.fix(41+273.15)
    # blk.fs.previous_process_outlet_temperature.fix(36+273.15)
    blk.fs.previous_grid_duty.fix(10)
    blk.fs.previous_acc_grid_duty.fix(0)

    blk.fs.fpc_outlet.properties[0].temperature.fix(blk.fs.previous_fpc_outlet_temperature())
    blk.fs.tes_hx_outlet.properties[0].temperature.fix(blk.fs.previous_tes_tank_temp())
    blk.fs.hx_solar_hot_outlet.properties[0].temperature.fix(blk.fs.previous_hx_solar_hot_outlet_temperature())
    blk.fs.hx_solar_cold_outlet.properties[0].temperature.fix(blk.fs.previous_hx_solar_cold_outlet_temperature())
    blk.fs.process_outlet.properties[0].temperature.fix(blk.fs.previous_process_outlet_temperature())
    blk.fs.tes_process_outlet.properties[0].temperature.fix(blk.fs.previous_tes_tank_temp())

    # Initializing and propagate the FPC
    # FPC
    # blk.fs.previous_hx_solar_hot_outlet_temperature.fix()
    blk.fs.hx_solar_hot_outlet.initialize(optarg=optarg)

    propagate_state(blk.fs.hx_solar_fpc)
    # m.fs.fpc.inlet.fix()
    blk.fs.fpc.total_irradiance.fix(GHI)
    blk.fs.fpc.collector_area.fix(2)
    blk.fs.fpc.outlet.pressure.fix(101325)
    blk.fs.fpc.initialize()

    # Initialize and propragate for solar HX
    # blk.fs.previous_fpc_outlet_temperature.fix() 
    blk.fs.fpc_outlet.initialize()
    propagate_state(blk.fs.fpc_hx_solar)

    # blk.fs.previous_tes_tank_temp.fix()
    blk.fs.tes_hx_outlet.initialize()
    propagate_state(blk.fs.tes_hx_solar)

    blk.fs.hx_solar.effectiveness.fix(0.7)
    blk.fs.hx_solar.hot_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    blk.fs.hx_solar.cold_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)

    blk.fs.hx_solar.initialize_build()

    # Initialize and propagate from TES
    # blk.fs.previous_hx_solar_cold_outlet_temperature.fix()
    blk.fs.hx_solar_cold_outlet.initialize()
    propagate_state(blk.fs.hx_solar_tes)

    # blk.fs.previous_process_outlet_temperature.fix()
    blk.fs.process_outlet.initialize()
    propagate_state(blk.fs.process_tes)

    blk.fs.tes.dt.fix(dt)
    # m.fs.tes.tes_hx_inlet.fix()
    blk.fs.tes.tes_initial_temperature.fix(blk.fs.previous_tes_tank_temp())
    blk.fs.tes.tes_volume.fix(10)
    blk.fs.tes.hours_storage.fix(6)
    # m.fs.tes.heat_load.fix(0.5)

    # m.fs.tes.tes_process_inlet.temperature.fix(20 + 273.15)
    # m.fs.tes.tes_process_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(mass_fr_tes_process)
    # m.fs.tes.tes_process_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    # m.fs.tes.tes_process_inlet.pressure.fix(101325)

    # Fix outlet vapor flow to be 0
    blk.fs.tes.tes_hx_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    blk.fs.tes.tes_hx_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(mass_fr_tes_hx_solar)
    blk.fs.tes.tes_hx_outlet.pressure.fix()
    # m.fs.tes.tes_hx_outlet.temperature.fix(m.fs.previous_tes_tank_temp())

    blk.fs.tes.tes_process_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    blk.fs.tes.tes_process_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(mass_fr_tes_process)
    blk.fs.tes.tes_process_outlet.pressure.fix()
    # m.fs.tes.tes_process_outlet.temperature.fix(m.fs.previous_tes_tank_temp())
 
    blk.fs.tes.initialize()

    # Grid Heater
    # m.fs.grid_heater.outlet.temperature[0].setub(99+273.15)
    propagate_state(blk.fs.tes_gridHtr)
    blk.fs.tes_process_outlet.initialize()
    blk.fs.grid_heater.outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    # m.fs.grid_heater.outlet.pressure.fix()
    blk.fs.grid_heater.outlet.temperature[0].fix(process_inlet_temp + 273.15)

    blk.fs.grid_heater.initialize()



def create_mp_steady_state(
        m= None,
        GHI = 900, 
        elec_price = 0.07,
        mass_fr_fpc = 0.05,
        mass_fr_tes_hx_solar = 0.1,
        mass_fr_tes_process = 0.05,
):
    
    if m is None:
        m = ConcreteModel()

    m = build_thermal_flowsheet(m)

    create_coupling_variables(m)
    create_feed_streams(m, mass_fr_fpc, mass_fr_tes_hx_solar, mass_fr_tes_process)
    create_arcs(m)

    # Constraints on the TES output so that the TES tank temperature is the same as the temperature to process and solar hx

    @m.Constraint(doc='Temperature to the process and solar HX are the same')
    def eq_hx_process(b,t):
        return b.fs.tes.tes_process_outlet.temperature[0] == b.fs.tes.tes_hx_outlet.temperature[0]

    @m.Constraint(doc='Calculated TES temperature is the same as the outlet to solar HX')
    def eq_temp_process(b,t):
        return b.fs.tes.tes_temperature[0] == b.fs.tes.tes_hx_outlet.temperature[0]
    
    @m.Constraint(doc='Calculating the accumulated grid heat duty')
    def eq_acc_grid_duty(b):
        return b.fs.acc_grid_duty == b.fs.previous_grid_duty + b.fs.previous_acc_grid_duty

    return m


def print_results(m):

    print('\nTank Variables')
    print('Tank temperature:', value(m.fs.tes.tes_temperature[0])-273.15)
    print('Temperature of stream exiting to grid heater:',m.fs.tes.tes_process_outlet.temperature[0]()-273.15)
    print('Temperature of stream exiting to solar hx:',m.fs.tes.tes_hx_outlet.temperature[0]()-273.15)

    # FPC
    print('\nFlat plate collector')
    print('FPC inlet temperature:',m.fs.fpc.inlet.temperature[0].value-273.15)
    print('FPC outlet temperature:',m.fs.fpc.outlet.temperature[0]()-273.15)

    # Solar HX
    print('\nSolar HX')
    print('Solar HX hot side inlet:',m.fs.hx_solar.hot_side_inlet.temperature[0].value - 273.15)
    print('Solar HX hot side outlet:',m.fs.hx_solar.hot_side_outlet.temperature[0].value - 273.15)
    print('Solar HX cold side inlet:',m.fs.hx_solar.cold_side_inlet.temperature[0].value - 273.15)
    print('Solar HX cold side outlet:',m.fs.hx_solar.cold_side_outlet.temperature[0].value - 273.15)

    # Grid inlet temperature
    print('\nGrid powered heater')
    print('Grid heater inlet temperature:',m.fs.grid_heater.inlet.temperature[0]() - 273.15)
    print('Grid heater outlet temperature:',m.fs.grid_heater.outlet.temperature[0]() - 273.15)

    # Grid heat duty
    print('Grid heater heat duty:', m.fs.grid_heater.heat_duty[0]())
    print('Accumulated grid heat duty:', m.fs.acc_grid_duty())
    print('Degrees of freedom:', degrees_of_freedom(m))

def main():
    m = create_mp_steady_state(
        GHI = 0, 
        elec_price = 0.07,
        mass_fr_fpc = 0.05,
        mass_fr_tes_hx_solar = 0.1,
        mass_fr_tes_process = 0.05,         
    )
    
    fix_dof_and_initialize(m, 
                           dt = 3600, 
                           GHI = 0 ,
                           mass_fr_tes_hx_solar = 0.1, 
                           mass_fr_tes_process = 0.05, 
                           process_inlet_temp=60)

    # print_results(m)
    print(degrees_of_freedom(m))

    solver = get_solver()
    results = solver.solve(m)

    print(degrees_of_freedom(m))

    return m
    

if __name__ == "__main__":    
    m = main()
    print_results(m)
