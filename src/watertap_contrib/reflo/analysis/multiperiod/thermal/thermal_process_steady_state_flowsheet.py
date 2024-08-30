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


# from thermal_storage import ThermalEnergyStorage
from watertap_contrib.reflo.unit_models.zero_dimensional.thermal_storage import ThermalEnergyStorage

# from flat_plate_physical import FlatPlatePhysical
from watertap_contrib.reflo.solar_models.zero_order.flat_plate_physical import FlatPlatePhysical


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

def build_thermal_flowsheet(m=None):

    '''
    This function adds the relevant components to the flowsheet
    '''
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic = False)
    m.fs.properties = WaterParameterBlock()
    
    # Include FPC model
    m.fs.fpc =  FlatPlatePhysical(
        property_package = m.fs.properties
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

    # HX between the TES and process
    m.fs.hx_process = HeatExchanger(
                hot_side = {"property_package": m.fs.properties},
                cold_side = {"property_package": m.fs.properties},
                delta_temperature_callback = delta_temperature_lmtd_callback,
                flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )

    # Grid heater to meet process inlet temperature set point
    m.fs.grid_heater = Heater(
        property_package = m.fs.properties)
    
    # Create coupling variables

    m.fs.previous_fpc_outlet_temperature = Var(
        domain = NonNegativeReals,
        initialize = 30 +273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc = 'Outlet temperature from FPC from the previous time step'
    )

    m.fs.previous_tes_tank_temp = Var(
        domain = NonNegativeReals,
        initialize = 25+273.15,
        bounds = (20+273.15, 99+273.15),
        units = pyunits.K,
        doc= 'Temperature of the thermal storage tank from the previous time step'
    )

    return m


def create_mp_steady_state(
        m = None,
        GHI = 900, 
        elec_price = 0.07,
        md_set_point = 50 + 273.15,
        md_outlet_temp = 55 + 273.15,
        mass_fr_md = 0.01,
        mass_fr_fpc = 0.05,
        mass_fr_tes_hx_solar = 0.1,
        mass_fr_tes_process = 0.05,
):

    m = build_thermal_flowsheet(m)

    _create_feed_streams(m,GHI,mass_fr_fpc,mass_fr_tes_hx_solar,mass_fr_tes_process,md_set_point,md_outlet_temp,mass_fr_md)
    _create_arcs(m)
    _create_constraints(m)
    initialize_feed_arcs(m)
    # initialize_system(m)
    m = set_model_inputs(m)

    return m

def _create_feed_streams(m,
                         GHI=500, 
                         mass_fr_fpc = 0.05,
                         mass_fr_tes_hx_solar=0.1,
                         mass_fr_tes_process=0.05,
                         md_set_point=50+273.15,
                         md_outlet_temp=35+273.15,
                         mass_fr_md = 0.01):

    m.fs.GHI = Var(initialize = 500, units=pyunits.W / pyunits.m**2)
    m.fs.GHI.fix(GHI)    
    # Creating a state block for the outlet stream from FPC
    m.fs.fpc_outlet = Feed(property_package = m.fs.properties)

    m.fs.fpc_outlet.properties[0].temperature.fix(m.fs.previous_fpc_outlet_temperature)
    m.fs.fpc_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_fpc)
    m.fs.fpc_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.fpc_outlet.properties[0].pressure.fix(101325)

    # Creating state blocks for the outlet stream from TES with the TES tank temperature from the previous step
    # Outlet stream from the TES going back to the solar HX
    m.fs.tes_hx_outlet = Feed(property_package = m.fs.properties)

    m.fs.tes_hx_outlet.properties[0].temperature.fix(m.fs.previous_tes_tank_temp)
    m.fs.tes_hx_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_hx_solar)
    m.fs.tes_hx_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.tes_hx_outlet.properties[0].pressure.fix(101325)


    # Outlet stream from the TES going to grid heater/process
    m.fs.tes_process_outlet = Feed(property_package = m.fs.properties)

    m.fs.tes_process_outlet.properties[0].temperature.fix(m.fs.previous_tes_tank_temp)
    m.fs.tes_process_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_tes_process)
    m.fs.tes_process_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.tes_process_outlet.properties[0].pressure.fix(101325)

    # In the future, an MD unit model should be included instead of these feed streams

    # Outlet stream going from process to process HX
    m.fs.process_outlet = Feed(property_package = m.fs.properties)

    m.fs.process_outlet.properties[0].temperature.fix(md_outlet_temp)
    m.fs.process_outlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_md)
    m.fs.process_outlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.process_outlet.properties[0].pressure.fix(101325)

    # Outlet stream going from process HX to process
    m.fs.process_inlet = Feed(property_package = m.fs.properties)

    m.fs.process_inlet.properties[0].temperature.fix(md_set_point)
    m.fs.process_inlet.properties[0].flow_mass_phase_comp['Liq','H2O'].fix(mass_fr_md)
    m.fs.process_inlet.properties[0].flow_mass_phase_comp['Vap','H2O'].fix(0)
    m.fs.process_inlet.properties[0].pressure.fix(101325)


def _create_arcs(m):

    m.fs.fpc_hx_solar = Arc(source = m.fs.fpc_outlet.outlet, 
                            destination = m.fs.hx_solar.hot_side_inlet,
                            doc = 'Connect FPC outlet to the solar HX hotside inlet')
    
    m.fs.hx_solar_fpc = Arc(source = m.fs.hx_solar.hot_side_outlet,
                            destination = m.fs.fpc.inlet,
                            doc = 'Connect solar hx hot side outlet back to fpc')
    
    m.fs.tes_hx_solar = Arc(source = m.fs.tes_hx_outlet.outlet,
                            destination = m.fs.hx_solar.cold_side_inlet,
                            doc = 'Connect TES from previous time step to the solar HX cold side inlet')
    
    m.fs.hx_solar_tes = Arc(source = m.fs.hx_solar.cold_side_outlet,
                            destination = m.fs.tes.tes_hx_inlet,
                            doc = 'Connect solar HX cold side outlet back to TES inlet')
    
    m.fs.tes_gridHtr = Arc(source = m.fs.tes_process_outlet.outlet,
                           destination = m.fs.grid_heater.inlet,
                           doc = 'Connect TES to the grid heater inlet')
    
    m.fs.gridHtr_hx_process = Arc(source = m.fs.grid_heater.outlet,
                                  destination = m.fs.hx_process.hot_side_inlet,
                                  doc = 'Connect grid heater to the process HX hot side inlet')
    
    m.fs.hx_process_tes = Arc(source = m.fs.hx_process.hot_side_outlet,
                              destination = m.fs.tes.tes_process_inlet,
                              doc = 'Connect the process HX hot side outlet back to TES' )
    
    m.fs.process_hx_process = Arc(source = m.fs.process_outlet.outlet,
                                  destination = m.fs.hx_process.cold_side_inlet,
                                  doc = 'Connect process outlet to the process HX cold side inlet')
        

def _create_constraints(m):

    @m.Constraint(doc='Set the grid heater outlet temperature based on heat balance in process HX')
    def eq_grid_hx_process(b,t):
        return b.fs.grid_heater.outlet.temperature[0] == b.fs.hx_process.hot_side_inlet.temperature[0]

    @m.Constraint()
    def process_temp(b):
        return b.fs.hx_process.cold_side_outlet.temperature[0] == b.fs.process_inlet.properties[0].temperature


def initialize_feed_arcs(m):

    # Initialize Feed streams
    solver = get_solver()
    optarg = solver.options

    m.fs.fpc_outlet.initialize(optarg=optarg)
    m.fs.tes_hx_outlet.initialize(optarg=optarg)
    m.fs.tes_process_outlet.initialize(optarg=optarg)
    m.fs.process_outlet.initialize(optarg=optarg)
    m.fs.process_inlet.initialize(optarg=optarg)


def initialize_system(m):
    blk = build_thermal_flowsheet()
    _create_feed_streams(blk)
    _create_arcs(blk)
    _create_constraints(blk)
    initialize_feed_arcs(blk)

    # Initialize unit models
    print('Initializing system')
    # Solar hx
    propagate_state(blk.fs.fpc_hx_solar)
    propagate_state(blk.fs.tes_hx_solar)

    blk.fs.hx_solar.hot_side_inlet.fix()
    blk.fs.hx_solar.cold_side_inlet.fix()
    blk.fs.hx_solar.hot_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    blk.fs.hx_solar.cold_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    blk.fs.hx_solar.effectiveness.fix(0.7)
    blk.fs.hx_solar.initialize_build()

    # FPC
    propagate_state(blk.fs.hx_solar_fpc)
    blk.fs.fpc.inlet.fix()
    blk.fs.fpc.G_trans.set_value(400)
    blk.fs.fpc.area_coll.fix(1)
    blk.fs.fpc.initialize()

    # Grid heater
    propagate_state(blk.fs.tes_gridHtr)
    blk.fs.grid_heater.inlet.fix()

    blk.fs.grid_heater.outlet.temperature[0].setub(99+273.15)
    blk.fs.grid_heater.outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)

    blk.fs.grid_heater.initialize()

    # Process HX
    propagate_state(blk.fs.process_hx_process)
    propagate_state(blk.fs.gridHtr_hx_process)
    propagate_state(blk.fs.process_hx_process)

    blk.fs.hx_process.cold_side_inlet.fix()
    blk.fs.hx_process.hot_side_inlet.fix()
    blk.fs.hx_process.hot_side_inlet.temperature.unfix()
    blk.fs.hx_process.overall_heat_transfer_coefficient[0].fix(100)
    blk.fs.hx_process.area.fix(1)

    blk.fs.hx_process.hot_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    blk.fs.hx_process.cold_side_outlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(0)
    blk.fs.hx_process.initialize_build()

    # TES
    propagate_state(blk.fs.hx_solar_tes)
    propagate_state(blk.fs.hx_process_tes)

    blk.fs.tes.tes_hx_inlet.fix()
    blk.fs.tes.tes_hx_inlet.temperature.unfix()
    blk.fs.tes.tes_process_inlet.fix()
    blk.fs.tes.tes_process_inlet.temperature.unfix()


    blk.fs.tes.dt.fix()
    blk.fs.tes.tes_initial_temp.fix(blk.fs.previous_fpc_outlet_temperature())
    blk.fs.tes.tes_diameter.fix(1.25)

    blk.fs.tes.initialize()

    solver = get_solver()
    results = solver.solve(m)

    # Unfix variables
    blk.fs.hx_solar.hot_side_inlet.temperature.unfix()
    blk.fs.hx_solar.cold_side_inlet.temperature.unfix()
    blk.fs.fpc.inlet.temperature.unfix()
    blk.fs.grid_heater.inlet.temperature.unfix()
    blk.fs.hx_process.cold_side_inlet.temperature.unfix()
    blk.fs.tes.tes_initial_temp.unfix()

    return

    
def set_model_inputs(m):
    # fix model inputs
    # Fixing the solar HX inlet stream
    print('Setting time step inputs')
    # solar hx
    propagate_state(m.fs.fpc_hx_solar)
    propagate_state(m.fs.tes_hx_solar)

    m.fs.hx_solar.hot_side_inlet.fix()
    m.fs.hx_solar.cold_side_inlet.fix()
    m.fs.hx_solar.initialize_build()

    # FPC
    propagate_state(m.fs.hx_solar_fpc)
    m.fs.fpc.inlet.fix()
    m.fs.fpc.G_trans.set_value(m.fs.GHI)
    m.fs.fpc.initialize()
    
    # Fixing the TES outlet to solar HX and grid heater
    m.fs.hx_solar.cold_side_inlet.fix()
    m.fs.grid_heater.inlet.fix()
    # Fixing process outlets stream
    m.fs.hx_process.cold_side_inlet.fix()

    # Grid heater
    propagate_state(m.fs.tes_gridHtr)
    m.fs.grid_heater.inlet.fix()
    m.fs.grid_heater.initialize()
    
    # Process HX
    propagate_state(m.fs.process_hx_process)
    propagate_state(m.fs.gridHtr_hx_process)
    propagate_state(m.fs.process_hx_process)

    m.fs.hx_process.cold_side_inlet.fix()
    m.fs.hx_process.hot_side_inlet.fix()
    m.fs.hx_process.hot_side_inlet.temperature.unfix()
    m.fs.hx_process.initialize_build()
    
    # TES
    propagate_state(m.fs.hx_solar_tes)
    propagate_state(m.fs.hx_process_tes)

    m.fs.tes.tes_hx_inlet.fix()
    m.fs.tes.tes_hx_inlet.temperature.unfix()
    m.fs.tes.tes_process_inlet.fix()
    m.fs.tes.tes_process_inlet.temperature.unfix()


    m.fs.tes.tes_initial_temp.fix(m.fs.previous_fpc_outlet_temperature())
    m.fs.tes.initialize()

    return m

def print_results(m):
    print('Tank temperature:',m.fs.tes.tes_temp[0]()-273.15)

    # FPC
    print('FPC outlet temperature:',m.fs.fpc.outlet.temperature[0]()-273.15)

    # Grid inlet temperature
    print('Grid heater inlet temperature:',m.fs.grid_heater.inlet.temperature[0]() - 273.15)

    # Grid outlet temperature
    print('Grid heater outlet temperature:',m.fs.grid_heater.outlet.temperature[0]() - 273.15)

    # Grid heat duty
    print('Grid heater heat duty:', m.fs.grid_heater.heat_duty[0]())

    # Process inlet temperature (to MD)
    print('Process heater - hot side inlet temperature:', m.fs.hx_process.hot_side_inlet.temperature[0]() - 273.15)

    # Process inlet temperature (to MD)
    print('Process heater - cold side outlet temperature:', m.fs.hx_process.cold_side_outlet.temperature[0]() - 273.15)


def main():
    m = build_thermal_flowsheet(
        GHI = 900, 
        elec_price = 0.07,
        md_set_point = 50 + 273.15,
        md_outlet_temp = 35 + 273.15,
        mass_fr_md = 0.01,
        mass_fr_fpc = 0.05,
        mass_fr_tes_hx_solar = 0.1,
        mass_fr_tes_process = 0.05,         
    )
    
    print_results(m)
    print(degrees_of_freedom(m))

    solver = get_solver()
    results = solver.solve(m)

    return m
    

if __name__ == "__main__":    
    m = main()
    print_results(m)