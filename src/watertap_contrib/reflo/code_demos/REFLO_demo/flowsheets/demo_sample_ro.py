# Import relevant libraries
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    Block,
    units as pyunits,
)

from pyomo.network import Arc
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

from idaes.core.util.initialization import propagate_state as _prop_state
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.core import REFLODatabase

from watertap_contrib.reflo.costing import (
    TreatmentCosting,
)

from watertap.unit_models.pressure_changer import Pump
from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.ro_system import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.components.deep_well_injection import *

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

from idaes.models.unit_models import Product, Feed

from idaes.core.util.model_statistics import *
import idaes.core.util.scaling as iscale

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *

from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import numpy as np
import pandas as pd
import matplotlib.image as mpimg
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils.results_dict import *
from watertap_contrib.reflo.code_demos.REFLO_demo.flowsheets.sweep_functions.case_study_plotting import *

"""
This module builds a sample reverse osmosis (RO) system using the REFLO framework.
It includes the following components:
- A flowsheet block to hold the system
- A NaCl parameter block for the RO properties
- A treatment block that includes feed, product, and waste streams
"""

# Build the model
def build_system(RE=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Get unit model input data from a database
    m.db = REFLODatabase()

    # Add uni
    m.fs.RO_properties = NaClParameterBlock()

    return m


# Build the treatment system
def build_treatment(m):
    # Create a treatment block
    m.fs.treatment = Block()
    # Add feed, product and waste streams
    m.fs.treatment.feed = Feed(property_package=m.fs.RO_properties)
    m.fs.treatment.product = Product(property_package=m.fs.RO_properties)

    # Add unit models for treatment and disposal
    m.fs.treatment.pump = Pump(property_package=m.fs.RO_properties)
    m.fs.treatment.RO = FlowsheetBlock(dynamic=False)
    m.fs.treatment.DWI = FlowsheetBlock(dynamic=False)

    # Build unit models
    build_ro(m, m.fs.treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)
    build_DWI(m, m.fs.treatment.DWI, prop_package=m.fs.RO_properties)


# Add connecions between unit models
def add_connections(m):
    # Connect feed to pump
    m.fs.treatment.feed_to_pump = Arc(
        source=m.fs.treatment.feed.outlet,
        destination=m.fs.treatment.pump.inlet,
    )

    # Connect pump to RO
    m.fs.treatment.pump_to_ro = Arc(
        source=m.fs.treatment.pump.outlet,
        destination=m.fs.treatment.RO.feed.inlet,
    )

    # Connect RO to product
    m.fs.treatment.ro_to_product = Arc(
        source=m.fs.treatment.RO.product.outlet,
        destination=m.fs.treatment.product.inlet,
    )

    # Connect RO to DWI
    m.fs.treatment.ro_to_dwi = Arc(
        source=m.fs.treatment.RO.disposal.outlet,
        destination=m.fs.treatment.DWI.unit.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

# Set operating conditions for the RO unit model
def set_operating_conditions(m):
    # Set feed flow rate
    m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(171.12)  
    m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.638) 
    m.fs.treatment.feed.properties[0].temperature.fix(298.15)  # 25 C
    m.fs.treatment.feed.properties[0].pressure.fix(101325)  # 101.325 kPa

def set_pump_operating_conditions(m, pump_pressure=30e5):
    # Set pump operating conditions
    m.fs.treatment.pump.efficiency_pump.fix(0.8)  # 80% efficiency
    m.fs.treatment.pump.control_volume.properties_out[0].pressure.fix(pump_pressure)

def set_ro_operating_conditions(m,ro_mem_area=10000):
    # Set RO operating conditions
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.95  # spacer porosity in membrane stage [-]
    length = 7  # effective membrane width [m]
    width = 20000  # effective membrane width [m]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    for idx, stage in m.fs.treatment.RO.stage.items():
        stage.module.A_comp.fix(mem_A)
        stage.module.B_comp.fix(mem_B)
        stage.module.area.fix(ro_mem_area / idx)
        stage.module.feed_side.velocity[0, 0].fix(0.35)
        # stage.module.length.fix(length)
        stage.module.width.setub(width)
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        stage.module.feed_side.channel_height.fix(height)
        stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

        stage.module.feed_side.friction_factor_darcy.setub(50)

        for e in stage.module.flux_mass_phase_comp:
            if e[-1] == "H2O":
                stage.module.flux_mass_phase_comp[e].setlb(1e-5)
                stage.module.flux_mass_phase_comp[e].setub(0.99)

    m.fs.treatment.RO.total_membrane_area = Var(
        initialize=10000,
        domain=NonNegativeReals,
        units=pyunits.m**2,
        doc="Total RO System Membrane Area",
    )

    m.fs.treatment.RO.eq_total_membrane_area = Constraint(
        expr=m.fs.treatment.RO.total_membrane_area
        == sum([stage.module.area for idx, stage in m.fs.treatment.RO.stage.items()])
    )

def add_ro_scaling(m, blk):

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    # Set scaling factors for module parameters in each stage in the RO unit model
    for idx, stage in blk.stage.items():
        module = stage.module
        iscale.set_scaling_factor(module.area, 1e5)
        iscale.set_scaling_factor(module.feed_side.area, 1)
        iscale.set_scaling_factor(module.width, 1e4)
        set_scaling_factor(module.length, 1e1)
        # set_scaling_factor(module.feed_side.velocity, 10)
        set_scaling_factor(module.feed_side.N_Sh_comp, 1e-4)

    # 
        for e in module.feed_side.properties:
            set_scaling_factor(
                module.feed_side.properties[e].flow_mass_phase_comp["Liq", "NaCl"], 1
            )
            set_scaling_factor(module.feed_side.properties[e].dens_mass_phase["Liq"], 1)
            # set_scaling_factor(module.feed_side.properties[e].dens_mass_phase["Liq"], 1e3)
            set_scaling_factor(module.feed_side.properties[e].mass_frac_phase_comp["Liq", "NaCl"], 1e1)

        for temp_stream in [
            module.eq_permeate_isothermal,
            module.feed_side.eq_equal_temp_interface,
            module.feed_side.eq_feed_isothermal,
            module.eq_permeate_outlet_isothermal,
        ]:
            for e in temp_stream:
                constraint_scaling_transform(temp_stream[e], 1e-2)
            for pressure_stream in [
                module.eq_permeate_outlet_isobaric,
                module.feed_side.eq_equal_pressure_interface,
            ]:
                for e in pressure_stream:
                    constraint_scaling_transform(pressure_stream[e], 1e-5)
            for e in module.eq_pressure_drop:
                constraint_scaling_transform(module.eq_pressure_drop[e], 1e-7)

        for e in module.feed_side.eq_N_Sh_comp:
            if e[-1] == "NaCl":
                constraint_scaling_transform(module.feed_side.eq_N_Sh_comp[e], 1e-3)

        for e in module.feed_side.eq_friction_factor:
            constraint_scaling_transform(module.feed_side.eq_friction_factor[e], 1e-2)
        for e in module.feed_side.eq_dP_dx:
            constraint_scaling_transform(module.feed_side.eq_dP_dx[e], 1e-2)

        set_scaling_factor(module.mixed_permeate[0.0].dens_mass_phase["Liq"], 1)
        set_scaling_factor(module.mixed_permeate[0.0].flow_vol_phase["Liq"], 100)
        constraint_scaling_transform(
            module.mixed_permeate[0.0].eq_flow_vol_phase["Liq"], 100
        )

        for e in module.recovery_mass_phase_comp:
            if e[-1] == "H2O":
                set_scaling_factor(module.recovery_mass_phase_comp, 1e1)

    calculate_scaling_factors(m)


def initialize_system(m):
    # Initialize the feed stream
    m.fs.treatment.feed.initialize()

    # Propagate state from feed to pump
    _prop_state(m.fs.treatment.feed_to_pump)

    # Initialize the pump
    m.fs.treatment.pump.initialize()

    # Propagate state from pump to RO
    _prop_state(m.fs.treatment.pump_to_ro)

    # Initialize the RO unit model
    init_ro_system(m, m.fs.treatment.RO)

    # Propagate state from RO to product
    _prop_state(m.fs.treatment.ro_to_product)

    # Initialize the product stream
    m.fs.treatment.product.initialize()

    # Propagate state from RO to DWI
    _prop_state(m.fs.treatment.ro_to_dwi)
    # Initialize the DWI unit model
    m.fs.treatment.DWI.unit.initialize()


def add_ro_recovery_constraint(m, blk, ro_recovery):
    m.fs.treatment.ro_water_recovery = Var(
        initialize=ro_recovery,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="RO Water Recovery",
    )

    blk.eq_water_recovery = Constraint(
        expr=blk.feed.properties[0].flow_vol * m.fs.treatment.ro_water_recovery
        == blk.product.properties[0].flow_vol
    )

def add_costing(m,electricity_price=0.058):
    # Add costing blocks

    # Add costing for treatment units
    m.fs.treatment.costing = TreatmentCosting()

    # Add pump costing
    m.fs.treatment.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing,
    )

    # Add RO costing
    for stage in m.fs.treatment.RO.stage.values():
        stage.module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.treatment.costing
        )
    
    # Add DWI costing
    add_DWI_costing(m.fs.treatment, m.fs.treatment.DWI, costing_blk=m.fs.treatment.costing)

    m.fs.treatment.costing.electricity_cost.fix(electricity_price)

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(m.fs.treatment.product.properties[0].flow_vol)

    product_m3h = pyunits.convert(
            m.fs.treatment.product.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.h,
        )
    m.fs.treatment.costing._add_flow_component_breakdowns(
        "electricity", "SEC_elec", product_m3h, period=pyunits.hr
    )

    m.fs.treatment.costing.initialize()


def set_ro_op_bounds(m, membrane_area=20000, lb=100, ub=900):
    
    m.fs.treatment.pump.control_volume.properties_out[0].pressure.setlb(
        lb * pyunits.psi
    )
    m.fs.treatment.pump.control_volume.properties_out[0].pressure.setub(
        ub * pyunits.psi
    )

    for _, stage in m.fs.treatment.RO.stage.items():
        stage.module.width.setub(5000)
        stage.module.feed_side.velocity[0, 0].unfix()
        stage.module.feed_side.velocity[0, 1].setlb(0.0)
        stage.module.feed_side.K.setlb(1e-6)
        stage.module.feed_side.friction_factor_darcy.setub(50)
        stage.module.flux_mass_phase_comp.setub(1)
        stage.module.feed_side.cp_modulus.setub(10)
        stage.module.rejection_phase_comp.setlb(1e-4)
        stage.module.feed_side.N_Re.setlb(1)
        stage.module.recovery_mass_phase_comp.setlb(1e-7)
 

def add_membrane_optimization(m):

    membrane_area=20000
    m.fs.treatment.pump.control_volume.properties_out[0].pressure.unfix()
    for _, stage in m.fs.treatment.RO.stage.items():
        stage.module.area.unfix()

    m.fs.treatment.RO.total_membrane_area.fix(membrane_area)
    m.fs.treatment.ro_water_recovery.fix(0.7)

    # Optimize the model
    set_ro_op_bounds(m)

    m.fs.membrane_area_objective = Objective(
    expr=m.fs.treatment.costing.LCOW, sense="minimize"
    )

    solver = get_solver()
    results = solver.solve(m, tee=False)

    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    return m


def run_recovery_sweep(sweep_type="ro_water_recovery"):

    sweep_dict = {
    'ro_water_recovery':np.linspace(0.3,0.8,5),
    }   
    
    input_dict = {
        'ro_water_recovery':0.8,
    }

    m = build_demo_sample_ro_sweep(
                ro_recovery=input_dict["ro_water_recovery"],
                )
    results_dict_test = build_results_dict(m)

    for i in sweep_dict[sweep_type]:
        input_dict[sweep_type] = i
        print(input_dict)
        m = build_demo_sample_ro_sweep(    
            ro_recovery=input_dict["ro_water_recovery"],
            )
        
        results_dict_test = results_dict_append(m, results_dict_test)

    df = pd.DataFrame.from_dict(results_dict_test)
    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/sample_ro_" + sweep_type + ".csv"
    df.to_csv(filename)


def plot_recovery_sweep(sweep_type="ro_water_recovery"):

    filename = "/Users/mhardika/Documents/watertap-seto/Mukta-Work/Demo/sample_ro_" + sweep_type + ".csv"
    df = pd.read_csv(filename).drop(columns="Unnamed: 0")
    df['ro_water_recovery'] = df["fs.treatment.ro_water_recovery"]*100

    xcol_dict = {
        "ro_water_recovery":"ro_water_recovery",
    }
    ax_dict = {
        "ro_water_recovery": "RO Water Recovery (%)",
    }

    unit_dict = {
        "Pump": "fs.treatment.pump.costing",
        "RO": "fs.treatment.RO.stage[1].module.costing",
        "DWI" : "fs.treatment.DWI.unit.costing",

    }
    agg_flows = {
        "Electricity":"electricity",
    }

    # xcol = "fs.water_recovery"
    flow_col = "fs.treatment.product.properties[0.0].flow_vol_phase[Liq]"
    ax_dict = dict(xlabel=ax_dict["ro_water_recovery"], ylabel="LCOW (\$/m$^3$)")


    # fig, ax = case_study_stacked_plot(
    #     df,
    #     unit_dict=unit_dict,
    #     costing_blk="fs.costing",
    #     agg_flows=agg_flows,
    #     xcol=xcol_dict['ro_water_recovery'],
    #     flow_col=flow_col,
    #     ax_dict=ax_dict,
    #     opex_hatch="\\\\\\",
    #     flow_hatch="..",
    # )

    # fig1 = plot_elec(df,flow_col,ax_dict)

    # Current file's directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # # Go one folder up
    parent_dir = os.path.dirname(current_dir)


    # fig.savefig(parent_dir+'/figures/results/sample-ro-recovery-sweep.png')
    # fig1.savefig(parent_dir+'/figures/results/sample-ro-recovery-sweep-sec.png')

    # Load the saved figures from disk
    img1 = mpimg.imread(parent_dir+'/figures/results/sample-ro-recovery-sweep.png')  # or .jpg, .pdf (with limitations), etc.
    img2 = mpimg.imread(parent_dir+'/figures/results/sample-ro-recovery-sweep-sec.png')
    # Now display both images in one new figure
    fig_combined, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.imshow(img1)
    ax1.axis('off')
 

    ax2.imshow(img2)
    ax2.axis('off')

    plt.tight_layout()
    plt.show()


def build_demo_sample_ro_sweep(ro_recovery=0.85):
    m = build_demo_sample_ro(ro_recovery)

    m.fs.treatment.ro_water_recovery.fix(ro_recovery)
    m.fs.treatment.pump.control_volume.properties_out[0].pressure.unfix()

    m.fs.membrane_area_objective = Objective(
    expr=m.fs.treatment.costing.LCOW, sense="minimize"
    )

    solver = get_solver()
    results = solver.solve(m, tee=False)

    return m


def build_demo_sample_ro(ro_recovery=0.85):
    # Build the system
    m = build_system()

    # Build the treatment system
    build_treatment(m)

    # Add connections between unit models
    add_connections(m)

    # Set operating conditions
    set_operating_conditions(m)

    # Set pump operating conditions
    set_pump_operating_conditions(m, pump_pressure=3e6)  # 300 kPa

    # Set RO operating conditions
    set_ro_operating_conditions(m, ro_mem_area=10000)  # 100

    # Add RO scaling factors
    add_ro_scaling(m, m.fs.treatment.RO)

    # Add RO recovery constraint
    add_ro_recovery_constraint(m, m.fs.treatment.RO, ro_recovery)

    # Initialize the system
    initialize_system(m)

    print("Degree of freedom:", degrees_of_freedom(m))

    solver = get_solver()

    # Add costing blocks
    add_costing(m, electricity_price=0.058)  # $0.058/kWh

    results = solver.solve(m, tee=False)

    return m

if __name__ == "__main__":
#     m=  build_demo_sample_ro(ro_recovery=0.7)

#     print("SEC electrical:", m.fs.treatment.costing.SEC_elec_component.display())
#     print("LCOW:", value(m.fs.treatment.costing.LCOW))
#     print("Membrane area:", value(m.fs.treatment.RO.total_membrane_area))
#     print("RO Water Recovery:", value(m.fs.treatment.ro_water_recovery))

    # m = add_membrane_optimization(m)

    # print("\nDegree of freedom:", degrees_of_freedom(m))
    # print("SEC electrical:", m.fs.treatment.costing.SEC_elec_component.display())
    # print("LCOW:", value(m.fs.treatment.costing.LCOW))
    # print("Membrane area:", value(m.fs.treatment.RO.total_membrane_area))
    # print("RO Water Recovery:", value(m.fs.treatment.ro_water_recovery))

    # run_recovery_sweep(sweep_type="ro_water_recovery")
    plot_recovery_sweep(sweep_type="ro_water_recovery")
