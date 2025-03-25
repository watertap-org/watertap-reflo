import os
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def save_water_dams(
        df,
        feed_conditions = None,
        global_costing_blk="fs.treatment.costing",
        xcol = None,
        treatment_unit_costing_dict = None,
        energy_unit_costing_dict = None,
        heat_gen_unit_dict = None,
        elec_gen_unit_dict = None,
        ):
    

    '''
    df : Dataframe that saves all relevant costing columns
    '''

    # Get costing params
    costing_params = dict()

    global_params = [
        "maintenance_labor_chemical_factor",
        "utilization_factor",
        "capital_recovery_factor",
        "total_investment_factor",
    ]

    for gp in global_params:
        # make sure all the global parameters have the same value
        if not len(df[f"{global_costing_blk}.{gp}"].unique()) == 1:
            print(df[f"{global_costing_blk}.{gp}"].unique())
            raise ValueError(f"Global parameter {gp} does not have uniform value.")
        costing_params[gp] = df[f"{global_costing_blk}.{gp}"].iloc[0]


    ### Start saving to WaterDAMS
    water_dams_df = pd.DataFrame()

    # Add LCOW
    water_dams_df['LCOW ($/m3)'] = df['fs.costing.LCOT']

    for u,k in xcol_dict.items():
        water_dams_df[u] = df[k]

    # Add Treatment heat requirement
    water_dams_df['Treatment Units Heat Consumed (MWh/year)'] = df['fs.treatment.costing.aggregate_flow_heat']*pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh/pyunits.year)()
    # Add treatment electricity requirement
    water_dams_df['Treatment Units Electricity Consumed (MWh/year)'] = df['fs.treatment.costing.aggregate_flow_electricity']*pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh/pyunits.year)()

    # Calculate total heat and electricity generated
    water_dams_df['Heat Generated (MWh/year)'] = 0 
    water_dams_df['Electricity Generated (MWh/year)'] = 0 
    
    if heat_gen_unit_dict!= None:
        # Add heat generating unit electricity requirement
        for key,value in heat_gen_unit_dict.items():
            water_dams_df['Heat Units Electricity Consumed (MWh/year)'] = df[value + '.electricity']*pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh/pyunits.year)()

            water_dams_df['Heat Generated (MWh/year)'] = water_dams_df['Heat Generated (MWh/year)'] - df[value + '.heat']*pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh/pyunits.year)()

    if elec_gen_unit_dict!= None:
        water_dams_df['Electricity Generated (MWh/year)'] = water_dams_df['Electricity Generated (MWh/year)'] - df[value + '.electricity']*pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh/pyunits.year)()

    # Add Treatment Capex and Opex
    total_treatment_capex = 0
    total_treatment_opex = 0

    for u, b in treatment_unit_costing_dict.items():
        ## CAPEX
        unit_capex = 0

        try:
            print(f"{b}.capital_cost")
            unit_capex += (
                df[f"{b}.capital_cost"]
                * costing_params["total_investment_factor"]
            )  # USD2023

        except KeyError:
            print(f"No CAPEX for {u} found.")
            pass

        water_dams_df[u +' Capex ($)'] = unit_capex

        total_treatment_capex += unit_capex  # $

        ### OPEX
        unit_opex_total = 0
        unit_opex_total += (
            unit_capex * costing_params["maintenance_labor_chemical_factor"]
        )  # $ / year

        try:
            unit_opex_total += df[f"{b}.fixed_operating_cost"]
        except KeyError:
            print(f"No Fixed OPEX for {u} found.")
            pass

        try:
            unit_opex_total += df[f"{b}.variable_operating_cost"]
        except KeyError:
            print(f"No Variable OPEX for {u} found.")
            pass

        water_dams_df[u +' Opex ($/year)'] = unit_opex_total

        total_treatment_opex += unit_opex_total  # $ / year

    # Add Energy Capex and Opex
    total_energy_capex = 0
    total_energy_opex = 0

    for u, b in energy_unit_costing_dict.items():
        ## CAPEX
        unit_capex = 0

        try:
            print(f"{b}.capital_cost")
            unit_capex += (
                df[f"{b}.capital_cost"]
                * costing_params["total_investment_factor"]
            )  # USD2023

        except KeyError:
            print(f"No CAPEX for {u} found.")
            pass

        water_dams_df[u +' Capex ($)'] = unit_capex

        total_energy_capex += unit_capex  # $

        ### OPEX
        unit_opex_total = 0
        unit_opex_total += (
            unit_capex * costing_params["maintenance_labor_chemical_factor"]
        )  # $ / year

        try:
            unit_opex_total += df[f"{b}.fixed_operating_cost"]
        except KeyError:
            print(f"No Fixed OPEX for {u} found.")
            pass

        try:
            unit_opex_total += df[f"{b}.variable_operating_cost"]
        except KeyError:
            print(f"No Variable OPEX for {u} found.")
            pass

        water_dams_df[u +' Opex ($/year)'] = unit_opex_total

        total_energy_opex += unit_opex_total  # $ / year

    water_dams_df.insert(1,'Total Capex ($)', total_treatment_capex + total_energy_capex)
    water_dams_df.insert(2,'Total Opex ($/year)', total_treatment_opex + total_energy_opex)

    if feed_conditions!=None:
        water_dams_df.insert(3,'Feed Flow Rate (kg/s)', df[feed_conditions["feed_flow_mass_water"]] + df[feed_conditions['feed_flow_mass_tds']])
        water_dams_df.insert(4,'TDS Flow Rate (kg/s)', df[feed_conditions['feed_flow_mass_tds']])
    else:
        water_dams_df.insert(3,'Feed Flow Rate (kg/s)', 'Unknown')
        water_dams_df.insert(4,'TDS Flow Rate (kg/s)', 'Unknown')

    return water_dams_df



if __name__ == "__main__":

    # TODO
    # Flow , salinity
   
    # File path
    filepath = os.path.abspath(__file__)
    util_dir = os.path.dirname(filepath)

    test_filename = util_dir + "/test_water_dams.csv" 
   
    df = pd.read_csv(test_filename).drop(columns="Unnamed: 0")

    feed_conditions = {
        'feed_flow_mass_water': 'fs.treatment.feed.properties[0.0].flow_mass_phase_comp[Liq,H2O]',
        'feed_flow_mass_tds':'fs.treatment.feed.properties[0.0].flow_mass_phase_comp[Liq,TDS]'
    }

    # Sweep Variables
    xcol_dict = {
    "Water Recovery": "fs.water_recovery",
    "Heat Price ($/kWh)": "fs.costing.heat_cost_buy",
    "Thermal Hours Storage (h)":"fs.energy.FPC.hours_storage",
    "Fraction of Heat from Grid":"fs.costing.frac_heat_from_grid",
    "Injecton Cost ($/m3)": "fs.treatment.costing.deep_well_injection.dwi_lcow",
    "FPC Collector Cost ($/m2)":"fs.energy.costing.flat_plate.cost_per_area_collector",
    "Thermal Storage Cost ($/m3)":"fs.energy.costing.flat_plate.cost_per_volume_storage",
    }

    # Costing for each treatment unit model
    treatment_unit_costing_dict = {
    "MD": "fs.treatment.md.unit",
    "DWI": "fs.treatment.dwi.unit.costing",
    }

    # Costing for each energy unit model
    energy_unit_costing_dict = {
        "FPC": "fs.energy.FPC.costing", 
    }

    # Heat generating unit models
    heat_gen_unit_dict = {
    "FPC": "fs.energy.FPC",
        }
    
    # Electricity generating unit models
    elec_gen_unit_dict = None


    results = save_water_dams(
    df=df,
    feed_conditions = feed_conditions,
    global_costing_blk="fs.treatment.costing",
    xcol = xcol_dict,
    treatment_unit_costing_dict = treatment_unit_costing_dict,
    energy_unit_costing_dict = energy_unit_costing_dict,
    heat_gen_unit_dict = heat_gen_unit_dict,
    elec_gen_unit_dict = elec_gen_unit_dict,
    )

    results.to_csv('/Users/mhardika/Documents/watertap-seto/Mukta-Work/kbhdp-case-study-md/water_dams/test.csv')