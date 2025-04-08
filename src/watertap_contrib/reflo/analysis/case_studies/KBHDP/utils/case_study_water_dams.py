import os
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def save_water_dams(
    df,
    feed_conditions=None,
    sys_costing_blk="fs.treatment.costing",
    treat_costing_blk=None,
    en_costing_blk=None,
    xcol_dict=None,
    treatment_unit_costing_dict=None,
    energy_unit_costing_dict=None,
    heat_gen_unit_dict=None,
    elec_gen_unit_dict=None,
    agg_dict=None,
    case_study_name="",
    scenario_name="", 
    costing_params_dict=None,
):
    """
    df : Dataframe that saves all relevant costing columns
    """

    global costing_params
    found_flows = list()
    print("found aggregated flows:")
    for c in df.columns:
        if "aggregate_flow_costs" in c:
            flow = c.split("[")[-1].replace("]", "")
            found_flows.append(flow)
    # for ff in found_flows:
            print(f"\t{c}\n\t{flow}")

    # Get costing params

    if costing_params_dict is not None:
        costing_params = costing_params_dict
    else:
        costing_params = dict()

        global_params = [
            "maintenance_labor_chemical_factor",
            "utilization_factor",
            "capital_recovery_factor",
            "total_investment_factor",
        ]

        for gp in global_params:
            # make sure all the global parameters have the same value
            if not len(df[f"{treat_costing_blk}.{gp}"].unique()) == 1:
                print(df[f"{treat_costing_blk}.{gp}"].unique())
                raise ValueError(f"Global parameter {gp} does not have uniform value.")
            costing_params[gp] = df[f"{treat_costing_blk}.{gp}"].iloc[0]
        print(costing_params)

    ### Start saving to WaterDAMS
    water_dams_df = pd.DataFrame()




    try:
        # Add Treatment heat requirement
        water_dams_df["Treatment Units Heat Consumed (MWh/year)"] = (
            df[f"{treat_costing_blk}.aggregate_flow_heat"]
            * pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh / pyunits.year)()
        )
    except:
        print(f"aggregate_flow_heat not found!")

    try:
        # Add treatment electricity requirement
        water_dams_df["Treatment Units Electricity Consumed (MWh/year)"] = (
            df[f"{treat_costing_blk}.aggregate_flow_electricity"]
            * pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh / pyunits.year)()
        )
    except:
        print(f"aggregate_flow_electricity not found!")

    # Calculate total heat and electricity generated
    water_dams_df["Heat Generated (MWh/year)"] = 0
    water_dams_df["Electricity Generated (MWh/year)"] = 0

    if heat_gen_unit_dict != None:
        # Add heat generating unit electricity requirement
        for k, v in heat_gen_unit_dict.items():
            water_dams_df["Heat Units Electricity Consumed (MWh/year)"] = (
                df[v + ".electricity"]
                * pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh / pyunits.year)()
            )

            water_dams_df["Heat Generated (MWh/year)"] = (
                water_dams_df["Heat Generated (MWh/year)"]
                - df[v + ".heat"]
                * pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh / pyunits.year)()
            )

    if elec_gen_unit_dict != None:
        for k, v in elec_gen_unit_dict.items():
            water_dams_df["Electricity Generated (MWh/year)"] = (
                water_dams_df["Electricity Generated (MWh/year)"]
                - df[v + ".electricity"]
                * pyunits.convert(1 * pyunits.kW, to_units=pyunits.MWh / pyunits.year)()
            )

    # Add Treatment Capex and Opex
    total_treatment_capex = 0
    total_treatment_opex = 0

    for u, b in treatment_unit_costing_dict.items():

        if u == "RO":
            for s, bb in b.items():

                print(s, bb)

                water_dams_df, total_treatment_capex, total_treatment_opex = (
                    agg_unit_costing(
                        df,
                        f"RO Stage {s}",
                        bb,
                        water_dams_df,
                        total_treatment_capex,
                        total_treatment_opex,
                    )
                )
        else:
            water_dams_df, total_treatment_capex, total_treatment_opex = (
                agg_unit_costing(
                    df, u, b, water_dams_df, total_treatment_capex, total_treatment_opex
                )
            )

    # Add Energy Capex and Opex
    total_energy_capex = 0
    total_energy_opex = 0

    for u, b in energy_unit_costing_dict.items():
        ## CAPEX
        # unit_capex = 0

        water_dams_df, total_energy_capex, total_energy_opex = agg_unit_costing(
            df, u, b, water_dams_df, total_energy_capex, total_energy_opex
        )

    # water_dams_df = pd.DataFrame()
    if agg_dict is not None:
        # print("agg_dict")
        for flow, funits in agg_dict.items():
            if funits in ("kg/day", "kg/d"):
                water_dams_df[
                    f"{flow.title().replace('_', ' ')} Annual Use (kg/yr)"
                ] = (df[f"{sys_costing_blk}.aggregate_flow_{flow}"] * 365)
            else:
                if flow in ["electricity", "heat"]:
                    water_dams_df[
                        f"Treatment System {flow.title().replace('_', ' ')} Power Required (kW)"
                    ] = df[f"{treat_costing_blk}.aggregate_flow_{flow}"]
                    if en_costing_blk is not None:
                        water_dams_df[
                            f"Energy System {flow.title().replace('_', ' ')} Power Required (kW)"
                        ] = df[f"{en_costing_blk}.aggregate_flow_{flow}"]
                    water_dams_df[
                        f"System Net Grid {flow.title().replace('_', ' ')} Power Required (kW)"
                    ] = df[f"{sys_costing_blk}.aggregate_flow_{flow}"]
                    if sys_costing_blk == treat_costing_blk:
                        water_dams_df[
                            f"System Net Grid {flow.title().replace('_', ' ')} Annual Cost ($/yr)"
                        ] = df[f"{sys_costing_blk}.aggregate_flow_costs[{flow}]"]
                    elif flow == "electricity":
                        water_dams_df[
                            f"System Net Grid Electricity Annual Cost ($/yr)"
                        ] = df[f"{sys_costing_blk}.total_electric_operating_cost"]
                    elif flow == "heat":
                        water_dams_df[
                            f"System Net Grid Heat Annual Cost ($/yr)"
                        ] = df[f"{sys_costing_blk}.total_heat_operating_cost"]
                    continue
                else:
                    water_dams_df[
                        f"{flow.title().replace('_', ' ')} Annual Use ({funits})"
                    ] = df[f"{treat_costing_blk}.aggregate_flow_{flow}"]

            water_dams_df[f"{flow.title().replace('_', ' ')} Annual Cost ($/yr)"] = df[
                f"{treat_costing_blk}.aggregate_flow_costs[{flow}]"
            ]


    if treat_costing_blk is not None:

        water_dams_df.insert(0, "Total Treatment System Capex ($)", df[f"{treat_costing_blk}.total_capital_cost"])
        water_dams_df.insert(
            1, "Total Treatment System Opex ($/year)", df[f"{treat_costing_blk}.total_operating_cost"]
        )

    if en_costing_blk is not None:

        water_dams_df.insert(2, "Total Energy System Capex ($)", df[f"{en_costing_blk}.total_capital_cost"])
        water_dams_df.insert(
            3, "Total Energy System Opex ($/year)", df[f"{en_costing_blk}.total_operating_cost"]
        )

    water_dams_df.insert(0, "Total System Capex ($)", df[f"{sys_costing_blk}.total_capital_cost"])
    water_dams_df.insert(
        1, "Total System Opex ($/year)", df[f"{sys_costing_blk}.total_operating_cost"]
    )
    # Add LCOW
    try:
        water_dams_df.insert(0, "LCOW ($/m3)", df[f"{sys_costing_blk}.LCOT"])
    except KeyError:
        print("LCOT not found; trying LCOW...")
        water_dams_df.insert(0, "LCOW ($/m3)", df[f"{sys_costing_blk}.LCOW"])

    water_dams_df.insert(0, "Inlet Flow Rate (MGD)", feed_conditions["feed_flow_mgd"])  # MGD
    water_dams_df.insert(1, "Inlet TDS (g/L)", feed_conditions["feed_tds"])  #g/L
    for u, k in xcol_dict.items():
        water_dams_df.insert(0, u, df[k]) 
    water_dams_df.insert(0, "Sweep Column", list(xcol_dict.keys())[0])  
    water_dams_df.insert(0, "Scenario Name", scenario_name)  
    water_dams_df.insert(0, "Case Study Name", case_study_name)  


    return water_dams_df


def agg_unit_costing(df, u, b, water_dams_df, total_capex, total_opex):
    print(u, b)

    unit_capex = 0

    try:
        # print(f"{b}.capital_cost")
        unit_capex += (
            df[f"{b}.capital_cost"] * costing_params["total_investment_factor"]
        )  # USD2023

    except KeyError:
        print(f"No CAPEX for {u} found.")
        pass

    water_dams_df[u + " Capex ($)"] = unit_capex

    total_capex += unit_capex  # $

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

    water_dams_df[u + " Opex ($/year)"] = unit_opex_total

    total_opex += unit_opex_total  # $ / year

    return water_dams_df, total_capex, total_opex


if __name__ == "__main__":

    # TODO
    # Flow , salinity

    # File path
    filepath = os.path.abspath(__file__)
    util_dir = os.path.dirname(filepath)
    project_folder = os.path.dirname(util_dir)
    water_dams_folder = os.path.join(project_folder, "water_dams")

    test_filename = util_dir + "/test_water_dams.csv"

    df = pd.read_csv(test_filename).drop(columns="Unnamed: 0")

    feed_conditions = {
        "feed_flow_mass_water": "fs.treatment.feed.properties[0.0].flow_mass_phase_comp[Liq,H2O]",
        "feed_flow_mass_tds": "fs.treatment.feed.properties[0.0].flow_mass_phase_comp[Liq,TDS]",
    }

    # Sweep Variables
    xcol_dict = {
        "Water Recovery": "fs.water_recovery",
        "Heat Price ($/kWh)": "fs.costing.heat_cost_buy",
        "Thermal Hours Storage (h)": "fs.energy.FPC.hours_storage",
        "Fraction of Heat from Grid": "fs.costing.frac_heat_from_grid",
        "Injecton Cost ($/m3)": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        "FPC Collector Cost ($/m2)": "fs.energy.costing.flat_plate.cost_per_area_collector",
        "Thermal Storage Cost ($/m3)": "fs.energy.costing.flat_plate.cost_per_volume_storage",
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
        feed_conditions=feed_conditions,
        sys_costing_blk="fs.treatment.costing",
        xcol=xcol_dict,
        treatment_unit_costing_dict=treatment_unit_costing_dict,
        energy_unit_costing_dict=energy_unit_costing_dict,
        heat_gen_unit_dict=heat_gen_unit_dict,
        elec_gen_unit_dict=elec_gen_unit_dict,
    )

    results.to_csv(os.path.join(water_dams_folder, "test.csv"))
