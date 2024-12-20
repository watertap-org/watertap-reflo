costing_data_keys = [
    {
        "filekey": "fs.treatment.costing.co2.cost",
        "return_key": "fs.treatment.costing.co2.cost",
        # "units": "%",
    },
    {
        "filekey": "fs.treatment.costing.lime.cost",
        "return_key": "fs.treatment.costing.lime.cost",
        # "units": "%",
    },
    {
        "filekey": "fs.treatment.costing.mgcl2.cost",
        "return_key": "fs.treatment.costing.mgcl2.cost",
        # "units": "%",
    },
    {
        "filekey": "fs.treatment.costing.soda_ash.cost",
        "return_key": "fs.treatment.costing.soda_ash.cost",
        # "units": "%",
    },
    {
        "filekey": "fs.costing.frac_heat_from_grid",
        "return_key": "Grid Frac Heat",
        # "units": "%",
    },
    {
        "filekey": "fs.costing.heat_cost_buy",
        "return_key": "fs.costing.heat_cost_buy",
        # "units": "USD/kWh",
    },
    {
        "filekey": "fs.energy.costing.flat_plate.cost_per_area_collector",
        "return_key": "fs.energy.costing.flat_plate.cost_per_area_collector",
        # "units": "USD/kWh",
    },
    {
        "filekey": "fs.costing.electricity_cost_buy",
        "return_key": "fs.costing.electricity_cost_buy",
        # "units": "USD/kWh",
    },
    {
        "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
        "return_key": "FPC Cost",
        "units": "USD/a/kW",
    },
    {
        "filekey": "fs.water_recovery",
        "return_key": "fs.water_recovery",
        # "units": "%",
    },
    {
        "filekey": "fs.treatment.costing.LCOW",
        "return_key": "LCOW",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.costing.LCOT",
        "return_key": "LCOT",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.energy.pv.annual_energy",
        "return_key": "fs.energy.pv.annual_energy",
        # "units": "kWh",
    },
    {
        "filekey": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        "return_key": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        # "units": "kWh",
    },
    {
        "filekey": "fs.treatment.costing.deep_well_injection.dwi_lcow",
        "return_key": "Disposal Cost",
        "units": "USD/m**3",
    },
    {
        "filekey": "fs.costing.frac_heat_from_grid",
        "return_key": "fs.costing.frac_heat_from_grid",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
        "return_key": "fs.energy.costing.flat_plate.fixed_operating_by_capacity",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.treatment.md.unit.md_costing.capital_cost",
        "return_key": "fs.treatment.md.unit.md_costing.capital_cost",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.treatment.md.unit.md_costing.fixed_operating_cost",
        "return_key": "fs.treatment.md.unit.md_costing.fixed_operating_cost",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.treatment.costing.aggregate_flow_costs[heat]",
        "return_key": "fs.treatment.costing.aggregate_flow_costs[heat]",
        # "units": "USD/m**3",
    },
    {
        "filekey": "fs.treatment.costing.aggregate_flow_costs[electricity]",
        "return_key": "fs.treatment.costing.aggregate_flow_costs[electricity]",
        # "units": "USD/m**3",
    },
]

default_device_groups = {
    "Heat": {
        "OPEX": {
            "units": {
                "fs.costing.total_heat_operating_cost",
            },
        },
    },
    "Electricity": {
        "OPEX": {
            "units": {
                "fs.costing.aggregate_flow_electricity_purchased",
            },
        },
    },
}

figure_device_groups = {
    "KBHDP_SOA_1": {
        "Heat": {
            "OPEX": {
                "units": {
                    "fs.costing.total_heat_operating_cost",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Electricity": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_electricity_purchased",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Chemicals": {
            "OPEX": {
                "units": {
                    "fs.treatment.costing.aggregate_flow_costs[lime]",
                    "fs.treatment.costing.aggregate_flow_costs[mgcl2]",
                    "fs.treatment.costing.aggregate_flow_costs[soda_ash]",
                    "fs.treatment.costing.aggregate_flow_costs[co2]",
                },
            },
        },
        "Injection": {
            "OPEX": {
                "units": {"fs.treatment.DWI.unit.costing.variable_operating_cost"},
            },
        },
        "Pumps": {
            "CAPEX": {
                "units": {"fs.treatment.pump.costing.capital_cost"},
            },
        },
        "UF": {
            "CAPEX": {
                "units": {
                    "fs.treatment.UF.unit.costing.capital_cost",
                },
            },
        },
        "RO": {
            "CAPEX": {
                "units": {
                    "fs.treatment.RO.stage[1].module.costing.capital_cost",
                    "fs.treatment.RO.stage[2].module.costing.capital_cost",
                }
            },
            "OPEX": {
                "units": {
                    "fs.treatment.RO.stage[1].module.costing.fixed_operating_cost",
                    "fs.treatment.RO.stage[2].module.costing.fixed_operating_cost",
                }
            },
        },
        "Softening": {
            "CAPEX": {
                "units": {
                    "fs.treatment.softener.unit.costing.capital_cost",
                }
            },
            "OPEX": {
                "units": {
                    "fs.treatment.softener.unit.costing.fixed_operating_cost",
                }
            },
        },
    },
    "KBHDP_RPT_1": {
        "Heat": {
            "OPEX": {
                "units": {
                    "fs.costing.total_heat_operating_cost",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Electricity": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_electricity_purchased",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Injection": {
            "OPEX": {
                "units": {"fs.treatment.DWI.unit.costing.variable_operating_cost"},
            },
        },
        "Pumps": {
            "CAPEX": {
                "units": {"fs.treatment.pump.costing.capital_cost"},
            },
        },
        "PV": {
            "CAPEX": {
                "units": {"fs.energy.pv.costing.capital_cost"},
            },
            "OPEX": {
                "units": {
                    "fs.energy.pv.costing.fixed_operating_cost",
                }
            },
        },
        "EC": {
            "CAPEX": {
                "units": {
                    "fs.treatment.EC.ec.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.treatment.EC.ec.costing.fixed_operating_cost",
                    "fs.treatment.costing.aggregate_flow_costs[aluminum]",
                },
            },
        },
        "RO": {
            "CAPEX": {
                "units": {
                    "fs.treatment.RO.stage[1].module.costing.capital_cost",
                    "fs.treatment.RO.stage[2].module.costing.capital_cost",
                }
            },
            "OPEX": {
                "units": {
                    "fs.treatment.RO.stage[1].module.costing.fixed_operating_cost",
                    "fs.treatment.RO.stage[2].module.costing.fixed_operating_cost",
                }
            },
        },
        "UF": {
            "CAPEX": {
                "units": {
                    "fs.treatment.UF.unit.costing.capital_cost",
                },
            },
        },
    },
    "KBHDP_RPT_2": {
        "FPC": {
            "CAPEX": {
                "units": {"fs.energy.FPC.costing.capital_cost"},
            },
            "OPEX": {
                "units": {
                    "fs.energy.FPC.costing.fixed_operating_cost",
                }
            },
        },
        "Heat": {
            "OPEX": {
                "units": {
                    "fs.costing.total_heat_operating_cost",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Electricity": {
            "OPEX": {
                "units": {
                    "fs.costing.aggregate_flow_electricity_purchased",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Injection": {
            "OPEX": {
                "units": {"fs.treatment.DWI.unit.costing.variable_operating_cost"},
            },
        },
        "Pumps": {
            "CAPEX": {
                "units": {"fs.treatment.pump.costing.capital_cost"},
            },
        },
        "EC": {
            "CAPEX": {
                "units": {
                    "fs.treatment.EC.ec.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.treatment.EC.ec.costing.fixed_operating_cost",
                    "fs.treatment.costing.aggregate_flow_costs[aluminum]",
                },
            },
        },
        "LTMED": {
            "CAPEX": {
                "units": {
                    "fs.treatment.LTMED.unit.costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.treatment.LTMED.unit.costing.fixed_operating_cost",
                },
            },
        },
        "UF": {
            "CAPEX": {
                "units": {
                    "fs.treatment.UF.unit.costing.capital_cost",
                },
            },
        },
    },
    "KBHDP_RPT_3": {
        "Heat": {
            "OPEX": {
                "units": {
                    "fs.costing.total_heat_operating_cost",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "Electricity": {
            "OPEX": {
                "units": {
                    "fs.costing.total_electric_operating_cost",
                    #   "fs.costing.aggregate_flow_electricity_sold"
                },
            },
        },
        "FPC": {
            "CAPEX": {
                "units": {"fs.energy.FPC.costing.capital_cost"},
            },
            "OPEX": {
                "units": {
                    "fs.energy.FPC.costing.fixed_operating_cost",
                }
            },
        },
        "Injection": {
            "OPEX": {
                "units": {"fs.treatment.DWI.unit.costing.variable_operating_cost"},
            },
        },
        "MD": {
            "CAPEX": {
                "units": {
                    "fs.treatment.md.unit.md_costing.capital_cost",
                },
            },
            "OPEX": {
                "units": {
                    "fs.treatment.md.unit.md_costing.fixed_operating_cost",
                    #fs.treatment.md.unit.md_costing.fixed_operating_cost
                    # "fs.treatment.costing.aggregate_flow_costs[heat]",
                    # "fs.treatment.costing.aggregate_flow_costs[electricity]",
                },
            },
        },
    },
}
