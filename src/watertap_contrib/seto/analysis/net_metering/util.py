from pyomo.environ import value, Var, units as pyunits
import os

absolute_path = os.path.dirname(__file__)

__all__ = ["display_ro_pv_results", "display_pv_results", "display_costing_breakdown"]

def display_ro_pv_results(m, sep="."):
    ro = m.fs.treatment.ro
    erd = m.fs.treatment.erd
    liq = "Liq"
    nacl = "NaCl"
    header = f'{"PARAM":<25s}{"VALUE":<25s}{"UNITS":<10s}'
    prop_in = ro.feed_side.properties_in[0]
    prop_out = ro.feed_side.properties_out[0]
    prop_perm = ro.mixed_permeate[0]
    # pv_cost = m.fs.energy.costing.pv_surrogate
    line = f'\n{f"{sep*60}":<60s}'
    flux_lmh = pyunits.convert(
        ro.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        / ro.feed_side.properties_in[0].dens_mass_phase["Liq"],
        to_units=pyunits.liter / pyunits.m**2 / pyunits.hr,
    )()
    flow_out_L = pyunits.convert(
        prop_perm.flow_vol, to_units=pyunits.liter / pyunits.hr
    )
    title = f'\n{"=======> INTEGRATED SYSTEM RESULTS <=======":^60}\n'
    print(line)
    print(title)
    print(header)
    if hasattr(m.fs, "sys_costing"):
        print(f'{"LCOW":<24s}{f"${m.fs.sys_costing.LCOW():<25.3f}"}{"$/m3":<15s}')
        if hasattr(m.fs.sys_costing, "LCOE"):
            print(f'{"LCOE":<24s}{f"${m.fs.sys_costing.LCOE():<25.3f}"}{"$/kWh":<15s}')
        print(
            f'{"Total Capital Cost":<24s}{f"${m.fs.sys_costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"Total Op. Cost":<24s}{f"${m.fs.sys_costing.total_operating_cost():<25,.0f}"}{"$/yr":<10s}'
        )
        print(
            f'{"SEC":<25s}{f"{m.fs.sys_costing.specific_electric_energy_consumption():<25.2f}"}{"kWh/m3":<10s}'
        )
        print(
            f'{"PV Avg. Elect. Gen":<25s}{f"{-1*m.fs.energy.pv.electricity():<25.1f}"}{"kW":<10s}'
        )
        print(
            f'{"Treatment Elect. Use":<25s}{f"{m.fs.treatment.costing.aggregate_flow_electricity():<25.1f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Electrical %":<25s}{f"{-100*(m.fs.energy.pv.electricity()/m.fs.treatment.costing.aggregate_flow_electricity()):<25.1f}"}{"%":<10s}'
        )
        print(
            f'{"Plant Lifetime":<25s}{f"{m.fs.sys_costing.plant_lifetime():<25.1f}"}{"yrs":<10s}'
        )
        title = f'\n{"=======> PV SYSTEM RESULTS <=======":^60}\n'
        print(title)
        print(header)
        print(
            f'{"PV Capital Cost":<24s}{f"${m.fs.energy.costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"PV Design Size":<25s}{f"{m.fs.energy.pv.design_size():<25,.0f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Annual Generation":<25s}{f"{pyunits.convert(-1*m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)():<25,.0f}"}{"kWh/yr":<10s}'
        )
        print(
            f'{"PV Nameplate Capacity":<25s}{f"{pyunits.convert(-1*m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25,.0f}"}{"MWh/yr":<10s}'
        )
        print(
            f'{"PV Avg. Gen":<25s}{f"{-1*m.fs.energy.pv.electricity():<25,.0f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Land Required":<25s}{f"{m.fs.energy.pv.land_req():<25.1f}"}{"acres":<10s}'
        )
        print(
            f'{"PV Land Cost":<24s}{f"${m.fs.energy.pv.costing.land_cost():<25,.0f}"}{"$":<10s}'
        )
        
        title = f'\n{"=======> RO SYSTEM RESULTS <=======":^60}\n'
        print(title)
        print(header)
        print(
            f'{"Tot. Treat. Capital":<24s}{f"${m.fs.treatment.costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"RO Capital Cost":<24s}{f"${m.fs.treatment.ro.costing.capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"Pump Capital Cost":<24s}{f"${m.fs.treatment.p1.costing.capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"ERD Capital Cost":<24s}{f"${m.fs.treatment.erd.costing.capital_cost():<25,.0f}"}{"$":<10s}'
        )

        print(
            f'{"RO Fixed Op. Cost":<24s}{f"${m.fs.treatment.ro.costing.fixed_operating_cost():<25,.0f}"}{"$/yr":<10s}'
        )
        print(
            f'{"Pumping Power":<25s}{f"{m.fs.treatment.p1.control_volume.work[0]()/1000:<25,.0f}"}{"kW":<10s}'
        )
        print(
            f'{"RO Capital Cost":<24s}{f"${m.fs.treatment.ro.costing.capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"RO Operating Cost":<24s}{f"${m.fs.treatment.costing.total_operating_cost():<25,.0f}"}{"$/yr":<10s}'
        )

    print(
        f'{"RO Pressure":<25s}{f"{pyunits.convert(ro.inlet.pressure[0], to_units=pyunits.psi)():<25.1f}"}{"psi":<10s}'
    )
    print(f'{"Membrane Area":<25s}{f"{ro.area():<25.1f}"}{"m2":<10s}')
    print(f'{"Flux":<25s}{f"{flux_lmh:<25.1f}"}{"LMH":<10s}')
    print(
        f'{"Vol. Recovery":<25s}{f"{100 * ro.recovery_vol_phase[0, liq]():<25.1f}"}{"%":<10s}'
    )
    print(f'{"Flow In":<25s}{f"{prop_in.flow_vol():<25.3f}"}{"m3/s":<10s}')
    print(
        f'{"Flow In [MGD]":<25s}{f"{pyunits.convert(prop_in.flow_vol, to_units=pyunits.Mgallons/pyunits.day)():<25.2f}"}{"MGD":<10s}'
    )
    print(f'{"Flow Out":<25s}{f"{prop_perm.flow_vol():<25.3f}"}{"m3/s":<10s}')
    print(
        f'{"Feed Velocity":<25s}{f"{(value(ro.feed_side.velocity[0, 0])):<25.3f}"}{"m3/s":<10s}'
    )
    print(
        f'{"Conc. In":<25s}{f"{pyunits.convert(prop_in.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25,.0f}"}{"mg/L":<10s}'
    )
    print(
        f'{"Conc. Reject":<25s}{f"{pyunits.convert(prop_out.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25,.0f}"}{"mg/L":<10s}'
    )
    print(
        f'{"Conc. Perm":<25s}{f"{pyunits.convert(prop_perm.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25.1f}"}{"mg/L":<10s}'
    )
    print(
        f'{"Pump Pressure":<25s}{f"{pyunits.convert(erd.inlet.pressure[0], to_units=pyunits.psi)():<25.1f}"}{"psi":<10s}'
    )
    print(
        f'{"ERD Pressure":<25s}{f"{pyunits.convert(m.fs.treatment.p1.outlet.pressure[0], to_units=pyunits.psi)():<25.1f}"}{"psi":<10s}'
    )
    print(
        f'{"ERD Power Recovered":<25s}{f"{-1 * erd.work_mechanical[0]() * 1e-3:<25.1f}"}{"kW":<10s}'
    )


def display_pv_results(m, sep="."):
    line = f'\n{f"{sep*60}":<60s}'
    header = f'{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<10s}'
    print(line)
    title = f'\n{"=======> m.fs.energy.pv.costing <=======":^60}\n'
    print(title)
    print(header)
    print(
        f'{"PV Total Capital Cost":<34s}{f"${m.fs.energy.costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
    )
    print(
        f'{"PV Agg Cap Cost":<34s}{f"${m.fs.energy.costing.aggregate_capital_cost():<25,.0f}"}{"$":<10s}'
    )
    print(
        f'{"PV Factor Tot Investment Cost":<34s}{f"${m.fs.energy.costing.factor_total_investment():<25.4f}"}{"dimless":<10s}'
    )
    print(
        f'{"PV Annual Gen":<35s}{f"{pyunits.convert(-1*m.fs.energy.pv.electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25,.0f}"}{"MWh/yr":<10s}'
    )
    title = f'\n{"=======> m.fs.energy.costing or m.fs.sys_costing <=======":^60}\n'
    print(title)
    print(header)
    
    print(
        f'{"PV Total Capital Cost":<34s}{f"${m.fs.energy.costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
    )
    print(
            f'{"PV Annual Gen":<25s}{f"{m.fs.energy.pv.costing.annual_generation():<25.0f}"}{"MWh/yr":<25s}'
        )
    print(
        f'{"PV Factor Cap Annualization":<35s}{f"{m.fs.sys_costing.factor_capital_annualization():<25.1f}"}{"unitless":<10s}'
    )
    print(
        f'{"PV Util Factor":<35s}{f"{m.fs.sys_costing.utilization_factor():<25.4f}"}{"dimless":<10s}'
    )
    print(
        f'{"LCOE":<34s}{f"${(m.fs.sys_costing.LCOE()):<25.3f}"}{"$/kWh":<10s}'
    )
    print("\n")

def display_costing_breakdown(m, sep="."):
    header = f'{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<10s}'
    title = f'\n{"=======> Electricity  <=======":^60}\n'
    print(title)
    print(header)

    for f in m.fs.energy.costing.used_flows:
        print(f'{f"En. Flow [{f}]":<35s}{value(getattr(m.fs.energy.costing, f"aggregate_flow_{f}")):<25.1f}{"kW":<10s}')
    for f in m.fs.treatment.costing.used_flows:
        print(f'{f"Treat Flow [{f}]":<35s}{value(getattr(m.fs.treatment.costing, f"aggregate_flow_{f}")):<25.1f}{"kW":<10s}')
    
    print(f'{f"System Flow [electricity]":<35s}{value(m.fs.sys_costing.aggregate_flow_electricity):<25.1f}{"kW":<10s}')
    
    print(f'{f"Treat Total [{f}]":<35s}{pyunits.convert(m.fs.treatment.costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25,.0f}{"MWh/yr":<10s}')
    print(f'{f"System Total [{f}]":<35s}{pyunits.convert(m.fs.sys_costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25,.0f}{"MWh/yr":<10s}')
    print(f'{"Energy Op Cost Electric":<34s}{f"${m.fs.energy.costing.total_electric_operating_cost():<25,.0f}"}{"$/yr":<10s}')
    print(f'{"Treat Op Cost Electric":<34s}{f"${m.fs.treatment.costing.total_electric_operating_cost():<25,.0f}"}{"$/yr":<10s}')
    print(f'{"System Op Cost Electric":<34s}{f"${m.fs.sys_costing.total_electric_operating_cost():<25,.0f}"}{"$/yr":<10s}')
    print(f'{"System Agg Flow Electric":<35s}{f"{m.fs.sys_costing.aggregate_flow_electricity():<25.2f}"}{"kW":<10s}')
    print(f'{"System Net Electric":<35s}{f"{m.fs.sys_costing.net_flow_electricity():<25.2f}"}{"kW":<10s}')
    print(f'{"System Net Electric":<35s}{f"{pyunits.convert(m.fs.sys_costing.net_flow_electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25.2f}"}{"MWh/yr":<10s}')
    print(f'{"System Agg Electric":<35s}{f"{pyunits.convert(m.fs.sys_costing.aggregate_flow_electricity, to_units=pyunits.kWh/pyunits.year)()/1000:<25.2f}"}{"MWh/yr":<10s}')
    print(f'{"Net Electric Cost":<34s}{f"${m.fs.sys_costing.net_flow_electricity_cost():<25,.0f}"}{"$/yr":<10s}')
    print(f'{"Total Op Cost":<34s}{f"${m.fs.sys_costing.total_operating_cost():<25,.0f}"}{"$/yr":<10s}')


def generate_costing_report(m, filepath = None):
    costing_data = {"Attribute":'Value'}
    for v in m.fs.energy.pv.costing.component_data_objects(Var, descend_into=True):
        costing_data[str(v)] = value(v)
    for v in m.fs.energy.costing.component_data_objects(Var, descend_into=True):
        costing_data[str(v)] = value(v)
    for v in m.fs.treatment.costing.component_data_objects(Var, descend_into=True):
        costing_data[str(v)] = value(v)
    for v in m.fs.sys_costing.component_data_objects(Var, descend_into=True):
        costing_data[str(v)] = value(v)

    if filepath != None:
        with open(filepath, 'w') as f:
            for key in costing_data.keys():
                f.write("%s,%s\n"%(key,costing_data[key]))

    return costing_data

def extract_values(block):
    costing_data = dict.fromkeys(['key_0','key_1','key_2','key_3'], None)
    for v in block.component_data_objects(Var, descend_into=False):
        keys = str(v).split('.')
        for idx, key in enumerate(keys[:-1]):
            costing_data['key_'+str(idx)] = key
        costing_data[keys[-1]] = value(v)
    return costing_data

def generate_detailed_costing_report(m, level, filepath = None):
    dict_track = []
    if level == 'Process':
        cost_blocks = [m.fs.energy.costing, m.fs.treatment.costing, m.fs.sys_costing]
    else:
        cost_blocks = [m.fs.energy.pv.costing, m.fs.energy.costing.pv_surrogate]
    for block in cost_blocks:
        dict_track.append(extract_values(block))
    return dict_track