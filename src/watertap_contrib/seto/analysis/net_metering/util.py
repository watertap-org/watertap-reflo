from pyomo.environ import value, Var, units as pyunits
import os

absolute_path = os.path.dirname(__file__)

__all__ = ["display_ro_pv_results", "display_pv_results"]

def display_ro_pv_results(m, sep="."):
    ro = m.fs.treatment.ro
    erd = m.fs.treatment.erd
    liq = "Liq"
    nacl = "NaCl"
    header = f'{"PARAM":<25s}{"VALUE":<25s}{"UNITS":<10s}'
    prop_in = ro.feed_side.properties_in[0]
    prop_out = ro.feed_side.properties_out[0]
    prop_perm = ro.mixed_permeate[0]
    pv_cost = m.fs.energy.costing.pv_surrogate
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
            f'{"RO Electricity Use":<25s}{f"{m.fs.treatment.costing.aggregate_flow_electricity():<25.1f}"}{"kW":<10s}'
        )
        print(
            f'{"Overall Elect. Use":<25s}{f"{m.fs.sys_costing.aggregate_flow_electricity():<25.1f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Electrical %":<25s}{f"{((-100 * m.fs.energy.pv.electricity())/m.fs.sys_costing.aggregate_flow_electricity()):<25.1f}"}{"%":<10s}'
        )
        print(
            f'{"Plant Lifetime":<25s}{f"{m.fs.sys_costing.plant_lifetime():<25.1f}"}{"yrs":<10s}'
        )
        title = f'\n{"=======> PV SYSTEM RESULTS <=======":^60}\n'
        print(title)
        print(header)
        print(
            f'{"PV Capital Cost":<24s}{f"${m.fs.energy.pv.costing.capital_cost():<25,.0f}"}{"$":<10s}'
        )
        print(
            f'{"PV Fixed Op.":<24s}{f"${m.fs.energy.pv.costing.fixed_operating_cost():<25.0f}"}{"$/yr":<10s}'
        )
        print(
            f'{"PV Fixed Op. by Capcity":<24s}{f"${pv_cost.fixed_operating_by_capacity():<25.2f}"}{"$/kW":<10s}'
        )
        print(
            f'{"PV Var Op.":<24s}{f"${m.fs.energy.pv.costing.variable_operating_cost():<25.0f}"}{"$/yr":<10s}'
        )
        print(
            f'{"PV Var Op. by Ann Gen.":<24s}{f"${pv_cost.variable_operating_by_generation():<25.0f}"}{"$/MWh":<10s}'
        )
        print(
            f'{"PV Annual Gen":<25s}{f"{m.fs.energy.pv.costing.annual_generation()/1000:<25.0f}"}{"MWh/yr":<10s}'
        )
        print(
            f'{"PV Nameplate Capacity":<25s}{f"{m.fs.energy.pv.costing.system_capacity()/1000:<25,.0f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Land Required":<25s}{f"{m.fs.energy.pv.costing.land_area():<25.4f}"}{"acres":<10s}'
        )

        print(
            f'{"PV Avg. Gen":<25s}{f"{-1*m.fs.energy.pv.electricity():<25.0f}"}{"kW":<10s}'
        )
        print(
            f'{"PV Design Size":<25s}{f"{m.fs.energy.pv.design_size():<25.1f}"}{"kW":<10s}'
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
            f'{"Pumping Power":<25s}{f"{m.fs.treatment.p1.control_volume.work[0]():<25,.0f}"}{"W":<10s}'
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
    # print(f'{"Flux Check":<25s}{f"{flow_out_L() / ro.area():<25.4f}"}{"LMH":<10s}')
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
    # print(f'{"Feed Velocity":<25s}{f"{prop_in.velocity:<25.4f}"}{"m/s":<10s}')
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
        f'{"PV Total Capital Cost":<34s}{f"${m.fs.energy.costing.total_capital_cost():<25,.2f}"}{"$":<10s}'
    )
    print(
        f'{"PV Total Operating Cost":<34s}{f"${m.fs.energy.costing.total_operating_cost():<25,.2f}"}{"$":<10s}'
    )
    print(
        f'{"PV Agg Cap Cost":<34s}{f"${m.fs.energy.costing.aggregate_capital_cost():<25.4f}"}{"$":<10s}'
    )
    print(
        f'{"PV Factor Tot Investment Cost":<34s}{f"${m.fs.energy.costing.factor_total_investment():<25.4f}"}{"dimless":<10s}'
    )
    print(
        f'{"PV Annual Gen":<35s}{f"{m.fs.energy.pv.costing.annual_generation()/1000:<25.1f}"}{"MWh/yr":<10s}'
    )
    # print(
    #     f'{"LCOE":<34s}{f"${((m.fs.energy.pv.costing.capital_cost()*m.fs.sys_costing.factor_capital_annualization())+(m.fs.energy.pv.costing.fixed_operating_cost()))/(1000*(m.fs.energy.pv.costing.annual_generation())):<25.4f}"}{"$/kWh":<10s}'
    # )
    title = f'\n{"=======> m.fs.energy.costing or m.fs.sys_costing <=======":^60}\n'
    print(title)
    print(header)
    
    print(
        f'{"PV Total Capital Cost":<34s}{f"${m.fs.energy.costing.total_capital_cost():<25,.0f}"}{"$":<10s}'
    )
    print(
        f'{"PV Total Op Cost":<34s}{f"${m.fs.energy.costing.total_operating_cost():<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"PV Annual Gen":<35s}{f"{m.fs.sys_costing.annual_energy_generated():<25,.0f}"}{"kWh/yr":<10s}'
    )
    print(
        f'{"PV Factor Cap Annualization":<35s}{f"{m.fs.sys_costing.factor_capital_annualization():<25.4f}"}{"unitless":<10s}'
    )
    print(
        f'{"PV Util Factor":<35s}{f"{m.fs.sys_costing.utilization_factor():<25.4f}"}{"dimless":<10s}'
    )
    print(
        f'{"LCOE":<34s}{f"${((m.fs.energy.costing.total_capital_cost()*m.fs.sys_costing.factor_capital_annualization())+(m.fs.energy.costing.total_operating_cost()))/(m.fs.sys_costing.annual_energy_generated()):<25.4f}"}{"$/kWh":<10s}'
    )
    print("\n")
    print(
        f'{"MLC Op Cost":<34s}{f"${m.fs.energy.costing.maintenance_labor_chemical_operating_cost():<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"Agg Fixed Op Cost":<34s}{f"${m.fs.energy.costing.aggregate_fixed_operating_cost():<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"Agg Var Op Cost":<34s}{f"${m.fs.energy.costing.aggregate_variable_operating_cost():<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"Treatment Agg Flow Costs [elec]":<34s}{f"${(value(list(m.fs.treatment.costing.aggregate_flow_costs.values())[0])):<25,.0f}"}{"$/yr":<10s}'
    )
    print(
        f'{"Energy Agg Flow Costs [elec]":<34s}{f"${(value(list(m.fs.energy.costing.aggregate_flow_costs.values())[0])):<25,.0f}"}{"$/yr":<10s}'
    )

    title = f'\n{"=======> m.fs.sys_costing.LCOE() <=======":^60}\n'
    print(title)
    print(header)
    print(f'{"LCOE":<34s}{f"${m.fs.sys_costing.LCOE():<25.4f}"}{"$/kWh":<10s}')
    print(f'{"Treatment Agg Flow Electricity":<35s}{f"{m.fs.treatment.costing.aggregate_flow_electricity():<25.1f}"}{"kW":<10s}')
    print(f'{"PV Agg Flow Electricity":<35s}{f"{m.fs.energy.costing.aggregate_flow_electricity():<25.1f}"}{"kW":<10s}')
    print("\n\n")
