from pyomo.environ import units as pyunits

__all__ = ["display_ro_pv_results"]


def display_ro_pv_results(m, sep="."):
    ro = m.fs.treatment.ro
    erd = m.fs.treatment.erd
    liq = "Liq"
    nacl = "NaCl"
    header = f'{"PARAM":<25s}{"VALUE":<25s}{"UNITS":<25s}'
    prop_in = ro.feed_side.properties_in[0]
    prop_out = ro.feed_side.properties_out[0]
    prop_perm = ro.mixed_permeate[0]
    pv_cost = m.fs.energy.costing.photovoltaic
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
        print(
            f'{"LCOW":<25s}{f"{m.fs.sys_costing.LCOW():<25.4f}"}{"$/m3":<25s}'
        )
        if hasattr(m.fs.sys_costing, "LCOE"):
            print(
                f'{"LCOE":<25s}{f"{m.fs.sys_costing.LCOE():<25.4f}"}{"$/kWh":<25s}'
            )
        print(
            f'{"Total Capital Cost":<25s}{f"{m.fs.sys_costing.total_capital_cost():<25.2f}"}{"$":<25s}'
        )
        print(
            f'{"Total Op. Cost":<25s}{f"{m.fs.sys_costing.total_operating_cost():<25.2f}"}{"$/yr":<25s}'
        )
        print(
            f'{"SEC":<25s}{f"{m.fs.sys_costing.specific_electric_energy_consumption():<25.4f}"}{"kWh/m3":<25s}'
        )
        print(
            f'{"PV Avg. Elect. Gen":<25s}{f"{m.fs.energy.pv.electricity():<25.4f}"}{"kW":<25s}'
        )
        print(
            f'{"RO Electricity Use":<25s}{f"{m.fs.treatment.costing.aggregate_flow_electricity():<25.4f}"}{"kW":<25s}'
        )
        print(
            f'{"Overall Elect. Use":<25s}{f"{m.fs.sys_costing.aggregate_flow_electricity():<25.4f}"}{"kW":<25s}'
        )
        title = f'\n{"=======> PV SYSTEM RESULTS <=======":^60}\n'
        print(title)
        print(header)
        print(
            f'{"PV Capital Cost":<25s}{f"{m.fs.energy.pv.costing.capital_cost():<25.2f}"}{"$":<25s}'
        )
        print(
            f'{"PV Fixed Op.":<25s}{f"{m.fs.energy.pv.costing.fixed_operating_cost():<25.2f}"}{"$/yr":<25s}'
        )
        print(
            f'{"PV Fixed Op. by Capcity":<25s}{f"{pv_cost.fixed_operating_by_capacity():<25.2f}"}{"$/kW":<25s}'
        )
        print(
            f'{"PV Var Op.":<25s}{f"{m.fs.energy.pv.costing.variable_operating_cost():<25.2f}"}{"$/yr":<25s}'
        )
        print(
            f'{"PV Var Op. by Ann Gen.":<25s}{f"{pv_cost.variable_operating_by_generation():<25.2f}"}{"$/MWh":<25s}'
        )
        print(
            f'{"PV Annual Gen":<25s}{f"{m.fs.energy.pv.costing.annual_generation():<25.4f}"}{"MWh/yr":<25s}'
        )
        print(
            f'{"PV Nameplate Capacity":<25s}{f"{m.fs.energy.pv.costing.system_capacity():<25.4f}"}{"W":<25s}'
        )
        print(
            f'{"PV Land Required":<25s}{f"{m.fs.energy.pv.costing.land_area():<25.4f}"}{"acres":<25s}'
        )

        print(
            f'{"PV Avg. Gen":<25s}{f"{-1 * m.fs.energy.pv.electricity():<25.4f}"}{"kW":<25s}'
        )
        title = f'\n{"=======> RO SYSTEM RESULTS <=======":^60}\n'
        print(title)
        print(header)
        print(
            f'{"Tot. Treat. Capital":<25s}{f"{m.fs.treatment.costing.total_capital_cost():<25.4f}"}{"$":<25s}'
        )
        print(
            f'{"RO Capital Cost":<25s}{f"{m.fs.treatment.ro.costing.capital_cost():<25.2f}"}{"$":<25s}'
        )
        print(
            f'{"Pump Capital Cost":<25s}{f"{m.fs.treatment.p1.costing.capital_cost():<25.2f}"}{"$":<25s}'
        )
        print(
            f'{"ERD Capital Cost":<25s}{f"{m.fs.treatment.erd.costing.capital_cost():<25.2f}"}{"$":<25s}'
        )

        print(
            f'{"RO Fixed Op. Cost":<25s}{f"{m.fs.treatment.ro.costing.fixed_operating_cost():<25.2f}"}{"$/yr":<25s}'
        )
        print(
            f'{"Pumping Power":<25s}{f"{m.fs.treatment.p1.control_volume.work[0]():<25.2f}"}{"W":<25s}'
        )
        print(
            f'{"RO Capital Cost":<25s}{f"{m.fs.treatment.ro.costing.capital_cost():<25.2f}"}{"$":<25s}'
        )
        print(
            f'{"RO Operating Cost":<25s}{f"{m.fs.treatment.costing.total_operating_cost():<25.4f}"}{"$/yr":<25s}'
        )

    print(
        f'{"RO Pressure":<25s}{f"{pyunits.convert(ro.inlet.pressure[0], to_units=pyunits.psi)():<25.4f}"}{"psi":<25s}'
    )
    print(f'{"Membrane Area":<25s}{f"{ro.area():<25.4f}"}{"m2":<25s}')
    print(f'{"Flux":<25s}{f"{flux_lmh:<25.4f}"}{"LMH":<25s}')
    # print(f'{"Flux Check":<25s}{f"{flow_out_L() / ro.area():<25.4f}"}{"LMH":<25s}')
    print(
        f'{"Vol. Recovery":<25s}{f"{100 * ro.recovery_vol_phase[0, liq]():<25.4f}"}{"%":<25s}'
    )
    print(f'{"Flow In":<25s}{f"{prop_in.flow_vol():<25.4f}"}{"m3/s":<25s}')
    print(
        f'{"Flow In [MGD]":<25s}{f"{pyunits.convert(prop_in.flow_vol, to_units=pyunits.Mgallons/pyunits.day)():<25.4f}"}{"MGD":<25s}'
    )
    print(
        f'{"Flow Out":<25s}{f"{prop_perm.flow_vol():<25.4f}"}{"m3/s":<25s}'
    )
    print(
        f'{"Conc. In":<25s}{f"{pyunits.convert(prop_in.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25.4f}"}{"mg/L":<25s}'
    )
    print(
        f'{"Conc. Reject":<25s}{f"{pyunits.convert(prop_out.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25.4f}"}{"mg/L":<25s}'
    )
    print(
        f'{"Conc. Perm":<25s}{f"{pyunits.convert(prop_perm.conc_mass_phase_comp[liq, nacl], to_units=pyunits.mg/pyunits.L)():<25.4f}"}{"mg/L":<25s}'
    )
    print(
        f'{"Pump Pressure":<25s}{f"{pyunits.convert(erd.inlet.pressure[0], to_units=pyunits.psi)():<25.4f}"}{"psi":<25s}'
    )
    print(
        f'{"ERD Pressure":<25s}{f"{pyunits.convert(m.fs.treatment.p1.outlet.pressure[0], to_units=pyunits.psi)():<25.4f}"}{"psi":<25s}'
    )
    print(
        f'{"ERD Power Recovered":<25s}{f"{-1 * erd.work_mechanical[0]() * 1e-3:<25.4f}"}{"kW":<25s}'
    )

