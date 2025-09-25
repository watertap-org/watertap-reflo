from watertap_contrib.reflo.flowsheets.KBHDP.KBHDP_ZLD import *

def main(Qin=4, ro_recovery=0.8, md_water_recovery=0.7):

    # Build the reverse osmosis flowsheet
    m = build_zld_ro(ro_recovery=ro_recovery, Qin=Qin)
    
    # Add the membrane distillation unit model
    add_zld_md(m, md_water_recovery=md_water_recovery)
    
    results = solve(m)
    assert_optimal_termination(results)
    
    # Add the multi-effect crystallizer unit model
    add_zld_mec(m)
    
    # Connect all product streams from the unit models
    add_product_stream(m)
    
    # Add energy unit models - PV and CST
    build_energy(m)

    results = solve(m.fs.treatment.mec)
    assert_optimal_termination(results)
    results = solve(m)
    assert_optimal_termination(results)

    # Add treatment costing
    add_treatment_costing(m.fs.treatment, heat_price=0, electricity_price=0)
    results = solve(m)
    assert_optimal_termination(results)

 
    # Add REFLOSystem costing
    add_system_costing(m)
    m.fs.energy.cst.unit.system_capacity.unfix()
    m.fs.energy.pv.unit.system_capacity.unfix()

    
    iscale.calculate_scaling_factors(m)

    results = solve(m)
    assert_optimal_termination(results)

    # Update grid fraction of electricity
    m.fs.costing.frac_heat_from_grid.fix(0.9)
    m.fs.costing.frac_elec_from_grid.fix(0.5)

    results = solve(m)
    assert_optimal_termination(results)

    # Update grid fraction of heat
    m.fs.costing.frac_heat_from_grid.fix(0.5)
    
    # Add objective function to optimize the LCOT
    m.fs.obj = Objective(expr=m.fs.costing.LCOT)

    # Solve!
    results = solve(m)
    assert_optimal_termination(results)

    return m


if __name__ == "__main__":
    m = main()
