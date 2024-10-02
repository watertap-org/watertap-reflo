from idaes.core.solvers import get_solver
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_SOA import (
    build_system,
    add_connections,
    add_constraints,
    set_operating_conditions,
    init_system,
    add_costing,
    optimize,
    solve,
    print_all_results,
    )

def build():

    m = build_system()
    add_connections(m)
    add_constraints(m)
    set_operating_conditions(m)
    init_system(m)
    add_costing(m)
    optimize(m, ro_mem_area=None)
    return m

def solve_system(m, solver=None, use_model_state_storage=False, **kwargs):
  
    results = solve(m, solver=solver, raise_on_failure=True)

    print("solved okay!")
    return results


if __name__ == "__main__":
    m = build()
    results = solve_system(m)
    print_all_results(m)
