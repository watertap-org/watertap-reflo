from idaes.core.solvers import get_solver
from pyomo.environ import SolverFactory
from watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_3 import (
    main,
    save_results,
)


def sweep(sweep_type, water_recovery=0.5, heat_price=0.07, frac_heat_from_grid=0.01):
    m = main(
        water_recovery=water_recovery,
        heat_price=heat_price,
        frac_heat_from_grid=frac_heat_from_grid,
    )
    save_results(m, sweep_type=sweep_type)


if __name__ == "__main__":
    water_recovery_sweep = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    heat_price_sweep = [0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1]
    for recovery in water_recovery_sweep:
        sweep(
            "water_recovery_1",
            water_recovery=recovery,
            heat_price=0.08,
            frac_heat_from_grid=0.01,
        )

    # for heat_price in heat_price_sweep:
    #     sweep('heat_price_1',water_recovery=0.5,heat_price = heat_price,frac_heat_from_grid=0.01)
