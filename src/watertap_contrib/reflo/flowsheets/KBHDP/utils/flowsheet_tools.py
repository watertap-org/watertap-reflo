import math
import numpy as np

from pyomo.environ import (
    Var,
    Constraint,
    assert_optimal_termination,
    check_optimal_termination,
    units as pyunits,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import *

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.solvers import get_solver

__all__ = [
    "solve",
    "check_jac",
    "calc_scale",
    "print_fixed_and_unfixed_vars",
    "breakdown_dof",
]


def solve(m, solver=None, tee=False, raise_on_failure=True, debug=False):

    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        if debug:
            print("\n--------- CHECKING JACOBIAN ---------\n")
            print("\n--------- WHOLE FLOWSHEET ---------\n")
            check_jac(m)
            if hasattr(m.fs, "treatment"):
                print("\n--------- TREATMENT ---------\n")
                check_jac(m.fs.treatment)
            print("\n--------- ENERGY ---------\n")
            if hasattr(m.fs, "energy"):
                check_jac(m.fs.energy)

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(m)
        print(f'\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print(f'\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print(f'\n{"=======> CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)

        raise RuntimeError(msg)
    else:
        print("\n--------- FAILED SOLVE!!! ---------\n")
        print(msg)
        assert False


def check_jac(m, print_extreme_jacobian_values=True):
    def calc_scale(value):
        if value != 0 or value != None:
            try:
                return round(-1 * math.log(abs(value), 10), 1)
            except:
                return 0
        else:
            return 0

    jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(m, min_scale=1e-8)
    try:
        cond_number = iscale.jacobian_cond(m, jac=jac_scaled) / 1e10
        print("--------------------------")
        print("COND NUMBER:", cond_number)
    except:
        print("Cond number failed")
        cond_number = None
    if print_extreme_jacobian_values:
        print("--------------------------")
        print("Extreme Jacobian entries:")
        extreme_entries = iscale.extreme_jacobian_entries(
            m, jac=jac_scaled, nlp=nlp, zero=1e-20, large=100
        )
        for val, con, var in extreme_entries:
            if val >= 100:
                print(val, con.name, var.name, value(var), calc_scale(value(var)))
        print("--------------------------")
        print("Extreme Jacobian columns:")
        extreme_cols = iscale.extreme_jacobian_columns(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, var in extreme_cols:
            if val >= 100:
                print(val, var.name)
        print("------------------------")
        print("Extreme Jacobian rows:")
        extreme_rows = iscale.extreme_jacobian_rows(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, con in extreme_rows:
            if val >= 100:
                print(val, con.name)
    return cond_number


def calc_scale(value):
    return math.floor(math.log(value, 10))


def print_fixed_and_unfixed_vars(blk):
    fixed_vars = [
        (v.name, v.value)
        for v in blk.component_data_objects(ctype=Var, active=True, descend_into=False)
        if v.fixed
    ]
    unfixed_vars = [
        (v.name, v.value)
        for v in blk.component_data_objects(ctype=Var, active=True, descend_into=False)
        if v.fixed is False
    ]
    constraints = [
        c.name
        for c in blk.component_data_objects(
            ctype=Constraint, active=True, descend_into=True
        )
    ]
    print(f"{blk.name} Degrees of Freedom: {degrees_of_freedom(blk)}")
    print(f"Fixed Vars: ({len(fixed_vars)})")
    for v in fixed_vars:
        print(f"   {v[0]}: {v[1]}")

    print(f"Unfixed Vars: ({len(unfixed_vars)})")
    for v in unfixed_vars:
        print(f"   {v[0]}: {v[1]}")
    print(f"Constraints: ({len(constraints)})")
    for c in constraints:
        print(f"   {c}")


def display_dof_breakdown(blk, decend=False):
    w = 30
    print(
        "\n\n-------------------- DEGREE OF FREEDOM BREAKDOWN --------------------\n\n"
    )
    print(f'{"BLOCK":<40s}{"DEGREES OF FREEDOM":<{w}s}')
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(f"{v.name:<40s}{degrees_of_freedom(v)}")


def breakdown_dof(blk, detailed=False):
    all_vars = [
        v for v in blk.component_data_objects(ctype=Var, active=True, descend_into=True)
    ]
    equalities = [c for c in activated_equalities_generator(blk)]
    active_vars = variables_in_activated_equalities_set(blk)
    fixed_active_vars = fixed_variables_in_activated_equalities_set(blk)
    unfixed_active_vars = unfixed_variables_in_activated_equalities_set(blk)
    print("\n ===============DOF Breakdown================\n")
    print(f"Degrees of Freedom: {degrees_of_freedom(blk)}")

    if detailed:
        print(f"Activated Variables: ({len(active_vars)})")
        for v in active_vars:
            print(f"   {v}")
        print(f"Activated Equalities: ({len(equalities)})")
        for c in equalities:
            print(f"   {c}")

        print(f"Fixed Active Vars: ({len(fixed_active_vars)})")
        for v in fixed_active_vars:
            print(f"   {v}")

        print(f"Unfixed Active Vars: ({len(unfixed_active_vars)})")
        for v in unfixed_active_vars:
            print(f"   {v}")
        print("\n")
    print(f" {f' Active Vars':<30s}{len(active_vars)}")
    print(f"{'-'}{f' Fixed Active Vars':<30s}{len(fixed_active_vars)}")
    print(f"{'-'}{f' Activated Equalities':<30s}{len(equalities)}")
    print(f"{'='}{f' Degrees of Freedom':<30s}{degrees_of_freedom(blk)}")


def report_MCAS_stream_conc(m, stream):
    solute_set = m.fs.MCAS_properties.solute_set
    print(f"\n\n-------------------- {stream} CONCENTRATIONS --------------------\n\n")
    print(f'{"Component":<15s}{"Conc.":<10s}{"Units":10s}')
    for i in solute_set:
        print(
            f"{i:<15s}: {stream.conc_mass_phase_comp['Liq', i].value:<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp['Liq', i])}"
        )
    print(
        f'{"Overall TDS":<15s}: {sum(value(stream.conc_mass_phase_comp["Liq", i]) for i in solute_set):<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp["Liq", "Ca_2+"])}'
    )
    print(
        f"{'Vol. Flow Rate':<15s}: {stream.flow_mass_phase_comp['Liq', 'H2O'].value:<10.3f}{pyunits.get_units(stream.flow_mass_phase_comp['Liq', 'H2O'])}"
    )


def get_poorly_scaled_vars(blk):
    for var, sv in iscale.badly_scaled_var_generator(blk):
        cur_scale = iscale.get_scaling_factor(var)
        try:
            need_scale = math.log(1 / value(var), 10)
        except ValueError:
            need_scale = 1 / value(var)
        if cur_scale is None:
            cur_scale = None
            delta = None
        else:
            cur_scale = math.log(cur_scale, 10)

            delta = cur_scale - need_scale
        print(
            "Poorly_sclaed_var",
            var,
            "current_scale",
            cur_scale,
            "need_scale",
            need_scale,
            "delta",
            delta,
        )


def auto_scale_bad_vars(
    m,
    blk,
    display=False,
    auto_scale_all=False,
    small_values=1e-10,
    upscale_small_values=None,
    update_scaling=True,
    display_with_out_scale=False,
    rescale_magnitude=0,
    loosen=0,
):
    if auto_scale_all is False:
        for var, sv in iscale.badly_scaled_var_generator(blk):
            # print("poor scaling:", var, sv)
            sc = calc_scale(sv) - loosen
            if update_scaling:
                iscale.set_scaling_factor(
                    var,
                    10**sc,
                )
            if display:
                print("auto_scaled_bar_var", var, sv, 10**sc)
    else:
        rescaled_var_count = 0
        for var in blk.component_data_objects(Var):
            # print(var)
            _add = 0
            if (
                (var.value) != 0
                and var.value != None
                and display_with_out_scale == False
            ):
                val = abs(var.value)
                cur_scale = iscale.get_scaling_factor(var)
                if cur_scale == None:
                    cur_scale = 1
                sc = calc_scale(val) + _add - loosen
                if cur_scale == None:
                    cur_scale = sc + rescale_magnitude * 2
                if abs(sc - np.log10(cur_scale)) > rescale_magnitude and update_scaling:
                    # print(sc, np.log10(cur_scale), abs(10**sc - cur_scale))
                    iscale.set_scaling_factor(
                        var,
                        10**sc,
                    )
                    rescaled_var_count += 1
                if display:
                    print("auto_scaled", var, var.value, 10**sc)
            elif iscale.get_scaling_factor(var) is None:
                print("Var no scale", var)

        print("auto scaled", rescaled_var_count)
    iscale.calculate_scaling_factors(m)


def print_system_scaling_report(m):
    badly_scaled_var_list = iscale.list_badly_scaled_variables(m)
    if len(badly_scaled_var_list) > 0:
        print("Variables are not scaled well")
        print(
            f'{"Variable":<83s}{"Val":<15s}{"Val Scale":<10s}{"SF":<10s}{"Diff":<10s}'
        )
        print("Treatment:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "treatment"
        ]
        print("Energy:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "energy"
        ]
        print("Costing:")
        [
            print(
                f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
            )
            for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
            if var.name.split(".")[1] == "costing"
        ]
    else:
        print("Variables are scaled well")


def display_unfixed_vars(blk, report=True):
    print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
    print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
    print(f"{blk.name:<40s}{number_unused_variables(blk)}")
    for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
        print(f"{v.name:<40s}{number_unused_variables(v)}")
        for v2 in unused_variables_set(v):
            print(f"\t{v2.name:<40s}")


def get_scaling_factors(m):
    for var in [m.fs.treatment.costing.aggregate_flow_electricity]:
        val = value(var)
        scale = calc_scale(val)
        sf = iscale.get_scaling_factor(var)
        if sf is None:
            sf = scale
            iscale.set_scaling_factor(var, sf)

        print(f"{var.name:<50s}{val:<20.3f}{scale:<20.3f}{sf:<20.3f}")
