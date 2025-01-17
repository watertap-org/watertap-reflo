import math
import idaes.core.util.scaling as iscale
from pyomo.environ import (
    Var,
    Constraint,
    assert_optimal_termination,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
import time
from watertap.core.util.model_diagnostics.infeasible import *
import numpy as np

__all__ = [
    "check_jac",
    "calc_scale",
    "print_fixed_and_unfixed_vars",
    "breakdown_dof",
    "get_poorly_scaled_vars",
]


def check_jac(m, print_extreme_jacobian_values=True):
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
                print(
                    val,
                    con.name,
                    var.name,
                    value(var),
                    calc_scale(value(var)),
                    iscale.get_scaling_factor(var),
                )
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
    if value != 0 or value != None:
        try:
            return round(-1 * math.log(abs(value), 10), 1)
        except:
            return 0
    else:
        return 0


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
            # subsystem = v.name.split(".")[1]
            # if subsystem == "treatment":
            #     component = v.name.split(".")[2]
            #     if (component != "RO") and (component != "EC"):
            # print(f"   {v}")
            print(f"   {v}")
        print("\n")
    print(f" {f' Active Vars':<30s}{len(active_vars)}")
    print(f"{'-'}{f' Fixed Active Vars':<30s}{len(fixed_active_vars)}")
    print(f"{'-'}{f' Activated Equalities':<30s}{len(equalities)}")
    print(f"{'='}{f' Degrees of Freedom':<30s}{degrees_of_freedom(blk)}")

    # if degrees_of_freedom != 0:
    #     print("\nSuggested Variables to Fix:")
    #     unfixed_vars_without_constraint = [
    #         v for v in all_vars if v not in unfixed_active_vars
    #     ]
    #     for v in unfixed_vars_without_constraint:
    #         if v.fixed:
    #             print(f"   {v}       Fixed={v.fixed}")


def standard_solve(
    m,
    solver=None,
    tee=False,
    check_termination=False,
    expected_DOFs=None,
    debug=False,
    check_close_to_bounds=False,
    check_var_scailing=False,
    auto_rescale=False,
    check_jacobian=False,
    get_stats=True,
    symbolic_solver_labels=False,
    print_level=5,
    check_dofs=False,
    solver_name=None,
):
    print("SOLVING!")
    if check_dofs:
        dofs = degrees_of_freedom(m)
        print("DOFs:", dofs)
    else:
        dofs = 1
    solve_start = time.time()
    if expected_DOFs != None:
        assert dofs == expected_DOFs
        print("DOFS ok", expected_DOFs)

    if dofs >= 0:
        # check_jac(m)
        if solver is None:
            if solver_name is not None:
                solver = get_solver(solver_name)
                if solver_name == "cyipopt-watertap":
                    # cy_solver.options["bound_relax_factor"] = 0.0  # 1e-8
                    # cy_solver.options["constr_viol_tol"] = 1e-8
                    # cy_solver.options["acceptable_constr_viol_tol"] = 1e-8
                    # cy_solver.options["honor_original_bounds"] = "no"
                    # cy_solver.options["tol"] = 1e-8
                    solver.options["hessian_approximation"] = "limited-memory"
                    solver.options["linear_solver"] = "ma27"
            else:
                print("getting solver")
                solver = get_solver()
        solver.options["max_iter"] = 2000
        fs = m.find_component("fs")
        solver.options["print_level"] = print_level
        results = solver.solve(m, tee=tee)  # , symbolic_solver_labels=tee)

        if get_stats and fs is not None and fs.find_component("num_DOF") is not None:
            with open(tempfile, "r") as f:
                iters = 0
                solve_time = 0
                for line in f:
                    # print(line)
                    if line.startswith("Number of Iterations....:"):
                        tokens = line.split()
                        iters = int(tokens[3])
                    # elif line.s?tartswith(
                    #     "Total CPU secs in IPOPT (w/o function evaluations)   ="
                    # ):
                    #     tokens = line.split()
                    #     solve_time += float(tokens[9])
                    # elif line.startswith(
                    #     "Total CPU secs in NLP function evaluations           ="
                    # ):
                    #     tokens = line.split()
                    #     solve_time += float(tokens[8])
                    elif line.startswith("Objective...............:"):
                        tokens = line.split()
                        # print(tokens)
                        m.fs.solver_result_unscaled["Objective"].fix(float(tokens[2]))
                        m.fs.solver_result_scaled["Objective"].fix(float(tokens[1]))
                    elif line.startswith("Dual infeasibility......:"):
                        tokens = line.split()
                        # print(tokens, tokens[1])
                        m.fs.solver_result_unscaled["Dual infeasibility"].fix(
                            float(tokens[3])
                        )
                        m.fs.solver_result_scaled["Dual infeasibility"].fix(
                            float(tokens[2])
                        )
                    elif line.startswith("Constraint violation....:"):
                        tokens = line.split()
                        # print(tokens)
                        m.fs.solver_result_unscaled["Constraint violation"].fix(
                            float(tokens[3])
                        )
                        m.fs.solver_result_scaled["Constraint violation"].fix(
                            float(tokens[2])
                        )
                    elif line.startswith("Complementarity.........:"):
                        tokens = line.split()
                        # print(tokens)
                        m.fs.solver_result_unscaled["Complementarity"].fix(
                            float(tokens[2])
                        )
                        m.fs.solver_result_scaled["Complementarity"].fix(
                            float(tokens[1])
                        )
                    elif line.startswith("Overall NLP error.......:"):
                        tokens = line.split()
                        # print(tokens[1], tokens[2])
                        m.fs.solver_result_unscaled["Overall NLP error"].fix(
                            float(tokens[4])
                        )
                        m.fs.solver_result_scaled["Overall NLP error"].fix(
                            float(tokens[3])
                        )
                    elif line.startswith("Total seconds in IPOPT"):
                        tokens = line.split("=")
                        # print(tokens)
                        solve_time = float(tokens[-1])

                # m.fs.num_DOF.fix(dofs)
                m.fs.solver_time.fix(solve_time)
                m.fs.iteration.fix(iters)
                if check_dofs:
                    m.fs.num_DOF.fix(dofs)
                    m.fs.number_unfixed_variables_in_activated_equalities.fix(
                        number_unfixed_variables_in_activated_equalities(m)
                    )
                    m.fs.number_activated_equalities.fix(number_activated_equalities(m))

                print("iters:", iters, "time:", solve_time)
                # # report_statistics(m)
                # print("number_activated_constraints", number_activated_constraints(m))
                # print("number_unfixed_variables", number_unfixed_variables(m))
                # print(
                #     "number_unfixed_variables_in_activated_equalities",
                #     number_unfixed_variables_in_activated_equalities(m),
                # )
                # print(
                #     "number_activated_equalities",
                #     number_activated_equalities(m),
                # )
        if check_termination:
            assert_optimal_termination(results)
        if check_close_to_bounds:
            print("------vars_close_to_bound-tests---------")
            print_variables_close_to_bounds(m)
            print("------constraints_close_to_bound-tests---------")

            print_constraints_close_to_bounds(m)
        if check_var_scailing:
            print("------poor_scaling_vars-tests---------")
            get_poorly_scaled_vars(m)
            if auto_rescale:
                print("\n\n------auto_scaling_vars-tests---------\n\n")
                auto_scale_bad_vars(m, m.fs, display=True, auto_scale_all=True)
        if check_jacobian:
            print("------jac-tests---------")
            check_jac(m)
        try:
            assert_optimal_termination(results)
            succes = "Optimal solution found!"
            print("--------------------------")
        except:
            # print("------ifeasible_constraints-test---------")
            # print_infeasible_constraints(m, tol=1e-8)
            succes = "Solution NOT optimal !!!!!!!!!!!!!!"
            print("!-!-!-!-!-!-!-!-!-!-!-!-!")
        if m.find_component("fs.costing.LCOW") is not None:
            print("SOLVED, LCOW", value(m.fs.costing.LCOW), "$/m3", succes)
        else:
            print("Solved model, not costing to report")
        print("--------------------------")
    else:
        print("TO FEW DOFS NOT SOLVING!")
        results = None
        if check_termination:
            assert_optimal_termination(results)
    print("---------solve took: {}--------".format(time.time() - solve_start))
    return results


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

### ALL FUNCTIONS MOVED HERE FROM KBHDP_ZLD 1/17/2025 -KAS


# def propagate_state(arc, detailed=False):
#     _prop_state(arc)

#     if detailed:
#         print(
#             f"\nPropogation of {arc.source.name} to {arc.destination.name} successful."
#         )
#         arc.source.display()
#         print(arc.destination.name)
#         arc.destination.display()
#         print("\n")
# def display_system_stream_table(m):
#     treatment = m.fs.treatment
#     print("\n\n-------------------- SYSTEM STREAM TABLE --------------------\n\n")
#     print(
#         f'{"NODE":<20s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (G/L)":<20s}'
#     )
#     print(
#         f'{"Feed":<20s}{treatment.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{value(pyunits.convert(treatment.feed.properties[0.0].pressure, to_units=pyunits.bar)):<30.1f}'
#     )
#     display_flow_table(treatment.RO)
#     print("\n\n")


# def display_system_build(blk):
#     blocks = []
#     for v in blk.component_data_objects(ctype=Block, active=True, descend_into=False):
#         print(v)



# def report_MCAS_stream_conc(m, stream):
#     solute_set = m.fs.MCAS_properties.solute_set
#     print(f"\n\n-------------------- {stream} CONCENTRATIONS --------------------\n\n")
#     print(f'{"Component":<15s}{"Conc.":<10s}{"Units":10s}')
#     for i in solute_set:
#         print(
#             f"{i:<15s}: {stream.conc_mass_phase_comp['Liq', i].value:<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp['Liq', i])}"
#         )
#     print(
#         f'{"Overall TDS":<15s}: {sum(value(stream.conc_mass_phase_comp["Liq", i]) for i in solute_set):<10.3f}{pyunits.get_units(stream.conc_mass_phase_comp["Liq", "Ca_2+"])}'
#     )
#     print(
#         f"{'Vol. Flow Rate':<15s}: {stream.flow_mass_phase_comp['Liq', 'H2O'].value:<10.3f}{pyunits.get_units(stream.flow_mass_phase_comp['Liq', 'H2O'])}"
#     )


# def report_stream_ion_conc(m, stream):
#     print(f"\n\n-------------------- {stream} CONCENTRATIONS --------------------\n\n")
#     for ion_conc in stream.conc_mass_phase_comp:
#         print(
#             f"{ion_conc[1]:<15s}: {stream.conc_mass_phase_comp[ion_conc].value:<10.3f}{str(pyunits.get_units(stream.conc_mass_phase_comp[ion_conc]))}"
#         )


# def report_pump(m, pump):
#     print(f"\n\n-------------------- SYSTEM PUMP --------------------\n\n")
#     print(
#         f'{"Pump Pressure":<30s}{value(pyunits.convert(pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)):<10.1f}{"bar"}'
#     )
#     print(
#         f'{"Pump Work":<30s}{value(pyunits.convert(pump.control_volume.work[0], to_units=pyunits.kW)):<10.3f}{"kW"}'
#     )


# def print_system_scaling_report(m):
#     badly_scaled_var_list = iscale.list_badly_scaled_variables(m)
#     if len(badly_scaled_var_list) > 0:
#         print("Variables are not scaled well")
#         print(
#             f'{"Variable":<83s}{"Val":<15s}{"Val Scale":<10s}{"SF":<10s}{"Diff":<10s}'
#         )
#         print("Treatment:")
#         [
#             print(
#                 f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
#             )
#             for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
#             if var.name.split(".")[1] == "treatment"
#         ]
#         print("Energy:")
#         [
#             print(
#                 f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
#             )
#             for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
#             if var.name.split(".")[1] == "energy"
#         ]
#         print("Costing:")
#         [
#             print(
#                 f"   {var.name:<80s}{val:<15.1f}{-1*calc_scale(val):<10.1f}{-1*calc_scale(iscale.get_scaling_factor(var)):<10.1f}"
#             )
#             for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)
#             if var.name.split(".")[1] == "costing"
#         ]

#         # for var in badly_scaled_var_list:
#         #     keys = var[0].name.split(".")
#         #     val_scale = -1 * calc_scale(var[1])
#         #     sf_scale = -1 * calc_scale(iscale.get_scaling_factor(var[0]))
#         #     scale_diff = val_scale - sf_scale
#         #     print(f"{var[0].name:<80s}{val_scale:<10.1f}{sf_scale:<10.1f}{scale_diff:<10.1f}")
#         #     print(f"{var[0].name:<80s}{var[1]}")
#         # [print(i[0].name, i[1]) for i in badly_scaled_var_list]
#         # [print(f'Variable: {var.name:<80s} Value Scale:{-1*calc_scale(val):<5.1f} Scale Factor:{-1*calc_scale(iscale.get_scaling_factor(var))}') for var, val in iscale.list_badly_scaled_variables(m, include_fixed=True)]
#     else:
#         print("Variables are scaled well")




# def display_costing_breakdown(m):
#     header = f'\n{"PARAM":<35s}{"VALUE":<25s}{"UNITS":<25s}'
#     print(header)
#     print(
#         f'{"Product Flow":<35s}{f"{value(pyunits.convert(m.fs.treatment.product.properties[0].flow_vol, to_units=pyunits.m **3 * pyunits.yr ** -1)):<25,.1f}"}{"m3/yr":<25s}'
#     )
#     print(f'{"LCOW":<34s}{f"${m.fs.treatment.costing.LCOW():<25.3f}"}{"$/m3":<25s}\n')

#     print_EC_costing_breakdown(m.fs.treatment.EC)
#     print_UF_costing_breakdown(m.fs.treatment.UF)
#     print_RO_costing_breakdown(m.fs.treatment.RO)
#     # print_DWI_costing_breakdown(m.fs.treatment.DWI)
#     # print_PV_costing_breakdown(m.fs.energy.pv)


# def set_prob_for_box_solve(m):
#     treatment = m.fs.treatment
#     # m.fs.water_recovery.unfix()
#     # m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(75e5)
#     treatment.pump.control_volume.properties_out[0].pressure.unfix()
#     for idx, stage in treatment.RO.stage.items():
#         stage.module.recovery_vol_phase[0.0, "Liq"].fix(0.5)


# def box_solve_problem(m):
#     set_prob_for_box_solve(m)
#     breakdown_dof(m, detailed=False)
#     fsTools.standard_solve(
#         m,
#         tee=True,
#         check_close_to_bounds=True,
#         check_var_scailing=True,
#         check_dofs=True,
#         expected_DOFs=0,
#     )



# def get_scaling_factors(m):
#     for var in [m.fs.treatment.costing.aggregate_flow_electricity]:
#         val = value(var)
#         scale = calc_scale(val)
#         sf = get_scaling_factor(var)
#         if sf is None:
#             sf = scale
#             set_scaling_factor(var, sf)

#         print(f"{var.name:<50s}{val:<20.3f}{scale:<20.3f}{sf:<20.3f}")



# def display_unfixed_vars(blk, report=True):
#     print("\n\n-------------------- UNFIXED VARIABLES --------------------\n\n")
#     print(f'{"BLOCK":<40s}{"UNFIXED VARIABLES":<30s}')
#     print(f"{blk.name:<40s}{number_unused_variables(blk)}")
#     for v in blk.component_data_objects(ctype=Block, active=True, descend_into=True):
#         print(f"{v.name:<40s}{number_unused_variables(v)}")
#         for v2 in unused_variables_set(v):
#             print(f"\t{v2.name:<40s}")
