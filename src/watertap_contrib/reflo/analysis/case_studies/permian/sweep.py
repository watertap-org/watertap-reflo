import os
import time
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
from watertap.core.solvers import get_solver

import watertap_contrib.reflo.analysis.case_studies.permian.permian_RPT1_MD as ST1
import watertap_contrib.reflo.analysis.case_studies.permian.permian_ZLD1_MD as ST2

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(parent_dir, "sweep_yamls")
save_dir = os.path.join(parent_dir, "sweep_results")

solver = get_solver()
solver.options["max_iter"] = 3000

def run_case_sweep(case, case_name=None, yaml_file=None):
    if yaml_file == None:
        map_yaml_file = os.path.join(sweep_yaml_dir, "Permian_ST_1.yaml")
    else:
        map_yaml_file = os.path.join(sweep_yaml_dir, yaml_file)

    save_name = case_name
    # print(case_name)
    for root, dirs, files in os.walk(os.path.join(save_dir, "output")):
        for file in files:
            file_id = file.split("_")
            case_id = case_name.split("_")
            if file_id[:3] == case_id[:3]:
                timestr = time.strftime("%Y%m%d-%H%M%S")
                print("Moving Prior Data to Archive")
                original_file = os.path.join(save_dir, "output", file)
                archive_file = os.path.join(
                    save_dir, "archive", file[:-3] + "_" + timestr + file[-3:]
                )
                os.rename(original_file, archive_file)

    lT = loopTool(
        map_yaml_file,
        build_function=case.build_sweep,
        optimize_function=case.solve,
        solver=solver,
        save_name=save_name,
        saving_dir=save_dir,
        execute_simulations=True,
        number_of_subprocesses=2,
    )

def run_diff_sweep(case, case_name=None, diff_yaml_file=None):
    diff_yaml_file = os.path.join(sweep_yaml_dir, diff_yaml_file)

    for root, dirs, files in os.walk(os.path.join(save_dir, "output")):
        for file in files:
            file_id = file.split("_")
            case_id = case_name.split("_")
            print(case_id)
            print(file_id)
            print(file_id[4])
            if (file_id[0] == 'diff') and (file_id[4] == case_id[1]):
                print(file_id, case_id, file)
                timestr = time.strftime("%Y%m%d-%H%M%S")
                print("Moving Prior Data to Archive")
                original_file = os.path.join(save_dir, "output", file)
                archive_file = os.path.join(
                    save_dir, "archive", file[:-3] + "_" + timestr + file[-3:]
                )
                os.rename(original_file, archive_file)

    start_time = time.time()
    print("Starting Diff Sweep")
    lT = loopTool(
        diff_yaml_file,
        build_function=case.build_sweep,
        optimize_function=case.solve,
        solver=solver,
        save_name="diff_sweep",
        saving_dir=save_dir,
        execute_simulations=True,
        number_of_subprocesses=4,
    )
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    print("Diff Sweep Complete")

def run_all_cases():
    cases = [
        {"case": ST1, "case_name": "Permian_ST_1", "yaml_file": "Permian_ST_1.yaml"},
        # {"case": RPT1, "case_name": "KBHDP_RPT_1", "yaml_file": "KBHDP_RPT_1.yaml"},
        # {"case": RPT2, "case_name": "KBHDP_RPT_2", "yaml_file": "KBHDP_RPT_2.yaml"},
    ]

    for case in cases:
        run_case_sweep(
            case["case"], case_name=case["case_name"], yaml_file=case["yaml_file"]
        )

def run_all_diff_cases():
    cases = [
        # {"case": ST1, "case_name": "Permian_ST1_diff", "yaml_file": "Permian_ST_1_diff.yaml"},
        {"case": ST2, "case_name": "Permian_ST2_diff", "yaml_file": "Permian_ST_2_diff.yaml"},
        # {"case": RPT3, "case_name": "KBHDP_RPT3_diff", "yaml_file": "KBHDP_RPT_3_diff.yaml"},
    ]

    for case in cases:
        run_diff_sweep(
            case["case"], case_name=case["case_name"], diff_yaml_file=case["yaml_file"]
        )

if __name__ == "__main__":
    run_all_diff_cases()
    # run_all_cases()
    # ST1.build_sweep()
    # ST1.solve()