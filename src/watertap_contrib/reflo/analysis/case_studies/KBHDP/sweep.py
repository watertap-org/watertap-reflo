import os
import time
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
from idaes.core.solvers import get_solver

import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_SOA as SOA
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_1 as RPT1
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_2 as RPT2

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(parent_dir, "sweep_yamls")
save_dir = os.path.join(parent_dir, "sweep_results")

solver = get_solver()
solver.options["max_iter"] = 3000


def run_case_sweep(case, case_name=None, yaml_file=None):
    if yaml_file == None:
        map_yaml_file = os.path.join(sweep_yaml_dir, "KBHDP.yaml")
    else:
        map_yaml_file = os.path.join(sweep_yaml_dir, yaml_file)

    save_name = case_name

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
        number_of_subprocesses=1,
    )


def run_all_cases():
    cases = [
        {"case": SOA, "case_name": "KBHDP_SOA_1", "yaml_file": "KBHDP_SOA_1.yaml"},
        # {"case": SOA, "case_name": "KBHDP_SOA_1", "yaml_file": "KBHDP_SOA_1_diff.yaml"},
        # {"case": RPT1, "case_name": "KBHDP_RPT_1", "yaml_file": "KBHDP_RPT_1.yaml"},
        # {"case": RPT2, "case_name": "KBHDP_RPT_2", "yaml_file": "KBHDP_RPT_2.yaml"},
    ]

    for case in cases:
        run_case_sweep(
            case["case"], case_name=case["case_name"], yaml_file=case["yaml_file"]
        )


if __name__ == "__main__":
    run_all_cases()

    # run_case_sweep(SOA, case_name="KBHDP_SOA_1", yaml_file= 'KBHDP_SOA_1.yaml')
    # run_case_sweep(RPT2, case_name="KBHDP_RPT_2", yaml_file='KBHDP_RPT_2.yaml')
