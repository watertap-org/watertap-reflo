import os
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir

from watertap_contrib.reflo.analysis.case_studies.KBHDP.sweep_setup import build, solve_system

filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(parent_dir, "sweep_yamls")
save_dir = os.path.join(parent_dir, "sweep_results")



def run_water_recovery_sweep():
    map_yaml_file = os.path.join(sweep_yaml_dir, "KBHDP_SOA_1.yaml")
    save_name = "KBHDP_SOA_1"
    output_file_path = os.path.join(
        save_dir + "/output", save_name + "_analysisType_KBHDP_SOA_1_sweep.h5"
    )

    print(map_yaml_file)
    if os.path.exists(output_file_path):
        print(f"Found {output_file_path}")
        print(f"Removing {output_file_path}")
        os.remove(output_file_path)
    else:
        print(f"Storing results in {output_file_path}")

    lT = loopTool(
        map_yaml_file,
        build_function=build,
        optimize_function=solve_system,
        save_name=save_name,
        saving_dir=save_dir,
        execute_simulations=True,
        number_of_subprocesses=8,
        # parallel_back_end="RayIo",
    )




if __name__ == "__main__":
    run_water_recovery_sweep()
