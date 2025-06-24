import yaml
import pprint
import os
import time
import h5py
import pandas as pd
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
from idaes.core.solvers import get_solver

from parameter_sweep import ParameterSweep, DifferentialParameterSweep
from parameter_sweep import ParameterSweepReader
from parameter_sweep.writer import ParameterSweepWriter
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_SOA as SOA
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_1 as RPT1
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_2 as RPT2
import watertap_contrib.reflo.analysis.case_studies.KBHDP.KBHDP_RPT_3 as RPT3
# import watertap_contrib.reflo.analysis.case_studies.KBHDP.pretreatment_2 as PT2

from parameter_sweep.parallel.parallel_manager_factory import create_parallel_manager
filepath = os.path.abspath(__file__)
parent_dir = os.path.dirname(filepath)
sweep_yaml_dir = os.path.join(parent_dir, "sweep_yamls")
save_dir = os.path.join(parent_dir, "sweep_results")

solver = get_solver()
solver.options["max_iter"] = 3000

skips = [
    "diffus_phase",
    "diffus_param",
    "dens_mass_param",
    "dh_vap_w_param",
    "cp_phase_param",
    "pressure_sat_param_psatw",
    "enth_mass_param",
    "osm_coeff_param",
    "visc_d_param",
    "therm_cond_phase_param",
    "pressure_sat_param",
    "bpe_",
    "TIC",
    "TPEC",
    "blocks[",
    "yearly_heat_production",
    "yearly_electricity_production",
    "cp_param_NaCl_liq",
    "_translator",
    "permeate_side",
    "properties_interface",
    "material_flow_dx",
    "._flow_terms",
    "pressure_dx",
    "MCAS_properties",
    "cp_param_NaCl_solid",
    "cp_vap_param",
    "temp_sat_solvent",
    "m_mec.fs.properties",
    "cp_mass_phase",
]


def run_diff_sweep(case, case_name=None, diff_yaml_file=None):
    diff_yaml_file = os.path.join(sweep_yaml_dir, diff_yaml_file)

    for root, dirs, files in os.walk(os.path.join(save_dir, "output")):
        for file in files:
            file_id = file.split("_")
            case_id = case_name.split("_")
            # print(case_id)
            if (file_id[0] == "diff") and (file_id[3] == case_id[1]):
                # print(file_id, case_id, file)
                timestr = time.strftime("%Y%m%d-%H%M%S")
                print("Moving Prior Data to Archive")
                original_file = os.path.join(save_dir, "output", file)
                archive_file = os.path.join(
                    save_dir, "archive", file[:-3] + "_" + timestr + file[-3:]
                )
                os.rename(original_file, archive_file)

    lT = loopTool(
        diff_yaml_file,
        build_function=case.build_sweep,
        optimize_function=case.solve,
        solver=solver,
        save_name="diff_sweep",
        saving_dir=save_dir,
        execute_simulations=True,
        number_of_subprocesses=8,
    )


def run_case_sweep(case, case_name=None, yaml_file=None):
    if yaml_file == None:
        map_yaml_file = os.path.join(sweep_yaml_dir, "KBHDP.yaml")
    else:
        map_yaml_file = os.path.join(sweep_yaml_dir, yaml_file)

    yamld = yaml.safe_load(open(map_yaml_file, "r"))

    h5s = list()
    for k, v in yamld.items():
        _h5 = f"{save_dir}/output/{case_name}_analysisType_{case_name}_{list(v['sweep_param_loop'].keys())[0]}_sweep.h5"
        h5s.append(_h5)
                
    save_name = case_name
    # for k, v in 
    # for root, dirs, files in os.walk(os.path.join(save_dir, "output")):
    #     # print(files)
    #     for file in files:
    #         file_id = file.split("_")
    #         case_id = case_name.split("_")
    #         new_file = f"{save_dir}/output/{file}\n"

    #         if file_id[:3] == case_id[:3]:

    #             timestr = time.strftime("%Y%m%d-%H%M%S")
    #             print("Moving Prior Data to Archive")
    #             original_file = os.path.join(save_dir, "output", file)
    #             archive_file = os.path.join(
    #                 save_dir, "archive", file[:-3] + "_" + timestr + file[-3:]
    #             )
    #             os.rename(original_file, archive_file)

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

    for _h5 in h5s:
        _h5_csv = _h5.replace(".h5", ".csv")
        d = h5py.File(_h5, "r")
        # only key is name of sweep
        d = d[list(d.keys())[0]]
        # only key is analysis type
        d = d[list(d.keys())[0]]
        # default keys are ['outputs', 'solve_successful', 'sweep_params'], we want "outputs"
        d = d[list(d.keys())[0]]
        # only key here is "value"
        x = dict()
        for k in d.keys():
            if any(s in k for s in skips):
                continue
            # print(k, d[k]["value"][:])
            x[k] = d[k]["value"][:]
        df = pd.DataFrame(x)
        df.to_csv(_h5_csv, index=False)
    
    # assert False


    # print(f"\n\n{f'Parameter':<40}{'Value':<10}{'Result':<10}")
    # for item in lT.model_manager.ps.writer.parallel_manager.results.results[
    #     "sweep_params"
    # ].items():
    #     for idx, value in enumerate(item[1]["value"]):
    #         result = lT.ps.model_manager.ps.writer.parallel_manager.results.results[
    #             "solve_successful"
    #         ][idx]
    #         print(f"{f'{item[0]}':<40}{value:<10.2f}{result}")


def run_all_cases():
    cases = [
        # {"case": SOA, "case_name": "KBHDP_SOA_1", "yaml_file": "KBHDP_SOA_1.yaml"},
        # {"case": RPT1, "case_name": "KBHDP_RPT_1", "yaml_file": "KBHDP_RPT_1.yaml"},
        # {"case": RPT2, "case_name": "KBHDP_RPT_2", "yaml_file": "KBHDP_RPT_2.yaml"},
        {"case": RPT3, "case_name": "KBHDP_RPT_3", "yaml_file": "KBHDP_RPT_3.yaml"},
    ]

    for case in cases:
        run_case_sweep(
            case["case"], case_name=case["case_name"], yaml_file=case["yaml_file"]
        )


def run_all_diff_cases():
    cases = [
        {
            "case": RPT1,
            "case_name": "KBHDP_RPT1_diff",
            "yaml_file": "KBHDP_RPT_1_diff.yaml",
        },
        # {"case": RPT2, "case_name": "KBHDP_RPT2_diff", "yaml_file": "KBHDP_RPT_2_diff.yaml"},
        # {"case": RPT3, "case_name": "KBHDP_RPT3_diff", "yaml_file": "KBHDP_RPT_3_diff.yaml"},
    ]

    for case in cases:
        run_diff_sweep(
            case["case"], case_name=case["case_name"], diff_yaml_file=case["yaml_file"]
        )


if __name__ == "__main__":
    # SOA.build_sweep()
    run_all_cases()
    # run_all_diff_cases()
    # run_diff_sweep(RPT1, case_name="KBHDP_RPT_1_diff", diff_yaml_file='KBHDP_RPT_1_diff.yaml')
    # run_diff_sweep(RPT2, case_name="KBHDP_RPT2_diff", diff_yaml_file='KBHDP_RPT_2_diff.yaml')
    # run_case_sweep(SOA, case_name="KBHDP_SOA_1", yaml_file= 'KBHDP_SOA_1.yaml')
    # run_case_sweep(RPT2, case_name="KBHDP_RPT_2", yaml_file='KBHDP_RPT_2.yaml')
