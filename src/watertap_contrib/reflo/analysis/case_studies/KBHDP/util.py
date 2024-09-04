from pyomo.environ import Var, value, units as pyunits
import pandas as pd
import idaes.core.util.scaling as iscale


def m3s_to_gpm(Q):
    return [q * 1.585032e4 for q in Q]


def m3s_to_mgd(Q):
    return [q * 22.827 for q in Q]


def mgd_to_m3s(Q):
    return [q * 0.0438 for q in Q]


def generate_costing_report(m, export=False, filepath=None):
    costing_data = {"Attribute": "Value"}
    costing_objs = []
    costing_vals = []
    costing_units = []
    # for v in m.fs.energy.pv.costing.component_data_objects(Var, descend_into=True):
    #     costing_data[str(v)] = value(v)
    # for v in m.fs.energy.costing.component_data_objects(Var, descend_into=True):
    #     costing_data[str(v)] = value(v)
    for v in m.fs.costing.component_data_objects(Var, descend_into=True):
        costing_data[str(v)] = value(v)
        costing_objs.append(str(v))
        costing_vals.append(value(v))
        costing_units.append(pyunits.get_units(v))
    # for v in m.fs.sys_costing.component_data_objects(Var, descend_into=True):
    #     costing_data[str(v)] = value(v)

    # if filepath != None:
    #     with open(filepath, "w") as f:
    #         for key in costing_data.keys():
    #             f.write("%s,%s\n" % (key, costing_data[key]))
    costing_df = pd.DataFrame(columns=["Attribute", "Value", "Units"])
    costing_df["Attribute"] = costing_objs
    costing_df["Value"] = costing_vals
    costing_df["Units"] = costing_units
    if export:
        costing_df.to_csv(filepath, index=False)

    return costing_data, costing_df


def extract_values(block):
    costing_data = dict.fromkeys(["key_0", "key_1", "key_2", "key_3"], None)
    for v in block.component_data_objects(Var, descend_into=False):
        keys = str(v).split(".")
        for idx, key in enumerate(keys[:-1]):
            costing_data["key_" + str(idx)] = key
        costing_data[keys[-1]] = value(v)
    return costing_data


def generate_detailed_costing_report(m, level, filepath=None):
    dict_track = []
    if level == "Process":
        cost_blocks = [m.fs.energy.costing, m.fs.treatment.costing, m.fs.sys_costing]
    else:
        cost_blocks = [m.fs.energy.pv.costing, m.fs.energy.costing.pv_surrogate]
    for block in cost_blocks:
        dict_track.append(extract_values(block))
    return dict_track


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
        for val, var, con in extreme_entries:
            print(val, var.name, con.name)
        print("--------------------------")
        print("Extreme Jacobian columns:")
        extreme_cols = iscale.extreme_jacobian_columns(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, var in extreme_cols:
            print(val, var.name)
        print("------------------------")
        print("Extreme Jacobian rows:")
        extreme_rows = iscale.extreme_jacobian_rows(
            m, jac=jac_scaled, nlp=nlp, small=1e-3
        )
        for val, con in extreme_rows:
            print(val, con.name)
    return cond_number
