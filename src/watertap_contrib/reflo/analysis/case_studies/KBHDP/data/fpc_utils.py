import os
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    check_optimal_termination,
    SolverFactory,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.analysis.case_studies.KBHDP.utils import *
from watertap_contrib.reflo.analysis.case_studies.KBHDP.data import *
from watertap_contrib.reflo.solar_models.surrogate.flat_plate.flat_plate_surrogate import (
    FlatPlateSurrogate,
)

__all__ = [
    "plot_contour",
    "run_pysam_fpc_scaled_draw_sweep",
    "plot_pysam_fpc_scaled_draw_sweep",
    "plot_pysam_fpc_time_series",
    "build_and_run_fpc_surrogate",
    "run_fpc_sweep",
    "run_reflo_pysam_fpc_comparison",
    "run_reflo_pysam_fpc_comparison_sweep",
]

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
weather_file = os.path.join(__location__, "el_paso_texas-KBHDP-weather.csv")
param_file = os.path.join(__location__, "fpc/solar_water_heating-kbhdp.json")

temperatures = {
    "T_cold": 20,
    "T_hot": 70,  # this will be overwritten by temperature_hot value
    "T_amb": 18,
}


def plot_contour(
    df,
    x="",
    y="",
    z="fs.costing.LCOW",
    min_x=None,
    min_y=None,
    max_x=None,
    max_y=None,
    approach="interpolate",
    cmap="turbo",
    interp_method="cubic",
    grid_len=100,
    levels=20,
    set_dict=dict(),
    cb_title="",
    contour_label_fmt="  %#.2e  ",
    add_contour_labels=False,
):
    """
    Generic function to create contour plot from
    three columns in DataFrame
    """

    fig, ax = plt.subplots()

    if min_x is not None:
        df = df[df[x] >= min_x].copy()
    if min_y is not None:
        df = df[df[y] >= min_y].copy()

    if max_x is not None:
        df = df[df[x] <= min_x].copy()
    if max_y is not None:
        df = df[df[y] <= min_y].copy()

    if approach == "pivot":
        pivot = df.pivot(index=y, columns=x, values=z)

        x = pivot.columns.values
        y = pivot.index.values
        z = pivot.values

        contourf = ax.contourf(x, y, z, levels=levels, cmap=cmap)
        cb = plt.colorbar(contourf, label=cb_title)

    if approach == "interpolate":

        x = df[x]
        y = df[y]
        print(z)
        z = df[z]

        xi = np.linspace(min(x), max(x), grid_len)
        yi = np.linspace(min(y), max(y), grid_len)
        xi, yi = np.meshgrid(xi, yi)

        zi = griddata((x, y), z, (xi, yi), method=interp_method)

        contourf = ax.contourf(xi, yi, zi, levels=levels, cmap=cmap)
        cb = plt.colorbar(contourf, label=cb_title)

    ax.set(**set_dict)

    if add_contour_labels:
        contour = ax.contour(
            contourf,
            levels,
            colors="k",
            linestyles="dashed",
        )
        ax.clabel(contour, colors="black", fmt=contour_label_fmt)


def run_pysam_fpc_scaled_draw_sweep(
    heat_load_mwt=None, hours_storage=None, temperature_hot=None, num_pts=25
):
    """
    Run PySAM FPC model across range of scaled_draw values.
    This was mostly used for testing the implementation.
    """

    config_data = read_module_datafile(param_file)

    if "solar_resource_file" in config_data:
        del config_data["solar_resource_file"]
    tech_model = setup_pysam_fpc_model(
        temperatures=temperatures,
        weather_file=weather_file,
        config_data=config_data,
    )
    # run with scaled_draw=None to find the current calculated value for scaled_draw as the upperbound for the sweep
    _, tech_model = run_pysam_fpc_model(
        tech_model,
        heat_load_mwt=heat_load_mwt,
        hours_storage=hours_storage,
        temperature_hot=temperature_hot,
        return_tech_model=True,
    )
    scaled_draws = np.linspace(1, tech_model.value("scaled_draw")[0], num_pts).tolist()

    q_annual = list()
    q_loss = list()
    q_aux = list()
    q_auxonly = list()
    temp_deliv = list()

    for sd in scaled_draws:
        config_data = read_module_datafile(param_file)

        if "solar_resource_file" in config_data:
            del config_data["solar_resource_file"]
        tech_model = setup_pysam_fpc_model(
            temperatures=temperatures,
            weather_file=weather_file,
            config_data=config_data,
        )
        _, tech_model = run_pysam_fpc_model(
            tech_model,
            heat_load_mwt=heat_load_mwt,
            hours_storage=hours_storage,
            temperature_hot=temperature_hot,
            scaled_draw=sd,
        )
        q_annual.append(tech_model.Outputs.annual_Q_deliv)
        q_auxonly.append(sum(tech_model.Outputs.Q_auxonly))
        q_aux.append(sum(tech_model.Outputs.Q_aux))
        q_loss.append(sum(tech_model.Outputs.Q_loss))
        temp_deliv.append(np.mean(tech_model.Outputs.T_deliv))

    df = pd.DataFrame.from_dict(
        {
            "q_annual": q_annual,
            "q_auxonly": q_auxonly,
            "q_aux": q_aux,
            "q_loss": q_loss,
            "temp_delivered_avg": temp_deliv,
        }
    )

    df["heat_load"] = heat_load_mwt
    df["hours_storage"] = hours_storage
    df["temperature_hot"] = temperature_hot
    df["scaled_draw"] = scaled_draws

    return df


def plot_pysam_fpc_scaled_draw_sweep(df, xcol="scaled_draw"):
    """
    Plot results of scaled draw sweep
    """

    heat_load = df.heat_load.unique()[0]
    hours_storage = df.hours_storage.unique()[0]
    temp_hot = df.temperature_hot.unique()[0]

    # plot heat cols
    fig, ax = plt.subplots()

    ax.plot(df[xcol], df.q_annual, label="Heat Delivered", marker=".")
    ax.plot(df[xcol], df.q_aux, label="Auxiliary Heat Required", marker=".")
    ax.plot(df[xcol], df.q_loss, label="Heat Loss", marker=".")
    ax.set(
        xlabel="Scaled Draw [kg/hr]",
        ylabel="Heat [W]",
        title=f"{heat_load} MW, {hours_storage} hr Storage, {temp_hot}C",
    )
    ax.legend()

    fig, ax = plt.subplots()
    ax.plot(
        df[xcol],
        df.temp_delivered_avg,
        marker=".",
        label="Temperature Delivered",
        color="r",
    )
    ax.plot(
        df[xcol],
        [temp_hot for _ in df[xcol]],
        ls=":",
        label="Temperature Desired",
        color="k",
        alpha=0.5,
    )
    ax.set(
        xlabel="Scaled Draw [kg/hr]",
        ylabel="Temperature [C]",
        title=f"{heat_load} MW, {hours_storage} hr Storage, {temp_hot}C",
    )
    ax.legend()


def plot_pysam_fpc_time_series(
    tech_model, surr_heat=None, treat_heat_req=None, starting_week=15, num_weeks=1
):
    """
    Plot time-series results from FPC PySAM for:
        1. Temperature
        2. Draw flow rate
        3. Heat
    """

    start = 168 * starting_week
    end = start + (168 * num_weeks)

    figs = list()
    axs = list()

    fig, ax = plt.subplots()

    ax.plot(
        tech_model.Outputs.T_deliv[start:end], label="Delivered Temperature", color="r"
    )
    ax.legend()
    ax.set(xlabel="Time [hr]", ylabel="Temperature [C]")

    figs.append(fig)
    axs.append(ax)

    fig, ax = plt.subplots()
    ax.plot(tech_model.Outputs.draw[start:end], label="Draw Flow")
    ax.legend()
    ax.set(xlabel="Time [hr]", ylabel="Draw Flow [kg/hr]")
    figs.append(fig)
    axs.append(ax)

    fig, ax = plt.subplots()
    ax.plot(tech_model.Outputs.Q_loss[start:end], label="Heat Loss")
    ax.plot(tech_model.Outputs.Q_deliv[start:end], label="Heat Delivered")
    ax.plot(tech_model.Outputs.Q_aux[start:end], label="Auxiliary Heat Required")
    
    if surr_heat is not None:
        ax.plot(
            [surr_heat for _ in tech_model.Outputs.Q_aux[start:end]],
            ls=":",
            color="k",
            label="Heat Delivered by Surrogate",
        )
    if treat_heat_req is not None:
        ax.plot(
            [treat_heat_req for _ in tech_model.Outputs.Q_aux[start:end]],
            ls=":",
            color="r",
            label="Heat Required for Treatment",
        )
    ax.legend()
    ax.set(xlabel="Time [hr]", ylabel="Heat [kW]")
    figs.append(fig)
    axs.append(ax)

    return figs, axs


def build_and_run_fpc_surrogate(
    heat_load=1,
    hours_storage=1,
    temperature_hot=80,
    input_variables=dict(),
    output_variables=dict(),
    dataset_filename=None,
    surrogate_model_file=None,
    surrogate_filename_save="test",
    heat_cost=0,
    solver=None,
    **kwargs,
):
    """
    Build and run an FPC flowsheet.
    """
    global m

    if solver is None:
        solver = SolverFactory("ipopt")

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = EnergyCosting()
    m.fs.FPC = FlatPlateSurrogate(
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
        surrogate_model_file=surrogate_model_file,
        surrogate_filename_save=surrogate_filename_save,
        scale_training_data=True,
    )
    m.fs.FPC.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.heat_cost.fix(heat_cost)
    print(f"dof = {degrees_of_freedom(m)}")

    m.fs.FPC.heat_load.fix(heat_load)
    m.fs.FPC.hours_storage.fix(hours_storage)
    m.fs.FPC.temperature_hot.fix(temperature_hot)

    print(f"dof = {degrees_of_freedom(m)}")

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOH()

    m.fs.FPC.initialize()
    m.fs.costing.initialize()

    for y, c in m.fs.costing.yearly_heat_production_constraint.items():
        cvc(m.fs.costing.yearly_heat_production[y], c)
    cvc(
        m.fs.costing.lifetime_heat_production,
        m.fs.costing.lifetime_heat_production_constraint,
    )
    results = solver.solve(m)
    if not check_optimal_termination(results):
        results = solver.solve(m)

    # assert_optimal_termination(results)

    return m


def run_fpc_sweep(
    m=None,
    heat_loads=list(),
    hours_storages=list(),
    temperature_hots=list(),
    surrogate_model_file=None,
    dataset_filename=None,
    save_sweep=None,
    **kwargs,
):
    # need an initial model to build_results_dict with
    if m is None:
        m = build_and_run_fpc_surrogate(
            heat_load=5,
            hours_storage=12,
            temperature_hot=80,
            surrogate_model_file=surrogate_model_file,
            dataset_filename=dataset_filename,
            **kwargs,
        )

    rdf = pd.DataFrame()

    for hl in heat_loads:
        rd = build_results_dict(m)
        for hs in hours_storages:
            for th in temperature_hots:
                print(f"Running FPC Surrogate:")
                print(f"\tHeat load: {hl:.2f} MW")
                print(f"\tHours storage: {hs:.2f} hr")
                print(f"\tHot temperature: {th:.2f} C")
                try:
                    m = build_and_run_fpc_surrogate(
                        heat_load=hl,
                        hours_storage=hs,
                        temperature_hot=th,
                        surrogate_model_file=surrogate_model_file,
                        dataset_filename=dataset_filename,
                        heat_cost=0,
                    )
                except:
                    continue
                rd = results_dict_append(m, rd)
        rdf = pd.concat([rdf, pd.DataFrame.from_dict(rd)])

    if save_sweep is not None:
        rdf.to_csv(save_sweep, index=False)

    return rdf


def run_reflo_pysam_fpc_comparison(
    heat_load=None,
    hours_storage=None,
    temperature_hot=None,
    surrogate_model_file=None,
    dataset_filename=None,
    **fpc_build_kwargs,
):
    config_data = read_module_datafile(param_file)

    if "solar_resource_file" in config_data:
        del config_data["solar_resource_file"]

    tech_model = setup_pysam_fpc_model(
        temperatures=temperatures,
        weather_file=weather_file,
        config_data=config_data,
    )
    result, tech_model = run_pysam_fpc_model(
        tech_model,
        heat_load_mwt=heat_load,
        hours_storage=hours_storage,
        temperature_hot=temperature_hot,
        return_tech_model=True,
    )

    m = build_and_run_fpc_surrogate(
        heat_load=heat_load,
        hours_storage=hours_storage,
        temperature_hot=temperature_hot,
        surrogate_model_file=surrogate_model_file,
        dataset_filename=dataset_filename,
        heat_cost=0,
        **fpc_build_kwargs,
    )
    heat_annual_diff_rel = (m.fs.FPC.heat_annual() - result["heat_annual"]) / result[
        "heat_annual"
    ]
    elect_annual_diff_rel = (
        m.fs.FPC.electricity_annual() - result["electricity_annual"]
    ) / result["electricity_annual"]

    print(
        f"\n{'Parameter':<30} {'PySAM Results':<30} {'REFLO Results':<30} {'Relative Difference':<30}"
    )
    print(
        f"{'Annual Heat Produced':<30} {result['heat_annual']:<30.2f} {m.fs.FPC.heat_annual():<30.2f} {heat_annual_diff_rel*100:.2f}%"
    )
    print(
        f"{'Annual Electricity Reqd.':<30} {result['electricity_annual']:<30.2f} {m.fs.FPC.electricity_annual():<30.2f} {elect_annual_diff_rel*100:.2f}%"
    )
    # print(f"{'Annual Heat Produced:':<40} {tech_model.value('annual_Q_deliv'):<40.2f} {} MW")


def run_reflo_pysam_fpc_comparison_sweep(pysam_df, fpc_build_kwargs=dict()):
    pysam_df.reset_index(inplace=True)

    heat_annual_error = list()
    heat_annual_model = list()
    elec_annual_error = list()
    elec_annual_model = list()

    for i, row in pysam_df.iterrows():
        print(f"Running FPC Surrogate:")
        print(f"\tHeat load: {row.heat_load:.2f} MW")
        print(f"\tHours storage: {row.hours_storage:.2f} hr")
        print(f"\tHot temperature: {row.temperature_hot:.2f} C")
        try:
            m = build_and_run_fpc_surrogate(
                heat_load=row.heat_load,
                hours_storage=row.hours_storage,
                temperature_hot=row.temperature_hot,
                **fpc_build_kwargs,
            )
            heat_annual_error.append(value(m.fs.FPC.heat_annual) - row.heat_annual)
            elec_annual_error.append(
                value(m.fs.FPC.electricity_annual) - row.electricity_annual
            )
            heat_annual_model.append(value(m.fs.FPC.heat_annual))
            elec_annual_model.append(value(m.fs.FPC.electricity_annual))
        except:
            heat_annual_error.append(None)
            elec_annual_error.append(None)
            heat_annual_model.append(None)
            elec_annual_model.append(None)

    pysam_df["heat_annual_model"] = heat_annual_model
    pysam_df["heat_annual_error"] = heat_annual_error
    pysam_df["heat_annual_error_rel"] = (
        pysam_df.heat_annual_model - pysam_df.heat_annual
    ) / pysam_df.heat_annual
    pysam_df["electricity_annual_model"] = elec_annual_model
    pysam_df["electricity_annual_error"] = elec_annual_error
    pysam_df["electricity_annual_error_rel"] = (
        pysam_df.electricity_annual_model - pysam_df.electricity_annual
    ) / pysam_df.electricity_annual

    return pysam_df


# if __name__ == "__main__":
#     surrogate_model_file = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/data/fpc/FPC_KBHDP_el_paso_HIGH_heat_load_1-50_hours_storage_1-24_temperature_hot_50-98-rerun.json"
#     dataset_filename = "/Users/ksitterl/Documents/Python/watertap-reflo/watertap-reflo/src/watertap_contrib/reflo/analysis/case_studies/KBHDP/data/fpc/FPC_KBHDP_el_paso_HIGH_heat_load_1-50_hours_storage_1-24_temperature_hot_50-98-rerun.pkl"
#     input_bounds = dict(
#         heat_load=[1, 50], hours_storage=[1, 24], temperature_hot=[50, 98]
#     )
#     input_units = dict(heat_load="MW", hours_storage="hour", temperature_hot="degK")
#     input_variables = {
#         "labels": ["heat_load", "hours_storage", "temperature_hot"],
#         "bounds": input_bounds,
#         "units": input_units,
#     }

#     output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
#     output_variables = {
#         "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
#         "units": output_units,
#     }
#     run_reflo_pysam_fpc_comparison(
#         heat_load=10,
#         hours_storage=12,
#         temperature_hot=80,
#         surrogate_model_file=surrogate_model_file,
#         dataset_filename=dataset_filename,
#         input_variables=input_variables,
#         output_variables=output_variables,
#     )

#     df = pd.read_pickle(dataset_filename)
#     df = df[df.temperature_hot == 80].copy()
#     df = df[df.heat_load <= 25].copy()

#     fpc_build_kwargs = dict(
#         input_variables=input_variables,
#         output_variables=output_variables,
#         dataset_filename=dataset_filename,
#         surrogate_model_file=surrogate_model_file,
#     )

#     df = run_reflo_pysam_fpc_comparison_sweep(df, fpc_build_kwargs=fpc_build_kwargs)

#     plot_contour(
#         df,
#         x="heat_load",
#         y="hours_storage",
#         z="heat_annual_error_rel",
#         approach="interpolate",
#         cmap="turbo",
#         interp_method="cubic",
#         grid_len=100,
#         levels=10,
#         set_dict=dict(),
#         cb_title="",
#         contour_label_fmt="  %#.2f  ",
#         add_contour_labels=True,
#     )
#     plt.show()
