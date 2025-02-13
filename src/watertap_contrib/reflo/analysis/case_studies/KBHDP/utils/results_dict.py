import itertools
import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    Expression,
    value,
    Var,
    Param,
    Constraint,
    assert_optimal_termination,
    check_optimal_termination,
    TransformationFactory,
    units as pyunits,
)
from idaes.core.solvers import get_solver

__all__ = [
    "build_results_dict",
    "results_dict_append",
    "plot_results_dict"
]


def build_results_dict(
    b,
    components=[Var, Expression, Param],
    skips=[],
):
    """
    Function to build a results dictionary with all components on block
    components: Pyomo objects to include in results_dict
    skips: component names (or parts of names) to skip

    returns: dictionary with desired model component results to store

    To be used with results_dict_append function
    """

    results_dict = {}

    for c in components:
        for v in b.component_objects(c):
            if any(s in v.name for s in skips):
                continue
            elif v.is_indexed():
                for i, vv in v.items():
                    results_dict[vv.name] = []
            else:
                results_dict[v.name] = []

    return results_dict


def results_dict_append(
    b,
    results_dict,
    components=[Var, Expression, Param],
    tmp_results_dict=None,
):
    """
    Add to results_dict after model solve

    results_dict: dictionary meant to store model results
    components: model components included in dictionary
    tmp_results_dict: temporary results dictionary for nested loop meant to add to larger results dictionary

    returns: the same results_dict, but with another "row" added
    ___________________________________________________________
    EXMPLE FOR 1 PARAM SWEEP:

    m = build_model_func()

    results_dict = build_results_dict(m)

    for x in x_params:
        m = build_model_func()

        m.fs.unit_model.x_var.fix(x)

        results = solver.solve(m)

        results_dict = results_dict_append(m, results_dict)
    ___________________________________________________________
    EXAMPLE FOR 2 PARAM SWEEP:

    m = build_model_func()

    results_dict = build_results_dict(m)

    for x in x_params:
        tmp_results_dict = build_results_dict(m)
        for y in y_params:

            m = build_model_func()

            m.fs.unit_model.x_var.fix(x)
            m.fs.unit_model.y_var.fix(y)

            results = solver.solve(m)

            tmp_results_dict = results_dict_append(m, tmp_results_dict)

        results_dict = results_dict_append(b, results_dict, tmp_results_dict=tmp_results_dict)

    """

    appended = list()

    if tmp_results_dict is None:
        for c in components:
            for v in b.component_objects(c):
                if v.is_indexed():
                    # print(v.name)
                    # idx = [*v._index_set]
                    for i, vv in v.items():
                        if vv.name in results_dict.keys():
                            if vv.name in appended:
                                continue
                            try:
                                results_dict[vv.name].append(value(vv))
                                appended.append(vv.name)
                            except:
                                pass
                else:
                    if v.name in results_dict.keys():
                        if v.name in appended:
                            continue
                        try:
                            results_dict[v.name].append(value(v))
                            appended.append(v.name)
                        except:
                            pass

    else:
        for k, v in tmp_results_dict.items():
            if k in appended or k not in results_dict.keys():
                continue
            results_dict[k].append(v)
            appended.append(k)

    return results_dict


def plot_results_dict(
    results_dict,
    xcol="",
    xlabel="",
    ycol=None,
    ylabel="",
    cols=[],
    cmap="hsv",
    mode="2d",
    num_pts=5,
    plot_type=None,
    plot_kwargs=dict(),
    leg_col=None,
    fig_title=None,
    plot_levels=False,
    cbar_label=None,
    add_levels_dict=None,
    sharex=False,
):
    """
    Plotting results_dict
    xcol: The column used for the x-axis. Required for both 2D and 3D plots.
    xlabel: str label for xcol
    ycol: The column used for the y-axis. Triggers 3D plot. Must be independent variable.
    ylabel: str label for ycol
    cols: Keys in results_dict to be plotted; for both 2D and 3D plots.
    cmap: colormap to generate colors for each plot
    mode: Category of plot (2D or 3D); probably redundant.
    num_pts: Number of contours for 3D plots.
    plot_type: For 2D plots, the type of plot to make (e.g., scatter, plot). Required to be str of method of Axes object (e.g. ax.scatter - "scatter" is the plot_type)
    plot_kwargs: kwargs for plot
    leg_col: Column to use as legend.
    """

    if plot_type is None:
        plot_type = "plot"
    if leg_col is not None:
        leg_labes = list()
        for x in results_dict[leg_col]:
            if isinstance(x, list):
                assert all(i == x[0] for i in x)
                leg_labes.append(x[0])
            else:
                leg_labes.append(x)

    if ycol is not None:
        mode = "3d"
        min_level = 100000000000
        max_level = -10000000000
        levels_dict = dict()
        for col in cols:
            z = results_dict[col]
            for zs in z:
                if min(zs) < min_level:
                    min_level = min(zs)
                if max(zs) > max_level:
                    max_level = max(zs)
            levels = np.linspace(min_level, max_level, num_pts)
            if add_levels_dict:
                levels = np.append(levels, add_levels_dict[col])
            levels_dict[col] = sorted(levels)

    cmap = plt.cm.get_cmap(cmap)
    rgbas = itertools.cycle([cmap(x) for x in np.linspace(0, 0.9, len(cols))])
    if mode == "3d":
        x = results_dict[xcol]
        y = results_dict[ycol]
        fig, axs = plt.subplots(
            nrows=len(cols), ncols=1, figsize=(8, 5 * len(cols)), sharex=False
        )
        # for ax, col in zip(axs, cols):
        for i, col in enumerate(cols):
            if len(cols) == 1:
                ax = axs
            else:
                ax = axs[i]
            z = results_dict[col]
            cs1 = ax.contourf(x, y, z, 100, cmap=cmap)
            if plot_levels:
                cs2 = ax.contour(
                    x,
                    y,
                    z,
                    levels_dict[col],
                    colors="k",
                    linestyles="dashed",
                )
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            title = col.split(".")[-1].replace("_", " ").upper()
            # ax.set_title(title)
            cbar = fig.colorbar(cs1, ax=ax)
            if cbar_label is None:
                cbar.set_label(title)
            else:
                cbar.set_label(cbar_label)

    if mode == "2d":

        fig, axs = plt.subplots(
            nrows=len(cols), ncols=1, figsize=(8, 5 * len(cols)), sharex=sharex
        )
        if len(cols) == 1:
            axs = [axs]
        x = results_dict[xcol]

        for j, (ax, col) in enumerate(zip(axs, cols)):
            plot_func = getattr(ax, plot_type)
            print(col)
            y = results_dict[col]
            if isinstance(y[0], list):
                rgbas2 = itertools.cycle([cmap(c) for c in np.linspace(0, 0.9, len(y))])
                for i, (xx, yy) in enumerate(zip(x, y)):
                    plot_kwargs["color"] = next(rgbas2)
                    if leg_col is not None:
                        plot_kwargs["label"] = (
                            leg_col.split(".")[-1].replace("_", " ").title()
                            + ": "
                            + str(leg_labes[i])
                        )
                    plot_func(xx, yy, **plot_kwargs)

            else:
                # ax.plot(x, y, color=next(rgbas))
                plot_kwargs["color"] = next(rgbas)
                if leg_col is not None:
                    plot_kwargs["label"] = (
                        leg_col.split(".")[-1].replace("_", " ").title()
                        + ": "
                        + str(leg_labes[j])
                    )
                plot_func(x, y, **plot_kwargs)
            ax.set_title(col.split(".")[-1].replace("_", " ").title())
            ax.set_xlabel(xlabel)
            if "label" in plot_kwargs:
                ax.legend()
    if fig_title is not None:
        fig.suptitle(fig_title, fontsize="xx-large", fontweight="bold", y=1.05, x=0.4)
    return fig, axs


if __name__ == "__main__":

    import numpy as np
    import pandas as pd
    from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
    from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
    from watertap_contrib.reflo.analysis.case_studies.KBHDP import *
    from watertap_contrib.reflo.costing import TreatmentCosting
    from idaes.core.util.model_statistics import *

    solver = get_solver()

    def build_model(flow_mgd=1):
        """
        Function to build the model
        """

        rho = 1000 * pyunits.kg / pyunits.m**3
        flow_mgd = flow_mgd * pyunits.Mgallons / pyunits.day

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = TreatmentCosting()
        inlet_conc = {"Na_+": 120, "Cl_-": 100}
        m.fs.flow_mgd = Var(initialize=flow_mgd)
        m.fs.flow_mgd.fix()

        m.fs.properties = MCASParameterBlock(
            solute_list=inlet_conc.keys(), material_flow_basis=MaterialFlowBasis.mass
        )

        m.fs.DWI = FlowsheetBlock(dynamic=False)
        build_DWI(m, m.fs.DWI, m.fs.properties)

        m.fs.DWI.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={
                "cost_method": "blm"
            },  # could be "as_capex" or "blm"
        )

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.DWI.feed.properties[0].flow_vol_phase["Liq"])

        flow_mass_phase_water = pyunits.convert(
            flow_mgd * rho, to_units=pyunits.kg / pyunits.s
        )

        prop = m.fs.DWI.unit.properties[0]
        prop.temperature.fix()
        prop.pressure.fix()

        for solute, conc in inlet_conc.items():
            mass_flow_solute = pyunits.convert(
                flow_mgd * conc * pyunits.kg / pyunits.m**3,
                to_units=pyunits.kg / pyunits.s,
            )
            prop.flow_mass_phase_comp["Liq", solute].fix(mass_flow_solute)
            prop.set_default_scaling(
                "flow_mass_phase_comp",
                value(1 / mass_flow_solute),
                index=("Liq", solute),
            )
        prop.set_default_scaling(
            "flow_mass_phase_comp",
            value(1 / flow_mass_phase_water),
            index=("Liq", "H2O"),
        )
        prop.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_phase_water)
        TransformationFactory("network.expand_arcs").apply_to(m)

        return m

    m = build_model()

    results_dict = build_results_dict(m, skips=["diffus_phase_comp"])

    flows_mgd = [1, 2, 3, 4, 5]
    dwi_lcows = np.linspace(1000, 10000, 11)

    for x in dwi_lcows:
        m = build_model()
        m.fs.DWI.unit.injection_well_depth.set_value(x)
        results = solver.solve(m)
        assert_optimal_termination(results)
        results_dict = results_dict_append(m, results_dict)
        print(x, m.fs.costing.LCOW())

    df = pd.DataFrame.from_dict(results_dict)
    plot_results_dict(
        results_dict,
        xcol="fs.DWI.unit.injection_well_depth",
        cols=[
            "fs.costing.LCOW",
            "fs.costing.total_capital_cost",
            "fs.costing.total_operating_cost",
        ],
    )
    # plt.show()

    results_dict = build_results_dict(m, skips=["diffus_phase_comp"])
    for x in dwi_lcows:
        tmp_results_dict = build_results_dict(m, skips=["diffus_phase_comp"])
        for y in flows_mgd:
            m = build_model(flow_mgd=y)
            m.fs.DWI.unit.injection_well_depth.set_value(x)
            results = solver.solve(m)
            assert_optimal_termination(results)
            tmp_results_dict = results_dict_append(m, tmp_results_dict)
        results_dict = results_dict_append(
            m, results_dict, tmp_results_dict=tmp_results_dict
        )

    plot_results_dict(
        results_dict,
        xcol="fs.flow_mgd",
        ycol="fs.DWI.unit.injection_well_depth",
        cols=[
            "fs.costing.LCOW",
            "fs.costing.total_capital_cost",
            "fs.costing.total_operating_cost",
        ],
        cmap="viridis",
        plot_levels=False,
    )
    plt.show()
