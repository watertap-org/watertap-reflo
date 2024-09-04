from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    value,
    Var,
    units as pyunits,
)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import idaes.core.util.scaling as iscale


def plot_flux_profile(m, blk):
    P0 = value(pyunits.convert(blk.feed_side.pressure[0, 0], to_units=pyunits.psi))
    P2 = -2 * value(
        pyunits.convert(
            blk.feed_side.deltaP[0, 0], to_units=pyunits.psi * pyunits.m**-1
        )
    )
    P3 = 2 * value(
        pyunits.convert(
            blk.feed_side.properties_interface[0, 1].pressure_osm_phase["Liq"],
            to_units=pyunits.psi,
        )
    )

    module_type = ["Flat-plate", "Sprial-wound"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax2 = ax.twinx()
    ax3 = ax.twinx()
    ax4 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.2))
    ax4.spines.right.set_position(("axes", 1.35))
    line_style = ["-", "--", "-.", ":"]
    color1 = "tab:blue"
    color2 = "tab:orange"
    color3 = "tab:green"
    color4 = "tab:red"
    # for idx, blk in enumerate([m.fs.unit1, m.fs.unit2]):
    # blk.feed_side.deltaP.display()
    x = list(blk.length_domain)
    y1 = [
        value(
            pyunits.convert(
                blk.permeate_side[0.0, i].flow_mass_phase_comp["Liq", "H2O"],
                to_units=pyunits.kg * pyunits.hr**-1,
            )
        )
        for i in x[1:]
    ]
    y2 = [
        value(pyunits.convert(blk.feed_side.pressure[0, i], to_units=pyunits.psi))
        for i in x[1:]
    ]
    y3 = [
        -1
        * value(
            pyunits.convert(
                blk.feed_side.deltaP[0, i], to_units=pyunits.psi * pyunits.m**-1
            )
        )
        for i in x[1:]
    ]
    y4 = [
        value(
            pyunits.convert(
                blk.feed_side.properties_interface[0, i].pressure_osm_phase["Liq"],
                to_units=pyunits.psi,
            )
        )
        for i in x[1:]
    ]
    ax.plot(
        x[1:],
        y1,
        lw=2,
        ls=line_style[0],
        color=color1,
        label=f"{module_type[0]} Water Flux",
    )
    ax2.plot(
        x[1:],
        y2,
        lw=2,
        ls=line_style[0],
        color=color2,
        label=f"{module_type[0]} Pressure",
    )
    ax3.plot(
        x[1:],
        y3,
        lw=2,
        ls=line_style[0],
        color=color3,
        label=f"{module_type[0]} Pressure Loss",
    )
    ax4.plot(
        x[1:],
        y4,
        lw=2,
        ls=line_style[0],
        color=color4,
        label=f"{module_type[0]} Osm Pressure",
    )

    ax.set_xlabel("Distance Along Membrane Length (m)")
    ax.set_ylabel("Flux (kg/m^2/hr)")
    ax2.set_ylabel("Pressure (psi)")
    ax3.set_ylabel("Pressure Loss (psi/m)")
    ax4.set_ylabel("Osmotic Pressure @ Interface (psi)")

    # custom_lines = [
    #     Line2D([0], [0], color="k", lw=2, linestyle=line_style[0]),
    #     Line2D([0], [4], color="k", lw=2, linestyle=line_style[1]),
    # ]
    # ax.legend(
    #     custom_lines,
    #     ["Flat-plate", "Spiral-wound"],
    #     loc="upper left",
    #     frameon=False,
    #     fontsize=8,
    # )
    # ax2.legend(loc='upper right', frameon=False, fontsize=8)

    ax2.spines["left"].set_color(color1)
    ax2.spines["right"].set_color(color2)
    ax3.spines["right"].set_color(color3)
    ax4.spines["right"].set_color(color4)
    ax.tick_params(axis="y", colors=color1)
    ax2.tick_params(axis="y", colors=color2)
    ax3.tick_params(axis="y", colors=color3)
    ax4.tick_params(axis="y", colors=color4)
    ax.yaxis.label.set_color(color1)
    ax2.yaxis.label.set_color(color2)
    ax3.yaxis.label.set_color(color3)
    ax4.yaxis.label.set_color(color4)

    plt.title("Flux Profile")
    # # Set x and y limits
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 50)
    ax2.set_ylim(0, P0)
    ax3.set_ylim(0, P2)
    ax4.set_ylim(0, P3)
    fig.tight_layout()
    plt.show()


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
