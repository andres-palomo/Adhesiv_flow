#!/usr/bin/env python3
"""
make_figures.py -- turn the finite_vol_*.c output files into PNG figures.

This is a deliberately SIMPLE version: everywhere there was a choice
between a short "clever" numpy trick and a longer, more obvious one
(a for-loop, an if/else, a plain list), the obvious one was picked.
Nothing here needs more than: reading a text file with np.loadtxt,
indexing a 2D array with [row, col], a for loop, and basic
matplotlib calls (plt.plot, plt.pcolormesh, plt.streamplot).

WHAT THIS DOES
--------------
Each finite_vol_*.c program writes plain-text .dat files prefixed with
its own name, e.g. finite_vol_4_seringe_SIMPLE_final_pressure.dat. This
script reads those files for the five milestone cases and writes:

    <case>_fig_pressure.png      -- pressure heatmap
    <case>_fig_velocity.png      -- speed heatmap + streamlines
    <case>_fig_flowrate.png      -- Q(x) along the duct
    <case>_fig_convergence.png   -- |v|, mdiv and Q(in/out) vs iteration

and, for the two shear-thinning cases (5 and 6) only:

    <case>_fig_viscosity_profile.png -- shear rate and eta_a across
                                         the channel at one x-slice

It also writes the four comparison figures the report (tutorial.tex)
actually uses:

    fig_flowrate_comparison.png          -- legacy pressure vs SIMPLE
    fig_viscosity_profile.png            -- Newtonian vs power-law eta_a
    fig_velocity_profile_comparison.png  -- Newtonian vs power-law v_x,
                                             plain duct (left) and
                                             syringe throat (right),
                                             both at x=90 (see the
                                             comment on
                                             X_SLICE_THROAT_COMPARISON
                                             below for why)
    fig_flowrate_4case.png               -- Q(x) for all four cases

HOW TO RUN
----------
    cd shear-thinning          # the folder with the .dat files in it
    python3 make_figures.py

Needs numpy and matplotlib (pip install numpy matplotlib --break-system-packages
if they are missing).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")   # draw to PNG files, do not try to open a window
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Constants. Every case uses the same grid, so we just hardcode it
# instead of figuring it out from the data.
# ---------------------------------------------------------------------

NX = 100          # number of cells across (must match Nx in the C code)
NY = 20            # number of cells tall  (must match Ny in the C code)
U_IN = 0.1         # inlet velocity (must match U_in in the C code)
EXPECTED_Q = U_IN * NY   # the flow rate Q(x) should have everywhere
X_SLICE = 90       # the x column used for the per-case eta/shear-rate
                   # profile plots (plot_viscosity_profile) and for
                   # make_report_fig_viscosity_profile

# x-column used specifically for the Newtonian-vs-power-law v_x
# comparison figure (make_report_fig_velocity_profile_comparison).
# This one is NOT the same as X_SLICE.
#
# Three things changed here after comparing the numbers directly:
#
# 1. The original version of this figure put finite_vol_2_duct.c
#    (plain duct, Newtonian) next to
#    finite_vol_5_duct_SIMPLE_shearthinning.c (plain duct, power-law).
#    That's not a fair rheology comparison: finite_vol_2_duct.c still
#    uses the *legacy* pressure formula, which Table tab:case-comparison
#    in tutorial.tex shows loses 58% of the flow rate along a plain
#    duct with no geometry change at all (Q_in=1.673, Q_out=0.704,
#    against a target of 2.0). By any x station downstream, that run's
#    velocities are suppressed by the solver bug, not by rheology, so
#    the "shear-thinning is faster at the core" story that comparison
#    told was really "the mass-conserving run is faster than the one
#    that already lost most of its flow rate".
#
# 2. finite_vol_duct_SIMPLE.c fills that gap: plain duct, Newtonian,
#    SIMPLE pressure coupling -- the duct-side counterpart to
#    finite_vol_4_seringe_SIMPLE.c. It replaces finite_vol_2_duct.c in
#    this figure, so the duct panel is now finite_vol_duct_SIMPLE.c
#    (Newtonian) against finite_vol_5_duct_SIMPLE_shearthinning.c
#    (power-law) -- both SIMPLE, both mass-conserving (Q(x) flat at
#    ~2.0 the whole way), an apples-to-apples, rheology-only
#    comparison, the same way the syringe panel already was.
#
# 3. x=70 (the old slice) sits right at the syringe throat's entrance.
#    Scanning peak/average v_x across x showed:
#      - on the duct, the profile is fully developed (peak/avg ~1.456,
#        the classic parabolic value) from about x=30 all the way to
#        x=90, so any x in that range would do:
#      - on the syringe, the wide upstream section never actually
#        reaches fully-developed (parabolic) flow on its own -- it
#        stays close to plug-like (peak/avg ~1.08-1.11) out to about
#        x=40, then the ratio jumps sharply toward ~1.5 right before
#        the contraction (x=45-49), because the pressure field "feels"
#        the upcoming constriction and pulls the core forward in
#        advance -- pressure information travels upstream in an
#        incompressible solve. The throat itself only settles into its
#        own steady, fully-developed shape from about x=55 onward, and
#        stays flat out to x=90+.
#    x=90 is comfortably inside the fully-developed region for BOTH
#    geometries, well past every entrance/transient effect, so it's
#    used for both panels below.
X_SLICE_THROAT_COMPARISON = 90

# The six cases this script knows how to plot. "has_solid" is True for
# the syringe-geometry cases, whose .dat files carry one extra is_solid
# column marking the blocked-off constriction cells. "shear_thinning"
# is True for the two power-law cases, which also write a
# *_shear_eta.dat file.
CASES = [
    {"name": "finite_vol_2_duct",                       "has_solid": False, "shear_thinning": False},
    {"name": "finite_vol_3_seringe",                     "has_solid": True,  "shear_thinning": False},
    {"name": "finite_vol_4_seringe_SIMPLE",              "has_solid": True,  "shear_thinning": False},
    {"name": "finite_vol_duct_SIMPLE",                   "has_solid": False, "shear_thinning": False},
    {"name": "finite_vol_5_duct_SIMPLE_shearthinning",   "has_solid": False, "shear_thinning": True},
    {"name": "finite_vol_6_seringe_SIMPLE_shearthinning","has_solid": True,  "shear_thinning": True},
]


# ---------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------

def load_flowrate_profile(filename):
    """flow_rate_profile.dat is always two columns: x, Q(x)."""
    data = np.loadtxt(filename, comments="#")
    x = data[:, 0]
    Q = data[:, 1]
    return x, Q


def load_convergence_log(filename):
    """convergence_log.dat has a '# it |v| check ... Q_in ... Q_out ...'
    header whose middle columns differ a bit between cases (a duct case
    has Q_mid, a syringe case has Q_pre and Q_throat instead). We only
    need a few columns, so we just look up their position by name in
    the header instead of hardcoding a column number for each case."""
    with open(filename) as f:
        header_line = f.readline()
    column_names = header_line.replace("#", "").split()

    data = np.loadtxt(filename, comments="#")
    if data.ndim == 1:
        # Only one data row in the file (the case converged almost
        # immediately) -- np.loadtxt gives back a 1D array in that
        # case, so we reshape it to a 1-row 2D array to keep the
        # column indexing below working the same way every time.
        data = data.reshape(1, -1)

    iteration = data[:, column_names.index("it")]
    v_norm    = data[:, column_names.index("|v|")]
    check     = data[:, column_names.index("check")]
    mdiv      = data[:, column_names.index("mdiv")]
    q_in      = data[:, column_names.index("Q_in")]
    q_out     = data[:, column_names.index("Q_out")]
    return iteration, v_norm, check, mdiv, q_in, q_out


def load_2d_field(filename, num_value_columns, has_solid_column):
    """Read a 2D field file. Every line is:

        x  y  value_1  [value_2 ...]  [is_solid]

    with blank lines separating x-blocks (np.loadtxt skips blank lines
    on its own). Returns:

        values -- a list of NY-by-NX arrays, one per value column
        solid  -- an NY-by-NX array of True/False (only if
                  has_solid_column is True, otherwise None)

    We build the grids with a plain for loop over the rows, instead of
    a vectorized numpy trick, so it is obvious what is happening: for
    every (x, y) in the file, put its value(s) into that cell of the
    grid.
    """
    data = np.loadtxt(filename, comments="#")

    values = []
    for i in range(num_value_columns):
        values.append(np.zeros((NY, NX)))

    solid = None
    if has_solid_column:
        solid = np.zeros((NY, NX), dtype=bool)

    number_of_rows = data.shape[0]
    for row in range(number_of_rows):
        x = int(data[row, 0])
        y = int(data[row, 1])
        for i in range(num_value_columns):
            values[i][y, x] = data[row, 2 + i]
        if has_solid_column:
            is_solid_value = data[row, 2 + num_value_columns]
            solid[y, x] = (is_solid_value > 0.5)

    return values, solid


# ---------------------------------------------------------------------
# Small helpers for dealing with the syringe's solid (blocked) cells.
# No masked arrays here -- just plain arrays with NaN or 0 written into
# the solid cells, using an ordinary for loop.
# ---------------------------------------------------------------------

def hide_solid_cells(field, solid):
    """Return a copy of field with NaN wherever solid is True. A
    pcolormesh plot leaves NaN cells unfilled, so they show through as
    whatever background color the axes has (we set that to light gray
    below) -- a cheap way to draw the constriction wall."""
    result = field.copy()
    for y in range(NY):
        for x in range(NX):
            if solid[y, x]:
                result[y, x] = np.nan
    return result


def zero_solid_cells(field, solid):
    """Return a copy of field with 0 wherever solid is True.
    streamplot cannot handle NaN values, so for the velocity used in
    streamlines we zero out the solid region instead of using NaN; it
    is already hidden visually by the gray speed plot underneath it."""
    result = field.copy()
    for y in range(NY):
        for x in range(NX):
            if solid[y, x]:
                result[y, x] = 0.0
    return result


def remove_solid_rows(y_values, field_column, solid_column):
    """Given a 1D slice through the grid (e.g. one x-column of vx) and
    the matching 1D slice of the solid mask, return only the entries
    that are NOT solid. Used for line plots (profiles), where a solid
    cell should simply not appear instead of being drawn as zero."""
    keep = ~solid_column   # True where the cell is fluid, not solid
    return y_values[keep], field_column[keep]


# ---------------------------------------------------------------------
# Per-case figures
# ---------------------------------------------------------------------

def plot_pressure(case_name, has_solid):
    filename = case_name + "_final_pressure.dat"
    values, solid = load_2d_field(filename, num_value_columns=1, has_solid_column=has_solid)
    pressure = values[0]
    if has_solid:
        pressure = hide_solid_cells(pressure, solid)

    plt.figure(figsize=(8, 8 * NY / NX + 1))
    axis = plt.gca()
    axis.set_facecolor("lightgray")
    mesh = plt.pcolormesh(pressure, cmap="viridis")
    plt.colorbar(mesh, label="pressure p")
    plt.xlabel("x (cell)")
    plt.ylabel("y (cell)")
    plt.title(case_name + ": final pressure field")
    axis.set_aspect("equal")
    plt.tight_layout()

    outname = case_name + "_fig_pressure.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print("wrote", outname)


def plot_velocity(case_name, has_solid):
    filename = case_name + "_final_velocity.dat"
    values, solid = load_2d_field(filename, num_value_columns=2, has_solid_column=has_solid)
    vx, vy = values[0], values[1]

    speed = np.sqrt(vx**2 + vy**2)

    if has_solid:
        speed_to_plot = hide_solid_cells(speed, solid)
        vx_for_lines = zero_solid_cells(vx, solid)
        vy_for_lines = zero_solid_cells(vy, solid)
    else:
        speed_to_plot = speed
        vx_for_lines = vx
        vy_for_lines = vy

    plt.figure(figsize=(8, 8 * NY / NX + 1))
    axis = plt.gca()
    axis.set_facecolor("lightgray")
    mesh = plt.pcolormesh(speed_to_plot, cmap="magma")
    plt.colorbar(mesh, label="speed |v|")

    # streamplot wants the coordinates of the CENTER of each cell
    x_centers = np.arange(NX) + 0.5
    y_centers = np.arange(NY) + 0.5
    plt.streamplot(x_centers, y_centers, vx_for_lines, vy_for_lines,
                    color="white", density=1.1, linewidth=0.7, arrowsize=0.8)

    plt.xlabel("x (cell)")
    plt.ylabel("y (cell)")
    axis.set_aspect("equal")
    axis.set_xlim(0, NX)
    axis.set_ylim(0, NY)
    plt.tight_layout()

    outname = case_name + "_fig_velocity.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print("wrote", outname)


def plot_flowrate(case_name, has_solid):
    x, Q = load_flowrate_profile(case_name + "_flow_rate_profile.dat")

    plt.figure(figsize=(7, 4))
    plt.plot(x, Q, linewidth=1.8, color="tab:blue")
    plt.axhline(EXPECTED_Q, linestyle=":", color="black",
                label=f"expected Q = U_in * Ny = {EXPECTED_Q:.2f}")
    if has_solid:
        plt.axvline(49, linestyle="--", color="gray", label="constriction starts")
    plt.xlabel("x (cell)")
    plt.ylabel("Q(x)")
    plt.title(case_name + ": flow-rate profile")
    plt.legend()
    plt.tight_layout()

    outname = case_name + "_fig_flowrate.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print("wrote", outname)


def plot_convergence(case_name):
    iteration, v_norm, check, mdiv, q_in, q_out = load_convergence_log(case_name + "_convergence_log.dat")

    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(11, 4))

    # marker="o" makes sure a case that converged in a single logged
    # iteration (only one row in the file) still shows up as a visible
    # dot instead of an invisible zero-length line.
    left_axis.plot(iteration, v_norm, marker="o", label="|v|")
    left_axis.plot(iteration, check, marker="o", label="relative change")
    left_axis.plot(iteration, mdiv, marker="o", label="max|div(v)|")
    left_axis.set_yscale("log")
    left_axis.set_xlabel("iteration")
    left_axis.set_title(case_name + ": convergence")
    left_axis.legend(fontsize=8)

    right_axis.plot(iteration, q_in, marker="o", label="Q_in")
    right_axis.plot(iteration, q_out, marker="o", label="Q_out")
    right_axis.axhline(EXPECTED_Q, linestyle=":", color="black", label="expected Q")
    right_axis.set_xlabel("iteration")
    right_axis.set_title(case_name + ": flow rate vs iteration")
    right_axis.legend(fontsize=8)

    figure.tight_layout()
    outname = case_name + "_fig_convergence.png"
    figure.savefig(outname, dpi=150)
    plt.close(figure)
    print("wrote", outname)


def plot_viscosity_profile(case_name, has_solid):
    filename = case_name + "_shear_eta.dat"
    values, solid = load_2d_field(filename, num_value_columns=2, has_solid_column=has_solid)
    gamma_dot, eta = values[0], values[1]

    y = np.arange(NY)
    gamma_column = gamma_dot[:, X_SLICE]
    eta_column = eta[:, X_SLICE]
    if has_solid:
        solid_column = solid[:, X_SLICE]
        y, gamma_column = remove_solid_rows(y, gamma_column, solid_column)
        _, eta_column = remove_solid_rows(np.arange(NY), eta_column, solid_column)

    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(11, 4.3))

    left_axis.plot(y, gamma_column, "o-", markersize=3)
    left_axis.set_xlabel("y (cell)")
    left_axis.set_ylabel("shear rate gamma_dot")
    left_axis.set_title(f"shear rate at x={X_SLICE}", fontsize=10)

    right_axis.plot(y, eta_column, "o-", markersize=3, color="tab:orange")
    right_axis.set_xlabel("y (cell)")
    right_axis.set_ylabel("eta_a (apparent viscosity)")
    right_axis.set_title(f"eta_a at x={X_SLICE}", fontsize=10)

    figure.suptitle(case_name, fontsize=11)
    figure.tight_layout()
    outname = case_name + "_fig_viscosity_profile.png"
    figure.savefig(outname, dpi=150)
    plt.close(figure)
    print("wrote", outname)


# ---------------------------------------------------------------------
# Report figures (the exact filenames tutorial.tex includes)
# ---------------------------------------------------------------------

def get_vx_profile(case_name, has_solid, x_slice=X_SLICE):
    """v_x(y) at the given x column for one case, with any solid cells
    dropped. x_slice defaults to X_SLICE (kept for anything else that
    still wants the old default), but callers pass the slice they
    actually want explicitly -- see X_SLICE_THROAT_COMPARISON above."""
    values, solid = load_2d_field(case_name + "_final_velocity.dat",
                                   num_value_columns=2, has_solid_column=has_solid)
    vx = values[0]
    y = np.arange(NY)
    vx_column = vx[:, x_slice]
    if has_solid:
        solid_column = solid[:, x_slice]
        y, vx_column = remove_solid_rows(y, vx_column, solid_column)
    return y, vx_column


def make_report_fig_flowrate_comparison():
    x_legacy, Q_legacy = load_flowrate_profile("finite_vol_3_seringe_flow_rate_profile.dat")
    x_simple, Q_simple = load_flowrate_profile("finite_vol_4_seringe_SIMPLE_flow_rate_profile.dat")

    plt.figure(figsize=(7.5, 4.5))
    plt.plot(x_legacy, Q_legacy, color="tab:red", linewidth=1.8,
             label="finite_vol_3_seringe (legacy pressure)")
    plt.plot(x_simple, Q_simple, color="tab:green", linewidth=1.8,
             label="finite_vol_4_seringe_SIMPLE (SIMPLE)")
    plt.axvline(49, linestyle="--", color="gray")
    plt.axhline(EXPECTED_Q, linestyle=":", color="black")
    plt.xlabel("x (cell)")
    plt.ylabel("Q(x)")
    plt.legend()
    plt.tight_layout()

    outname = "fig_flowrate_comparison.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print("wrote", outname)


def make_report_fig_viscosity_profile(mu_const=1.0):
    values, solid = load_2d_field("finite_vol_5_duct_SIMPLE_shearthinning_shear_eta.dat",
                                   num_value_columns=2, has_solid_column=False)
    eta_powerlaw = values[1]
    eta_column = eta_powerlaw[:, X_SLICE]

    y = np.arange(NY)
    eta_newtonian_column = np.full(NY, mu_const)   # constant viscosity, same value everywhere

    plt.figure(figsize=(6, 4.5))
    plt.plot(y, eta_newtonian_column, "o-", markersize=3, color="tab:blue",
             label="finite_vol_2_duct (Newtonian)")
    plt.plot(y, eta_column, "o-", markersize=3, color="tab:orange",
             label="finite_vol_5_duct_SIMPLE_shearthinning (power-law)")
    plt.xlabel("y (cell)")
    plt.ylabel("eta_a")
    plt.title(f"Apparent viscosity at x={X_SLICE}")
    # The power-law curve peaks in the middle of the plot, so a legend
    # placed "best" (the default) tends to land right on top of it.
    # Put it below the plot instead, where there is always empty space.
    plt.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=1)
    plt.tight_layout()

    outname = "fig_viscosity_profile.png"
    plt.savefig(outname, dpi=150, bbox_inches="tight")   # keep the legend below the axes from being cut off
    plt.close()
    print("wrote", outname)


def make_report_fig_velocity_profile_comparison():
    """Newtonian vs power-law v_x, both at x=X_SLICE_THROAT_COMPARISON
    (90): plain duct (left) and syringe throat (right). Both panels now
    compare a SIMPLE/SIMPLE pair (finite_vol_duct_SIMPLE vs
    finite_vol_5_duct_SIMPLE_shearthinning on the left,
    finite_vol_4_seringe_SIMPLE vs
    finite_vol_6_seringe_SIMPLE_shearthinning on the right), so each
    panel isolates the rheology difference instead of mixing it with a
    mass-conservation difference -- see the comment on
    X_SLICE_THROAT_COMPARISON above for why finite_vol_2_duct.c was
    dropped from the duct panel and why x=90 was picked over the
    original x=70."""
    x_slice = X_SLICE_THROAT_COMPARISON

    y_duct_newt, vx_duct_newt = get_vx_profile("finite_vol_duct_SIMPLE", has_solid=False, x_slice=x_slice)
    y_duct_pl,   vx_duct_pl   = get_vx_profile("finite_vol_5_duct_SIMPLE_shearthinning", has_solid=False, x_slice=x_slice)
    y_ser_newt,  vx_ser_newt  = get_vx_profile("finite_vol_4_seringe_SIMPLE", has_solid=True, x_slice=x_slice)
    y_ser_pl,    vx_ser_pl    = get_vx_profile("finite_vol_6_seringe_SIMPLE_shearthinning", has_solid=True, x_slice=x_slice)

    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    left_axis.plot(vx_duct_newt, y_duct_newt, "o-", markersize=3, color="tab:blue",
                    label="finite_vol_duct_SIMPLE (Newtonian)")
    left_axis.plot(vx_duct_pl, y_duct_pl, "o-", markersize=3, color="tab:orange",
                    label="finite_vol_5_duct_SIMPLE_shearthinning (power-law)")
    left_axis.set_xlabel(r"$v_x$")
    left_axis.set_ylabel("y (cell)")
    left_axis.set_title(f"Plain duct, x={x_slice}")
    left_axis.legend(fontsize=7)

    right_axis.plot(vx_ser_newt, y_ser_newt, "o-", markersize=3, color="tab:blue",
                     label="finite_vol_4_seringe_SIMPLE (Newtonian)")
    right_axis.plot(vx_ser_pl, y_ser_pl, "o-", markersize=3, color="tab:orange",
                     label="finite_vol_6_seringe_SIMPLE_shearthinning (power-law)")
    right_axis.set_xlabel(r"$v_x$")
    right_axis.set_title(f"Syringe throat, x={x_slice}")
    right_axis.legend(fontsize=7)

    figure.tight_layout()
    outname = "fig_velocity_profile_comparison.png"
    figure.savefig(outname, dpi=150)
    plt.close(figure)
    print("wrote", outname)


def make_report_fig_flowrate_4case():
    x1, Q1 = load_flowrate_profile("finite_vol_2_duct_flow_rate_profile.dat")
    x2, Q2 = load_flowrate_profile("finite_vol_5_duct_SIMPLE_shearthinning_flow_rate_profile.dat")
    x3, Q3 = load_flowrate_profile("finite_vol_4_seringe_SIMPLE_flow_rate_profile.dat")
    x4, Q4 = load_flowrate_profile("finite_vol_6_seringe_SIMPLE_shearthinning_flow_rate_profile.dat")

    plt.figure(figsize=(8, 4.8))
    plt.plot(x1, Q1, color="tab:blue",   linestyle="-",  linewidth=1.8, label="duct, Newtonian")
    plt.plot(x2, Q2, color="tab:orange", linestyle="-",  linewidth=1.8, label="duct, power-law")
    plt.plot(x3, Q3, color="tab:blue",   linestyle="--", linewidth=1.8, label="syringe, Newtonian")
    plt.plot(x4, Q4, color="tab:orange", linestyle="--", linewidth=1.8, label="syringe, power-law")
    plt.axvline(49, linestyle=":", color="gray", label="constriction starts")
    plt.axhline(EXPECTED_Q, linestyle=":", color="black", label=f"expected Q = {EXPECTED_Q:.2f}")
    plt.xlabel("x (cell)")
    plt.ylabel("Q(x)")
    plt.legend(fontsize=8)
    plt.tight_layout()

    outname = "fig_flowrate_4case.png"
    plt.savefig(outname, dpi=150)
    plt.close()
    print("wrote", outname)


# ---------------------------------------------------------------------
# Main: just call every figure function, one after another.
# ---------------------------------------------------------------------

def main():
    for case in CASES:
        name = case["name"]
        has_solid = case["has_solid"]

        plot_pressure(name, has_solid)
        plot_velocity(name, has_solid)
        plot_flowrate(name, has_solid)
        plot_convergence(name)
        if case["shear_thinning"]:
            plot_viscosity_profile(name, has_solid)

    make_report_fig_flowrate_comparison()
    make_report_fig_viscosity_profile()
    make_report_fig_velocity_profile_comparison()
    make_report_fig_flowrate_4case()

    print("\nDone.")


if __name__ == "__main__":
    main()