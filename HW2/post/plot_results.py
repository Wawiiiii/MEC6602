# plot_results.py
# Post-processing for the MEC6602E HW2 2D Euler solver (NACA0012, O-grid).
#
# Reads the three kinds of output files written by euler_main:
#   - history file  : iteration, ||R_rho||, ||R_rho||/||R_rho||_0, CL, CD, CM, seconds
#   - surface file  : wall Cp (Tecplot POINT zone), with CL/CD/CM in a comment line
#   - field file    : Tecplot ordered BLOCK zone (X, Y, Rho, U, V, P, Mach, Cp, Entropy)
# and writes PNG figures into an output folder.
#
# Every input option can be repeated to overlay several runs on the same figure:
#   python plot_results.py --history h33.dat --history h65.dat --surface s33.dat --surface s65.dat
#
# Usage:
#   python plot_results.py --history h.dat --surface s.dat --field f.dat --out plots
#   (optional) --case "M = 0.5, alpha = 1.25 deg"  to force the text used in the titles

import argparse
import os
import re

import numpy as np

# Use the non-interactive Agg backend: figures are only saved to PNG files
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Window used for the field contour plots (in chords)
X_MIN = -0.5
X_MAX = 1.5
Y_MIN = -0.75
Y_MAX = 0.75


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def label_from_path(path):
    # Legend label = file name without folder and without extension
    base = os.path.basename(path)
    label = os.path.splitext(base)[0]
    return label


def find_number(text, name):
    # Look for "name = <number>" in a line of text and return the number (or None)
    pattern = name + r"\s*=\s*([-+0-9.eE]+)"
    match = re.search(pattern, text)
    if match is None:
        return None
    return float(match.group(1))


def case_string(mach, alpha):
    # Build a short "M = ..., alpha = ..." string for the titles
    if mach is None or alpha is None:
        return ""
    return "M = %.3g, $\\alpha$ = %.3g$^\\circ$" % (mach, alpha)


def save_figure(fig, out_dir, name, written):
    # Save one figure as PNG, close it, and remember the file name
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    written.append(path)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_history(path):
    # History file: one header line starting with '#', then 7 columns per row
    data = np.loadtxt(path, comments="#", ndmin=2)
    hist = {}
    hist["iter"] = data[:, 0]
    hist["res"] = data[:, 1]
    hist["res_rel"] = data[:, 2]
    hist["CL"] = data[:, 3]
    hist["CD"] = data[:, 4]
    hist["CM"] = data[:, 5]
    hist["time"] = data[:, 6]
    hist["label"] = label_from_path(path)
    return hist


def read_surface(path):
    # Surface file: '#' comment lines, Tecplot header lines, then rows "x/c y/c Cp p"
    rows = []
    mach = None
    alpha = None
    cl = None
    cd = None
    cm = None
    with open(path, "r") as f:
        for line in f:
            text = line.strip()
            if text == "":
                continue
            # Comment lines may contain the case parameters and the force coefficients
            if text.startswith("#"):
                if mach is None:
                    mach = find_number(text, "M")
                if alpha is None:
                    alpha = find_number(text, "alpha")
                if cl is None:
                    cl = find_number(text, "CL")
                if cd is None:
                    cd = find_number(text, "CD")
                if cm is None:
                    cm = find_number(text, "CM")
                continue
            # Try to read the line as numbers; Tecplot header lines (TITLE, ZONE...) fail and are skipped
            try:
                values = [float(v) for v in text.split()]
            except ValueError:
                continue
            if len(values) >= 3:
                rows.append(values[:4])
    data = np.array(rows)
    surf = {}
    surf["x"] = data[:, 0]
    surf["y"] = data[:, 1]
    surf["Cp"] = data[:, 2]
    surf["p"] = data[:, 3] if data.shape[1] > 3 else None
    surf["M"] = mach
    surf["alpha"] = alpha
    surf["CL"] = cl
    surf["CD"] = cd
    surf["CM"] = cm
    surf["label"] = label_from_path(path)
    return surf


def read_field(path):
    # Field file: Tecplot ASCII ordered zone with DATAPACKING=BLOCK
    with open(path, "r") as f:
        lines = f.readlines()

    # Line 1: TITLE (may contain M and alpha)
    title_line = lines[0]
    mach = find_number(title_line, "M")
    alpha = find_number(title_line, "alpha")

    # Line 2: VARIABLES = "X" "Y" ...  -> list of names between quotes
    names = re.findall(r'"([^"]*)"', lines[1])

    # Line 3: ZONE ... I=.., J=..
    zone_line = lines[2]
    ni = int(re.search(r"\bI\s*=\s*(\d+)", zone_line).group(1))
    nj = int(re.search(r"\bJ\s*=\s*(\d+)", zone_line).group(1))

    # All remaining numbers, in block order: variable 1 for all nodes, then variable 2, ...
    numbers = np.array(" ".join(lines[3:]).split(), dtype=float)
    n_nodes = ni * nj
    expected = n_nodes * len(names)
    if numbers.size != expected:
        raise ValueError("%s: found %d numbers, expected %d" % (path, numbers.size, expected))

    # Reshape each variable to (J, I): I is the fastest index in Tecplot ordering
    field = {}
    for k, name in enumerate(names):
        block = numbers[k * n_nodes:(k + 1) * n_nodes]
        field[name] = block.reshape(nj, ni)
    field["names"] = names
    field["ni"] = ni
    field["nj"] = nj
    field["M"] = mach
    field["alpha"] = alpha
    field["label"] = label_from_path(path)
    return field


# ---------------------------------------------------------------------------
# Quick sanity checks printed to the terminal
# ---------------------------------------------------------------------------

def check_history(hist):
    # Print final residual and final coefficients
    n = len(hist["iter"])
    orders = -np.log10(hist["res_rel"][-1]) if hist["res_rel"][-1] > 0 else np.inf
    print("[history %s] %d iterations, final ||R||/||R0|| = %.3e (%.1f orders), %.2f s"
          % (hist["label"], n, hist["res_rel"][-1], orders, hist["time"][-1]))
    print("    final CL = %.6f  CD = %.6f  CM = %.6f"
          % (hist["CL"][-1], hist["CD"][-1], hist["CM"][-1]))
    bad = np.sum(~np.isfinite(np.column_stack([hist["res"], hist["CL"], hist["CD"], hist["CM"]])))
    if bad > 0:
        print("    WARNING: %d non-finite values in the history" % bad)


def check_surface(surf):
    # Print number of points, Cp range and the trailing edge values
    n = len(surf["x"])
    bad = np.sum(~np.isfinite(surf["Cp"]))
    print("[surface %s] %d points, Cp in [%.4f, %.4f], CL = %s CD = %s CM = %s"
          % (surf["label"], n, np.nanmin(surf["Cp"]), np.nanmax(surf["Cp"]),
             surf["CL"], surf["CD"], surf["CM"]))
    print("    trailing edge Cp: first point %.4f, last point %.4f"
          % (surf["Cp"][0], surf["Cp"][-1]))
    if bad > 0:
        print("    WARNING: %d non-finite Cp values" % bad)


def check_field(field):
    # Print non-finite counts, min density/pressure, and the mismatch across the i=0 / i=I-1 cut
    print("[field %s] I = %d, J = %d, variables: %s"
          % (field["label"], field["ni"], field["nj"], " ".join(field["names"])))
    for name in field["names"]:
        values = field[name]
        bad = np.sum(~np.isfinite(values))
        cut = np.nanmax(np.abs(values[:, 0] - values[:, -1]))
        print("    %-8s min %+.4e  max %+.4e  |cut mismatch| %.2e%s"
              % (name, np.nanmin(values), np.nanmax(values), cut,
                 ("  NON-FINITE: %d" % bad) if bad > 0 else ""))


# ---------------------------------------------------------------------------
# History plots
# ---------------------------------------------------------------------------

def plot_residual(histories, case, out_dir, written):
    # Figure 1: relative density residual vs iteration
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for hist in histories:
        ax.semilogy(hist["iter"], hist["res_rel"], label=hist["label"])
    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$\|R_\rho\| \,/\, \|R_\rho\|_0$")
    ax.set_title("Density residual history" + ("  (" + case + ")" if case else ""))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save_figure(fig, out_dir, "residual_iter.png", written)

    # Figure 2: relative density residual vs wall-clock time
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for hist in histories:
        ax.semilogy(hist["time"], hist["res_rel"], label=hist["label"])
    ax.set_xlabel("Wall-clock time [s]")
    ax.set_ylabel(r"$\|R_\rho\| \,/\, \|R_\rho\|_0$")
    ax.set_title("Density residual vs run time" + ("  (" + case + ")" if case else ""))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save_figure(fig, out_dir, "residual_time.png", written)


def plot_forces(histories, case, out_dir, written):
    # Figure: CL, CD, CM vs iteration (3 stacked subplots sharing the x axis)
    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    keys = ["CL", "CD", "CM"]
    for ax, key in zip(axes, keys):
        for hist in histories:
            ax.plot(hist["iter"], hist[key], label="%s (final %.5f)" % (hist["label"], hist[key][-1]))
        ax.set_ylabel(key)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

        # The first iterations can be far from the converged value: zoom the y axis
        # on the second half of the run so the converged level is readable
        tail_values = np.concatenate([h[key][len(h[key]) // 2:] for h in histories])
        all_values = np.concatenate([h[key] for h in histories])
        tail_low = np.min(tail_values)
        tail_high = np.max(tail_values)
        full_low = np.min(all_values)
        full_high = np.max(all_values)
        span = max(tail_high - tail_low, 0.05 * (full_high - full_low), 1e-6)
        ax.set_ylim(max(full_low, tail_low - 10 * span) - 0.02 * span,
                    min(full_high, tail_high + 10 * span) + 0.02 * span)
    axes[-1].set_xlabel("Iteration")
    axes[0].set_title("Force coefficient history" + ("  (" + case + ")" if case else ""))
    fig.tight_layout()
    save_figure(fig, out_dir, "forces_iter.png", written)


# ---------------------------------------------------------------------------
# Surface Cp plot
# ---------------------------------------------------------------------------

def split_surface(surf):
    # Points go TE -> lower -> LE -> upper -> TE, so the leading edge is the minimum x/c.
    # The LE point is kept in both halves so the two curves join.
    i_le = int(np.argmin(surf["x"]))
    lower = slice(0, i_le + 1)
    upper = slice(i_le, len(surf["x"]))
    return lower, upper


def plot_cp(surfaces, case, out_dir, written):
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for k, surf in enumerate(surfaces):
        lower, upper = split_surface(surf)
        if len(surfaces) == 1:
            # One file: two colours for upper / lower surface
            ax.plot(surf["x"][upper], surf["Cp"][upper], "o-", ms=3, color=colors[0], label="upper surface")
            ax.plot(surf["x"][lower], surf["Cp"][lower], "s-", ms=3, color=colors[1], label="lower surface")
        else:
            # Several files: one colour per file, solid = upper, dashed = lower
            color = colors[k % len(colors)]
            name = surf["label"]
            if surf["CL"] is not None:
                name = name + " (CL=%.4f, CD=%.4f)" % (surf["CL"], surf["CD"])
            ax.plot(surf["x"][upper], surf["Cp"][upper], "o-", ms=3, color=color, label=name + " upper")
            ax.plot(surf["x"][lower], surf["Cp"][lower], "s--", ms=3, color=color, label=name + " lower")

    # Suction side on top: negative Cp upward
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel(r"$C_p$")
    ax.set_title("Surface pressure coefficient" + ("  (" + case + ")" if case else ""))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    # Text box with the force coefficients when only one file is plotted
    if len(surfaces) == 1 and surfaces[0]["CL"] is not None:
        s = surfaces[0]
        text = "CL = %.5f\nCD = %.5f\nCM = %.5f" % (s["CL"], s["CD"], s["CM"])
        ax.text(0.97, 0.97, text, transform=ax.transAxes, ha="right", va="top",
                fontsize=9, family="monospace",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    save_figure(fig, out_dir, "cp_surface.png", written)


# ---------------------------------------------------------------------------
# Field plots
# ---------------------------------------------------------------------------

def values_in_window(field, values):
    # Keep only the nodes inside the plotting window (used to choose colour ranges)
    x = field["X"]
    y = field["Y"]
    inside = (x >= X_MIN) & (x <= X_MAX) & (y >= Y_MIN) & (y <= Y_MAX)
    return values[inside & np.isfinite(values)]


def draw_airfoil(ax, field):
    # The j = 0 line is the airfoil wall: fill it so the contours stop at the body
    xw = field["X"][0, :]
    yw = field["Y"][0, :]
    ax.fill(xw, yw, facecolor="lightgrey", edgecolor="black", linewidth=1.0, zorder=3)


def setup_window(ax, title):
    # Common axis settings for the near-field plots
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect("equal")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.set_title(title)


def plot_contour(field, name, levels, cmap, cbar_label, title, extend, out_dir, file_name, written):
    # Filled contours on the structured (J, I) arrays + a few black iso-lines
    fig, ax = plt.subplots(figsize=(8, 5.5))
    x = field["X"]
    y = field["Y"]
    q = field[name]
    cs = ax.contourf(x, y, q, levels=levels, cmap=cmap, extend=extend)
    ax.contour(x, y, q, levels=levels[::4], colors="k", linewidths=0.4, alpha=0.5)
    draw_airfoil(ax, field)
    setup_window(ax, title)
    cbar = fig.colorbar(cs, ax=ax, shrink=0.85)
    cbar.set_label(cbar_label)
    save_figure(fig, out_dir, file_name, written)


def plot_field(field, case, multi, out_dir, written):
    # Suffix with the file name only when several field files are given
    suffix = ("_" + field["label"]) if multi else ""
    head = field["label"] + ": " if multi else ""
    tail = ("  (" + case + ")") if case else ""

    # Mach number: colour range from the values seen near the airfoil
    mach = values_in_window(field, field["Mach"])
    levels = np.linspace(np.min(mach), np.max(mach), 33)
    plot_contour(field, "Mach", levels, "jet", "Mach", head + "Mach number" + tail,
                 "neither", out_dir, "field_mach" + suffix + ".png", written)

    # Pressure coefficient: reversed colour map so suction (negative Cp) is red
    cp = values_in_window(field, field["Cp"])
    levels = np.linspace(np.min(cp), np.max(cp), 33)
    plot_contour(field, "Cp", levels, "jet_r", r"$C_p$", head + "Pressure coefficient" + tail,
                 "neither", out_dir, "field_cp" + suffix + ".png", written)

    # Static pressure (non-dimensional, p_inf = 1/gamma)
    p = values_in_window(field, field["P"])
    levels = np.linspace(np.min(p), np.max(p), 33)
    plot_contour(field, "P", levels, "jet", r"$p$ (p$_\infty$ = 1/$\gamma$)", head + "Pressure" + tail,
                 "neither", out_dir, "field_pressure" + suffix + ".png", written)

    # Entropy error (p/rho^gamma relative to free stream - 1): values are tiny and can have
    # both signs, so use a symmetric range set by the 99th percentile of |s| near the body
    s = values_in_window(field, field["Entropy"])
    smax = np.percentile(np.abs(s), 99.0)
    if smax <= 0.0:
        smax = 1e-12
    levels = np.linspace(-smax, smax, 33)
    plot_contour(field, "Entropy", levels, "RdBu_r",
                 r"$(p/\rho^\gamma)/(p_\infty/\rho_\infty^\gamma) - 1$",
                 head + "Entropy error" + tail, "both", out_dir,
                 "field_entropy" + suffix + ".png", written)

    # Mesh: near-field view and full-domain view side by side
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    x = field["X"]
    y = field["Y"]
    for ax in axes:
        # Lines of constant j (around the airfoil) and constant i (wall to farfield)
        for j in range(field["nj"]):
            ax.plot(x[j, :], y[j, :], color="k", linewidth=0.3, alpha=0.6)
        for i in range(field["ni"]):
            ax.plot(x[:, i], y[:, i], color="k", linewidth=0.3, alpha=0.6)
        draw_airfoil(ax, field)
        ax.set_aspect("equal")
        ax.set_xlabel("x/c")
        ax.set_ylabel("y/c")
    axes[0].set_xlim(X_MIN, X_MAX)
    axes[0].set_ylim(Y_MIN, Y_MAX)
    axes[0].set_title(head + "Mesh near the airfoil (%d x %d)" % (field["ni"], field["nj"]))
    axes[1].set_title(head + "Full domain")
    save_figure(fig, out_dir, "mesh" + suffix + ".png", written)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot HW2 Euler solver results (history, surface Cp, field).")
    parser.add_argument("--history", action="append", default=[], help="history file (repeat to overlay)")
    parser.add_argument("--surface", action="append", default=[], help="surface Cp file (repeat to overlay)")
    parser.add_argument("--field", action="append", default=[], help="Tecplot field file (repeat for several)")
    parser.add_argument("--out", default="plots", help="output folder for the PNG files (default: plots)")
    parser.add_argument("--case", default=None, help="text for the titles, e.g. 'M = 0.5, alpha = 1.25 deg'")
    args = parser.parse_args()

    if not args.history and not args.surface and not args.field:
        parser.error("give at least one of --history, --surface, --field")

    # Create the output folder if needed
    os.makedirs(args.out, exist_ok=True)

    # Read everything first
    histories = [read_history(p) for p in args.history]
    surfaces = [read_surface(p) for p in args.surface]
    fields = [read_field(p) for p in args.field]

    # Case text for the history titles: --case if given, otherwise (only when there is a
    # single history file) M/alpha from the first surface or field file that has them,
    # because the history file itself does not store M and alpha
    case = args.case
    if case is None:
        case = ""
        candidates = surfaces + fields if len(histories) <= 1 else []
        for item in candidates:
            case = case_string(item["M"], item["alpha"])
            if case:
                break

    # Sanity checks in the terminal
    for hist in histories:
        check_history(hist)
    for surf in surfaces:
        check_surface(surf)
    for field in fields:
        check_field(field)

    # Plots
    written = []
    if histories:
        plot_residual(histories, case, args.out, written)
        plot_forces(histories, case, args.out, written)
    if surfaces:
        # With several surface files the M/alpha may differ: only show it for one file
        if args.case:
            surf_case = args.case
        elif len(surfaces) == 1:
            surf_case = case_string(surfaces[0]["M"], surfaces[0]["alpha"])
        else:
            surf_case = ""
        plot_cp(surfaces, surf_case, args.out, written)
    for field in fields:
        field_case = args.case if args.case else case_string(field["M"], field["alpha"])
        plot_field(field, field_case, len(fields) > 1, args.out, written)

    print("PNG files written:")
    for path in written:
        print("  " + path)


if __name__ == "__main__":
    main()
