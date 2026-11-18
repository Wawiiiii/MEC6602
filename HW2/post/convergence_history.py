"""Convergence histories of the report: density residual drop against the iteration, one figure
per case shown in the report body, one curve per grid (rewritten 2026-11-18: one panel per case,
600 DPI; the user found the multi-panel figures unreadable).

Inputs (read only): results/<folder>/<grid>_terminal*.txt, the per-iteration rows printed by the
solver (fine grid only; the coarse levels of the full-multigrid start are not printed), and
results/<folder>/<grid>_summary.txt. When a grid has several attempts (_terminal_a2.txt, ...) the
attempt used is the one whose last row (iteration count and time) matches the summary.
M 0.5 alpha 0 kept no terminal output (2026-11-13 runs) and has no figure (user 2026-11-18: the
alpha 1.25 figure is enough); a partial rerun (129 to 513) is in results/Mach_05_alpha_0_history.

Figures, report/figures/conv_<case>.png, 600 DPI, 3.4 x 2.6 in (as the Mach and Cp figures):
  conv_M05_a125, conv_M08_a0, conv_M08_a125, conv_M125_a0, conv_M125_a125 (grids
  129 to 2049), conv_M08_CL04 (fixed CL 0.4, grids 129 to 1025; the jumps are the alpha changes)
  conv_forces_M05_a125, conv_forces_M08_a125, conv_forces_M125_a125, conv_forces_M08_CL04: |C^n - C|
  of CL, CD, CM on the finest grid (2049, 1025 for fixed CL) against the iteration (C = converged
  value of the summary file)
Also prints, per case: grid, threads, iterations, fine-grid time, final drop, status, file used.
Run from the HW2 folder:  python3 post/convergence_history.py [case ...]
"""
import glob
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator

RESULTS = "results"
FIGDIR = os.path.join("report", "figures")
DPI = 600

# (figure tag, folder, grids)
GRIDS = [129, 257, 513, 1025, 2049]
CASES = [("M05_a125", "Mach_05_alpha_125", GRIDS),
         ("M08_a0", "Mach_08_alpha_0", GRIDS),
         ("M08_a125", "Mach_08_alpha_125", GRIDS),
         ("M125_a0", "Mach_125_alpha_0", GRIDS),
         ("M125_a125", "Mach_125_alpha_125", GRIDS),
         ("M08_CL04", "Mach_08_CL_04", GRIDS[:-1])]

# ordered colours light (coarse) -> dark (fine) that also change hue (viridis samples), so
# neighbouring grids stay distinguishable; markers repeat the identity for greyscale print
COLOURS = {129: "#7ad151", 257: "#22a884", 513: "#2a788e", 1025: "#414487", 2049: "#440154"}
MARKERS = {65: "o", 129: "s", 257: "^", 513: "D", 1025: "v", 2049: "P"}

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "font.size": 9,
})

NUM = r"([-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?)"
ROW = re.compile(r"^\s*(\d+)\s+" + r"\s+".join([NUM] * 6) + r"(\s|$)")  # iter + 6 numbers
FMG = re.compile(r"FMG start:\s*(\d+) x \d+, (\d+) iterations, ([0-9.]+) s, drop ([0-9.eE+-]+)"
                 r"(?:, alpha ([0-9.+-]+) deg)?")
CHANGE = re.compile(r"fixed CL: CL [0-9.eE+-]+ -> alpha ([0-9.+-]+) -> ([0-9.+-]+) deg")
HEADER_ALPHA = re.compile(r"^\s+alpha\s+([0-9.+-]+) deg")


def parse_terminal(path):
    """Rows and side information of one terminal file.
    Returns a dict: it, res, drop, cl, time (lists, one entry per printed row), alpha (alpha in
    effect after each row, i.e. including a change printed right after it), fmg (list of
    (grid, iterations, seconds, drop)), threads, status."""
    d = {"it": [], "res": [], "drop": [], "cl": [], "cd": [], "cm": [], "time": [], "alpha": [], "fmg": [],
         "threads": None, "status": None}
    head_alpha = fmg_alpha = None
    changes = []  # (index of the row it follows, alpha before, alpha after)
    for line in open(path):
        m = ROW.match(line)
        if m:
            d["it"].append(int(m.group(1)))
            d["res"].append(float(m.group(2)))
            d["drop"].append(float(m.group(3)))
            d["cl"].append(float(m.group(4)))
            d["cd"].append(float(m.group(5)))
            d["cm"].append(float(m.group(6)))
            d["time"].append(float(m.group(7)))
            continue
        m = CHANGE.search(line)
        if m:
            changes.append((len(d["it"]) - 1, float(m.group(1)), float(m.group(2))))
            continue
        m = FMG.search(line)
        if m:
            d["fmg"].append((int(m.group(1)), int(m.group(2)), float(m.group(3)),
                             float(m.group(4))))
            if m.group(5):
                fmg_alpha = float(m.group(5))
            continue
        m = HEADER_ALPHA.match(line)
        if m and head_alpha is None and not d["it"]:
            head_alpha = float(m.group(1))
        if line.strip().startswith("threads") and d["threads"] is None:
            d["threads"] = int(line.split()[1])
        if line.strip().startswith("status"):
            d["status"] = line.split(None, 1)[1].strip()
    # starting alpha of the fine grid: first change "from", else last FMG alpha, else header
    a = changes[0][1] if changes else (fmg_alpha if fmg_alpha is not None else head_alpha)
    k = 0
    for i in range(len(d["it"])):
        while k < len(changes) and changes[k][0] <= i:
            a = changes[k][2]
            k += 1
        d["alpha"].append(a)
    d["n_changes"] = len(changes)
    return d


def summary_last(folder, g):
    # (iterations, seconds, status, drop) from the last line of the summary file
    path = os.path.join(RESULTS, folder, "%dx%d_summary.txt" % (g, g))
    w = [l.split() for l in open(path) if l.strip() and not l.startswith("#")][-1]
    return int(w[6]), float(w[7]), w[8], float(w[9])


def summary_forces(folder, g):
    # converged CL, CD, CM (10 digits) from the last line of the summary file
    path = os.path.join(RESULTS, folder, "%dx%d_summary.txt" % (g, g))
    w = [l.split() for l in open(path) if l.strip() and not l.startswith("#")][-1]
    return float(w[3]), float(w[4]), float(w[5])


def load_grid(folder, g):
    """Terminal data of the attempt that produced the final result of grid g: the file whose
    last row matches the summary (iteration count and time)."""
    it_s, t_s, status, drop_s = summary_last(folder, g)
    files = sorted(glob.glob(os.path.join(RESULTS, folder, "%dx%d_terminal*.txt" % (g, g))))
    match = []
    for f in files:
        d = parse_terminal(f)
        if d["it"] and d["it"][-1] == it_s and abs(d["time"][-1] - t_s) < 0.01:
            match.append((f, d))
    if len(match) != 1:
        raise RuntimeError("%s %d: %d terminal files match the summary" % (folder, g, len(match)))
    f, d = match[0]
    d["file"] = os.path.basename(f)
    d["others"] = [os.path.basename(x) for x in files if x != f]
    d["summary"] = (it_s, t_s, status, drop_s)
    return d


def style_log_axis(ax, step, ymin):
    # log residual axis down to the decade below ymin: a labelled tick every <step> decades,
    # minor ticks at the other decades
    kmax = int(math.ceil(-math.log10(ymin)))
    ax.set_yscale("log")
    ax.set_ylim(10.0 ** -kmax, 3.0)
    ax.set_yticks([10.0 ** -k for k in range(0, kmax + 1, step)])
    ax.set_yticks([10.0 ** -k for k in range(0, kmax + 1)], minor=True)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major", color="0.85", lw=0.5)


def plot_kw(g, n):
    # line + about ten markers per curve, always one on the last point (identity also in
    # greyscale print)
    every = max(1, round(n / 10.0))
    marks = sorted(set(range(0, n, every)) | {n - 1})
    return dict(color=COLOURS[g], lw=1.0, marker=MARKERS[g], ms=3.2, mew=0.7, mfc="white",
                mec=COLOURS[g], markevery=marks, label="%d $\\times$ %d" % (g, g))


def print_table(title, folder, grids, data):
    print("\n" + title)
    print("  %-10s %4s %10s %6s %11s %10s  %-10s %s" % ("grid", "thr", "FMG_s", "iter", "fine_s",
                                                         "drop", "status", "terminal file"))
    for g in grids:
        d = data[g]
        it_s, t_s, status, drop_s = d["summary"]
        fmg_t = sum(x[2] for x in d["fmg"])
        print("  %-10s %4d %10.2f %6d %11.2f %10.2e  %-10s %s%s" % (
            "%dx%d" % (g, g), d["threads"], fmg_t, it_s, t_s, drop_s, status, d["file"],
            ("   (not used: %s)" % ", ".join(d["others"])) if d["others"] else ""))


def case_figure(tag, folder, grids):
    data = {g: load_grid(folder, g) for g in grids}
    print_table("%s (%s)" % (tag, folder), folder, grids, data)
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for g in grids:
        d = data[g]
        ax.plot(d["it"], d["drop"], **plot_kw(g, len(d["it"])))
    ymin = min(min(data[g]["drop"]) for g in grids)
    style_log_axis(ax, 2, ymin)
    ax.set_ylabel(r"$\|R_\rho\|_2 \,/\, \|R_\rho\|_{2,1}$")
    ax.set_xlabel("Iteration")
    ax.set_xlim(0, max(data[g]["it"][-1] for g in grids) + 1)
    # legend in a block above the axes: never on top of a curve
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=7, frameon=False,
              handlelength=2.2, columnspacing=1.0, handletextpad=0.4, borderaxespad=0.0)
    fig.tight_layout()
    fname = "conv_%s.png" % tag
    fig.savefig(os.path.join(FIGDIR, fname), dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  -> %s" % fname)


# force convergence on the finest grid, alpha = 1.25 deg (Section 3 cases): (figure tag, folder)
# and the fixed-lift case on its finest grid (1025)
FORCE_CASES = [("M05_a125", "Mach_05_alpha_125", 2049),
               ("M08_a125", "Mach_08_alpha_125", 2049),
               ("M125_a125", "Mach_125_alpha_125", 2049),
               ("M08_CL04", "Mach_08_CL_04", 1025)]
# CL, CD, CM: distinct hues (colour-blind safe) + distinct markers and dashes
FORCE_STYLE = {"cl": ("#0072B2", "o", "-", r"$C_L$"),
               "cd": ("#D55E00", "s", "--", r"$C_D$"),
               "cm": ("#009E73", "^", ":", r"$C_M$")}


def force_figure(tag, folder, grid):
    """|C^n - C| of CL, CD, CM against the iteration on the finest grid, C = converged value of the
    summary file (10 digits); the iteration rows print 7 significant digits, which sets a floor
    of about 1e-8 (rows equal to C at that precision are not drawn on the log axis)."""
    d = load_grid(folder, grid)
    final = dict(zip(("cl", "cd", "cm"), summary_forces(folder, grid)))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    lo = 1.0
    for k in ("cl", "cd", "cm"):
        col, mk, ls, lab = FORCE_STYLE[k]
        pts = [(i, abs(v - final[k])) for i, v in zip(d["it"], d[k]) if abs(v - final[k]) > 0.0]
        lo = min(lo, min(e for _, e in pts))
        n = len(pts)
        every = max(1, round(n / 10.0))
        ax.semilogy([p[0] for p in pts], [p[1] for p in pts], color=col, lw=1.0, ls=ls, marker=mk,
                    ms=3.2, mew=0.7, mfc="white", mec=col,
                    markevery=sorted(set(range(0, n, every)) | {n - 1}), label=lab)
    kmax = int(math.ceil(-math.log10(lo)))
    kmin = int(math.floor(-math.log10(max(max(abs(v - final[k]) for v in d[k]) for k in final))))
    ax.set_ylim(10.0 ** -kmax, 10.0 ** -kmin * 1.5)
    ax.set_yticks([10.0 ** -k for k in range(kmin, kmax + 1)])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.grid(True, which="major", color="0.85", lw=0.5)
    ax.set_ylabel(r"$|C^{\,n}-C|$")
    ax.set_xlabel("Iteration (%d $\\times$ %d)" % (grid, grid))
    ax.set_xlim(0, d["it"][-1] + 1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=7, frameon=False,
              handlelength=2.6, columnspacing=1.2, handletextpad=0.4, borderaxespad=0.0)
    fig.tight_layout()
    fname = "conv_forces_%s.png" % tag
    fig.savefig(os.path.join(FIGDIR, fname), dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  -> %s (%s %dx%d: %d iterations, status %s)" % (fname, folder, grid, grid,
                                                            d["it"][-1], d["summary"][2]))


if __name__ == "__main__":
    import sys
    os.makedirs(FIGDIR, exist_ok=True)
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    want = sys.argv[1:]
    for tag, folder, grids in CASES:
        if not want or tag in want:
            case_figure(tag, folder, grids)
    for tag, folder, grid in FORCE_CASES:
        if not want or "forces_" + tag in want:
            force_figure(tag, folder, grid)
