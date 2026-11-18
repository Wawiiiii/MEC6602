"""Grid-convergence figures and table rows of the report (Section 2 and its appendices).

For each case (Mach number, limiter) reads results/<folder>/<grid>_summary.txt (last line of each
file) at alpha = 0 and alpha = 1.25 deg and writes, in the style of Vassberg & Jameson (2010),
Figs. 17-23, four panels report/figures/gc_<tag>_{cd_a0,cl_a125,cd_a125,cm_a125}.png:
  Log(abs(F - F*)) against Log(1/NC), open symbols F > F*, filled symbols F < F*,
  dotted line O(p) through the finest grid of the triplet used, grid sizes written vertically.
F* and p come from the finest triplet that gives an order (monotone with p > 0); the triplet used
is printed for the captions. The script also prints the LaTeX rows of the tables:
  force table with run times (M 0.5, Hemker-Koren) and, per case, the order table of the two
  finest triplets.
Run from the HW2 folder:  python3 post/grid_convergence_plots.py
"""
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

RESULTS = "results"
FIGDIR = os.path.join("report", "figures")
GRIDS = [9, 17, 33, 65, 129, 257, 513, 1025, 2049]  # nodes per direction; NC = nodes - 1 cells
DPI = 300

# tag, Mach, limiter, folder alpha = 0, folder alpha = 1.25
CASES = [
    ("M05HK", 0.5, "Hemker--Koren", "Mach_05_alpha_0", "Mach_05_alpha_125"),
    ("M05VA", 0.5, "Van Albada", "Mach_05_alpha_0_vanalbada", "Mach_05_alpha_125_vanalbada"),
    ("M08HK", 0.8, "Hemker--Koren", "Mach_08_alpha_0", "Mach_08_alpha_125"),
    ("M08VA", 0.8, "Van Albada", "Mach_08_alpha_0_vanalbada", "Mach_08_alpha_125_vanalbada"),
    ("M125HK", 1.25, "Hemker--Koren", "Mach_125_alpha_0", "Mach_125_alpha_125"),
]
# panels: (alpha 0 or 1.25, coefficient index in (CL, CD, CM), name, file suffix)
PANELS = [(0, 1, "Cd", "cd_a0"), (1, 0, "Cl", "cl_a125"), (1, 1, "Cd", "cd_a125"), (1, 2, "Cm", "cm_a125")]

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "font.size": 10,
})


def read_case(folder):
    # {NC: (CL, CD, CM, iterations, seconds)} from the last line of every summary file
    out = {}
    for g in GRIDS:
        path = os.path.join(RESULTS, folder, "%dx%d_summary.txt" % (g, g))
        rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
        w = rows[-1]
        out[g - 1] = (float(w[3]), float(w[4]), float(w[5]), int(w[6]), float(w[7]))
    return out


def richardson(fc, fm, ff):
    """Slides 26-28, refinement ratio 2: R = (Ff - Fm)/(Fm - Fc), p = -log2 R,
    F* = Ff + (Ff - Fm)/(2^p - 1). Returns (None, None) without an order: R <= 0 (not monotone)
    or R >= 1 (diverging, p <= 0)."""
    if fm == fc:
        return None, None
    r = (ff - fm) / (fm - fc)
    if r <= 0.0 or r >= 1.0:
        return None, None
    p = -math.log2(r)
    return p, ff + (ff - fm) / (2.0 ** p - 1.0)


def finest_order(nc, f):
    # finest triplet with an order: (index of its finest grid, p, F*)
    for i in range(len(f) - 1, 1, -1):
        p, fs = richardson(f[i - 2], f[i - 1], f[i])
        if p is not None:
            return i, p, fs
    raise ValueError("no triplet gives an order")


def fmt_star(v):
    return "%+.9f" % v


def error_plot(data, k, name, fname):
    nc = sorted(data)
    f = [data[n][k] for n in nc]
    i, p, fstar = finest_order(nc, f)
    x = [math.log10(1.0 / n) for n in nc]
    y = [math.log10(abs(v - fstar)) for v in f]

    fig, ax = plt.subplots(figsize=(4.6, 3.7))
    # axes and dotted grid as in the paper
    xmin, xmax = -3.5, -0.5
    ymin = math.floor(min(y)) - 0.5
    ymax = math.ceil(max(y)) + 0.5
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_locator(MultipleLocator(0.5 if a is ax.xaxis else 1.0))
        a.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.grid(True, which="major", color="black", ls=":", lw=0.6)
    ax.tick_params(direction="out", length=3, width=0.6)
    for s in ax.spines.values():
        s.set_linewidth(0.8)

    # O(p) through the finest grid of the triplet used
    xr = [xmin, xmax]
    ax.plot(xr, [y[i] + p * (xx - x[i]) for xx in xr], ":", color="black", lw=0.8, zorder=1)
    # solution curve, open F > F*, filled F < F*
    ax.plot(x, y, "-", color="black", lw=1.0, zorder=2)
    for xi, yi, v in zip(x, y, f):
        filled = v < fstar
        ax.plot(xi, yi, "o", ms=5, mec="black", mfc="black" if filled else "white", mew=0.9, zorder=3)
    # grid sizes, vertical, above local maxima and below the other points
    for j, (xi, yi, n) in enumerate(zip(x, y, nc)):
        up = (j == 0 or yi >= y[j - 1]) and (j == len(y) - 1 or yi >= y[j + 1])
        ax.annotate(str(n), (xi, yi), textcoords="offset points", xytext=(0, 7 if up else -7),
                    rotation=90, ha="center", va="bottom" if up else "top", fontsize=7)

    handles = [
        Line2D([], [], color="black", lw=1.0, marker="o", ms=4, mfc="white", mec="black",
               label="%s > %s*" % (name, name)),
        Line2D([], [], color="black", lw=1.0, marker="o", ms=4, mfc="black", mec="black",
               label="%s < %s*" % (name, name)),
        Line2D([], [], color="black", lw=0.8, ls=":", marker="+", ms=5, label="O( %.3f )" % p),
    ]
    leg = ax.legend(handles=handles, loc="upper left", fontsize=7, frameon=True, fancybox=False,
                    edgecolor="black", handlelength=3.0, borderpad=0.4)
    leg.get_frame().set_linewidth(0.6)

    ax.set_title("%s* = %s" % (name, fmt_star(fstar)), fontsize=10)
    ax.set_xlabel("Log ( 1 / NC )", fontweight="bold")
    ax.set_ylabel("Log ( abs( %s - %s* ) )" % (name, name), fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, fname), dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return nc[i - 2], nc[i - 1], nc[i], p, fstar


def order_rows(a0, a125):
    # LaTeX rows of an order table: two finest triplets, p and F*; "--" when there is no order
    rows = []
    for lab, d, k, sym in (("$0^\\circ$", a0, 1, "C_D"), ("$1.25^\\circ$", a125, 0, "C_L"),
                           ("", a125, 1, "C_D"), ("", a125, 2, "C_M")):
        cells = []
        for t in ((256, 512, 1024), (512, 1024, 2048)):
            p, fs = richardson(d[t[0]][k], d[t[1]][k], d[t[2]][k])
            cells.append("-- & --" if p is None else "%.2f & $%s$" % (p, fmt_star(fs)))
        rows.append("%s & $%s$ & %s \\\\" % (lab, sym, " & ".join(cells)))
    return rows


def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.makedirs(FIGDIR, exist_ok=True)
    for tag, mach, lim, f0, f125 in CASES:
        a0, a125 = read_case(f0), read_case(f125)
        print("\n%%%% %s: M %g, %s" % (tag, mach, lim))
        if tag == "M05HK":
            print("% force table rows: NC & NC^2 & CD(a0) & t(a0) & CL & CD & CM & t(a1.25)")
            for n in sorted(a0):
                print("%s & %s & $%+.9f$ & %s & $%+.9f$ & $%+.9f$ & $%+.9f$ & %s \\\\" % (
                    "{:,}".format(n).replace(",", "{,}"), "{:,}".format(n * n).replace(",", "{,}"),
                    a0[n][1], tsec(a0[n][4]), a125[n][0], a125[n][1], a125[n][2], tsec(a125[n][4])))
        print("% order table rows")
        for r in order_rows(a0, a125):
            print(r)
        print("% figures (triplet used, p, F*)")
        for ia, k, name, suf in PANELS:
            fname = "gc_%s_%s.png" % (tag, suf)
            c, m, ff, p, fs = error_plot(a125 if ia else a0, k, name, fname)
            print("%s: %d-%d-%d, p %.3f, F* %s" % (fname, c, m, ff, p, fmt_star(fs)))


def tsec(t):
    # run time in seconds, thousands separated for LaTeX
    if t < 0.005:
        return "$<$0.01"
    if t < 10:
        return "%.2f" % t
    if t < 1000:
        return "%.1f" % t
    return "{:,}".format(int(round(t))).replace(",", "{,}")


if __name__ == "__main__":
    main()
