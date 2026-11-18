"""Flow-field and surface-pressure figures and force tables of the report (Section 3 and Appendix J).

Figures (finest grid, 2049 x 2049, Hemker-Koren limiter), written to report/figures/ at 300 DPI,
for each Mach number <tag> = 05, 08, 125 and alpha <a> = 0, 125 (1.25 deg):
  mach_M<tag>_a<a>.png    : Mach number around the airfoil (colour map + iso-lines, sonic line in bold)
  cp_M<tag>_a<a>.png      : surface Cp (small black open circles)
Tables (printed, LaTeX rows like Table 1): CD at alpha = 0, CL CD CM at alpha = 1.25, run times
(fine-grid solve, seconds), dagger = stopped by the Cauchy criterion.
Also prints the numbers quoted in the text (peak Mach, shock positions, bow-shock stand-off).

The 2049 field files are ~600 MB of ASCII; X, Y, Mach, p near the airfoil are cached in
post/.cache/ (npz) after the first read.
Run from the HW2 folder:  python3 post/flow_results.py
"""
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable

RESULTS = "results"
FIGDIR = os.path.join("report", "figures")
CACHE = os.path.join("post", ".cache")
GRIDS = [9, 17, 33, 65, 129, 257, 513, 1025, 2049]
GAMMA = 1.4
DPI = 300
RMAX = 4.0  # field rows kept: nodes within this distance (chords) of the mid-chord

MACHS = [("05", 0.5), ("08", 0.8), ("125", 1.25)]

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "font.size": 9,
})


# ---------------------------------------------------------------- tables

def read_case(folder):
    # {NC: (CL, CD, CM, iterations, seconds, status)} from the last line of every summary file
    out = {}
    for g in GRIDS:
        path = os.path.join(RESULTS, folder, "%dx%d_summary.txt" % (g, g))
        w = [l.split() for l in open(path) if l.strip() and not l.startswith("#")][-1]
        out[g - 1] = (float(w[3]), float(w[4]), float(w[5]), int(w[6]), float(w[7]), w[8])
    return out


def tsec(t, status):
    if t < 0.005:
        s = "$<$0.01"
    elif t < 10:
        s = "%.2f" % t
    elif t < 1000:
        s = "%.1f" % t
    else:
        s = "{:,}".format(int(round(t))).replace(",", "{,}")
    return s + ("$^\\dagger$" if status == "cauchy" else "")


def table_rows(f0, f125):
    a0, a125 = read_case(f0), read_case(f125)
    rows = []
    for n in sorted(a0):
        rows.append("%s & %s & $%+.9f$ & %s & $%+.9f$ & $%+.9f$ & $%+.9f$ & %s \\\\" % (
            "{:,}".format(n).replace(",", "{,}"), "{:,}".format(n * n).replace(",", "{,}"),
            a0[n][1], tsec(a0[n][4], a0[n][5]), a125[n][0], a125[n][1], a125[n][2],
            tsec(a125[n][4], a125[n][5])))
    return rows


# ---------------------------------------------------------------- field

def load_field(folder):
    # X, Y, Mach, p/p_inf of the 2049 field near the airfoil (rows j with r < RMAX), cached
    os.makedirs(CACHE, exist_ok=True)
    cache = os.path.join(CACHE, folder + "_f2049.npz")
    if os.path.exists(cache):
        d = np.load(cache)
        return d["X"], d["Y"], d["M"], d["P"]
    path = os.path.join(RESULTS, folder, "f2049x2049.dat")
    with open(path, "rb") as fh:
        head = [fh.readline() for _ in range(3)]
        zone = head[2].decode()
        ni = int(zone.split("I=")[1].split(",")[0])
        nj = int(zone.split("J=")[1].split(",")[0])
        data = np.array(fh.read().split(), dtype=float)
    n = ni * nj
    # variables: X Y Rho U V P Mach Cp Entropy
    X, Y = data[0:n].reshape(nj, ni), data[n:2 * n].reshape(nj, ni)
    M = data[6 * n:7 * n].reshape(nj, ni)
    P = data[5 * n:6 * n].reshape(nj, ni) * GAMMA  # p / p_inf (p_inf = 1/gamma)
    r = np.sqrt((X - 0.5) ** 2 + Y ** 2).max(axis=1)
    jmax = int(np.argmax(r > RMAX)) + 1
    X, Y, M, P = X[:jmax], Y[:jmax], M[:jmax], P[:jmax]
    np.savez(cache, X=X, Y=Y, M=M, P=P)
    return X, Y, M, P


# colour range and iso-line step of each field, per free-stream Mach number (same for both alphas)
RANGES = {
    "Mach": {0.5: (0.2, 0.70, 0.025), 0.8: (0.3, 1.40, 0.05), 1.25: (0.3, 1.70, 0.05)},
    "p": {0.5: (0.85, 1.20, 0.01), 0.8: (0.45, 1.55, 0.05), 1.25: (0.30, 2.70, 0.10)},
}


def field_plot(F, X, Y, label, mach, fname):
    # colour map + iso-lines (offset by half a step so none falls on the uniform free stream),
    # sonic line in bold on the Mach field
    vmin, vmax, dm = RANGES[label][mach]
    fig, ax = plt.subplots(figsize=(3.4, 2.55))
    pc = ax.pcolormesh(X, Y, F, shading="gouraud", cmap="jet", vmin=vmin, vmax=vmax, rasterized=True)
    ax.contour(X, Y, F, levels=np.arange(vmin + dm / 2, vmax, dm), colors="black", linewidths=0.25)
    if label == "Mach" and F.max() > 1.0:
        ax.contour(X, Y, F, levels=[1.0], colors="black", linewidths=1.0)
    ax.add_patch(Polygon(np.c_[X[0], Y[0]], closed=True, fc="white", ec="black", lw=0.6, zorder=3))
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.75, 0.75)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.06)
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label("$p/p_\\infty$" if label == "p" else "Mach")
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, fname), dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def stagnation_line(X, Y, M):
    # bow-shock stand-off: along the grid line through the leading edge, the outermost sonic crossing
    i_le = int(np.argmin(X[0]))
    x, m = X[:, i_le], M[:, i_le]
    for j in range(len(m) - 1, 0, -1):  # from outside inwards
        if m[j] > 1.0 >= m[j - 1]:
            t = (m[j] - 1.0) / (m[j] - m[j - 1])
            xs = x[j] + t * (x[j - 1] - x[j])
            return X[0, i_le] - xs
    return float("nan")


# ---------------------------------------------------------------- surface Cp

def load_surface(folder):
    rows = []
    for l in open(os.path.join(RESULTS, folder, "s2049x2049.dat")):
        w = l.split()
        if not w or l.startswith("#") or not (w[0][0].isdigit() or w[0][0] in "-."):
            continue
        rows.append([float(v) for v in w[:3]])
    return np.array(rows)  # x/c, y/c, Cp


def cp_star(mach):
    return 2.0 / (GAMMA * mach ** 2) * (((2.0 + (GAMMA - 1.0) * mach ** 2) / (GAMMA + 1.0))
                                         ** (GAMMA / (GAMMA - 1.0)) - 1.0)


def cp_plot(mach, srf, fname):
    # surface Cp of one case, small black open circles (no Cp* line, user 2026-11-18)
    fig, ax = plt.subplots(figsize=(3.4, 1.9))  # shorter than the Mach plots (report space)
    ax.plot(srf[:, 0], srf[:, 2], "o", ms=1.6, mfc="none", mec="black", mew=0.4)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.04)  # the closed trailing edge of the Vassberg-Jameson NACA0012 is at x/c = 1.00893
    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$C_p$")
    ax.grid(True, color="0.85", lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, fname), dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def shocks(s):
    # strongest Cp rise (compression) on the upper (y > 0) and lower surfaces, x/c of the jump
    out = {}
    for side, sel in (("upper", s[:, 1] > 0), ("lower", s[:, 1] < 0)):
        p = s[sel]
        p = p[np.argsort(p[:, 0])]
        p = p[(p[:, 0] > 0.05) & (p[:, 0] < 0.98)]
        dcp = np.diff(p[:, 2]) / np.diff(p[:, 0])
        k = int(np.argmax(dcp))
        out[side] = (0.5 * (p[k, 0] + p[k + 1, 0]), dcp[k])
    return out


# ---------------------------------------------------------------- main

def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.makedirs(FIGDIR, exist_ok=True)
    for tag, folder0, folder125 in (("M08HK", "Mach_08_alpha_0", "Mach_08_alpha_125"),
                                    ("M125HK", "Mach_125_alpha_0", "Mach_125_alpha_125"),
                                    ("M05VA", "Mach_05_alpha_0_vanalbada", "Mach_05_alpha_125_vanalbada"),
                                    ("M08VA", "Mach_08_alpha_0_vanalbada", "Mach_08_alpha_125_vanalbada")):
        print("%% table rows %s" % tag)
        for r in table_rows(folder0, folder125):
            print(r)
    for tag, mach in MACHS:
        fields = {}
        for atag in ("0", "125"):
            folder = "Mach_%s_alpha_%s" % (tag, atag)
            X, Y, M, P = load_field(folder)
            fields[atag] = (X, Y, M)
            field_plot(M, X, Y, "Mach", mach, "mach_M%s_a%s.png" % (tag, atag))
            srf = load_surface(folder)
            cp_plot(mach, srf, "cp_M%s_a%s.png" % (tag, atag))
            near = np.sqrt((X - 0.5) ** 2 + Y ** 2) < 1.0
            print("%% M %g alpha %s: peak Mach near airfoil %.3f, min Cp %.3f, Cp* %.3f"
                  % (mach, "1.25" if atag == "125" else "0", M[near].max(), srf[:, 2].min(), cp_star(mach)))
            if mach >= 0.8 and mach < 1.0:
                sh = shocks(srf)
                print("   upper shock x/c %.3f, lower shock x/c %.3f" % (sh["upper"][0], sh["lower"][0]))
            if mach > 1.0:
                print("   bow-shock stand-off %.4f c" % stagnation_line(X, Y, M))


if __name__ == "__main__":
    main()
