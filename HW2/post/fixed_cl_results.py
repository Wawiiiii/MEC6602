"""Fixed-lift results of the report (Section 4): M 0.8, CL 0.4, Hemker-Koren limiter.

Figures (report/figures/, 300 DPI, same style as post/flow_results.py), from the finest fixed-CL
grid (1025 x 1025):
  mach_M08_CL04.png : Mach number around the airfoil
  cp_M08_CL04.png   : surface Cp at CL 0.4
(convergence histories against the iteration: post/convergence_history.py)
Printed: LaTeX rows of the grid table (alpha, CD, CM, iterations, alpha changes, fine-grid time),
Richardson orders of the two finest triplets, numbers quoted in the text.
Run from the HW2 folder:  python3 post/fixed_cl_results.py
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flow_results as fr  # noqa: E402  (plot style, field / Cp plots, Cp*, shock positions)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

FOLDER = "Mach_08_CL_04"         # fixed-CL results
REF = "Mach_08_alpha_125"        # fixed alpha 1.25 deg, same Mach number
GRIDS = [9, 17, 33, 65, 129, 257, 513, 1025]
FINE = 1025                      # finest fixed-CL grid of the report
CL_T = 0.4                       # target lift coefficient

ROW = re.compile(r"^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([-0-9.eE+]+)\s+([-0-9.eE+]+)\s+([-0-9.eE+]+)\s+([0-9.]+)")
CHANGE = re.compile(r"fixed CL: CL ([0-9.]+) -> alpha ([-0-9.]+) -> ([-0-9.]+) deg")
FINAL = re.compile(r"^\s*alpha\s+([-0-9.]+) deg\s+\(fixed CL .*?(\d+) alpha changes")


def summary(folder, g):
    # last summary line: alpha, CL, CD, CM, iterations, seconds, status
    w = [l.split() for l in open(os.path.join(fr.RESULTS, folder, "%dx%d_summary.txt" % (g, g)))
         if l.strip() and not l.startswith("#")][-1]
    return float(w[2]), float(w[3]), float(w[4]), float(w[5]), int(w[6]), float(w[7]), w[8]


def terminal(g):
    # fine-grid iteration rows [(it, drop, CL)], alpha changes [(it, alpha_new)], first alpha,
    # final alpha (10 digits) and number of changes
    rows, changes, a0, final = [], [], None, None
    for l in open(os.path.join(fr.RESULTS, FOLDER, "%dx%d_terminal.txt" % (g, g))):
        m = ROW.match(l)
        if m:
            rows.append((int(m.group(1)), float(m.group(3)), float(m.group(4))))
            continue
        m = CHANGE.search(l)
        if m:
            if a0 is None:
                a0 = float(m.group(2))
            changes.append((rows[-1][0], float(m.group(3))))
            continue
        m = FINAL.match(l)
        if m:
            final = (float(m.group(1)), int(m.group(2)))
    return rows, changes, a0, final


def load_field(g):
    # X, Y, Mach of grid g near the airfoil (same file layout as the 2049 fields of flow_results.py)
    with open(os.path.join(fr.RESULTS, FOLDER, "f%dx%d.dat" % (g, g)), "rb") as fh:
        head = [fh.readline() for _ in range(3)]
        zone = head[2].decode()
        ni = int(zone.split("I=")[1].split(",")[0])
        nj = int(zone.split("J=")[1].split(",")[0])
        data = np.array(fh.read().split(), dtype=float)
    n = ni * nj
    X, Y = data[0:n].reshape(nj, ni), data[n:2 * n].reshape(nj, ni)
    M = data[6 * n:7 * n].reshape(nj, ni)
    r = np.sqrt((X - 0.5) ** 2 + Y ** 2).max(axis=1)
    jmax = int(np.argmax(r > fr.RMAX)) + 1
    return X[:jmax], Y[:jmax], M[:jmax]


def load_surface(folder, g):
    rows = []
    for l in open(os.path.join(fr.RESULTS, folder, "s%dx%d.dat" % (g, g))):
        w = l.split()
        if not w or l.startswith("#") or not (w[0][0].isdigit() or w[0][0] in "-."):
            continue
        rows.append([float(v) for v in w[:3]])
    return np.array(rows)  # x/c, y/c, Cp


def richardson(c, m, f):
    # order p and continuum value F* of one triplet (refinement ratio 2); None if not monotone
    r = (f - m) / (m - c)
    if not 0.0 < r < 1.0:
        return None
    p = math.log2(1.0 / r)
    return p, f + (f - m) / (2.0 ** p - 1.0)


def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.makedirs(fr.FIGDIR, exist_ok=True)

    # table rows: NC, alpha (final, 10 digits from the terminal), CD, CM, iterations, changes, time
    data = {}
    print("% table rows: NC & NC^2 & alpha & CD & CM & t")
    for g in GRIDS:
        a, cl, cd, cm, its, t, st = summary(FOLDER, g)
        _, _, _, final = terminal(g)
        a = final[0]
        data[g - 1] = (a, cl, cd, cm)
        n = g - 1
        print("%s & %s & $%.6f$ & $%+.9f$ & $%+.9f$ & %s \\\\ %% CL %.11f, %d iterations, %d alpha changes, %s" % (
            "{:,}".format(n).replace(",", "{,}"), "{:,}".format(n * n).replace(",", "{,}"),
            a, cd, cm, fr.tsec(t, st), cl, its, final[1], st))

    # Richardson on the two finest triplets
    for name, k in (("alpha", 0), ("CD", 2), ("CM", 3)):
        for c, m, f in ((128, 256, 512), (256, 512, 1024)):
            res = richardson(data[c][k], data[m][k], data[f][k])
            print("%% %-5s %d-%d-%d: %s" % (name, c, m, f, "not monotone" if res is None else
                                             "p %.2f  F* %+.9f" % res))

    # comparison with fixed alpha 1.25 on the same grid
    a, cl, cd, cm = data[FINE - 1]
    ra, rcl, rcd, rcm, _, _, _ = summary(REF, FINE)
    print("%% fixed alpha %.2f on %d: CL %.9f CD %.9f CM %.9f" % (ra, FINE, rcl, rcd, rcm))
    print("%% dCL/dalpha = %.4f /deg, CD +%.1f %%, CM %.6f -> %.6f"
          % ((cl - rcl) / (a - ra), 100 * (cd / rcd - 1), rcm, cm))

    # figures
    rows, changes, a0, final = terminal(FINE)
    print("%% 1025 history: %d iterations, alpha changes at iterations %s, start alpha %.6f"
          % (rows[-1][0], [k for k, _ in changes], a0))
    s_cl, s_ref = load_surface(FOLDER, FINE), load_surface(REF, FINE)
    fr.cp_plot(0.8, s_cl, "cp_M08_CL04.png")  # same style as the Section 3 Cp plots
    for tag, s in (("CL 0.4", s_cl), ("alpha 1.25", s_ref)):
        sh = fr.shocks(s)
        print("%% %s: upper shock x/c %.3f, lower shock x/c %.3f" % (tag, sh["upper"][0], sh["lower"][0]))
    X, Y, M = load_field(FINE)
    fr.field_plot(M, X, Y, "Mach", 0.8, "mach_M08_CL04.png")
    near = np.sqrt((X - 0.5) ** 2 + Y ** 2) < 1.0
    print("%% CL 0.4: peak Mach near the airfoil %.3f" % M[near].max())


if __name__ == "__main__":
    main()
