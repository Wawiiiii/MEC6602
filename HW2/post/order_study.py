# Grid convergence study: observed order of accuracy p and extrapolated value (post.pdf slides 26-28,
# Vassberg & Jameson 2010).
#
# Input: the surface Cp files of the solver (one per grid, coarse to fine, each grid twice as fine as
# the previous one). Their comment line "# CL = ... CD = ... CM = ..." holds the force coefficients.
#
# For three grids c (coarse), m (medium), f (fine), with refinement ratio r = 2 (slides 26-27):
#   R = (F_f - F_m) / (F_m - F_c)
#   G(p) = (4^p - 2^p) R - 2^p + 1 = 0,   G'(p) = ln(4) 4^p R - ln(2) 2^p (R + 1)
#   p is found by Newton iteration p <- p - G(p) / G'(p)
#   F_inf = F_f + (F_f - F_m) / (2^p - 1)       (slide 28, extrapolation to h -> 0)
# Subsonic drag: the exact value is 0, so CD itself is the error and p = log2(CD_m / CD_f) (slide 28).
#
# Usage:
#   python order_study.py s65.dat s129.dat s257.dat s513.dat [--cells 64 128 256 512] [--out table.txt]

import argparse
import math
import re
import sys


def read_forces(filename):
    # Return (CL, CD, CM) from the "# CL = ... CD = ... CM = ..." comment line of a surface file.
    pattern = re.compile(r"CL\s*=\s*([-+0-9.eE]+)\s+CD\s*=\s*([-+0-9.eE]+)\s+CM\s*=\s*([-+0-9.eE]+)")
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                match = pattern.search(line)
                if match:
                    return tuple(float(v) for v in match.groups())
    raise ValueError(f"no 'CL = ... CD = ... CM = ...' line in {filename}")


def order_newton(fc, fm, ff, p0=2.0, tol=1e-12, max_iter=100):
    # Observed order p from three grids (r = 2) by Newton iteration on G(p) (slide 27).
    # Returns None when the three values are not monotone (R <= 0: no asymptotic range yet).
    if fm == fc:
        return None
    R = (ff - fm) / (fm - fc)        # ratio of the changes between the grids
    if R <= 0.0:
        return None                  # oscillating convergence: p is not defined
    p = p0                           # starting guess
    for _ in range(max_iter):
        G = (4.0**p - 2.0**p) * R - 2.0**p + 1.0                              # G(p)
        Gp = math.log(4.0) * 4.0**p * R - math.log(2.0) * 2.0**p * (R + 1.0)  # G'(p)
        if Gp == 0.0:
            return None
        step = G / Gp                # Newton step
        p -= step
        if abs(step) < tol:
            return p
    return p


def main():
    parser = argparse.ArgumentParser(description="Observed order of accuracy (post.pdf slides 26-28)")
    parser.add_argument("files", nargs="+", help="surface files, coarse to fine (refinement ratio 2)")
    parser.add_argument("--cells", nargs="*", type=int, help="cells in i of each grid (labels only)")
    parser.add_argument("--out", default="", help="also write the table to this file")
    args = parser.parse_args()

    if len(args.files) < 3:
        sys.exit("need at least 3 grids")
    labels = args.cells if args.cells and len(args.cells) == len(args.files) else [f"grid {k}" for k in range(len(args.files))]
    forces = [read_forces(name) for name in args.files]   # (CL, CD, CM) of every grid

    lines = []                                             # output lines
    lines.append(f"{'grid':>10} {'CL':>14} {'CD':>14} {'CM':>14}")
    for lab, (cl, cd, cm) in zip(labels, forces):
        lines.append(f"{str(lab):>10} {cl:14.8f} {cd:14.8f} {cm:14.8f}")
    lines.append("")

    # Order p and extrapolated value for every set of three consecutive grids.
    lines.append("Observed order p (Newton, slides 26-27) and extrapolated value F_inf (slide 28)")
    lines.append(f"{'grids':>22} {'p CL':>8} {'CL_inf':>12} {'p CD':>8} {'CD_inf':>12} {'p CM':>8} {'CM_inf':>12}"
                 f" {'p CD (CD = error)':>18}")
    for k in range(len(forces) - 2):
        c, m, f = forces[k], forces[k + 1], forces[k + 2]
        row = f"{str(labels[k]) + '-' + str(labels[k + 1]) + '-' + str(labels[k + 2]):>22}"
        for idx in (0, 1, 2):      # CL, CD (three-grid formula, limit not assumed 0) and CM
            p = order_newton(c[idx], m[idx], f[idx])
            if p is None:
                row += f" {'n/a':>8} {'n/a':>12}"
            else:
                finf = f[idx] + (f[idx] - m[idx]) / (2.0**p - 1.0)
                row += f" {p:8.3f} {finf:12.8f}"
        # Drag with the exact value 0 assumed (subsonic, slide 28): p = log2(CD_m / CD_f), two finer grids.
        if m[1] > 0.0 and f[1] > 0.0:
            row += f" {math.log2(m[1] / f[1]):18.3f}"
        else:
            row += f" {'n/a':>18}"
        lines.append(row)

    text = "\n".join(lines)
    print(text)
    if args.out:
        with open(args.out, "w") as out:
            out.write(text + "\n")
        print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
