"""
Drag convergence plot in the style of Vassberg & Jameson (2010), Fig. 17:
log10|CD - CD*| against log10(1/NC), NC = cells per direction.

Two panels:
  left  : CD* from Richardson extrapolation of the three finest grids (post.pdf slides 26-28)
  right : CD* = 0, the exact value for subsonic inviscid flow
Open circle: CD > CD*, filled circle: CD < CD* (as in V&J).

Usage:
  python post/drag_convergence.py results/Mach_05_alpha_0 --out results/Mach_05_alpha_0/cd_convergence.png
"""
import argparse
import glob
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_summaries(folder):
    """Last line of every *_summary.txt in folder -> sorted lists (NC, CD)."""
    rows = {}
    for path in glob.glob(os.path.join(folder, "*summary*.txt")):
        with open(path) as f:
            lines = [l for l in f if l.strip() and not l.startswith("#")]
        for l in lines:
            tok = l.split()
            ni = int(tok[0].split("x")[0])  # grid nodes per direction
            rows[ni - 1] = float(tok[4])    # NC = cells per direction -> CD (last run of that grid wins)
    nc = sorted(rows)
    return nc, [rows[n] for n in nc]


def richardson(fc, fm, ff):
    """Observed order p and extrapolated value for refinement ratio 2."""
    R = (ff - fm) / (fm - fc)  # ratio of successive differences
    p = math.log2(1.0 / R)     # root of G(p) = (4^p - 2^p) R - 2^p + 1
    return p, ff + (ff - fm) / (2.0**p - 1.0)


def panel(ax, nc, cd, cdstar, title, ref_slope=None):
    x = np.log10(1.0 / np.array(nc, dtype=float))
    err = np.abs(np.array(cd) - cdstar)
    y = np.log10(err)
    ax.plot(x, y, "-", color="0.15", lw=1.2, zorder=1)
    for xi, yi, c, n in zip(x, y, cd, nc):
        filled = c < cdstar
        ax.plot(xi, yi, "o", ms=8, mec="0.15", mfc="0.15" if filled else "white", mew=1.4, zorder=2)
        ax.annotate(str(n), (xi, yi), xytext=(0, -10), textcoords="offset points", rotation=90,
                    ha="center", va="top", fontsize=8)
    xl = np.array([x.min() - 0.3, x.max() + 0.3])
    if ref_slope is None:
        # Least-squares slope on the three finest grids, drawn across the plot (V&J "O(p)" line).
        s, b = np.polyfit(x[-3:], y[-3:], 1)
        ax.plot(xl, s * xl + b, ":", color="0.3", lw=1.2, label=f"O({s:.3f}), 3 finest grids")
    else:
        # Reference slope through the second-finest grid: shows where the finest grid departs from it.
        b = y[-2] - ref_slope * x[-2]
        ax.plot(xl, ref_slope * xl + b, ":", color="0.3", lw=1.2, label=f"O({ref_slope:g}) through NC = {nc[-2]}")
    ax.set_xlim(xl[0], xl[1])
    ax.set_ylim(y.min() - 1.3, y.max() + 1.0)
    ax.plot([], [], "o", mec="0.15", mfc="white", label="CD > CD*")
    ax.plot([], [], "o", mec="0.15", mfc="0.15", label="CD < CD*")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("log10( 1 / NC )")
    ax.set_ylabel("log10( |CD - CD*| )")
    ax.grid(True, ls=":", color="0.75")
    ax.legend(loc="upper left", fontsize=8, frameon=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="folder with the *_summary.txt files of one case")
    ap.add_argument("--out", default="cd_convergence.png")
    ap.add_argument("--title", default="NACA0012, M = 0.50, alpha = 0.0 deg, BDF1 + GMRES")
    a = ap.parse_args()

    nc, cd = read_summaries(a.folder)
    p, cdstar = richardson(cd[-3], cd[-2], cd[-1])
    print(f"grids (NC): {nc}")
    print(f"Richardson {nc[-3]}-{nc[-2]}-{nc[-1]}: p = {p:.3f}, CD* = {cdstar:.4e}")

    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    panel(axs[0], nc, cd, cdstar, f"CD* = {cdstar:+.3e}  (Richardson, p = {p:.2f})")
    p_pair = math.log2(cd[-3] / cd[-2]) # p_bar of the two grids before the finest
    panel(axs[1], nc, cd, 0.0, "CD* = 0  (exact)", ref_slope=round(p_pair, 2))
    axs[1].yaxis.set_tick_params(labelleft=True)
    fig.suptitle(a.title, fontsize=12)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
