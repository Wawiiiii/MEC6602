"""All-schemes profile figure at a single CFL for the report appendix.

Reads HW1/results/<scheme>_CFL_<cfl>.dat (written by main_cfl_study) and saves
HW1/schemes_plots/wave_scheme_comparison.png by default. Schemes that diverge at
every CFL (forward upwind, hybrid theta=0) are not plotted.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PANELS = [
    ("explicit_backward", "Backward upwind"),
    ("lax", "Lax"),
    ("lax_wendroff", "Lax–Wendroff"),
    ("leap_frog", "Leapfrog"),
    ("theta_0_5", r"Hybrid $\theta=\frac{1}{2}$"),
    ("theta_1_0", r"Hybrid $\theta=1$"),
    ("scheme_2space_4time", r"$(2,4)$ scheme"),
    ("scheme_4space_2time", r"$(4,2)$ scheme"),
]
EXACT_SUPPORT = (2.25, 2.75)  # pulse support at t_f = 3.5 (c = 0.5)
XLIM = (1.9, 3.05)


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=hw1_dir / "results")
    parser.add_argument("--cfl", type=float, default=0.50)
    parser.add_argument("--output", type=Path,
                        default=hw1_dir / "schemes_plots" / "wave_scheme_comparison.png")
    args = parser.parse_args()

    fig, axes = plt.subplots(2, 4, figsize=(11, 5), sharex=True, sharey=True)
    xa, xb = EXACT_SUPPORT
    for ax, (scheme, title) in zip(axes.flat, PANELS):
        data = np.loadtxt(args.results_dir / f"{scheme}_CFL_{args.cfl:.2f}.dat", comments="#")
        ax.plot([XLIM[0], xa, xa, xb, xb, XLIM[1]], [0, 0, 1, 1, 0, 0],
                "k--", linewidth=1.0, label="exact")
        ax.plot(data[:, 0], data[:, 1], linewidth=1.2, label=f"CFL = {args.cfl:.2f}")
        ax.set_title(title, fontsize=10)
        ax.set_xlim(*XLIM)
        ax.set_ylim(-0.4, 1.4)
        ax.grid(True, alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("x")
    for ax in axes[:, 0]:
        ax.set_ylabel("u")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.02))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {args.output}")


if __name__ == "__main__":
    main()
