"""Three-panel CFL comparison figure for the report.

Panels: backward upwind, Lax-Wendroff and hybrid theta=1, each showing the
profile at t_f for several CFL numbers against the exact pulse.

Reads HW1/results/<scheme>_CFL_<cfl>.dat (written by main_cfl_study) and saves
HW1/schemes_plots/wave_cfl_comparison.png by default.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PANELS = [
    ("explicit_backward", "Backward upwind"),
    ("lax_wendroff", "Lax–Wendroff"),
    ("theta_1_0", r"Hybrid $\theta=1$"),
]
CFLS = [0.25, 0.50, 0.75, 1.00]
EXACT_SUPPORT = (2.25, 2.75)  # pulse support at t_f = 3.5 (c = 0.5)
XLIM = (1.9, 3.05)


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=hw1_dir / "results")
    parser.add_argument("--output", type=Path,
                        default=hw1_dir / "schemes_plots" / "wave_cfl_comparison.png")
    args = parser.parse_args()

    fig, axes = plt.subplots(1, len(PANELS), figsize=(11, 3.2), sharey=True)
    xa, xb = EXACT_SUPPORT
    for ax, (scheme, title) in zip(axes, PANELS):
        ax.plot([XLIM[0], xa, xa, xb, xb, XLIM[1]], [0, 0, 1, 1, 0, 0],
                "k--", linewidth=1.0, label="exact")
        for cfl in CFLS:
            path = args.results_dir / f"{scheme}_CFL_{cfl:.2f}.dat"
            data = np.loadtxt(path, comments="#")
            ax.plot(data[:, 0], data[:, 1], linewidth=1.2, label=f"CFL = {cfl:.2f}")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_xlim(*XLIM)
        ax.set_ylim(-0.4, 1.4)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("u")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.08))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {args.output}")


if __name__ == "__main__":
    main()
