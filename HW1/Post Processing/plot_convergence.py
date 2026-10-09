"""Plot the space_*.dat files produced by main_convergence.

Run from any directory. By default, read HW1/results_convergence and save PNGs
in HW1/convergence_plots.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Nominal spatial orders for schemes with a useful reference slope.
NOMINAL_ORDER = {
    "explicit_backward": 1,
    "lax": 1,
    "lax_wendroff": 2,
    "leap_frog": 2,
}


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=hw1_dir / "results_convergence")
    parser.add_argument("--output-dir", type=Path, default=hw1_dir / "convergence_plots")
    args = parser.parse_args()

    files = sorted(args.results_dir.glob("space_*.dat"))
    if not files:
        parser.error(f"No space_*.dat files in {args.results_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filepath in files:
        scheme = filepath.stem[len("space_"):]
        data = np.loadtxt(filepath, comments="#", ndmin=2)
        if data.shape[1] != 3:
            raise ValueError(f"Expected dx, error, order columns in {filepath}")

        dx, error = data[:, 0], data[:, 1]
        # Keep only positive finite values for the logarithmic comparison.
        valid = np.isfinite(dx) & (dx > 0) & np.isfinite(error) & (error > 0)
        if not np.any(valid):
            print(f"{scheme}: no positive finite error values, skipping")
            continue

        log_dx = np.log10(dx[valid])
        log_error = np.log10(error[valid])

        fig, ax = plt.subplots()
        ax.plot(log_dx, log_error, "o-", label=scheme)

        order = NOMINAL_ORDER.get(scheme)
        if order is not None:
            # Anchor the nominal-order reference at the finest valid mesh.
            finest = np.argmin(log_dx)
            reference = log_error[finest] + order * (log_dx - log_dx[finest])
            ax.plot(log_dx, reference, "k--", linewidth=1,
                    label=f"order {order} reference")

        ax.set(xlabel="log10(dx)", ylabel="log10(L2 error)",
               title=f"{scheme}: spatial convergence")
        ax.grid(True)
        ax.legend()
        output = args.output_dir / f"{scheme}.png"
        fig.savefig(output, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Figure saved: {output}")


if __name__ == "__main__":
    main()
