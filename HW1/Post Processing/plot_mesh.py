"""Plot mesh points from a main_cfl_study result file (column 1 is x)."""

import argparse
from pathlib import Path

import numpy as np


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="A results/*_CFL_*.dat file")
    parser.add_argument("--output", type=Path, default=hw1_dir / "mesh_plot.png")
    parser.add_argument("--no-show", action="store_true", help="Save without opening a window")
    args = parser.parse_args()

    # Use the first available result file when no input is specified.
    if args.input is None:
        files = sorted((hw1_dir / "results").glob("*_CFL_*.dat"))
        if not files:
            parser.error(f"No *_CFL_*.dat files in {hw1_dir / 'results'}")
        args.input = files[0]

    data = np.loadtxt(args.input, comments="#", ndmin=2)
    if data.shape[1] < 2:
        raise ValueError(f"Expected x and u columns in {args.input}")
    x = data[:, 0]
    if not np.all(np.isfinite(x)):
        raise ValueError(f"Non-finite x values in {args.input}")

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Display mesh coordinates along a single horizontal axis.
    fig, ax = plt.subplots()
    ax.plot(x, np.zeros_like(x), "o", markersize=3)
    ax.set(xlabel="x", title=f"Mesh points (n = {len(x)})")
    ax.set_yticks([])
    ax.grid(True, axis="x")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {args.output}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
