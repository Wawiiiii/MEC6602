"""Plot results from main_euler1d_MacCormack, run from HW1.

Save each field to a separate PNG. Use --no-show to skip opening windows.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true", help="Save PNG files without opening windows.")
    args = parser.parse_args()

    results_dir = Path(__file__).resolve().parent.parent / "Euler_1D_results"
    cases = [
        ("supersonic_MacCormack.dat", "Supersonic outlet", "-"),
        ("subsonic_MacCormack.dat", "Subsonic outlet", "--"),
    ]
    fields = [
        (1, "area", "Cross-sectional area", r"$A$ (m$^2$)"),
        (2, "density", "Density", r"$\rho$ (kg/m$^3$)"),
        (3, "velocity", "Velocity", r"$u$ (m/s)"),
        (4, "pressure", "Pressure", r"$p$ (Pa)"),
        (5, "temperature", "Temperature", r"$T$ (K)"),
        (6, "mach", "Mach number", r"$M$"),
    ]
    loaded_cases = []

    for filename, label, linestyle in cases:
        filepath = results_dir / filename
        if not filepath.is_file():
            print(f"Missing file: {filepath}")
            continue
        data = np.loadtxt(filepath, comments="#", ndmin=2)
        if data.shape[1] != 10 or not np.all(np.isfinite(data)):
            raise ValueError(f"Invalid data in {filepath}: expected 10 columns of finite values.")
        loaded_cases.append((data, label, linestyle))

    if not loaded_cases:
        raise SystemExit("No results to plot. Run main_euler1d_MacCormack from HW1.")

    figures = []
    for column, name, title, ylabel in fields:
        fig, ax = plt.subplots(figsize=(8, 5))
        figures.append(fig)
        for data, label, linestyle in loaded_cases:
            ax.plot(data[:, 0], data[:, column], linestyle, label=label)
        ax.set_title(f"1D Euler — MacCormack scheme\n{title}")
        ax.set_xlabel("x (m)")
        ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.legend()
        if name == "mach":
            ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        fig.tight_layout()
        output = results_dir / f"euler1d_{name}_MacCormack.png"
        fig.savefig(output, dpi=200)
        print(f"Figure saved: {output}")
    if not args.no_show:
        plt.show()
    for fig in figures:
        plt.close(fig)


if __name__ == "__main__":
    main()
