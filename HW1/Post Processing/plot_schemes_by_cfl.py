"""Trace une figure par schéma avec les courbes CFL disponibles superposées.

Utilisation : python "HW1/Post Processing/plot_schemes_by_cfl.py"
Les données sont lues dans HW1/results et les PNG sont créés dans
HW1/schemes_cfl_plots. Les deux dossiers peuvent être changés en ligne de commande.
"""

import argparse
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")  # Save figures without a graphical interface.
import matplotlib.pyplot as plt
import numpy as np


FILENAME = re.compile(r"(.+)_CFL_(\d+\.\d{2})\.dat$")
DIVERGENCE_LIMIT = 1e6


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=hw1_dir / "results")
    parser.add_argument("--output-dir", type=Path, default=hw1_dir / "schemes_cfl_plots")
    args = parser.parse_args()

    if not args.results_dir.is_dir():
        parser.error(f"Dossier de données introuvable : {args.results_dir}")

    # Group result files by scheme for the CFL comparisons.
    files_by_scheme = {}
    for path in args.results_dir.glob("*_CFL_*.dat"):
        match = FILENAME.fullmatch(path.name)
        if match:
            files_by_scheme.setdefault(match.group(1), []).append(
                (float(match.group(2)), path))
    if not files_by_scheme:
        parser.error(f"Aucun fichier *_CFL_*.dat dans {args.results_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scheme, files in sorted(files_by_scheme.items()):
        curves = []
        for cfl, path in sorted(files):
            data = np.loadtxt(path, comments="#", ndmin=2)
            if data.shape[1] < 2:
                raise ValueError(f"Deux colonnes x et u attendues dans {path}")
            x, u = data[:, 0], data[:, 1]
            finite = np.isfinite(u)
            diverged = not finite.all() or np.any(np.abs(u[finite]) > DIVERGENCE_LIMIT)
            curves.append((cfl, x, np.where(finite, u, np.nan), diverged))

        if not curves:
            continue

        # Use stable curves to set a readable vertical scale.
        stable = [u for _, _, u, diverged in curves if not diverged]
        if stable:
            values = np.concatenate(stable)
            lower = min(0.0, float(np.nanmin(values)))
            upper = max(1.0, float(np.nanmax(values)))
            margin = 0.08 * (upper - lower)
            y_limits = (lower - margin, upper + margin)
        else:
            y_limits = (-0.1, 1.1)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.set_xlim(min(float(np.min(x)) for _, x, _, _ in curves),
                    max(float(np.max(x)) for _, x, _, _ in curves))
        ax.set_ylim(*y_limits)
        ax.autoscale(enable=False, axis="y")
        for cfl, x, u, diverged in curves:
            label = f"CFL = {cfl:g}"
            if diverged:
                label += " (divergent)"
                # Hide out-of-range segments so they do not obscure other curves.
                u = np.where((u >= y_limits[0]) & (u <= y_limits[1]), u, np.nan)
            ax.plot(x, u, label=label, linewidth=1.5,
                    zorder=2 if diverged else 3)

        ax.set(xlabel="x", ylabel="u", title=f"{scheme} — CFL comparison")
        ax.grid(True, alpha=0.3)
        ax.legend()
        if any(diverged for _, _, _, diverged in curves):
            ax.text(0.01, 0.01, "Out of scale values are not displayed",
                    transform=ax.transAxes, fontsize=9, va="bottom")
        fig.tight_layout()
        output = args.output_dir / f"{scheme}_cfl_comparison.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
        print(f"Figure enregistrée : {output}")


if __name__ == "__main__":
    main()
