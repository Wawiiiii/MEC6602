import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt

CFL = 1.25 # Change this to select the CFL to plot

results_dir = os.path.join(os.path.dirname(__file__), "..", "results")
output_dir = os.path.join(os.path.dirname(__file__), "..", "schemes_plots")
os.makedirs(output_dir, exist_ok=True)
pattern = re.compile(r"(.+)_CFL_([0-9.]+)\.dat$")

cfl_text = f"{CFL:.2f}"

# Save one profile figure per scheme at the selected CFL.
for filepath in sorted(glob.glob(os.path.join(results_dir, f"*_CFL_{cfl_text}.dat"))):
    match = pattern.match(os.path.basename(filepath))
    if not match:
        continue
    scheme, _ = match.groups()

    data = np.loadtxt(filepath, comments="#")
    x, u = data[:, 0], data[:, 1]

    # Flag non-finite or excessively large values as divergence.
    diverged = not np.all(np.isfinite(u)) or np.max(np.abs(u[np.isfinite(u)]), initial=0.0) > 1e6

    plt.figure()
    if not diverged:
        plt.plot(x, u)
    else:
        finite_u = u[np.isfinite(u)]
        max_abs = np.max(np.abs(finite_u)) if finite_u.size else float("nan")
        print(f"{scheme}: diverged, max|u| = {max_abs:.3e}")
        plt.text(0.5, 0.5, f"diverged (max|u| = {max_abs:.1e})",
                 ha="center", va="center", transform=plt.gca().transAxes)

    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(f"{scheme}, CFL = {cfl_text}")
    plt.grid(True)

    output = os.path.join(output_dir, f"{scheme}_CFL_{cfl_text}.png")
    plt.savefig(output, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {output}")
