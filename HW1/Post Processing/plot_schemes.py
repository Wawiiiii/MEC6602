import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt

CFL = 1.0 # Change this to select the CFL to plot

results_dir = os.path.join(os.path.dirname(__file__), "..", "results")
pattern = re.compile(r"(.+)_CFL_([0-9.]+)\.dat$")

cfl_text = f"{CFL:.2f}"

for filepath in sorted(glob.glob(os.path.join(results_dir, f"*_CFL_{cfl_text}.dat"))):
    match = pattern.match(os.path.basename(filepath))
    if not match:
        continue
    scheme, _ = match.groups()

    data = np.loadtxt(filepath, comments="#")
    x, u = data[:, 0], data[:, 1]

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

plt.show()
