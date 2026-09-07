import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt

results_dir = os.path.join(os.path.dirname(__file__), "..", "results_convergence")
pattern = re.compile(r"space_(.+)\.dat$")

# Nominal spatial order per scheme (see the homework's theory table).
# explicit_forward is unconditionally unstable and has no order to reference.
nominal_order = {
    "explicit_backward": 1,
    "lax": 1,
    "lax_wendroff": 2,
    "leap_frog": 2,
}

for filepath in sorted(glob.glob(os.path.join(results_dir, "space_*.dat"))):
    match = pattern.match(os.path.basename(filepath))
    if not match:
        continue
    scheme = match.group(1)

    data = np.loadtxt(filepath, comments="#")
    dx, error = data[:, 0], data[:, 1]

    finite = np.isfinite(error)
    if not np.any(finite):
        print(f"{scheme}: no finite error values, skipping (diverged at every resolution)")
        continue

    log_dx = np.log10(dx[finite])
    log_error = np.log10(error[finite])

    plt.figure()
    plt.plot(log_dx, log_error, "o-", label=scheme)

    # Reference slope line matching this scheme's own nominal order,
    # anchored at the finest point
    order = nominal_order.get(scheme)
    if order is not None:
        i_min = np.argmin(log_dx)
        ref = log_error[i_min] + order * (log_dx - log_dx[i_min])
        plt.plot(log_dx, ref, "k--", linewidth=1, label=f"order {order} reference")

    plt.xlabel("log10(dx)")
    plt.ylabel("log10(L2 error)")
    plt.title(f"{scheme}: spatial convergence")
    plt.grid(True)
    plt.legend()

plt.show()
