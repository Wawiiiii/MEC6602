"""Compare 1D Euler results with isentropic and normal-shock nozzle solutions.

Run from any directory. Use --no-show to save figures without opening windows.
The inlet conditions and back-pressure ratio must match the Euler mains.
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


GAMMA = 1.4
GAS_CONSTANT = 287.0
INLET_TEMPERATURE = 300.0
INLET_PRESSURE = 101325.0


def area_mach_ratio(mach):
    """Return A/A* for an isentropic quasi-1D flow."""
    factor = 2.0 / (GAMMA + 1.0) * (1.0 + 0.5 * (GAMMA - 1.0) * mach**2)
    exponent = (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0))
    return factor**exponent / mach


def mach_from_area(area_ratio, branch):
    """Invert A/A* by bisection on the requested Mach branch."""
    if area_ratio < 1.0 - 1e-12:
        raise ValueError(f"A/A* must be at least one; got {area_ratio}")
    if abs(area_ratio - 1.0) < 1e-12:
        return 1.0

    if branch == "supersonic":
        lower, upper = 1.0, 2.0
        while area_mach_ratio(upper) < area_ratio:
            upper *= 2.0
    elif branch == "subsonic":
        lower, upper = 1e-12, 1.0
    else:
        raise ValueError(f"Unknown Mach branch: {branch}")

    for _ in range(70):
        middle = 0.5 * (lower + upper)
        if (area_mach_ratio(middle) < area_ratio) == (branch == "supersonic"):
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def isentropic_profile(x, area, critical_area, total_pressure, total_temperature, branch):
    """Return columns x, A, rho, u, p, T, M on one isentropic branch."""
    mach = np.array([mach_from_area(a / critical_area, branch) for a in area])
    factor = 1.0 + 0.5 * (GAMMA - 1.0) * mach**2
    temperature = total_temperature / factor
    pressure = total_pressure / factor ** (GAMMA / (GAMMA - 1.0))
    density = pressure / (GAS_CONSTANT * temperature)
    velocity = mach * np.sqrt(GAMMA * GAS_CONSTANT * temperature)
    return np.column_stack((x, area, density, velocity, pressure, temperature, mach))


def normal_shock_profile(x, area, critical_area, total_pressure, total_temperature,
                         back_pressure):
    """Join two isentropic branches across a normal shock chosen by outlet pressure."""
    exit_area = area[-1]

    def shock_state(shock_x):
        shock_area = float(np.interp(shock_x, x, area))
        mach_before = mach_from_area(shock_area / critical_area, "supersonic")
        factor_before = 1.0 + 0.5 * (GAMMA - 1.0) * mach_before**2
        pressure_before = total_pressure / factor_before ** (GAMMA / (GAMMA - 1.0))

        mach_after_sq = factor_before / (GAMMA * mach_before**2 - 0.5 * (GAMMA - 1.0))
        mach_after = math.sqrt(mach_after_sq)
        pressure_after = pressure_before * (
            1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (mach_before**2 - 1.0))
        factor_after = 1.0 + 0.5 * (GAMMA - 1.0) * mach_after**2
        total_pressure_after = pressure_after * factor_after ** (GAMMA / (GAMMA - 1.0))
        critical_area_after = shock_area / area_mach_ratio(mach_after)
        exit_mach = mach_from_area(exit_area / critical_area_after, "subsonic")
        exit_factor = 1.0 + 0.5 * (GAMMA - 1.0) * exit_mach**2
        exit_pressure = total_pressure_after / exit_factor ** (GAMMA / (GAMMA - 1.0))
        return exit_pressure, shock_area, critical_area_after, total_pressure_after

    lower, upper = float(x[0]), float(x[-1])
    pressure_at_lower = shock_state(lower)[0]
    pressure_at_upper = shock_state(upper)[0]
    if not min(pressure_at_lower, pressure_at_upper) <= back_pressure <= max(
        pressure_at_lower, pressure_at_upper
    ):
        raise ValueError(
            f"Back pressure {back_pressure:g} Pa does not place a normal shock "
            f"inside the nozzle (exit-pressure range: {pressure_at_lower:g} to "
            f"{pressure_at_upper:g} Pa)."
        )

    for _ in range(70):
        middle = 0.5 * (lower + upper)
        pressure_at_middle = shock_state(middle)[0]
        if (pressure_at_middle - back_pressure) * (pressure_at_lower - back_pressure) > 0.0:
            lower, pressure_at_lower = middle, pressure_at_middle
        else:
            upper = middle
    shock_x = 0.5 * (lower + upper)
    _, shock_area, critical_area_after, total_pressure_after = shock_state(shock_x)

    # Duplicate x at the shock so plotted reference curves show a vertical jump.
    before = x < shock_x
    after = x > shock_x
    upstream = isentropic_profile(
        np.r_[x[before], shock_x], np.r_[area[before], shock_area],
        critical_area, total_pressure, total_temperature, "supersonic"
    )
    downstream = isentropic_profile(
        np.r_[shock_x, x[after]], np.r_[shock_area, area[after]],
        critical_area_after, total_pressure_after, total_temperature, "subsonic"
    )
    return np.vstack((upstream, downstream)), shock_x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true", help="Save PNG files without opening windows.")
    parser.add_argument("--mach-in", type=float, default=1.25, help="Inlet Mach number used by the Euler mains.")
    parser.add_argument("--back-pressure-ratio", type=float, default=1.9,
                        help="Subsonic back pressure divided by static inlet pressure.")
    args = parser.parse_args()
    if args.mach_in <= 1.0 or args.back_pressure_ratio <= 0.0:
        parser.error("The reference requires inlet Mach > 1 and back-pressure ratio > 0.")

    results_dir = Path(__file__).resolve().parent.parent / "Euler_1D_results"
    cases = [
        ("supersonic_MacCormack.dat", "supersonic", "MacCormack"),
        ("subsonic_MacCormack.dat", "subsonic", "MacCormack"),
        ("supersonic_Implicit.dat", "supersonic", "Beam-Warming"),
        ("subsonic_Implicit.dat", "subsonic", "Beam-Warming"),
    ]
    method_styles = {
        "MacCormack": {"color": "tab:blue", "linestyle": "-", "marker": "o"},
        "Beam-Warming": {"color": "tab:orange", "linestyle": "--", "marker": "s"},
    }
    fields = [
        (1, "area", "Cross-sectional area", r"$A$ (m$^2$)"),
        (2, "density", "Density", r"$\rho$ (kg/m$^3$)"),
        (3, "velocity", "Velocity", r"$u$ (m/s)"),
        (4, "pressure", "Pressure", r"$p$ (Pa)"),
        (5, "temperature", "Temperature", r"$T$ (K)"),
        (6, "mach", "Mach number", r"$M$"),
    ]
    loaded_cases = {"supersonic": [], "subsonic": []}
    meshes = {}

    for filename, outlet, method in cases:
        filepath = results_dir / filename
        if not filepath.is_file():
            print(f"Missing file: {filepath}")
            continue
        data = np.loadtxt(filepath, comments="#", ndmin=2)
        if data.shape[1] != 10 or not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
            raise ValueError(f"Invalid data in {filepath}: expected 10 finite columns and increasing x.")
        loaded_cases[outlet].append((data, method))
        meshes.setdefault(outlet, data)

    if not any(loaded_cases.values()):
        raise SystemExit("No results to plot. Run main_euler1d_MacCormack or main_euler1d_implicit from HW1.")

    inlet_factor = 1.0 + 0.5 * (GAMMA - 1.0) * args.mach_in**2
    total_temperature = INLET_TEMPERATURE * inlet_factor
    total_pressure = INLET_PRESSURE * inlet_factor ** (GAMMA / (GAMMA - 1.0))
    references = {}
    if "supersonic" in meshes:
        data = meshes["supersonic"]
        critical_area = data[0, 1] / area_mach_ratio(args.mach_in)
        reference = isentropic_profile(
            data[:, 0], data[:, 1], critical_area, total_pressure,
            total_temperature, "supersonic"
        )
        references["supersonic"] = (reference, "Isentropic")
    if "subsonic" in meshes:
        data = meshes["subsonic"]
        critical_area = data[0, 1] / area_mach_ratio(args.mach_in)
        reference, shock_x = normal_shock_profile(
            data[:, 0], data[:, 1], critical_area, total_pressure,
            total_temperature, args.back_pressure_ratio * INLET_PRESSURE
        )
        references["subsonic"] = (reference, "Normal shock")
        print(f"Analytical normal shock at x = {shock_x:.4f} m")

    figures = []
    for column, name, title, ylabel in fields:
        if name == "area":
            # All cases use the same nozzle geometry; draw it only once.
            data = next(iter(meshes.values()))
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(data[:, 0], data[:, column], color="black", linewidth=1.6)
            ax.set_title("Nozzle cross-sectional area")
            ax.set_xlabel("x (m)")
            ax.set_ylabel(ylabel)
            ax.grid(True)
            fig.tight_layout()
        else:
            outlets = [outlet for outlet in ("supersonic", "subsonic") if loaded_cases[outlet]]
            fig, axes = plt.subplots(1, len(outlets), figsize=(12, 5) if len(outlets) == 2 else (7, 5),
                                     squeeze=False)
            for ax, outlet in zip(axes[0], outlets):
                for data, method in loaded_cases[outlet]:
                    ax.plot(data[:, 0], data[:, column], label=method, linewidth=1.8,
                            markersize=3, markevery=20, markerfacecolor="white",
                            **method_styles[method])
                reference, label = references[outlet]
                ax.plot(reference[:, 0], reference[:, column], color="black",
                        linewidth=1.0, zorder=4, label=label)
                ax.set_title(f"{outlet.capitalize()} outlet")
                ax.set_xlabel("x (m)")
                ax.set_ylabel(ylabel)
                ax.grid(True)
                if name == "mach":
                    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
                ax.legend(fontsize="small")
            fig.suptitle(f"1D Euler — {title}")
            fig.tight_layout()
        figures.append(fig)
        output = results_dir / f"euler1d_{name}_comparison.png"
        fig.savefig(output, dpi=200)
        print(f"Figure saved: {output}")
    if not args.no_show:
        plt.show()
    for fig in figures:
        plt.close(fig)


if __name__ == "__main__":
    main()
