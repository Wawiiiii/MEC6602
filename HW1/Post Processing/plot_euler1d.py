"""Compare 1D Euler CFL sweeps with nozzle solutions and plot residual histories.

Run from any directory. Figures are saved without opening windows by default.
Use --show to also display them.
The inlet conditions and back-pressure ratio must match the Euler mains.
MacCormack records relative state updates; Beam-Warming records the relative
discrete equation residual used for its convergence test.
"""

import argparse
import math
from pathlib import Path
import re

import matplotlib
import numpy as np


GAMMA = 1.4
GAS_CONSTANT = 287.0
INLET_TEMPERATURE = 300.0
INLET_PRESSURE = 101325.0
PROFILE_NAME = re.compile(r"(supersonic|subsonic)_(MacCormack|Implicit)_CFL_(\d+(?:\.\d+)?)\.dat")
RESIDUAL_NAME = re.compile(r"residual_(supersonic|subsonic)_(MacCormack|Implicit)_CFL_(\d+(?:\.\d+)?)\.dat")


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


def read_euler_input(filename):
    values = {}
    for raw_line in filename.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0]
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    try:
        return {
            "n": int(values["n"]),
            "mach_in": float(values["Mach_in"]),
            "back_pressure_ratio": float(values["back_pressure_ratio"]),
            "MacCormack": [float(value) for value in values["CFL_MacCormack"].split()],
            "Implicit": [float(value) for value in values["CFL_Implicit"].split()],
        }
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Invalid Euler input in {filename}: {exc}") from exc


def main():
    hw1_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    display = parser.add_mutually_exclusive_group()
    display.add_argument("--show", dest="no_show", action="store_false",
                         help="Display figures as well as saving PNG files.")
    display.add_argument("--no-show", dest="no_show", action="store_true",
                         help="Save PNG files without opening windows (default).")
    parser.set_defaults(no_show=True)
    parser.add_argument("--mach-in", type=float,
                        help="Override inlet Mach from euler1d_input.txt.")
    parser.add_argument("--back-pressure-ratio", type=float,
                        help="Override subsonic back pressure from euler1d_input.txt.")
    parser.add_argument("--results-dir", type=Path, default=hw1_dir / "Euler_1D_results")
    parser.add_argument("--output-dir", type=Path,
                        help="Figure root (default: Euler_1D_results/Figures).")
    parser.add_argument("--cfl", type=float, nargs="+",
                        help="Override both CFL lists from euler1d_input.txt.")
    parser.add_argument("--compare-cfl", type=float, default=1.0,
                        help="CFL for comparing both schemes (default: 1.0).")
    args = parser.parse_args()
    study = read_euler_input(hw1_dir / "euler1d_input.txt")
    mach_in = args.mach_in if args.mach_in is not None else study["mach_in"]
    back_pressure_ratio = (args.back_pressure_ratio if args.back_pressure_ratio is not None
                           else study["back_pressure_ratio"])
    if mach_in <= 1.0 or back_pressure_ratio <= 0.0:
        parser.error("The reference requires inlet Mach > 1 and back-pressure ratio > 0.")
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    method_cfls = {
        method: sorted(set(args.cfl if args.cfl is not None else study[method]))
        for method in ("MacCormack", "Implicit")
    }
    cfl_values = sorted(set(method_cfls["MacCormack"] + method_cfls["Implicit"]))
    if not cfl_values:
        parser.error("At least one CFL value is required.")
    load_cfl_values = sorted(set(cfl_values + [args.compare_cfl]))
    figure_root = args.output_dir or args.results_dir / "Figures"
    figure_root.mkdir(parents=True, exist_ok=True)

    def selected_cfl(value):
        return next((cfl for cfl in load_cfl_values if abs(cfl - value) < 1e-6), None)

    method_labels = {"MacCormack": "MacCormack", "Implicit": "Beam-Warming"}
    method_colors = {"MacCormack": "tab:blue", "Implicit": "tab:orange"}
    method_markers = {"MacCormack": "o", "Implicit": "s"}
    cfl_colors = {cfl: plt.get_cmap("tab10")(i % 10)
                  for i, cfl in enumerate(load_cfl_values)}
    profiles = {"supersonic": [], "subsonic": []}
    residuals = {"supersonic": [], "subsonic": []}
    meshes = {}

    for filepath in sorted(args.results_dir.rglob("*.dat")):
        match = PROFILE_NAME.fullmatch(filepath.name)
        if match:
            outlet, method, cfl_text = match.groups()
            cfl = selected_cfl(float(cfl_text))
            if cfl is None:
                continue
            data = np.loadtxt(filepath, comments="#", ndmin=2)
            if data.shape[1] != 10 or not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
                raise ValueError(f"Invalid data in {filepath}: expected 10 finite columns and increasing x.")
            profiles[outlet].append((data, method, cfl))
            meshes.setdefault(outlet, data)
            continue

        match = RESIDUAL_NAME.fullmatch(filepath.name)
        if match:
            outlet, method, cfl_text = match.groups()
            cfl = selected_cfl(float(cfl_text))
            if cfl is None:
                continue
            data = np.loadtxt(filepath, comments="#", ndmin=2)
            if (data.shape[1] != 2 or not np.all(np.isfinite(data))
                    or np.any(data[:, 0] < 1) or np.any(np.diff(data[:, 0]) <= 0)
                    or np.any(data[:, 1] < 0)):
                raise ValueError(f"Invalid residual data in {filepath}: expected iteration and nonnegative residual.")
            residuals[outlet].append((data, method, cfl))

    # Older single-CFL results can still be plotted from any results subfolder.
    if selected_cfl(0.5) is not None:
        for outlet in profiles:
            for method in method_labels:
                if any(item[1] == method for item in profiles[outlet] + residuals[outlet]):
                    continue
                candidates = sorted(args.results_dir.rglob(f"{outlet}_{method}.dat"))
                if candidates:
                    filepath = candidates[0]
                    data = np.loadtxt(filepath, comments="#", ndmin=2)
                    if data.shape[1] != 10 or not np.all(np.isfinite(data)):
                        raise ValueError(f"Invalid data in {filepath}")
                    profiles[outlet].append((data, method, 0.5))
                    meshes.setdefault(outlet, data)

    if not any(profiles.values()) and not any(residuals.values()):
        raise SystemExit("No Euler profiles or residual histories for the selected CFL values.")

    inlet_factor = 1.0 + 0.5 * (GAMMA - 1.0) * mach_in**2
    total_temperature = INLET_TEMPERATURE * inlet_factor
    total_pressure = INLET_PRESSURE * inlet_factor ** (GAMMA / (GAMMA - 1.0))
    references = {}
    if "supersonic" in meshes:
        data = meshes["supersonic"]
        critical_area = data[0, 1] / area_mach_ratio(mach_in)
        reference = isentropic_profile(
            data[:, 0], data[:, 1], critical_area, total_pressure,
            total_temperature, "supersonic"
        )
        references["supersonic"] = (reference, "Isentropic")
    if "subsonic" in meshes:
        data = meshes["subsonic"]
        critical_area = data[0, 1] / area_mach_ratio(mach_in)
        reference, shock_x = normal_shock_profile(
            data[:, 0], data[:, 1], critical_area, total_pressure,
            total_temperature, back_pressure_ratio * INLET_PRESSURE
        )
        references["subsonic"] = (reference, "Normal shock")
        print(f"Analytical normal shock at x = {shock_x:.4f} m")

    fields = [
        (1, "area", "Cross-sectional area", r"$A$ (m$^2$)"),
        (2, "density", "Density", r"$\rho$ (kg/m$^3$)"),
        (3, "velocity", "Velocity", r"$u$ (m/s)"),
        (4, "pressure", "Pressure", r"$p$ (Pa)"),
        (5, "temperature", "Temperature", r"$T$ (K)"),
        (6, "mach", "Mach number", r"$M$"),
    ]
    figures = []
    def save_figure(fig, filename):
        fig.tight_layout()
        if filename.startswith("euler1d_area_"):
            output = figure_root / "geometry" / filename
        elif "_comparison_CFL_" in filename:
            outlet = next(name for name in profiles if f"_{name}_" in filename)
            output = (figure_root / "scheme_comparison" / f"CFL_{args.compare_cfl:g}"
                      / outlet / filename)
        else:
            method = next((name for name in method_labels if f"_{name}_cfl" in filename), None)
            outlet = next((name for name in profiles if f"_{name}_" in filename), None)
            if method is not None and outlet is not None:
                folder = "Beam-Warming" if method == "Implicit" else method
                output = figure_root / folder / outlet / filename
            else:
                output = figure_root / "other" / filename
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200)
        if args.no_show:
            plt.close(fig)
        else:
            figures.append(fig)
        print(f"Figure saved: {output}")

    def finish_field_axis(ax, outlet, column, name, ylabel):
        reference, label = references[outlet]
        ax.plot(reference[:, 0], reference[:, column], color="black",
                linewidth=1.0, zorder=4, label=label)
        ax.set(xlabel="x (m)", ylabel=ylabel, title=f"{outlet.capitalize()} outlet")
        ax.grid(True)
        if name == "mach":
            ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        ax.legend(fontsize="small")

    if meshes:
        data = next(iter(meshes.values()))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(data[:, 0], data[:, 1], color="black", linewidth=1.6)
        ax.set(xlabel="x (m)", ylabel=r"$A$ (m$^2$)",
               title="Nozzle cross-sectional area")
        ax.grid(True)
        save_figure(fig, "euler1d_area_comparison.png")

    for column, name, title, ylabel in (fields[1:] if meshes else []):
        # One scheme and one outlet per figure, with all requested CFL values.
        for outlet in profiles:
            for method in method_labels:
                curves = sorted((item for item in profiles[outlet]
                                 if item[1] == method and item[2] in method_cfls[method]),
                                key=lambda item: item[2])
                if not curves:
                    continue
                fig, ax = plt.subplots(figsize=(8, 5))
                for data, _, cfl in curves:
                    ax.plot(data[:, 0], data[:, column], color=cfl_colors[cfl],
                            linewidth=1.6, marker=method_markers[method],
                            markersize=3, markevery=max(1, len(data) // 12),
                            markerfacecolor="white", label=f"CFL={cfl:g}")
                for cfl in method_cfls[method]:
                    if (not any(item[2] == cfl for item in curves)
                            and any(item[1] == method and item[2] == cfl
                                    for item in residuals[outlet])):
                        ax.plot([], [], color=cfl_colors[cfl], linestyle=":",
                                label=f"CFL={cfl:g} (diverged)")
                finish_field_axis(ax, outlet, column, name, ylabel)
                fig.suptitle(f"{method_labels[method]} - {outlet} - {title}")
                save_figure(fig, f"euler1d_{name}_{outlet}_{method}_cfl.png")

        # Two schemes at one CFL, with a separate figure for each outlet.
        comparison = {}
        for outlet in profiles:
            pair = {method: data for data, method, cfl in profiles[outlet]
                    if abs(cfl - args.compare_cfl) < 1e-6}
            if len(pair) == 2:
                comparison[outlet] = pair
        for outlet, pair in comparison.items():
            fig, ax = plt.subplots(figsize=(8, 5))
            for method, data in pair.items():
                ax.plot(data[:, 0], data[:, column],
                        color=method_colors[method], linewidth=1.8, linestyle="--",
                        marker=method_markers[method], markersize=3,
                        markevery=max(1, len(data) // 12),
                        markerfacecolor="white", label=method_labels[method])
            finish_field_axis(ax, outlet, column, name, ylabel)
            fig.suptitle(f"1D Euler - {outlet} - {title} - CFL={args.compare_cfl:g}")
            save_figure(fig,
                        f"euler1d_{name}_{outlet}_comparison_CFL_{args.compare_cfl:g}.png")

    residual_labels = {
        "MacCormack": "Relative state update",
        "Implicit": "Relative discrete equation residual",
    }
    for outlet in residuals:
        for method in method_labels:
            curves = sorted((item for item in residuals[outlet]
                             if item[1] == method and item[2] in method_cfls[method]),
                            key=lambda item: item[2])
            if not curves:
                continue
            fig, ax = plt.subplots(figsize=(8, 5))
            for data, _, cfl in curves:
                diverged = not any(item[1] == method and item[2] == cfl
                                   for item in profiles[outlet])
                ax.semilogy(data[:, 0], np.maximum(data[:, 1], np.finfo(float).tiny),
                            color=cfl_colors[cfl], linewidth=1.4,
                            label=f"CFL={cfl:g}" + (" (diverged)" if diverged else ""))
            ax.set(xlabel="Iteration", ylabel=residual_labels[method],
                   title=f"{method_labels[method]} - {outlet} - residual history")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend()
            save_figure(fig, f"euler1d_residual_{outlet}_{method}_cfl.png")

    residual_comparison = {}
    for outlet in residuals:
        pair = {method: data for data, method, cfl in residuals[outlet]
                if abs(cfl - args.compare_cfl) < 1e-6}
        if len(pair) == 2:
            residual_comparison[outlet] = pair
    for outlet, pair in residual_comparison.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        for method, data in pair.items():
            ax.semilogy(data[:, 0], np.maximum(data[:, 1], np.finfo(float).tiny),
                        color=method_colors[method], linewidth=1.4, linestyle="--",
                        label=f"{method_labels[method]} ({residual_labels[method]})")
        ax.set(xlabel="Iteration", ylabel="Solver convergence metric",
               title=f"{outlet.capitalize()} outlet")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize="small")
        fig.suptitle(f"Residual histories - {outlet} - CFL={args.compare_cfl:g}")
        save_figure(fig, f"euler1d_residual_{outlet}_comparison_CFL_{args.compare_cfl:g}.png")

    if not args.no_show:
        plt.show()
    for fig in figures:
        plt.close(fig)


if __name__ == "__main__":
    main()
