# Run every HW2 test program with its grid arguments and print one summary line per test.
# Usage: python tests/run_all_tests.py [build_dir]   (from HW2/, default build dir: build/Release, or build/ on Linux)
import os
import re
import subprocess
import sys
import tempfile

HW2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # HW2 folder


def grid(n):
    return os.path.join(HW2, "grids", f"{n}x{n}.x")


TESTS = [
    ("test_mesh", [grid(65)]),
    ("test_metrics", [grid(65), os.path.join(tempfile.gettempdir(), "hw2_metrics.dat")]),
    ("test_roe", [grid(33), grid(65), grid(129)]),          # nested grids, coarse to fine (order checks)
    ("test_bc", [grid(33), grid(65), grid(129)]),
    ("test_explicit", [grid(33)]),
    ("test_output", [grid(33)]),
    ("test_jacobian", [grid(33)]),
    ("test_linsolve", [grid(33)]),
    ("test_implicit", [grid(33)]),
    ("test_multigrid", [grid(33), grid(65)]),
    ("test_fixedcl", [grid(33)]),
    ("test_muscl_order", []),
]


def main():
    default = os.path.join(HW2, "build", "Release") if os.name == "nt" else os.path.join(HW2, "build")
    bdir = sys.argv[1] if len(sys.argv) > 1 else default
    ext = ".exe" if os.name == "nt" else ""
    bad = 0
    for name, args in TESTS:
        exe = os.path.join(bdir, name + ext)
        r = subprocess.run([exe] + args, capture_output=True, text=True)
        out = r.stdout + r.stderr
        m = re.findall(r"(\d+) (?:checks )?passed, (\d+) failed", out)  # "N passed, M failed" or "N checks passed, M failed"
        summary = f"{m[-1][0]} passed, {m[-1][1]} failed" if m else "(no summary line)"
        ok = r.returncode == 0 and (not m or m[-1][1] == "0")
        bad += 0 if ok else 1
        print(f"{name:16s} rc={r.returncode}  {summary}  {'OK' if ok else 'FAIL'}", flush=True)
        if not ok:
            for line in out.splitlines():
                if "FAIL" in line:
                    print("    " + line)
    print("all tests OK" if bad == 0 else f"{bad} test program(s) FAILED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
