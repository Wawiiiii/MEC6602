"""Grid-study report: reads every results/<case>/<grid>_summary.txt and writes results/MATRIX_REPORT.md.

Per case: CL, CD, CM, iterations, time and status on each grid, then the observed order.
- p_bar = log2(f_coarse / f_fine), only for CD of subsonic cases (exact CD = 0).
- Three-grid Richardson (post.pdf slides 26-28, refinement ratio 2):
    r = (f3 - f2) / (f2 - f1), p = log2(r), f* = f1 - (f2 - f1) / (2^p - 1)
  f3 coarse, f1 fine. r <= 0 means the three values are not monotone: no order.
Run from the HW2 folder:  python3 post/matrix_report.py
"""
import glob
import math
import os
import re
import datetime

RESULTS = "results"
OUT = os.path.join(RESULTS, "MATRIX_REPORT.md")
GRIDS = ["9x9", "17x17", "33x33", "65x65", "129x129", "257x257", "513x513", "1025x1025", "2049x2049"]


def read_summary(path):
    # last non-comment line of a summary file (new runs are appended)
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    if not rows:
        return None
    w = rows[-1]
    return dict(grid=w[0], mach=float(w[1]), alpha=float(w[2]), cl=float(w[3]), cd=float(w[4]),
                cm=float(w[5]), it=int(w[6]), time=float(w[7]), status=w[8], drop=w[9],
                limiter=None)


def last_input(case_dir, grid):
    # input file of the last attempt (run_queue.py retries: <grid>_input_a2.txt, _a3.txt) and its number
    for n in (4, 3, 2):
        p = os.path.join(case_dir, "%s_input_a%d.txt" % (grid, n))
        if os.path.exists(p):
            return p, n
    return os.path.join(case_dir, grid + "_input.txt"), 1


def limiter_of(case_dir, grid):
    # limiter on/off (and kappa study) from the saved input file (no file: the original M 0.5 alpha 0 study, limiter on)
    inp, _ = last_input(case_dir, grid)
    if not os.path.exists(inp):
        return "on"
    lim, kap = "on", ""
    for l in open(inp):
        m = re.match(r"\s*limiter_k\s*=\s*([^\s#]+)", l)
        if m:
            lim = "off" if float(m.group(1)) > 1e3 else "on (K %s)" % m.group(1)
        m = re.match(r"\s*muscl_kappa\s*=\s*([^\s#]+)", l)
        if m and m.group(1) != "off":
            kap = ", unlimited kappa = %s" % m.group(1)
        m = re.match(r"\s*limiter\s*=\s*van_albada", l)
        if m:
            kap = ", Van Albada (kappa 0)"
    return lim + kap


def attempt_note(case_dir, grid):
    # "" for a first-attempt run, else the retry settings of the last attempt
    inp, n = last_input(case_dir, grid)
    if n == 1:
        return ""
    keys = {}
    for l in open(inp):
        m = re.match(r"\s*(cfl_max|max_iter|fmg_levels)\s*=\s*([^\s#]+)", l)
        if m:
            keys[m.group(1)] = m.group(2)
    return "attempt %d (%s)" % (n, ", ".join("%s %s" % kv for kv in keys.items()))


def triplet(f3, f2, f1):
    d32, d21 = f3 - f2, f2 - f1
    if d21 == 0 or d32 / d21 <= 0:
        return None, None
    p = math.log2(d32 / d21)
    if abs(2 ** p - 1) < 1e-12:
        return p, None
    return p, f1 - d21 / (2 ** p - 1)


def case_section(case_dir):
    name = os.path.basename(case_dir)
    rows, missing = [], []
    for g in GRIDS:
        s = os.path.join(case_dir, g + "_summary.txt")
        r = read_summary(s) if os.path.exists(s) else None
        if r:
            r["limiter"] = limiter_of(case_dir, g)
            rows.append(r)
        elif os.path.exists(os.path.join(case_dir, g + "_input.txt")):
            missing.append(g)
    if not rows:
        return ""
    M, a = rows[0]["mach"], rows[0]["alpha"]
    lims = sorted(set(r["limiter"] for r in rows))
    out = ["## %s (M %g, alpha %g, limiter %s)\n" % (name, M, a, ", ".join(lims))]
    out.append("| grid | CL | CD | CM | iter | time (s) | status | drop | retry |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        out.append("| %s | %.10e | %.10e | %.10e | %d | %.1f | %s | %s | %s |" % (
            r["grid"], r["cl"], r["cd"], r["cm"], r["it"], r["time"], r["status"], r["drop"],
            attempt_note(case_dir, r["grid"])))
    for g in missing:
        out.append("| %s | (running, or failed: see %s_terminal*.txt) | | | | | | | %s |" % (g, g, attempt_note(case_dir, g)))
    bad = [r["grid"] for r in rows if r["status"] != "converged"]
    if bad:
        out.append("\n**Not converged:** " + ", ".join(bad) + " (check before using these values).")
    out.append("")

    cells = [int(r["grid"].split("x")[0]) - 1 for r in rows]
    if M < 1 and len(rows) >= 2:  # subsonic / transonic: p_bar of CD against 0 only makes sense subsonic
        if M <= 0.5:
            out.append("CD order against the exact value 0, p_bar = log2(CD_coarse / CD_fine):\n")
            out.append("| refinement (cells) | p_bar |")
            out.append("|---|---|")
            for i in range(1, len(rows)):
                c0, c1 = rows[i - 1]["cd"], rows[i]["cd"]
                pb = "%.2f" % math.log2(c0 / c1) if c0 > 0 and c1 > 0 else "n/a"
                out.append("| %d -> %d | %s |" % (cells[i - 1], cells[i], pb))
            out.append("")
    if len(rows) >= 3:
        keys = [k for k in ("cl", "cd", "cm") if max(abs(r[k]) for r in rows) > 1e-9]  # skip CL, CM at alpha 0
        out.append("Three-grid Richardson (r = (f3 - f2)/(f2 - f1), p = log2 r, f* = f1 - (f2 - f1)/(2^p - 1)):\n")
        out.append("| triplet (cells) | " + " | ".join("p %s | %s*" % (k.upper(), k.upper()) for k in keys) + " |")
        out.append("|---|" + "---|---|" * len(keys))
        for i in range(2, len(rows)):
            cellsr = "%d-%d-%d" % (cells[i - 2], cells[i - 1], cells[i])
            cols = []
            for k in keys:
                p, fs = triplet(rows[i - 2][k], rows[i - 1][k], rows[i][k])
                if p is None:
                    cols.append("not monotone | -")
                else:
                    cols.append("%.2f | %s" % (p, "%.6e" % fs if fs is not None else "n/a"))
            out.append("| %s | %s |" % (cellsr, " | ".join(cols)))
        out.append("")
    return "\n".join(out)


def main():
    cases = sorted(d for d in glob.glob(os.path.join(RESULTS, "*")) if os.path.isdir(d))
    log = os.path.join(RESULTS, "matrix.log")
    last = open(log).read().strip().splitlines()[-1] if os.path.exists(log) else ""
    head = ["# Grid-study matrix report (generated, do not edit)\n",
            "Generated %s by `post/matrix_report.py` from the `*_summary.txt` files." %
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "Driver: `run_matrix.sh`, log `results/matrix.log`. Last log line:\n",
            "    " + last + "\n"]
    body = [case_section(c) for c in cases]
    open(OUT, "w").write("\n".join(head + [b for b in body if b]) + "\n")


if __name__ == "__main__":
    main()
