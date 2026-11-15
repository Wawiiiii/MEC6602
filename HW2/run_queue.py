"""Grid-study queue with stall handling (replaces run_matrix.sh, 2026-11-14).

Runs every (case, grid) below one after the other, all cores each. Settings = euler_inputs.txt
(the M 0.5 alpha 0 study), only mesh, mach, alpha, fmg_levels, output paths and the case's
limiter / kappa keys change. A grid whose summary already says "converged" is skipped, so the
queue can be stopped and restarted at any time.

Stall handling (no one may be watching):
- the solver's terminal output is read every 30 s; the run is a STALL if the best residual drop
  has not improved by 10x during the last STALL_WINDOW iterations (reset at each FMG level);
- a stalled, max_iter or diverged run is retried with safer settings (ATTEMPTS below); each
  attempt keeps its own <grid>_input_aN.txt / <grid>_terminal_aN.txt (attempt 1: no suffix);
- if every attempt fails the grid is marked FAILED and the queue goes on.
Events go to results/matrix.log and results/ALERTS.txt. results/MATRIX_REPORT.md is rewritten
after every run and every 10 minutes (post/matrix_report.py).

Run from the HW2 folder:  nohup setsid python3 run_queue.py > /dev/null 2>&1 &
"""
import os
import re
import subprocess
import sys
import time
import datetime

sys.path.insert(0, "post")
import matrix_report  # noqa: E402

TEMPLATE = "euler_inputs.txt"
LOG = "results/matrix.log"
ALERTS = "results/ALERTS.txt"
ALL = ["9x9", "17x17", "33x33", "65x65", "129x129", "257x257", "513x513", "1025x1025", "2049x2049"]

# folder, Mach, alpha, limiter on/off, extra case-file keys, grids, executable
VA = {"limiter": "van_albada"}  # slide 16, kappa^ = 0 (DEVIATIONS A28); default = Hemker-Koren, kappa^ = 1/3
QUEUES = {
    "main": [
        # done 2026-11-14 (skipped): limiter-off studies at M 0.5
        ("Mach_05_alpha_0_nolimiter_kappa0", 0.5, 0.0, "off", {"muscl_kappa": "0"}, ALL[:8], "build/kappa/euler"),  # to 1025
        ("Mach_05_alpha_0_nolimiter", 0.5, 0.0, "off", {}, ALL, "build/euler"),
        ("Mach_05_alpha_125_nolimiter", 0.5, 1.25, "off", {}, ALL, "build/euler"),
        # limiter on, both slide 16 limiters (user's order, 2026-11-15). M 0.5 alpha 0 Hemker-Koren
        # exists already (results/Mach_05_alpha_0); M 0.8 alpha 0 Hemker-Koren 2049 was running.
        ("Mach_08_alpha_0", 0.8, 0.0, "on", {}, ALL, "build/euler"),
        ("Mach_05_alpha_0_vanalbada", 0.5, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_05_alpha_125", 0.5, 1.25, "on", {}, ALL, "build/euler"),
        ("Mach_05_alpha_125_vanalbada", 0.5, 1.25, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_08_alpha_0_vanalbada", 0.8, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_08_alpha_125", 0.8, 1.25, "on", {}, ALL, "build/euler"),
        ("Mach_08_alpha_125_vanalbada", 0.8, 1.25, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {}, ALL, "build/euler"),
        ("Mach_125_alpha_0_vanalbada", 1.25, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_125_alpha_125", 1.25, 1.25, "on", {}, ALL, "build/euler"),
        ("Mach_125_alpha_125_vanalbada", 1.25, 1.25, "on", VA, ALL, "build/kappa/euler"),
    ],
}

# attempt 1 = the user's settings; 2 = same + FMG start from 257 up (the shock is then roughly in
# place before the CFL ramps up; this is how the previous campaign converged M 0.8 on 513);
# 3 and 4 = FMG + CFL cap (B20: residual jumps when the CFL passes ~100-800 on the shock cases).
# Note (2026-11-14, M 0.8 alpha 0 513): a CFL cap alone stalls at drop ~5e-4 (no Newton
# convergence at CFL 50), so the cap comes after the FMG attempt. (key changes, stall window)
ATTEMPTS = [
    ({}, 50),
    ({"fmg": True}, 50),
    ({"cfl_max": "50", "max_iter": "1000", "fmg": True}, 80),
    ({"cfl_max": "10", "max_iter": "3000", "fmg": True}, 300),
]

ITER_LINE = re.compile(r"^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+[-0-9.eE+]+\s+[-0-9.eE+]+\s+[-0-9.eE+]+\s+[0-9.]+")


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg, alert=False):
    with open(LOG, "a") as f:
        f.write(msg + "\n")
    if alert:
        with open(ALERTS, "a") as f:
            f.write("%s  %s\n" % (now(), msg.strip()))


def report():
    try:
        matrix_report.main()
    except Exception as e:  # the report must never stop the queue
        log("    (report failed: %s)" % e)


def last_status(summary):
    if not os.path.exists(summary):
        return None
    rows = [l.split() for l in open(summary) if l.strip() and not l.startswith("#")]
    return rows[-1][8] if rows else None


def write_input(path, out, grid, mach, alpha, lim, extra, changes):
    fmg = 4 if grid == "2049x2049" or (changes.get("fmg") and int(grid.split("x")[0]) >= 257) else 1
    sub = {
        "mesh": "mesh  = grids/%s.x" % grid,
        "mach": "mach  = %g" % mach,
        "alpha": "alpha = %g" % alpha,
        "fmg_levels": "fmg_levels = %d" % fmg,
        "surface": "surface = %s/s%s.dat" % (out, grid),
        "field": "field   = %s/f%s.dat" % (out, grid),
        "summary": "summary = %s/%s_summary.txt" % (out, grid),
    }
    for k in ("cfl_max", "max_iter"):
        if k in changes:
            sub[k] = "%s = %s   # retry setting (run_queue.py)" % (k, changes[k])
    lines = []
    for l in open(TEMPLATE):
        m = re.match(r"^(\w+)\s*=", l)
        lines.append(sub[m.group(1)] + "\n" if m and m.group(1) in sub else l)
    if lim == "off":
        lines.append("\nlimiter_k = 1e6   # no limiting (smooth subsonic case)\n")
    for k, v in extra.items():
        lines.append("%s = %s   # case option (run_queue.py)\n" % (k, v))
    open(path, "w").writelines(lines)


def run_attempt(exe, inp, term, window):
    """Runs one solver; returns 'done' (solver exited by itself) or 'stall' (killed)."""
    with open(term, "w") as tf:
        # stdbuf -oL: line-buffered output, otherwise the file fills in 4 kB blocks and a stall is seen hours late
        p = subprocess.Popen(["stdbuf", "-oL", exe, inp], stdout=tf, stderr=subprocess.STDOUT)
    best, ref, it_ref, last_it, last_rep = float("inf"), float("inf"), 0, -1, time.time()
    pos = 0
    while p.poll() is None:
        time.sleep(30)
        with open(term) as tf:
            tf.seek(pos)
            chunk = tf.read()
            pos = tf.tell()
        for l in chunk.splitlines():
            m = ITER_LINE.match(l)
            if not m:
                continue
            it, drop = int(m.group(1)), float(m.group(3))
            if it <= last_it:  # new FMG level: restart the stall test
                best, ref, it_ref = float("inf"), float("inf"), it
            last_it = it
            best = min(best, drop)
            if best < ref / 10 or ref == float("inf"):
                ref, it_ref = best, it
        if last_it - it_ref > window:
            p.terminate()
            try:
                p.wait(60)
            except subprocess.TimeoutExpired:
                p.kill()
            return "stall at iteration %d (best drop %.2e, no 10x gain in %d iterations)" % (last_it, best, window)
        if time.time() - last_rep > 600:
            report()
            last_rep = time.time()
    return "done"


def other_queue_running():
    # another run_queue.py process (a second queue waits for the first one to finish)
    r = subprocess.run(["pgrep", "-f", "^python3 run_queue.py"], capture_output=True, text=True)
    return any(int(p) != os.getpid() for p in r.stdout.split())


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    name = sys.argv[1] if len(sys.argv) > 1 else "main"
    queue = QUEUES[name]
    while other_queue_running():
        time.sleep(60)
    # let a solver started outside the queue finish first (its summary is then 'converged')
    while subprocess.run(["pgrep", "-f", "build/euler results/"], capture_output=True).returncode == 0:
        time.sleep(30)
    log("=== %s  run_queue.py started (queue %s)" % (now(), name))
    for folder, mach, alpha, lim, extra, grids, exe in queue:
        out = "results/" + folder
        os.makedirs(out, exist_ok=True)
        for g in grids:
            summary = "%s/%s_summary.txt" % (out, g)
            if last_status(summary) == "converged":
                rows = [l for l in open(summary) if not l.startswith("#")]
                log("    skip %s %s (already converged): %s" % (folder, g, rows[-1].strip()))
                continue
            ok = False
            # after a restart, go on after the last attempt that already failed (from the log)
            done = [int(m) for m in re.findall(r"NOT CONVERGED %s %s attempt (\d+)" % (folder, g), open(LOG).read())]
            first = max(done) + 1 if done else 1
            if first > 1:
                log("    resume %s %s at attempt %d (attempts 1-%d failed before)" % (folder, g, first, first - 1))
            for n, (changes, window) in enumerate(ATTEMPTS, start=1):
                if n < first:
                    continue
                if g == "2049x2049" and changes == {"fmg": True}:
                    continue  # 2049 already starts with FMG: attempt 2 would repeat attempt 1
                sfx = "" if n == 1 else "_a%d" % n
                inp, term = "%s/%s_input%s.txt" % (out, g, sfx), "%s/%s_terminal%s.txt" % (out, g, sfx)
                write_input(inp, out, g, mach, alpha, lim, extra, changes)
                log("=== %s  %s  %s  attempt %d  (M %g, alpha %g, limiter %s%s%s)" % (
                    now(), folder, g, n, mach, alpha, lim,
                    "".join(", %s %s" % kv for kv in extra.items()),
                    "" if not changes else ", " + ", ".join("%s %s" % kv for kv in changes.items())))
                res = run_attempt(exe, inp, term, window)
                st = last_status(summary) if res == "done" else None
                if res == "done":
                    rows = [l for l in open(summary) if not l.startswith("#")] if os.path.exists(summary) else []
                    log("    " + (rows[-1].rstrip() if rows else "no summary line (crash? see %s)" % term))
                if res == "done" and st == "converged":
                    ok = True
                    if n > 1:
                        log("    RECOVERED %s %s on attempt %d" % (folder, g, n), alert=True)
                    break
                why = res if res != "done" else "status %s" % st
                log("    NOT CONVERGED %s %s attempt %d: %s" % (folder, g, n, why), alert=True)
            if not ok:
                log("    FAILED %s %s after %d attempts, moving on" % (folder, g, len(ATTEMPTS)), alert=True)
            report()
    log("=== %s  queue done" % now(), alert=True)
    report()


if __name__ == "__main__":
    main()
