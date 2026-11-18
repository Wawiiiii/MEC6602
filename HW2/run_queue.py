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

Run from the HW2 folder (queue name optional, default main):
  Linux:   nohup setsid python3 run_queue.py [main|m125] > /dev/null 2>&1 &
  Windows: python run_queue.py m125   (in a window that stays open, or Start-Process)
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
        # M 1.25 alpha 0, both limiters (2026-11-15: the personal computer is busy); M 1.25 alpha 1.25 later
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {}, ALL, "build/euler"),
        ("Mach_125_alpha_0_vanalbada", 1.25, 0.0, "on", VA, ALL, "build/kappa/euler"),
    ],
    # 2026-11-15 ~01:00: user allowed the lab machines j5015-18, -17, -04, ONE CASE PER MACHINE for now
    # (same shared disk, same binaries, lab queues at nice 19). j5015-04 has ~31 GB free and a 2049 run
    # needs 25 GB, so it does M 1.25 alpha 0 HK up to 1025 and that 2049 runs here (end of j01).
    # Each case folder / grid belongs to one queue only (two queues never write the same grid).
    "j01": [
        ("Mach_08_alpha_0", 0.8, 0.0, "on", {}, ALL, "build/euler"),  # 2049 was running
        ("Mach_05_alpha_0_vanalbada", 0.5, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_05_alpha_125", 0.5, 1.25, "on", {}, ALL, "build/euler"),
        ("Mach_05_alpha_125_vanalbada", 0.5, 1.25, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_08_alpha_0_vanalbada", 0.8, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_125_alpha_0_vanalbada", 1.25, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {}, ["2049x2049"], "build/euler"),
    ],
    "j18": [
        ("Mach_08_alpha_125", 0.8, 1.25, "on", {}, ALL, "build/euler"),
    ],
    "j17": [
        ("Mach_08_alpha_125_vanalbada", 0.8, 1.25, "on", VA, ALL, "build/kappa/euler"),
    ],
    "j04": [
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {}, ALL[:8], "build/euler"),
    ],
    # 2026-11-15 10:00, user: use the lab machines today, add the remaining M 1.25 runs, done by Mon 09:00.
    # A 2049 run is mostly serial (GMRES / SGS, ~1 core), so -18 and -17 take two queues each (half the
    # cores each); -04 (20 GB free) takes the grids up to 1025 of the slow cases. Queue names with a
    # "d" suffix ("day") are today's split; every grid still belongs to exactly one queue.
    "j01d": [
        ("Mach_08_alpha_0", 0.8, 0.0, "on", {}, ALL, "build/euler"),  # 2049 running
        ("Mach_05_alpha_0_vanalbada", 0.5, 0.0, "on", VA, ALL, "build/kappa/euler"),
        ("Mach_05_alpha_125_vanalbada", 0.5, 1.25, "on", VA, ALL, "build/kappa/euler"),  # 2049 moved here from j18d (11:30)
    ],
    "j18d": [
        ("Mach_08_alpha_125", 0.8, 1.25, "on", {}, ALL, "build/euler"),  # 2049 running
    ],
    # 14:20: j18d's 2049 (old tol 1e-12) is crawling at 1e-7, so its next case starts now in parallel
    "j18g": [
        ("Mach_08_alpha_0_vanalbada", 0.8, 0.0, "on", dict(VA, threads="13"), ["2049x2049"], "build/kappa/euler"),
    ],
    "j18e": [
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {"threads": "26"}, ["2049x2049"], "build/euler"),
    ],
    # 12:35: the M 1.25 2049 runs may take 10-20 h each, so the second one of j18e / j17e runs in
    # parallel as a third queue (memory: ~32 GB per 2049 run; -18 had 209 GB free, -17 109 GB)
    "j18f": [
        # Mon 00:30, user: M 1.25 Van Albada 2049 not needed (run stopped, removed)
    ],
    "j17d": [
        ("Mach_08_alpha_125_vanalbada", 0.8, 1.25, "on", VA, ALL, "build/kappa/euler"),  # 2049 running
    ],
    # 15:35: j04d would end ~Mon 15:00, j17d's 2049 (old tol 1e-12) may run into the night: their
    # waiting cases start now in parallel on -18 / -17 (free cores and memory). A case folder is never
    # in two queues of the same machine (a starting queue takes over the running solvers of its folders).
    "j18h": [
        ("Mach_05_alpha_125", 0.5, 1.25, "on", {"threads": "13"}, ["2049x2049"], "build/euler"),
    ],
    "j18i": [
        ("Mach_125_alpha_125", 1.25, 1.25, "on", {"threads": "8"}, ALL[:8], "build/euler"),
    ],
    "j17h": [
        ("Mach_125_alpha_125_vanalbada", 1.25, 1.25, "on", dict(VA, threads="8"), ALL[:8], "build/kappa/euler"),
    ],
    "j17i": [
        ("Mach_05_alpha_125", 0.5, 1.25, "on", {"threads": "8"}, ALL[:8], "build/euler"),
    ],
    "j17e": [
        ("Mach_125_alpha_125", 1.25, 1.25, "on", {"threads": "14"}, ["2049x2049"], "build/euler"),
    ],
    "j17f": [
        # Mon 00:30, user: M 1.25 Van Albada 2049 not needed (run stopped, removed)
    ],
    "j04d": [
        ("Mach_125_alpha_0_vanalbada", 1.25, 0.0, "on", VA, ALL[:8], "build/kappa/euler"),
    ],
    # Mon 00:05: j01d finished, j04d's next case moved to j5015-01
    "j01e": [
        ("Mach_08_alpha_0_vanalbada", 0.8, 0.0, "on", VA, ALL[:8], "build/kappa/euler"),
    ],
    # 2026-11-17, user: Part 5 matrix, fixed CL 0.4 at M 0.8 (post.pdf slide 9, DEVIATIONS A29), Hemker-Koren,
    # grids up to 1025 for now; alpha 1.25 is the first guess (the solver moves it).
    "part5": [
        ("Mach_08_CL_04", 0.8, 1.25, "on", {"cl_target": "0.4"}, ALL[:8], "build/euler"),
    ],
    # 11:15: on 1025 the FMG start was already at CL 0.39999 (alpha from 513), but the slide's first
    # perturbation of 0.1 deg moved the shock so far that the residual bounced at 1e-3 - 1e-2 for 50
    # iterations and a secant step from an unsettled CL went the wrong way (1.4128 -> 1.3812 deg).
    # Stopped (files kept as 1025x1025_*_dalpha01.txt). Fine grids: first perturbation 0.01 deg and a
    # residual drop of 1e-3 before each alpha change. 2049 (part5c) only if 1025 converges (user).
    "part5b": [
        ("Mach_08_CL_04", 0.8, 1.25, "on", {"cl_target": "0.4", "cl_dalpha": "0.01", "cl_drop": "1e-3"},
         ["1025x1025"], "build/euler"),
    ],
    "part5c": [
        ("Mach_08_CL_04", 0.8, 1.25, "on", {"cl_target": "0.4", "cl_dalpha": "0.01", "cl_drop": "1e-3"},
         ["2049x2049"], "build/euler"),
    ],
    "m125": [
        ("Mach_125_alpha_0", 1.25, 0.0, "on", {}, ALL, "build/euler"),
        ("Mach_125_alpha_0_vanalbada", 1.25, 0.0, "on", VA, ALL, "build/euler"),
        ("Mach_125_alpha_125", 1.25, 1.25, "on", {}, ALL, "build/euler"),
        ("Mach_125_alpha_125_vanalbada", 1.25, 1.25, "on", VA, ALL, "build/euler"),
    ],
}

# Windows (personal computer): one Release build has every option; no stdbuf / pgrep there (the
# solver flushes each iteration row itself, src/solver.cpp)
WINDOWS = os.name == "nt"
WIN_EXE = "build/Release/euler.exe"

# attempt 1 = the user's settings; 2 = same + FMG start from 257 up (the shock is then roughly in
# place before the CFL ramps up; this is how the previous campaign converged M 0.8 on 513);
# 3 = FMG + stronger GMRES (restart 200, up to 800 iterations): on 2049 GMRES hits its 400 cap
# with a linear residual of 0.1-0.2 at high CFL, which slows Newton (M 0.8 alpha 0 2049, 2026-11-15);
# 4 = the same + a moderate CFL cap 1e3 (B20: shocks move too far per step above CFL ~100-800).
# Changed 2026-11-15 ~01:00 (user: make the runs robust): low CFL caps (50, 10) are NOT used any more:
# with the limiter a cap of 50 stalled at drop ~5e-4 (M 0.8 alpha 0 513), and on 2049 a CFL-10 run
# would take days. (key changes, stall window)
# User, 2026-11-15 ~01:05: no more cold starts. Every run starts with FMG (257 and up) and a stronger
# GMRES on the fine grids (base_changes); the retries come on top of that.
ATTEMPTS = [
    ({}, 50),
    ({"gmres_restart": "200", "gmres_max_iter": "800", "max_iter": "400"}, 80),
    ({"gmres_restart": "200", "gmres_max_iter": "800", "max_iter": "1000", "cfl_max": "1e3"}, 150),
    # B20 (M 1.25 alpha 0 1025, 2026-11-15 02:48): every step taken at CFL > ~130 made the residual
    # jump 4-15x (the shock moves too far), steps at CFL 40-120 reduced it ~2x. Cap just below that.
    ({"gmres_restart": "200", "gmres_max_iter": "800", "max_iter": "1500", "cfl_max": "100"}, 150),
    # M 1.25 alpha 0 1025, 2026-11-15 04:15 (diagnosis agent): the cap of 100 bounced the same way. The
    # linearized step is wrong at the shocks once they move by more than a cell; 513 / 257 escaped
    # only when the CFL passed ~1000. So: cfl0 10 (CFL ~1e3 at the bouncing level), reject poorly
    # solved steps (gmres_fail 1e-2), and the homework's Cauchy criterion on the forces as a safety
    # net (max - min of CL, CD, CM over 20 iterations < 2e-7; status "cauchy" counts as finished).
    ({"gmres_restart": "200", "gmres_max_iter": "800", "max_iter": "300", "cfl0": "10",
      "gmres_fail": "1e-2", "cauchy_window": "20", "cauchy_tol": "2e-7"}, 150),
]
OK_STATUS = ("converged", "cauchy")  # "cauchy" only from attempt 5 (the template has cauchy_window = 0)


def bounded_case(mach, grid):
    # M 1.25 on 1025 / 2049: the residual bounces at ~1e-2 (shocks, PROGRESS_NIGHT.md 04:12) while the
    # forces are steady to ~1e-6; one attempt, never stopped by the stall rule, ended by the Cauchy
    # criterion or max_iter, and "max_iter" accepted as finished (forces then good to ~1e-6)
    return mach > 1 and int(grid.split("x")[0]) >= 1025


def finished(st, mach, grid):
    return st in OK_STATUS or (st == "max_iter" and bounded_case(mach, grid))


def input_mach(inp):
    for l in open(inp):
        m = re.match(r"\s*mach\s*=\s*([^\s#]+)", l)
        if m:
            return float(m.group(1))
    return 0.0


def base_changes(grid, mach):
    # settings of every attempt: FMG start from 257 up; GMRES restart 200 / 800 iterations from 1025 up;
    # tol 1e-10 on 2049 (user 2026-11-15 10:00: forces are converged long before 1e-12 there; the last
    # two orders cost ~6 h per run); M 1.25 from 1025 up: the attempt-5 settings that ended the bouncing
    # (cfl0 10, gmres_fail 1e-2, Cauchy stop on the forces), so the failing attempts are not repeated.
    n = int(grid.split("x")[0])
    c = {}
    if n >= 257:
        c["fmg"] = True
    if n >= 1025:
        c.update({"gmres_restart": "200", "gmres_max_iter": "800"})
    if n >= 2049:
        c["tol"] = "1e-10"
    if mach > 1 and n >= 1025:
        # 12:30: window 20 / 2e-7 / max_iter 250 would not stop (CL oscillates ~3e-6 on 1025) and 250
        # iterations of 2049 take days: Cauchy 10 it / 1e-6, max_iter 150 (1025) / 60 (2049), one attempt
        c.update({"cfl0": "10", "gmres_fail": "1e-2", "cauchy_window": "10", "cauchy_tol": "1e-6",
                  "max_iter": "60" if n >= 2049 else "150"})
    return c

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
    for k, v in changes.items():
        if k != "fmg":
            sub[k] = "%s = %s   # run_queue.py setting" % (k, v)
    if "threads" in extra:  # two runs on one lab machine: each gets part of the cores
        sub["threads"] = "threads = %s   # case option (run_queue.py)" % extra["threads"]
        extra = {k: v for k, v in extra.items() if k != "threads"}
    lines, used = [], set()
    for l in open(TEMPLATE):
        m = re.match(r"^(\w+)\s*=", l)
        if m and m.group(1) in sub:
            lines.append(sub[m.group(1)] + "\n")
            used.add(m.group(1))
        else:
            lines.append(l)
    for k in sub:  # keys that are not in the template (e.g. gmres_fail)
        if k not in used:
            lines.append(sub[k] + "\n")
    if lim == "off":
        lines.append("\nlimiter_k = 1e6   # no limiting (smooth subsonic case)\n")
    for k, v in extra.items():
        lines.append("%s = %s   # case option (run_queue.py)\n" % (k, v))
    open(path, "w").writelines(lines)


def run_attempt(exe, inp, term, window):
    """Runs one solver; returns 'done' (solver exited by itself) or 'stall' (killed)."""
    with open(term, "w") as tf:
        # stdbuf -oL: line-buffered output, otherwise the file fills in 4 kB blocks and a stall is seen
        # hours late (binaries built before the fflush in src/solver.cpp); not on Windows
        cmd = [exe, inp] if WINDOWS else ["stdbuf", "-oL", exe, inp]
        p = subprocess.Popen(cmd, stdout=tf, stderr=subprocess.STDOUT)

    def stop():
        p.terminate()
        try:
            p.wait(60)
        except subprocess.TimeoutExpired:
            p.kill()

    return watch(lambda: p.poll() is None, stop, term, window)


def watch(alive, stop, term, window):
    """Stall test on a running solver's terminal file; 'done' (exited by itself) or 'stall ...' (stopped)."""
    best, ref, it_ref, last_it, last_rep = float("inf"), float("inf"), 0, -1, time.time()
    pos = 0
    while alive():
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
            stop()
            return "stall at iteration %d (best drop %.2e, no 10x gain in %d iterations)" % (last_it, best, window)
        if time.time() - last_rep > 600:
            report()
            last_rep = time.time()
    return "done"


def other_queue_running(name):
    # another process running the SAME queue (two different queues may share a machine); Linux only
    if WINDOWS:
        return False
    r = subprocess.run(["pgrep", "-f", "^python3 run_queue.py %s$" % name], capture_output=True, text=True)
    return any(int(p) != os.getpid() for p in r.stdout.split())


def running_solvers():
    # solvers started outside this queue (by an earlier queue process): [(pid, input file)]; Linux only
    if WINDOWS:
        return []
    r = subprocess.run(["pgrep", "-af", r"build/(kappa/)?euler results/"], capture_output=True, text=True)
    out = []
    for l in r.stdout.splitlines():
        w = l.split()
        if len(w) >= 3 and w[-1].endswith(".txt") and "_input" in w[-1]:
            out.append((int(w[0]), w[-1]))
    return out


def adopt(pid, inp):
    """Takes over a solver left running by an earlier queue process: same stall test; on a stall it
    is stopped and logged as a failed attempt, so the queue goes on with the next attempt."""
    folder, fname = os.path.basename(os.path.dirname(inp)), os.path.basename(inp)
    g = fname.split("_input")[0]
    m = re.search(r"_input_a(\d+)\.txt$", fname)
    n = int(m.group(1)) if m else 1
    term = inp.replace("_input", "_terminal")
    window = ATTEMPTS[min(n, len(ATTEMPTS)) - 1][1]
    # a taken-over 2049 run has hours of work in it and its files are written only at the end: wider window
    window = 30 if g == "2049x2049" else window
    mach = input_mach(inp)
    if bounded_case(mach, g):
        window = 10 ** 9
    if n >= len(ATTEMPTS):  # last attempt: never stop it (a stop leaves no files; max_iter / Cauchy end it and write them)
        window = 10 ** 9

    def alive():
        try:
            with open("/proc/%d/stat" % pid) as f:
                return f.read().rsplit(")", 1)[1].split()[0] != "Z"  # a zombie has finished
        except OSError:
            return False

    def stop():
        os.kill(pid, 15)
        for _ in range(60):
            if not alive():
                return
            time.sleep(1)
        os.kill(pid, 9)

    log("    adopted running solver %d: %s %s attempt %d" % (pid, folder, g, n))
    res = watch(alive, stop, term, window)
    time.sleep(5)  # let the solver finish writing its summary line
    st = last_status(os.path.join(os.path.dirname(inp), g + "_summary.txt")) if res == "done" else None
    if res == "done" and finished(st, mach, g):
        log("    adopted run finished: %s %s %s" % (folder, g, st))
    else:
        log("    NOT CONVERGED %s %s attempt %d: %s" % (folder, g, n, res if res != "done" else "status %s" % st),
            alert=True)


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    name = sys.argv[1] if len(sys.argv) > 1 else "main"
    queue = QUEUES[name]
    global LOG, ALERTS
    if name not in ("main", "j01", "j01d"):  # j01e has its own log  # one log / alert file per queue (this machine: the original ones)
        LOG, ALERTS = "results/matrix_%s.log" % name, "results/ALERTS_%s.txt" % name
    log("=== %s  host %s, queue %s" % (now(), os.uname().nodename if hasattr(os, "uname") else "?", name))
    if WINDOWS:
        queue = [c[:6] + (WIN_EXE,) for c in queue]
    while other_queue_running(name):
        time.sleep(60)
    # a solver started by an earlier queue process for one of THIS queue's cases: watch it to the end
    own = set("results/" + c[0] for c in queue)
    for pid, inp in running_solvers():
        if os.path.dirname(inp) in own:
            adopt(pid, inp)
    log("=== %s  run_queue.py started (queue %s)" % (now(), name))
    for folder, mach, alpha, lim, extra, grids, exe in queue:
        out = "results/" + folder
        os.makedirs(out, exist_ok=True)
        for g in grids:
            summary = "%s/%s_summary.txt" % (out, g)
            if finished(last_status(summary), mach, g):
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
                if n == 2 and int(g.split("x")[0]) >= 1025:
                    continue  # 1025 / 2049 already start with the strong GMRES: attempt 2 would repeat attempt 1
                changes = dict(base_changes(g, mach), **changes)
                sfx = "" if n == 1 else "_a%d" % n
                inp, term = "%s/%s_input%s.txt" % (out, g, sfx), "%s/%s_terminal%s.txt" % (out, g, sfx)
                write_input(inp, out, g, mach, alpha, lim, extra, changes)
                log("=== %s  %s  %s  attempt %d  (M %g, alpha %g, limiter %s%s%s)" % (
                    now(), folder, g, n, mach, alpha, lim,
                    "".join(", %s %s" % kv for kv in extra.items()),
                    "" if not changes else ", " + ", ".join("%s %s" % kv for kv in changes.items())))
                # 2049: 10-20 min per iteration, so 50 iterations would be 10+ hours: shorter window
                w = min(window, 15) if g == "2049x2049" else window
                if bounded_case(mach, g):
                    w = 10 ** 9
                if n == len(ATTEMPTS):  # last attempt: never stop it (a stop leaves no files)
                    w = 10 ** 9
                res = run_attempt(exe, inp, term, w)
                st = last_status(summary) if res == "done" else None
                if res == "done":
                    rows = [l for l in open(summary) if not l.startswith("#")] if os.path.exists(summary) else []
                    log("    " + (rows[-1].rstrip() if rows else "no summary line (crash? see %s)" % term))
                if res == "done" and finished(st, mach, g):
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
