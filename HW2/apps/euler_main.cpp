/*
   2D Euler solver, MEC6602E HW2.

   Usage:
     euler <case.txt>                         all settings read from a case file (see cases/example.txt)
     euler <mesh.x> <mach> <alpha_deg> [--key value ...]   same settings given on the command line

   Case file format: one "key = value" per line, '#' starts a comment, blank lines ignored.
   Keys (default in brackets):
     mesh        path of the plot3D grid file, relative to the folder the program is run from  [required]
     mach        free-stream Mach number                                                         [required]
     alpha       angle of attack in degrees                                                      [required]
     order       spatial order: 1 or 2                                                           [2]
     scheme      time scheme: rk (explicit), bdf1 (implicit Euler), bdf2, row (Rosenbrock-Wanner),
                 mg (FAS multigrid cycle, slides 41-53)                                         [rk]
     --- multigrid (scheme = mg) ---
     mg_cycle    v or w                                                                          [w]
     mg_levels   number of grids (0 = as many as possible, down to 9x9)                          [0]
     mg_pre, mg_post, mg_coarse   smoothing iterations before / after / on the coarsest grid      [1, 1, 2]
     mg_smoother rk (multistage, slide 18) or sgs (implicit Euler + block SGS)                   [rk]
     mg_sgs_cfl, mg_sgs_sweeps    CFL and sweeps of the SGS smoother                             [20, 2]
     mg_coarse_order   spatial order on the coarse grids                                         [1]
     mg_fmg      full multigrid start inside the cycle: cycles on each coarse level (0 = off)     [0]
     --- full multigrid start for any scheme (slide 42) ---
     fmg_levels  total number of grids of the FMG start (1 = off): solve coarsest first, inject   [1]
     fmg_tol     residual drop to reach on each FMG coarse grid                                  [1e-6]
     --- explicit (scheme = rk) ---
     stages      1 = explicit Euler, 3/4/5 = multistage of basic_algo slide 18                   [5]
     cfl         CFL number; 0 = the slide 18 value of the chosen scheme                         [0]
     --- implicit (scheme = bdf1 / bdf2 / row) ---
     linear_solver  gmres (Jacobian-free Newton-Krylov) or sgs (block SGS on the 1st-order Jacobian) [gmres]
     cfl0, cfl_max  CFL ramp sigma = cfl0 ||R0||/||R||, clamped to [cfl0, cfl_max] (slide 40)  [1, 1e8]
     gmres_restart, gmres_max_iter, gmres_rtol   GMRES settings                              [30, 60, 1e-3]
     sgs_sweeps, sgs_omega   SGS sweeps (preconditioner or solver) and relaxation             [2, 1.0]
     max_change     update limiter: largest relative change of rho and p per iteration      [0.5]
     limiter_k   limiter constant K, eps = (K h)^3 (muscl.h); very large = no limiting (study only) [5]
     threads     OpenMP threads (0 = all logical cores; results do not depend on it)            [0]
     wall_pressure  wall ghost pressure: extrapolated (2 p0 - p1, 2nd order) or mirror (slide 14) [extrapolated]
     wall_extrap_switch  extrapolated only after the residual dropped by this (1 = from the start;
                    ~1e-3 for supersonic cold starts); mirror before                           [1]
     max_iter    iteration limit                                                                 [100000]
     tol         stop when ||R_rho|| has dropped by this factor                                  [1e-12]
     cauchy_window  Cauchy criterion on CL, CD, CM: window in iterations (0 = off)               [0]
     cauchy_tol     Cauchy criterion: stop when max - min of CL, CD, CM over the window < this   [1e-6]
     chord       reference chord for CL, CD, CM                                                 [1.0]
     xref, yref  moment reference point (quarter chord)                                         [0.25 * chord, 0]
     history     residual / force history file (iteration, ||R_rho||, drop, CL, CD, CM, s)      [history.dat]
     field       Tecplot field file (node values), "none" to skip                               [field.dat]
     surface     surface Cp file, "none" to skip                                                 [surface_cp.dat]
     print_every print one line every N iterations (0 = never)                                   [500]
   On the command line the same keys are given as --key value (e.g. --order 1 --tol 1e-10).
   Plots: python post/plot_results.py <history> <surface> <field>  (post/plot_results.py).
*/

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "implicit.h"
#include "mesh.h"
#include "multigrid.h"
#include "muscl.h"
#include "output.h"
#include "solver.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Remove spaces and tabs at both ends.
static std::string trim(const std::string& s)
{
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) {
        return "";
    }
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Read "key = value" lines; '#' starts a comment.
static std::map<std::string, std::string> read_case_file(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("cannot open case file: " + filename);
    }
    std::map<std::string, std::string> kv; // key -> value
    std::string line;
    int lineno = 0; // line number, for error messages
    while (std::getline(in, line)) {
        ++lineno;
        const size_t hash = line.find('#');
        if (hash != std::string::npos) {
            line = line.substr(0, hash);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        const size_t eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error(filename + " line " + std::to_string(lineno) + ": expected key = value");
        }
        const std::string key = trim(line.substr(0, eq));
        const std::string val = trim(line.substr(eq + 1));
        if (key.empty() || val.empty()) {
            throw std::runtime_error(filename + " line " + std::to_string(lineno) + ": empty key or value");
        }
        kv[key] = val;
    }
    return kv;
}

// Take a key out of the map (so that unknown keys can be detected afterwards).
static std::string take(std::map<std::string, std::string>& kv, const std::string& key, const std::string& def,
                        bool required = false)
{
    const auto it = kv.find(key);
    if (it == kv.end()) {
        if (required) {
            throw std::runtime_error("missing required setting: " + key);
        }
        return def;
    }
    const std::string v = it->second;
    kv.erase(it);
    return v;
}

int main(int argc, char** argv)
{
    if (argc != 2 && argc < 4) {
        std::fprintf(stderr, "usage: %s <case.txt>\n"
                             "   or: %s <mesh.x> <mach> <alpha_deg> [--order 1|2] [--stages 1|3|4|5] [--cfl x] "
                             "[--max_iter N] [--tol x] [--history file] [--print_every N]\n",
                     argv[0], argv[0]);
        return 1;
    }

    try {
        // Collect the settings, from the case file or from the command line.
        std::map<std::string, std::string> kv; // key -> value
        if (argc == 2) {
            kv = read_case_file(argv[1]);
        } else {
            kv["mesh"] = argv[1];
            kv["mach"] = argv[2];
            kv["alpha"] = argv[3];
            for (int a = 4; a < argc; a += 2) {
                std::string key = argv[a];
                if (key.rfind("--", 0) != 0 || a + 1 >= argc) {
                    throw std::runtime_error("bad option " + key + " (expected --key value)");
                }
                kv[key.substr(2)] = argv[a + 1];
            }
        }

        const std::string mesh_file = take(kv, "mesh", "", true);           // grid file
        const double mach = std::stod(take(kv, "mach", "", true));          // free-stream Mach number
        const double alpha = std::stod(take(kv, "alpha", "", true));        // angle of attack, degrees
        const int order = std::stoi(take(kv, "order", "2"));                // spatial order
        const std::string scheme = take(kv, "scheme", "rk");                // time scheme
        const int stages = std::stoi(take(kv, "stages", "5"));              // number of RK stages
        const double cfl = std::stod(take(kv, "cfl", "0"));                 // CFL (0: the table value)
        euler::ImplicitSettings is;                                         // implicit settings
        const std::string lin = take(kv, "linear_solver", "gmres");         // implicit linear solver
        is.cfl0 = std::stod(take(kv, "cfl0", "1"));
        is.cfl_max = std::stod(take(kv, "cfl_max", "1e8"));
        is.gmres_restart = std::stoi(take(kv, "gmres_restart", "30"));
        is.gmres_max_iter = std::stoi(take(kv, "gmres_max_iter", "60"));
        is.gmres_rtol = std::stod(take(kv, "gmres_rtol", "1e-3"));
        is.sgs_sweeps = std::stoi(take(kv, "sgs_sweeps", "2"));
        is.sgs_omega = std::stod(take(kv, "sgs_omega", "1.0"));
        is.max_change = std::stod(take(kv, "max_change", "0.5"));
        euler::LIMITER_K = std::stod(take(kv, "limiter_k", "5"));            // limiter eps = (K h)^3, muscl.h
        euler::set_threads(std::stoi(take(kv, "threads", "0")));           // OpenMP threads (0 = all)
        const std::string wallp = take(kv, "wall_pressure", "extrapolated"); // wall ghost pressure, bc.h
        if (wallp != "mirror" && wallp != "extrapolated") {
            throw std::runtime_error("wall_pressure must be mirror or extrapolated");
        }
        euler::WALL_PRESSURE_EXTRAPOLATE = wallp == "extrapolated";
        euler::WALL_EXTRAP_SWITCH = std::stod(take(kv, "wall_extrap_switch", "1")); // start-up of the extrapolation
        euler::ExplicitSettings s;
        s.max_iter = std::stoi(take(kv, "max_iter", "100000"));
        s.tol_drop = std::stod(take(kv, "tol", "1e-12"));
        s.history_file = take(kv, "history", "history.dat");
        if (s.history_file == "none") {
            s.history_file.clear(); // "none" = no history file, like field / surface
        }
        s.print_every = std::stoi(take(kv, "print_every", "500"));
        s.cauchy_window = std::stoi(take(kv, "cauchy_window", "0"));
        s.cauchy_tol = std::stod(take(kv, "cauchy_tol", "1e-6"));
        s.ref.chord = std::stod(take(kv, "chord", "1.0"));
        s.ref.xref = std::stod(take(kv, "xref", std::to_string(0.25 * s.ref.chord)));
        s.ref.yref = std::stod(take(kv, "yref", "0"));
        const std::string field_file = take(kv, "field", "field.dat");         // Tecplot field output
        const std::string surface_file = take(kv, "surface", "surface_cp.dat"); // surface Cp output

        // Multigrid (scheme = mg) and full multigrid start (any scheme).
        euler::MultigridSettings ms;                                          // multigrid settings
        ms.levels = std::stoi(take(kv, "mg_levels", "0"));
        const std::string mg_cycle = take(kv, "mg_cycle", "w");              // v or w
        ms.n_pre = std::stoi(take(kv, "mg_pre", "1"));
        ms.n_post = std::stoi(take(kv, "mg_post", "1"));
        ms.n_coarse = std::stoi(take(kv, "mg_coarse", "2"));
        ms.fmg_cycles = std::stoi(take(kv, "mg_fmg", "0"));
        const std::string mg_smoother = take(kv, "mg_smoother", "rk");       // rk or sgs
        ms.sgs_cfl = std::stod(take(kv, "mg_sgs_cfl", "20"));
        ms.sgs_sweeps = std::stoi(take(kv, "mg_sgs_sweeps", "2"));
        const int mg_coarse_order = std::stoi(take(kv, "mg_coarse_order", "1")); // spatial order of coarse levels
        const int fmg_levels = std::stoi(take(kv, "fmg_levels", "1"));       // FMG start for rk / implicit (1 = off)
        const double fmg_tol = std::stod(take(kv, "fmg_tol", "1e-6"));       // residual drop on each FMG coarse grid
        if (!kv.empty()) {
            throw std::runtime_error("unknown setting: " + kv.begin()->first);
        }
        if (order != 1 && order != 2) {
            throw std::runtime_error("order must be 1 or 2");
        }
        s.order = order == 1 ? euler::Order::First : euler::Order::Second;
        const bool multigrid = scheme == "mg";                // FAS multigrid cycle
        const bool implicit = scheme != "rk" && !multigrid;   // one of the implicit schemes
        if (multigrid) {
            ms.gamma = mg_cycle == "v" ? 1 : 2;
            if (mg_cycle != "v" && mg_cycle != "w") {
                throw std::runtime_error("mg_cycle must be v or w");
            }
            if (mg_smoother != "rk" && mg_smoother != "sgs") {
                throw std::runtime_error("mg_smoother must be rk or sgs");
            }
            ms.smoother = mg_smoother == "rk" ? euler::MgSmoother::RK : euler::MgSmoother::SGS;
            ms.fine_order = s.order;
            ms.coarse_order = mg_coarse_order == 2 ? euler::Order::Second : euler::Order::First;
            ms.fine_rk = euler::rk_scheme(ms.fine_order, stages);
            ms.coarse_rk = euler::rk_scheme(ms.coarse_order, stages);
            if (cfl > 0.0) {
                ms.fine_rk.cfl = cfl;
            }
            ms.max_change = is.max_change;
            ms.max_iter = s.max_iter;
            ms.tol_drop = s.tol_drop;
            ms.cauchy_window = s.cauchy_window;
            ms.cauchy_tol = s.cauchy_tol;
            ms.ref = s.ref;
            ms.print_every = s.print_every;
            ms.history_file = s.history_file;
        } else if (!implicit) {
            s.rk = euler::rk_scheme(s.order, stages);
            if (cfl > 0.0) {
                s.rk.cfl = cfl;
            }
        } else {
            if (scheme == "bdf1") {
                is.scheme = euler::TimeScheme::BDF1;
            } else if (scheme == "bdf2") {
                is.scheme = euler::TimeScheme::BDF2;
            } else if (scheme == "row") {
                is.scheme = euler::TimeScheme::ROW;
            } else {
                throw std::runtime_error("scheme must be rk, bdf1, bdf2 or row");
            }
            if (lin == "gmres") {
                is.linsolve = euler::LinearSolver::GMRES;
            } else if (lin == "sgs") {
                is.linsolve = euler::LinearSolver::SGS;
            } else {
                throw std::runtime_error("linear_solver must be gmres or sgs");
            }
            is.order = s.order;
            is.max_iter = s.max_iter;
            is.tol_drop = s.tol_drop;
            is.cauchy_window = s.cauchy_window;
            is.cauchy_tol = s.cauchy_tol;
            is.ref = s.ref;
            is.print_every = s.print_every;
            is.history_file = s.history_file;
        }

        euler::Mesh m = euler::read_plot3d(mesh_file);
        euler::compute_metrics(m);
        euler::Flow f;
        euler::allocate_flow(m, f);
        const euler::Freestream fs = euler::make_freestream(mach, alpha);
        euler::init_freestream(f, fs);

        std::printf("OpenMP threads: %d\n", euler::get_threads());

        // One run of the chosen solver on a given mesh (used for the FMG coarse grids too).
        auto run_solver = [&](const euler::Mesh& mesh, euler::Flow& flow, double tol, const std::string& hist,
                              int print) {
            if (multigrid) {
                euler::MultigridSettings m2 = ms;
                m2.tol_drop = tol;
                m2.history_file = hist;
                m2.print_every = print;
                return euler::run_multigrid(mesh, flow, fs, m2);
            }
            if (implicit) {
                euler::ImplicitSettings i2 = is;
                i2.tol_drop = tol;
                i2.history_file = hist;
                i2.print_every = print;
                return euler::run_implicit(mesh, flow, fs, i2);
            }
            euler::ExplicitSettings s2 = s;
            s2.tol_drop = tol;
            s2.history_file = hist;
            s2.print_every = print;
            return euler::run_explicit(mesh, flow, fs, s2);
        };

        if (multigrid) {
            std::printf("mesh %s (%d x %d), M = %.4f, alpha = %.4f deg, order %d, multigrid %s-cycle, %d levels, "
                        "smoother %s, fine CFL %.3f\n",
                        mesh_file.c_str(), m.ni, m.nj, mach, alpha, order, ms.gamma == 1 ? "V" : "W",
                        euler::multigrid_levels(m, ms), mg_smoother.c_str(),
                        ms.smoother == euler::MgSmoother::RK ? ms.fine_rk.cfl : ms.sgs_cfl);
        } else if (!implicit) {
            std::printf("mesh %s (%d x %d), M = %.4f, alpha = %.4f deg, order %d, explicit RK %d stage(s), CFL %.3f\n",
                        mesh_file.c_str(), m.ni, m.nj, mach, alpha, order, s.rk.stages, s.rk.cfl);
        } else {
            std::printf("mesh %s (%d x %d), M = %.4f, alpha = %.4f deg, order %d, implicit %s + %s, CFL %.3g -> %.3g\n",
                        mesh_file.c_str(), m.ni, m.nj, mach, alpha, order, scheme.c_str(), lin.c_str(), is.cfl0,
                        is.cfl_max);
        }

        /*
        Full multigrid start for any solver (slide 42): solve on the coarse grids first (coarsest
        first, residual drop fmg_tol), inject each solution into the next finer grid, and use the
        result as the initial condition of the fine run. fmg_levels = total number of grids (1 = off).
        */
        if (fmg_levels > 1) {
            const int nl = std::min(fmg_levels, euler::multigrid_levels(m, euler::MultigridSettings{}));
            const std::vector<euler::Mesh> coarse = euler::coarse_hierarchy(m, nl);
            euler::Flow prev; // solution on the previous (coarser) grid
            for (size_t l = 0; l < coarse.size(); ++l) {
                euler::Flow cf; // flow on this coarse grid
                euler::allocate_flow(coarse[l], cf);
                euler::init_freestream(cf, fs);
                if (l > 0) {
                    euler::inject_piecewise_constant(coarse[l], coarse[l - 1], prev, cf);
                }
                const euler::RunResult rc = run_solver(coarse[l], cf, fmg_tol, "", 0);
                std::printf("  FMG start: %d x %d, %d iterations, %.2f s, drop %.1e\n", coarse[l].ni, coarse[l].nj,
                            rc.iterations, rc.seconds, rc.res / rc.res0);
                prev = cf;
            }
            euler::inject_piecewise_constant(m, coarse.back(), prev, f);
        }

        const euler::RunResult r = run_solver(m, f, s.tol_drop, s.history_file, s.print_every);
        const char* status = r.diverged ? "DIVERGED"
                           : r.converged ? "converged (residual)"
                           : r.cauchy ? "converged (Cauchy criterion on the forces)"
                                      : "not converged (iteration limit)";
        std::printf("%s after %d iterations, %.2f s: ||R_rho|| %.3e -> %.3e (drop %.2e)\n", status, r.iterations,
                    r.seconds, r.res0, r.res, r.res / r.res0);

        // Outputs: forces, surface Cp, field (run_explicit leaves the ghosts up to date).
        std::vector<double> pw; // wall pressure
        euler::wall_pressure(m, f, pw);
        const euler::Forces F = euler::integrate_forces(m, pw, fs, s.ref);
        std::printf("CL = %.8f   CD = %.8f   CM = %.8f   (chord %.5f, moment about (%.4f, %.4f))\n", F.cl, F.cd, F.cm,
                    s.ref.chord, s.ref.xref, s.ref.yref);
        std::vector<double> pwx; // wall pressure extrapolated from the first two cells (DEVIATIONS.md B6)
        euler::wall_pressure_extrapolated(m, f, pwx);
        const euler::Forces Fx = euler::integrate_forces(m, pwx, fs, s.ref);
        std::printf("CL = %.8f   CD = %.8f   CM = %.8f   (wall pressure extrapolated from the first two cells)\n",
                    Fx.cl, Fx.cd, Fx.cm);
        if (!s.history_file.empty()) {
            std::printf("history  -> %s\n", s.history_file.c_str());
        }
        if (surface_file != "none") {
            euler::write_surface_cp(m, pw, fs, s.ref, F, surface_file);
            std::printf("surface  -> %s\n", surface_file.c_str());
        }
        if (field_file != "none") {
            euler::write_tecplot_field(m, f, fs, field_file);
            std::printf("field    -> %s\n", field_file.c_str());
        }
        return r.diverged ? 2 : 0;
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
