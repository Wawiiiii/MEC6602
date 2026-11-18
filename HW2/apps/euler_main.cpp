
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

// String without leading and trailing blanks.
static std::string trim(const std::string& s)
{
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) {
        return "";
    }
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// fichier txt de controle
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

// Value of a setting (or its default), removed from the map so that unknown keys can be reported.
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
        // With linear_solver = sgs the CFL is capped at 10 by default: above ~20 it is unstable at M 0.8.
        is.cfl_max = std::stod(take(kv, "cfl_max", lin == "sgs" ? "10" : "1e8"));
        is.gmres_restart = std::stoi(take(kv, "gmres_restart", "30"));
        is.gmres_max_iter = std::stoi(take(kv, "gmres_max_iter", "60"));
        is.gmres_rtol = std::stod(take(kv, "gmres_rtol", "1e-3"));
        is.gmres_fail = std::stod(take(kv, "gmres_fail", "0.5")); // linear solve failed above this: CFL backtracking
        is.sgs_sweeps = std::stoi(take(kv, "sgs_sweeps", "2"));
        is.sgs_omega = std::stod(take(kv, "sgs_omega", "1.0"));
        is.max_change = std::stod(take(kv, "max_change", "0.5"));
        const std::string cl_target = take(kv, "cl_target", "off");        // fixed-CL mode
        is.fixed_cl.on = cl_target != "off";
        if (is.fixed_cl.on) {
            is.fixed_cl.cl_target = std::stod(cl_target);
        }
        is.fixed_cl.dalpha0 = std::stod(take(kv, "cl_dalpha", "0.1"));
        is.fixed_cl.drop = std::stod(take(kv, "cl_drop", "1e-2"));
        is.fixed_cl.cl_tol = std::stod(take(kv, "cl_tol", "1e-6"));
        is.fixed_cl.max_dalpha = std::stod(take(kv, "cl_max_dalpha", "2"));
        euler::LIMITER_K = std::stod(take(kv, "limiter_k", "5"));            // limiter eps = (K h)^3, muscl.h
        const std::string kappa = take(kv, "muscl_kappa", "off");            // study only: unlimited kappa scheme, muscl.h
        euler::MUSCL_KAPPA_ON = kappa != "off";
        if (euler::MUSCL_KAPPA_ON) euler::MUSCL_KAPPA = std::stod(kappa);
        const std::string limiter = take(kv, "limiter", "hemker_koren");    // kappa^ 1/3 or 0, muscl.h
        if (limiter != "hemker_koren" && limiter != "van_albada") {
            throw std::runtime_error("limiter must be hemker_koren or van_albada");
        }
        euler::VAN_ALBADA = limiter == "van_albada";
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
        const std::string vortex = take(kv, "farfield_vortex", "off"); // farfield vortex correction, bc.h
        if (vortex != "on" && vortex != "off") {
            throw std::runtime_error("farfield_vortex must be on or off");
        }
        euler::FARFIELD_VORTEX = vortex == "on";
        euler::VORTEX_CHORD = s.ref.chord; // CL of the circulation uses the same reference as the output
        euler::VORTEX_XREF = s.ref.xref;   // vortex at the moment reference point (quarter chord)
        euler::VORTEX_YREF = s.ref.yref;
        const std::string field_file = take(kv, "field", "field.dat");         // Tecplot field output
        const std::string surface_file = take(kv, "surface", "surface_cp.dat"); // surface Cp output
        const std::string summary_file = take(kv, "summary", "summary.txt");    // one line per run, appended

        // Multigrid (scheme = mg) and full multigrid start (any scheme).
        euler::MultigridSettings ms;                                          // multigrid settings
        ms.levels = std::stoi(take(kv, "mg_levels", "0"));
        const std::string mg_cycle = take(kv, "mg_cycle", "w");              // v or w
        ms.n_pre = std::stoi(take(kv, "mg_pre", "1"));
        ms.n_post = std::stoi(take(kv, "mg_post", "1"));
        ms.n_coarse = std::stoi(take(kv, "mg_coarse", "2"));
        ms.fmg_cycles = std::stoi(take(kv, "mg_fmg", "0"));
        const std::string mg_smoother = take(kv, "mg_smoother", "rk");       // rk or sgs
        ms.sgs_cfl = std::stod(take(kv, "mg_sgs_cfl", "10"));   // 20 is unstable at M 0.8
        ms.sgs_sweeps = std::stoi(take(kv, "mg_sgs_sweeps", "2"));
        ms.relax = std::stod(take(kv, "mg_relax", "0")); // correction relaxation, 0 = auto
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
        if (is.fixed_cl.on && !implicit) {
            throw std::runtime_error("cl_target (fixed CL) needs an implicit scheme: bdf1, bdf2 or row");
        }
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
        euler::Freestream fs = euler::make_freestream(mach, alpha); // fixed CL: alpha is updated after each run
        euler::init_freestream(f, fs);

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

        // Short solver label (summary file) and a longer description (terminal).
        std::string label;  // e.g. bdf1-gmres, rk5, mg-w-sgs
        char solver[160];   // e.g. implicit bdf1 + gmres, CFL 1 -> 1e+08
        if (multigrid) {
            label = "mg-" + mg_cycle + "-" + mg_smoother;
            std::snprintf(solver, sizeof(solver), "multigrid %s-cycle, %d levels, smoother %s, fine CFL %.3g",
                          ms.gamma == 1 ? "V" : "W", euler::multigrid_levels(m, ms), mg_smoother.c_str(),
                          ms.smoother == euler::MgSmoother::RK ? ms.fine_rk.cfl : ms.sgs_cfl);
        } else if (implicit) {
            label = scheme + "-" + lin + (is.fixed_cl.on ? "-cl" : "");
            std::snprintf(solver, sizeof(solver), "implicit %s + %s, CFL %.3g -> %.3g", scheme.c_str(), lin.c_str(),
                          is.cfl0, is.cfl_max);
        } else {
            label = "rk" + std::to_string(stages);
            std::snprintf(solver, sizeof(solver), "explicit RK, %d stage(s), CFL %.3g", s.rk.stages, s.rk.cfl);
        }

        // Case block.
        std::printf("Case\n");
        std::printf("  mesh          %s  (%d x %d nodes)\n", mesh_file.c_str(), m.ni, m.nj);
        std::printf("  Mach          %.4f\n", mach);
        if (is.fixed_cl.on) {
            std::printf("  alpha         %.4f deg (first guess)\n", alpha);
            std::printf("  fixed CL      CL target %.6f, tol %.1e, first dalpha %.3g deg, alpha change at a drop of %.1e\n",
                        is.fixed_cl.cl_target, is.fixed_cl.cl_tol, is.fixed_cl.dalpha0, is.fixed_cl.drop);
        } else {
            std::printf("  alpha         %.4f deg\n", alpha);
        }
        std::printf("  order         %d\n", order);
        std::printf("  solver        %s\n", solver);
        std::printf("  threads       %d\n", euler::get_threads());

        /*
        Full multigrid start for any solver: solve on the coarse grids first (coarsest
        first, residual drop fmg_tol), inject each solution into the next finer grid, and use the
        result as the initial condition of the fine run. fmg_levels = total number of grids (1 = off).
        */
        const int fmg_nl = std::min(fmg_levels, euler::multigrid_levels(m, euler::MultigridSettings{})); // grids the mesh allows
        if (fmg_levels > 1 && fmg_nl < 2) {
            std::printf("  FMG start skipped: %d x %d cannot be coarsened\n", m.ni, m.nj);
        }
        if (fmg_nl > 1) {
            const std::vector<euler::Mesh> coarse = euler::coarse_hierarchy(m, fmg_nl);
            euler::Flow prev; // solution on the previous (coarser) grid
            for (size_t l = 0; l < coarse.size(); ++l) {
                euler::Flow cf; // flow on this coarse grid
                euler::allocate_flow(coarse[l], cf);
                euler::init_freestream(cf, fs);
                if (l > 0) {
                    euler::inject_piecewise_constant(coarse[l], coarse[l - 1], prev, cf);
                }
                const euler::RunResult rc = run_solver(coarse[l], cf, fmg_tol, "", 0);
                std::printf("  FMG start: %d x %d, %d iterations, %.2f s, drop %.1e", coarse[l].ni, coarse[l].nj,
                            rc.iterations, rc.seconds, rc.res / rc.res0);
                if (is.fixed_cl.on) {
                    fs = euler::make_freestream(mach, rc.alpha_deg); // the next grid starts from this alpha
                    std::printf(", alpha %.6f deg (CL %.6f)", rc.alpha_deg, rc.forces.cl);
                }
                std::printf("\n");
                prev = cf;
            }
            euler::inject_piecewise_constant(m, coarse.back(), prev, f);
        }

        const euler::RunResult r = run_solver(m, f, s.tol_drop, s.history_file, s.print_every);
        if (is.fixed_cl.on) {
            fs = euler::make_freestream(mach, r.alpha_deg); // outputs (force directions, title) at the final alpha
        }
        const char* status = r.diverged ? "DIVERGED"
                           : r.converged ? "converged (residual drop reached)"
                           : r.cauchy ? "converged (Cauchy criterion on the forces)"
                                      : "NOT converged (iteration limit)";

        // Forces from the wall pressure of the solver (the momentum flux through the wall).
        std::vector<double> pw; // wall pressure
        euler::wall_pressure(m, f, pw);
        const euler::Forces F = euler::integrate_forces(m, pw, fs, s.ref);
        std::vector<double> pwx; // wall pressure extrapolated from the first two cells (diagnostic only)
        euler::wall_pressure_extrapolated(m, f, pwx);
        const euler::Forces Fx = euler::integrate_forces(m, pwx, fs, s.ref);

        // Result block.
        std::printf("\nResult\n");
        std::printf("  status        %s\n", status);
        std::printf("  iterations    %d\n", r.iterations);
        std::printf("  time          %.2f s\n", r.seconds);
        std::printf("  residual      %.3e -> %.3e   (drop %.2e)\n", r.res0, r.res, r.res / r.res0);
        if (is.fixed_cl.on) {
            std::printf("  alpha         %.10f deg  (fixed CL %.6f, %d alpha changes on this grid)\n", fs.alpha_deg,
                        is.fixed_cl.cl_target, r.alpha_changes);
        }
        std::printf("  CL            %17.10e\n", F.cl);
        std::printf("  CD            %17.10e\n", F.cd);
        std::printf("  CM            %17.10e\n", F.cm);
        std::printf("  reference     chord %.5f, moment about (%.4f, %.4f)\n", s.ref.chord, s.ref.xref, s.ref.yref);
        std::printf("  diagnostic    CL %.8e  CD %.8e  CM %.8e  (p extrapolated to the wall in post-processing)\n",
                    Fx.cl, Fx.cd, Fx.cm);

        // Output files.
        std::printf("\nOutput\n");
        if (!s.history_file.empty()) {
            std::printf("  history       %s\n", s.history_file.c_str());
        }
        if (surface_file != "none") {
            euler::write_surface_cp(m, pw, fs, s.ref, F, surface_file);
            std::printf("  surface       %s\n", surface_file.c_str());
        }
        if (field_file != "none") {
            euler::write_tecplot_field(m, f, fs, field_file);
            std::printf("  field         %s\n", field_file.c_str());
        }
        if (summary_file != "none") {
            euler::RunSummary rs; // this run, one line of the summary file
            rs.ni = m.ni;
            rs.nj = m.nj;
            rs.mach = mach;
            rs.alpha_deg = fs.alpha_deg; // the final alpha in fixed-CL mode
            rs.order = order;
            rs.scheme = label;
            rs.status = r.diverged ? "diverged" : r.converged ? "converged" : r.cauchy ? "cauchy" : "max_iter";
            rs.iterations = r.iterations;
            rs.seconds = r.seconds;
            rs.drop = r.res / r.res0;
            rs.forces = F;
            euler::append_run_summary(summary_file, rs, s.ref);
            std::printf("  summary       %s  (line appended)\n", summary_file.c_str());
        }
        return r.diverged ? 2 : 0;
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
