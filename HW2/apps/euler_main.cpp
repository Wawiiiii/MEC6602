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
     stages      time scheme: 1 = explicit Euler, 3/4/5 = multistage of basic_algo slide 18      [5]
     cfl         CFL number; 0 = the slide 18 value of the chosen scheme                         [0]
     max_iter    iteration limit                                                                 [100000]
     tol         stop when ||R_rho|| has dropped by this factor                                  [1e-12]
     history     residual history output file                                                    [history.dat]
     print_every print one line every N iterations (0 = never)                                   [500]
   On the command line the same keys are given as --key value (e.g. --order 1 --tol 1e-10).

   Step 5: explicit solution only; forces, Cp and Tecplot output come with step 6.
*/

#include "bc.h"
#include "flow.h"
#include "mesh.h"
#include "solver.h"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

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
        const int stages = std::stoi(take(kv, "stages", "5"));              // number of RK stages
        const double cfl = std::stod(take(kv, "cfl", "0"));                 // CFL (0: the table value)
        euler::ExplicitSettings s;
        s.max_iter = std::stoi(take(kv, "max_iter", "100000"));
        s.tol_drop = std::stod(take(kv, "tol", "1e-12"));
        s.history_file = take(kv, "history", "history.dat");
        s.print_every = std::stoi(take(kv, "print_every", "500"));
        if (!kv.empty()) {
            throw std::runtime_error("unknown setting: " + kv.begin()->first);
        }
        if (order != 1 && order != 2) {
            throw std::runtime_error("order must be 1 or 2");
        }
        s.order = order == 1 ? euler::Order::First : euler::Order::Second;
        s.rk = euler::rk_scheme(s.order, stages);
        if (cfl > 0.0) {
            s.rk.cfl = cfl;
        }

        euler::Mesh m = euler::read_plot3d(mesh_file);
        euler::compute_metrics(m);
        euler::Flow f;
        euler::allocate_flow(m, f);
        const euler::Freestream fs = euler::make_freestream(mach, alpha);
        euler::init_freestream(f, fs);

        std::printf("mesh %s (%d x %d), M = %.4f, alpha = %.4f deg, order %d, %d stage(s), CFL %.3f\n",
                    mesh_file.c_str(), m.ni, m.nj, mach, alpha, order, s.rk.stages, s.rk.cfl);
        const euler::ExplicitResult r = euler::run_explicit(m, f, fs, s);
        std::printf("%s after %d iterations, %.2f s: ||R_rho|| %.3e -> %.3e (drop %.2e)\n",
                    r.diverged ? "DIVERGED" : (r.converged ? "converged" : "not converged"), r.iterations, r.seconds,
                    r.res0, r.res, r.res / r.res0);
        if (!s.history_file.empty()) {
            std::printf("residual history written to %s\n", s.history_file.c_str());
        }
        return r.diverged ? 2 : 0;
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
