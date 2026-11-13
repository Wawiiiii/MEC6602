/*
   Step 6 check: forces, moment, surface Cp, Tecplot field, Cauchy criterion
   (post.pdf slides 4-5, 16-18; Euler2D slide 20).

   H1 uniform pressure, fluid at rest: no force, no moment (closed contour)
   H2 linear wall pressure p = p_inf + k y: F = (0, -k A_airfoil), Mz = -k (x_centroid - xref) A
      (divergence theorem; A and the centroid from the wall polygon) -> sign and lever arm of the moment
   H3 CL, CD are the force turned by alpha: CL^2 + CD^2 = |F|^2 / (q c)^2, CL(alpha = 0) = Fy / (q c)
   H4 the forces are exactly the momentum fluxes through the wall faces of the residual
   H5 corner ghost cells are filled by apply_bc (copies of the wall/farfield ghosts across the cut)
   H6 free stream with all BCs farfield: field file = free stream at every node (Cp = 0, entropy = 0),
      surface Cp = 0; node i = 0 and node i = ni-1 (same point) carry the same values
   H7 converged alpha = 0: CL = CM = 0
   H8 Cauchy criterion stops the run early, with CL close to the fully converged value

   Usage: test_output <mesh.x> [<mesh.x> ...]  (H7/H8 run on the first grid only)
*/

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"
#include "output.h"
#include "residual.h"
#include "solver.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using euler::GAMMA;
using euler::NEQ;
using euler::NG;

constexpr double PI = 3.14159265358979323846; // MSVC has no M_PI by default

static int g_pass = 0; // number of passed checks
static int g_fail = 0; // number of failed checks

// Print one check line and count it.
static void check(bool ok, const char* fmt, ...)
{
    std::printf(ok ? "  [PASS] " : "  [FAIL] ");
    va_list args;
    va_start(args, fmt);
    std::vprintf(fmt, args);
    va_end(args);
    std::printf("\n");
    if (ok) {
        ++g_pass;
    } else {
        ++g_fail;
    }
}

// Signed area and centroid of the wall polygon (node line j = 0, node ni-1 = node 0).
static void wall_polygon(const euler::Mesh& m, double& area, double& xc, double& yc)
{
    double a = 0.0;  // twice the signed area
    double cx = 0.0; // centroid sums
    double cy = 0.0;
    for (int i = 0; i < m.ni - 1; ++i) {
        const int p = m.node(i, 0);
        const int q = m.node(i + 1, 0);
        const double cr = m.x[p] * m.y[q] - m.x[q] * m.y[p];
        a += cr;
        cx += (m.x[p] + m.x[q]) * cr;
        cy += (m.y[p] + m.y[q]) * cr;
    }
    xc = cx / (3.0 * a);
    yc = cy / (3.0 * a);
    area = std::abs(0.5 * a);
}

// Read all numbers of a text file after the first nskip lines and after '#' lines.
static std::vector<double> read_numbers(const std::string& file, int nskip)
{
    std::ifstream in(file);
    std::vector<double> v;
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        if (n++ < nskip || line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream ss(line);
        double x;
        while (ss >> x) {
            v.push_back(x);
        }
    }
    return v;
}

static void test_grid(const euler::Mesh& m, euler::Flow& f, const std::string& tmpdir)
{
    euler::ForceRef ref;      // chord 1, quarter chord
    const euler::Freestream fs = euler::make_freestream(0.5, 3.0);
    const double q = 0.5 * fs.mach * fs.mach; // dynamic pressure

    // H1: fluid at rest, uniform pressure.
    {
        const double qr[NEQ] = {1.0, 0.0, 0.0, fs.q[3] + 0.3};
        double w[NEQ];
        euler::prim_to_cons(qr, w);
        for (int k = 0; k < NEQ; ++k) {
            std::fill(f.W[k].begin(), f.W[k].end(), w[k]);
        }
        euler::apply_bc(m, f, fs);
        const euler::Forces F = euler::compute_forces(m, f, fs, ref);
        check(std::abs(F.fx) < 1e-14 && std::abs(F.fy) < 1e-14 && std::abs(F.mz) < 1e-14,
              "H1 uniform pressure at rest: F = (%.1e, %.1e), Mz = %.1e (expect 0)", F.fx, F.fy, F.mz);
    }

    // H2: p_w = p_inf + k y at the face midpoints.
    {
        const double k = 0.7; // pressure gradient
        std::vector<double> pw(m.nic);
        for (int i = 0; i < m.nic; ++i) {
            pw[i] = fs.q[3] + k * 0.5 * (m.y[m.node(i, 0)] + m.y[m.node(i + 1, 0)]);
        }
        double A, xc, yc;
        wall_polygon(m, A, xc, yc);
        const euler::Forces F = euler::integrate_forces(m, pw, fs, ref);
        const double fy_exact = -k * A;                    // buoyancy-like force (divergence theorem)
        const double mz_exact = -k * (xc - ref.xref) * A;  // its moment about the quarter chord
        const double efy = std::abs(F.fy - fy_exact) / std::abs(fy_exact);
        const double emz = std::abs(F.mz - mz_exact) / std::abs(mz_exact);
        // Mz: midpoint rule on a quadratic integrand -> O(h^2); measured err * nic^2 ~ 3.2 on 9..257.
        check(std::abs(F.fx) < 1e-14 && efy < 1e-12 && emz * m.nic * m.nic < 5.0,
              "H2 p = p_inf + k y: Fx %.1e, Fy rel. err %.1e (exact on a polygon), Mz rel. err %.2e "
              "(midpoint rule, err * nic^2 = %.2f, O(h^2))",
              F.fx, efy, emz, emz * m.nic * m.nic);
    }

    // H3: CL, CD from the force.
    {
        std::vector<double> pw(m.nic);
        for (int i = 0; i < m.nic; ++i) {
            const double xm = 0.5 * (m.x[m.node(i, 0)] + m.x[m.node(i + 1, 0)]);
            const double ym = 0.5 * (m.y[m.node(i, 0)] + m.y[m.node(i + 1, 0)]);
            pw[i] = fs.q[3] + 0.3 * xm * xm - 0.2 * ym; // any pressure
        }
        const euler::Forces F = euler::integrate_forces(m, pw, fs, ref);
        const double e1 = std::abs(F.cl * F.cl + F.cd * F.cd - (F.fx * F.fx + F.fy * F.fy) / (q * q));
        const euler::Freestream fs0 = euler::make_freestream(0.5, 0.0);
        const euler::Forces F0 = euler::integrate_forces(m, pw, fs0, ref);
        const double e2 = std::abs(F0.cl - F0.fy / q) + std::abs(F0.cd - F0.fx / q);
        check(e1 < 1e-12 && e2 < 1e-14, "H3 CL^2 + CD^2 = |F|^2/q^2 (%.1e), alpha = 0: CL = Fy/q, CD = Fx/q (%.1e)", e1, e2);
    }

    // H4 / H5 on a random field with the real wall.
    {
        std::vector<double> rnd(m.vol.size());
        for (size_t c = 0; c < rnd.size(); ++c) {
            rnd[c] = std::sin(1.7 * c) * 0.1; // deterministic "random" perturbation
        }
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const double qq[NEQ] = {1.0 + rnd[c], fs.q[1] + rnd[c], fs.q[2] - rnd[c], fs.q[3] * (1.0 + rnd[c])};
                double w[NEQ];
                euler::prim_to_cons(qq, w);
                for (int k = 0; k < NEQ; ++k) {
                    f.W[k][c] = w[k];
                }
            }
        }
        euler::apply_bc(m, f, fs);
        euler::compute_primitives(f);
        const euler::Forces F = euler::compute_forces(m, f, fs, ref);
        double fx = 0.0; // force from the residual's wall fluxes
        double fy = 0.0;
        for (int i = 0; i < m.nic; ++i) {
            double Fl[NEQ];
            euler::j_face_flux(m, f, euler::Order::Second, i, 0, Fl);
            const int fc = m.node(i, 0);
            fx -= Fl[1] - fs.q[3] * m.njx[fc] * m.sj[fc];
            fy -= Fl[2] - fs.q[3] * m.njy[fc] * m.sj[fc];
        }
        check(std::abs(F.fx - fx) < 1e-14 && std::abs(F.fy - fy) < 1e-14,
              "H4 forces = momentum fluxes through the wall faces of the residual: diff %.1e %.1e", F.fx - fx, F.fy - fy);

        bool finite = true; // corners finite
        double diff = 0.0;  // corner vs ghost across the cut
        for (int k = 0; k < NEQ; ++k) {
            for (int g = 1; g <= NG; ++g) {
                for (int h = 1; h <= NG; ++h) {
                    const double a = f.W[k][m.cell(-g, -h)];            // corner, lower left
                    const double b = f.W[k][m.cell(m.nic - 1 + g, -h)]; // corner, lower right
                    const double c = f.W[k][m.cell(-g, m.njc - 1 + h)]; // corner, upper left
                    finite = finite && std::isfinite(a) && std::isfinite(b) && std::isfinite(c);
                    diff = std::max(diff, std::abs(a - f.W[k][m.cell(m.nic - g, -h)]));
                    diff = std::max(diff, std::abs(b - f.W[k][m.cell(g - 1, -h)]));
                }
            }
        }
        check(finite && diff == 0.0, "H5 corner ghosts filled (finite) and = ghosts across the cut (diff %.1e)", diff);

        // Cut consistency of the field file on this non-uniform field: node (0, j) = node (ni-1, j).
        const std::string ff = tmpdir + "/_test_field.dat";
        euler::write_tecplot_field(m, f, fs, ff);
        const std::vector<double> v = read_numbers(ff, 3);
        const size_t nn = static_cast<size_t>(m.ni) * m.nj; // nodes
        double cut = 0.0;                                   // largest difference across the cut
        bool fin = v.size() == 9 * nn;                      // right count and all finite
        for (size_t k = 0; k < v.size(); ++k) {
            fin = fin && std::isfinite(v[k]);
        }
        if (v.size() == 9 * nn) {
            for (int var = 0; var < 9; ++var) {
                for (int j = 0; j < m.nj; ++j) {
                    const size_t a = var * nn + static_cast<size_t>(j) * m.ni;
                    cut = std::max(cut, std::abs(v[a] - v[a + m.ni - 1]));
                }
            }
        }
        check(fin && cut < 1e-8, "H6 field file: %zu values (expect %zu), all finite, node i=0 vs i=ni-1: %.1e",
              v.size(), 9 * nn, cut);
        std::remove(ff.c_str());
    }

    // H6: free stream with all BCs farfield -> field = free stream, Cp = 0.
    {
        euler::init_freestream(f, fs);
        euler::apply_bc(m, f, fs, euler::JminBc::Farfield);
        const std::string ff = tmpdir + "/_test_field.dat";
        const std::string sf = tmpdir + "/_test_surface.dat";
        euler::write_tecplot_field(m, f, fs, ff);
        std::vector<double> pw;
        euler::wall_pressure(m, f, pw);
        const euler::Forces F = euler::integrate_forces(m, pw, fs, ref);
        euler::write_surface_cp(m, pw, fs, ref, F, sf);
        const std::vector<double> v = read_numbers(ff, 3);
        const size_t nn = static_cast<size_t>(m.ni) * m.nj;
        const double expect[7] = {1.0, fs.q[1], fs.q[2], fs.q[3], fs.mach, 0.0, 0.0}; // Rho U V P Mach Cp Entropy
        double worst = 0.0;
        for (int var = 2; var < 9; ++var) {
            for (size_t n = 0; n < nn; ++n) {
                worst = std::max(worst, std::abs(v[var * nn + n] - expect[var - 2]));
            }
        }
        const std::vector<double> s = read_numbers(sf, 3); // after the 3 '#' lines: TITLE, VARIABLES, ZONE are text
        double cpmax = 0.0;
        int rows = 0;
        for (size_t r = 0; r + 3 < s.size(); r += 4) {
            cpmax = std::max(cpmax, std::abs(s[r + 2]));
            ++rows;
        }
        // The field file is printed with 10 significant digits (%.9e): read-back agrees to ~1e-10.
        check(worst < 1e-9 && rows == m.nic && cpmax < 1e-13 && std::abs(F.cl) < 1e-13,
              "H6 free stream: field = free stream at every node (%.1e, print precision), %d surface rows (expect %d), "
              "max |Cp| %.1e",
              worst, rows, m.nic, cpmax);
        std::remove(ff.c_str());
        std::remove(sf.c_str());
    }
}

// H7 / H8 on one grid (explicit runs).
static void test_runs(const euler::Mesh& m, euler::Flow& f)
{
    // H7: alpha = 0.
    {
        const euler::Freestream fs = euler::make_freestream(0.5, 0.0);
        euler::init_freestream(f, fs);
        euler::ExplicitSettings s;
        s.rk = euler::rk_scheme(euler::Order::Second, 5);
        s.tol_drop = 1e-11;
        s.max_iter = 100000;
        const euler::RunResult r = euler::run_explicit(m, f, fs, s);
        const euler::Forces F = euler::compute_forces(m, f, fs, s.ref);
        check(r.converged && std::abs(F.cl) < 1e-11 && std::abs(F.cm) < 1e-11 && F.cd > 0.0,
              "H7 alpha = 0 converged: CL %.1e, CM %.1e (expect 0), CD %.6f", F.cl, F.cm, F.cd);
    }

    // H8: Cauchy criterion vs full convergence.
    {
        const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
        euler::ExplicitSettings s;
        s.rk = euler::rk_scheme(euler::Order::Second, 5);
        s.max_iter = 100000;

        euler::init_freestream(f, fs);
        s.tol_drop = 1e-12;
        const euler::RunResult full = euler::run_explicit(m, f, fs, s);

        euler::init_freestream(f, fs);
        s.tol_drop = 1e-30;     // residual criterion off
        s.cauchy_window = 500;  // forces flat over 500 iterations
        s.cauchy_tol = 1e-6;
        const euler::RunResult cau = euler::run_explicit(m, f, fs, s);
        const double dcl = std::abs(cau.forces.cl - full.forces.cl);
        const double dcd = std::abs(cau.forces.cd - full.forces.cd);
        check(cau.cauchy && cau.iterations < full.iterations && dcl < 1e-5 && dcd < 1e-5,
              "H8 Cauchy stop after %d it (full: %d it, drop %.1e): |dCL| %.1e, |dCD| %.1e", cau.iterations,
              full.iterations, cau.res / cau.res0, dcl, dcd);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x> [<mesh.x> ...]\n", argv[0]);
        return 1;
    }
    try {
        const std::string tmpdir = "."; // temporary files go in the current folder and are removed
        for (int a = 1; a < argc; ++a) {
            euler::Mesh m = euler::read_plot3d(argv[a]);
            euler::compute_metrics(m);
            euler::Flow f;
            euler::allocate_flow(m, f);
            std::printf("\n=== %s (%d x %d nodes)\n", argv[a], m.ni, m.nj);
            test_grid(m, f, tmpdir);
            if (a == 1) {
                test_runs(m, f);
            }
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
