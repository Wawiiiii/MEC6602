/*
   Step 4 check: boundary conditions with ghost cells (Euler2D_basic_algo slides 7, 13-16, 20).

   D. Boundary states alone (no mesh)
      D1 wall ghost (slide 14): same rho and p, normal velocity reversed, tangential kept;
         averaged state has V.n = 0, so the wall flux is (0, p n, 0)
      D2 farfield (slides 15-16): free stream in -> free stream out, in all four regimes
      D3 farfield subsonic inflow/outflow obey the characteristic relations they come from
      D4 farfield supersonic inflow/outflow copy the free stream / the interior
   E. Boundary conditions on the airfoil grids
      E1 slide 20 sanity check: all BCs farfield + free stream -> residual 0 ("nothing moves")
      E2 every ghost cell the residual needs is filled (NaN start -> no NaN in R), both jmin types
      E3 i-cut: ghosts are exact copies, fluxes on the two sides of the cut identical
      E4 wall: no mass and no energy through any wall face, momentum flux = p n dS
      E5 farfield: 1/2 (W_ghost + W_interior) = Wb exactly
      E6 free stream + real wall: only the wall rows move (1 row at first order, 2 at second)
      E7 up/down symmetry (the grids are exactly symmetric): a symmetric field gives a
         mirror-symmetric residual, both orders
      info: which farfield regime each face takes at M = 0.5 and M = 1.5

   Usage: test_bc <mesh.x> [<mesh.x> ...]
   Exit code = number of failed checks.
*/

#include "bc.h"
#include "flow.h"
#include "mesh.h"
#include "residual.h"
#include "roe.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <limits>
#include <random>
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

// Largest |a_k - b_k| over the 4 components.
static double max_abs_diff(const double a[NEQ], const double b[NEQ])
{
    double d = 0.0;
    for (int k = 0; k < NEQ; ++k) {
        d = std::max(d, std::abs(a[k] - b[k]));
    }
    return d;
}

// Speed of sound of a primitive state.
static double sound(const double q[NEQ])
{
    return std::sqrt(GAMMA * q[3] / q[0]);
}

// ===================================================================== D. boundary states

static void test_states()
{
    std::printf("\n=== D. Boundary states (slides 14-16)\n");
    std::mt19937 gen(99);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    // D1: wall ghost.
    {
        double worst = 0.0; // largest error over all the properties checked
        for (int n = 0; n < 10000; ++n) {
            const double a = 2.0 * PI * U(gen); // wall normal direction
            const double nx = std::cos(a);
            const double ny = std::sin(a);
            const double qi[NEQ] = {0.5 + U(gen), 2.0 * U(gen) - 1.0, 2.0 * U(gen) - 1.0, 0.3 + U(gen)}; // interior
            double qg[NEQ];
            euler::wall_ghost(qi, nx, ny, qg);

            const double Vn_i = qi[1] * nx + qi[2] * ny;  // interior normal velocity
            const double Vn_g = qg[1] * nx + qg[2] * ny;  // ghost normal velocity
            const double Vt_i = -qi[1] * ny + qi[2] * nx; // interior tangential velocity
            const double Vt_g = -qg[1] * ny + qg[2] * nx; // ghost tangential velocity
            worst = std::max({worst, std::abs(qg[0] - qi[0]), std::abs(qg[3] - qi[3]), std::abs(Vn_g + Vn_i),
                              std::abs(Vt_g - Vt_i)});

            // Averaged state (slide 7 #2) and its flux: (0, p01 nx, p01 ny, 0).
            double wi[NEQ];
            double wg[NEQ];
            double w01[NEQ];
            double q01[NEQ];
            double F[NEQ];
            euler::prim_to_cons(qi, wi);
            euler::prim_to_cons(qg, wg);
            for (int k = 0; k < NEQ; ++k) {
                w01[k] = 0.5 * (wi[k] + wg[k]);
            }
            euler::cons_to_prim(w01, q01);
            euler::physical_flux(q01, nx, ny, F);
            const double Fexp[NEQ] = {0.0, q01[3] * nx, q01[3] * ny, 0.0};
            worst = std::max(worst, max_abs_diff(F, Fexp));
        }
        check(worst < 1e-14, "D1 wall ghost: rho, p kept, V.n reversed, V.t kept, wall flux (0, p n, 0): max err %.2e",
              worst);
    }

    // D2: interior = free stream -> boundary = free stream, for every Mach number and face direction.
    {
        double worst = 0.0; // largest difference
        for (double M : {0.3, 0.8, 0.99, 1.01, 1.5, 3.0}) {
            for (int n = 0; n < 360; ++n) {
                const double a = 2.0 * PI * n / 360.0; // face normal direction
                const euler::Freestream fs = euler::make_freestream(M, 7.0);
                double qb[NEQ];
                euler::farfield_state(fs.q, fs.q, std::cos(a), std::sin(a), qb);
                worst = std::max(worst, max_abs_diff(qb, fs.q));
            }
        }
        check(worst < 1e-15, "D2 farfield keeps the free stream (M 0.3..3, all directions): max diff %.2e", worst);
    }

    /*
       D3: the slide formulas come from the 1D characteristic relations along n (n outward),
       with rho0 c0 frozen at the interior state:
         along V + c (leaves the domain, from d) : p + rho0 c0 V = const
         along V - c (enters, from a)            : p - rho0 c0 V = const
         along V (entropy, tangential velocity)  : rho - p / c0^2 = const, V_t = const
       Subsonic inflow : the V + c relation from d, the other three from a.
       Subsonic outflow: p from a, the other three from d.
       Checked here as relations, independently of how the formulas are written.
    */
    {
        double worst_in = 0.0;  // subsonic inflow, largest error
        double worst_out = 0.0; // subsonic outflow, largest error
        int n_in = 0;           // number of inflow cases
        int n_out = 0;          // number of outflow cases
        const euler::Freestream fs = euler::make_freestream(0.6, 3.0);
        const double* qa = fs.q;
        for (int n = 0; n < 20000; ++n) {
            const double a = 2.0 * PI * U(gen);
            const double nx = std::cos(a);
            const double ny = std::sin(a);
            // Interior state: free stream with a random disturbance (stays subsonic along n).
            const double qd[NEQ] = {qa[0] * (1.0 + 0.2 * (U(gen) - 0.5)), qa[1] + 0.2 * (U(gen) - 0.5),
                                    qa[2] + 0.2 * (U(gen) - 0.5), qa[3] * (1.0 + 0.2 * (U(gen) - 0.5))};
            const double c0 = sound(qd);
            const double rc = qd[0] * c0;               // rho0 c0
            const double Vd = qd[1] * nx + qd[2] * ny;  // interior normal velocity
            if (std::abs(Vd) >= c0) {
                continue;
            }
            double qb[NEQ];
            euler::farfield_state(qd, qa, nx, ny, qb);
            const double Va = qa[1] * nx + qa[2] * ny;  // free-stream normal velocity
            const double Vb = qb[1] * nx + qb[2] * ny;  // boundary normal velocity
            const double Vta = -qa[1] * ny + qa[2] * nx; // tangential velocities
            const double Vtd = -qd[1] * ny + qd[2] * nx;
            const double Vtb = -qb[1] * ny + qb[2] * nx;
            if (Vd < 0.0) {
                ++n_in;
                worst_in = std::max({worst_in,
                                     std::abs((qb[3] + rc * Vb) - (qd[3] + rc * Vd)),                           // V + c from d
                                     std::abs((qb[3] - rc * Vb) - (qa[3] - rc * Va)),                           // V - c from a
                                     std::abs((qb[0] - qb[3] / (c0 * c0)) - (qa[0] - qa[3] / (c0 * c0))),       // entropy from a
                                     std::abs(Vtb - Vta)});                                                     // shear from a
            } else {
                ++n_out;
                worst_out = std::max({worst_out,
                                      std::abs(qb[3] - qa[3]),                                                  // p from a
                                      std::abs((qb[3] + rc * Vb) - (qd[3] + rc * Vd)),                          // V + c from d
                                      std::abs((qb[0] - qb[3] / (c0 * c0)) - (qd[0] - qd[3] / (c0 * c0))),      // entropy from d
                                      std::abs(Vtb - Vtd)});                                                    // shear from d
            }
        }
        check(n_in > 1000 && worst_in < 1e-14, "D3 subsonic inflow = characteristic relations, %d cases: max err %.2e",
              n_in, worst_in);
        check(n_out > 1000 && worst_out < 1e-14, "D3 subsonic outflow = characteristic relations, %d cases: max err %.2e",
              n_out, worst_out);
    }

    // D4: supersonic along n -> copy.
    {
        const euler::Freestream fs = euler::make_freestream(2.0, 0.0);
        const double qd[NEQ] = {1.1, 1.9, 0.1, 0.75}; // interior, normal Mach ~ 1.7 along +x
        double qb_out[NEQ];
        double qb_in[NEQ];
        euler::farfield_state(qd, fs.q, 1.0, 0.0, qb_out);  // V.n > 0: outflow
        euler::farfield_state(qd, fs.q, -1.0, 0.0, qb_in);  // V.n < 0: inflow
        check(max_abs_diff(qb_out, qd) == 0.0 && max_abs_diff(qb_in, fs.q) == 0.0,
              "D4 supersonic outflow -> Wd, supersonic inflow -> Wa (exact copies)");
    }
}

// ===================================================================== E. on the grids

// Largest |R| over the physical cells, divided by the cell perimeter; flags NaN.
static double max_residual(const euler::Mesh& m, const euler::Flow& f, bool& has_nan, int j_from = 0)
{
    double worst = 0.0;
    has_nan = false;
    for (int j = j_from; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double per = m.si[m.node(i, j)] + m.si[m.node(i + 1, j)] + m.sj[m.node(i, j)] + m.sj[m.node(i, j + 1)];
            for (int k = 0; k < NEQ; ++k) {
                const double r = f.R[k][m.cell(i, j)];
                if (std::isnan(r)) {
                    has_nan = true;
                }
                worst = std::max(worst, std::abs(r) / per);
            }
        }
    }
    return worst;
}

// Random positive field in the physical cells only; every ghost reset to NaN first.
static void fill_random_physical(const euler::Mesh& m, euler::Flow& f, std::mt19937& gen, double mach)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    for (int k = 0; k < NEQ; ++k) {
        std::fill(f.W[k].begin(), f.W[k].end(), nan);
    }
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double q[NEQ] = {1.0 + 0.2 * U(gen), mach + 0.2 * U(gen), 0.2 * U(gen), 1.0 / GAMMA + 0.1 * U(gen)};
            double w[NEQ];
            euler::prim_to_cons(q, w);
            for (int k = 0; k < NEQ; ++k) {
                f.W[k][m.cell(i, j)] = w[k];
            }
        }
    }
}

static void test_grid(const euler::Mesh& m, euler::Flow& f)
{
    const euler::Order orders[2] = {euler::Order::First, euler::Order::Second};
    const char* oname[2] = {"1st", "2nd"};

    // E1: slide 20 sanity check.
    {
        double worst = 0.0;
        bool any_nan = false;
        const double cases[4][2] = {{0.5, 1.25}, {0.8, 1.25}, {1.5, 0.0}, {0.3, -20.0}}; // (Mach, alpha in degrees)
        for (const auto& c : cases) {
            const euler::Freestream fs = euler::make_freestream(c[0], c[1]);
            euler::init_freestream(f, fs);
            euler::apply_bc(m, f, fs, euler::JminBc::Farfield);
            for (euler::Order o : orders) {
                euler::compute_residual(m, f, o);
                bool nan = false;
                worst = std::max(worst, max_residual(m, f, nan));
                any_nan = any_nan || nan;
            }
        }
        check(!any_nan && worst < 1e-13,
              "E1 all BCs farfield + free stream (M 0.5/0.8/1.5/0.3, 1st+2nd): max |R|/perimeter = %.2e (slide 20: nothing moves)",
              worst);
    }

    // E2: every ghost the residual reads is filled by apply_bc.
    {
        std::mt19937 gen(5);
        const euler::Freestream fs = euler::make_freestream(0.8, 1.25);
        bool any_nan = false;
        for (euler::JminBc jb : {euler::JminBc::Wall, euler::JminBc::Farfield}) {
            fill_random_physical(m, f, gen, 0.8);
            euler::apply_bc(m, f, fs, jb);
            for (euler::Order o : orders) {
                euler::compute_residual(m, f, o);
                bool nan = false;
                max_residual(m, f, nan);
                any_nan = any_nan || nan;
            }
        }
        check(!any_nan, "E2 ghosts start NaN, physical cells random: no NaN in R after apply_bc (wall/farfield jmin, 1st/2nd)");
    }

    // E3 - E5 on a random field with the real wall.
    std::mt19937 gen(6);
    const euler::Freestream fs = euler::make_freestream(0.8, 1.25);
    fill_random_physical(m, f, gen, 0.8);
    euler::apply_bc(m, f, fs, euler::JminBc::Wall);
    euler::compute_primitives(f);

    // E3: i-cut.
    {
        double copy = 0.0; // largest ghost/source difference
        for (int k = 0; k < NEQ; ++k) {
            for (int j = 0; j < m.njc; ++j) {
                for (int g = 1; g <= NG; ++g) {
                    copy = std::max(copy, std::abs(f.W[k][m.cell(-g, j)] - f.W[k][m.cell(m.nic - g, j)]));
                    copy = std::max(copy, std::abs(f.W[k][m.cell(m.nic - 1 + g, j)] - f.W[k][m.cell(g - 1, j)]));
                }
            }
        }
        double cut = 0.0; // largest flux difference across the cut
        for (euler::Order o : orders) {
            for (int j = 0; j < m.njc; ++j) {
                double Fa[NEQ];
                double Fb[NEQ];
                euler::i_face_flux(m, f, o, 0, j, Fa);
                euler::i_face_flux(m, f, o, m.ni - 1, j, Fb);
                cut = std::max(cut, max_abs_diff(Fa, Fb));
            }
        }
        check(copy == 0.0 && cut < 1e-14, "E3 i-cut: ghosts exact copies (%.1e), flux i=0 vs i=ni-1: %.2e", copy, cut);
    }

    // E4: wall faces.
    {
        double mass_energy = 0.0; // largest |mass| or |energy| flux through a wall face, per unit length
        double momentum = 0.0;    // largest error of the momentum flux vs p01 n dS
        for (euler::Order o : orders) {
            for (int i = 0; i < m.nic; ++i) {
                const int fc = m.node(i, 0);
                double F[NEQ];
                euler::j_face_flux(m, f, o, i, 0, F);
                double w01[NEQ];
                double q01[NEQ];
                for (int k = 0; k < NEQ; ++k) {
                    w01[k] = 0.5 * (f.W[k][m.cell(i, -1)] + f.W[k][m.cell(i, 0)]);
                }
                euler::cons_to_prim(w01, q01);
                mass_energy = std::max({mass_energy, std::abs(F[0]) / m.sj[fc], std::abs(F[3]) / m.sj[fc]});
                momentum = std::max({momentum, std::abs(F[1] - q01[3] * m.njx[fc] * m.sj[fc]) / m.sj[fc],
                                     std::abs(F[2] - q01[3] * m.njy[fc] * m.sj[fc]) / m.sj[fc]});
            }
        }
        check(mass_energy < 1e-15 && momentum < 1e-15,
              "E4 wall faces: mass/energy flux %.2e (expect 0), momentum - p n dS %.2e", mass_energy, momentum);
    }

    // E5: farfield face average = boundary state.
    {
        double worst = 0.0; // largest difference
        for (int i = 0; i < m.nic; ++i) {
            const int fc = m.node(i, m.nj - 1);
            const int cd = m.cell(i, m.njc - 1); // interior cell
            const int cg = m.cell(i, m.njc);     // first ghost
            double wd[NEQ];
            double qd[NEQ];
            double qb[NEQ];
            double wb[NEQ];
            double w01[NEQ];
            for (int k = 0; k < NEQ; ++k) {
                wd[k] = f.W[k][cd];
                w01[k] = 0.5 * (f.W[k][cg] + f.W[k][cd]);
            }
            euler::cons_to_prim(wd, qd);
            euler::farfield_state(qd, fs.q, m.njx[fc], m.njy[fc], qb);
            euler::prim_to_cons(qb, wb);
            worst = std::max(worst, max_abs_diff(w01, wb));
        }
        check(worst < 1e-14, "E5 farfield: 1/2 (W_ghost + W_interior) = Wb: max diff %.2e", worst);
    }

    /*
       E6: free stream everywhere + real wall. The wall ghosts differ from the free stream, so the
       first row (j = 0) sees it; at second order the MUSCL stencil of face j = 1 also reads the
       ghost, so row j = 1 sees it too; every other row must stay exactly at R = 0.
    */
    {
        const euler::Freestream fs0 = euler::make_freestream(0.5, 3.0);
        for (int o = 0; o < 2; ++o) {
            euler::init_freestream(f, fs0);
            euler::apply_bc(m, f, fs0, euler::JminBc::Wall);
            euler::compute_residual(m, f, orders[o]);
            bool nan = false;
            const int rows = o == 0 ? 1 : 2; // rows that may move
            const double rest = max_residual(m, f, nan, rows);
            double first = 0.0; // largest |R| in the rows that may move
            for (int j = 0; j < rows; ++j) {
                for (int i = 0; i < m.nic; ++i) {
                    first = std::max(first, std::abs(f.R[1][m.cell(i, j)]) + std::abs(f.R[2][m.cell(i, j)]));
                }
            }
            check(!nan && rest < 1e-13 && first > 1e-6,
                  "E6 %s order, free stream + wall: |R| in rows j<%d up to %.2e, all other rows %.2e (expect 0)",
                  oname[o], rows, first, rest);
        }
    }

    /*
       E7: up/down symmetry. Cell (i, j) mirrors cell (nic-1-i, j) (checked in the grid files:
       x(i) = x(ni-1-i), y(i) = -y(ni-1-i) exactly). Field: random in the lower half, mirrored
       into the upper half (rho, p same, u same, v opposite), alpha = 0 free stream.
       Then R0, R1, R3 must be equal and R2 opposite in mirrored cells.
    */
    {
        std::mt19937 g7(7);
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        const euler::Freestream fs0 = euler::make_freestream(0.8, 0.0);
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic / 2; ++i) {
                const double q[NEQ] = {1.0 + 0.2 * U(g7), 0.8 + 0.2 * U(g7), 0.2 * U(g7), 1.0 / GAMMA + 0.1 * U(g7)};
                const double qm[NEQ] = {q[0], q[1], -q[2], q[3]}; // mirror image
                double w[NEQ];
                double wm[NEQ];
                euler::prim_to_cons(q, w);
                euler::prim_to_cons(qm, wm);
                for (int k = 0; k < NEQ; ++k) {
                    f.W[k][m.cell(i, j)] = w[k];
                    f.W[k][m.cell(m.nic - 1 - i, j)] = wm[k];
                }
            }
        }
        euler::apply_bc(m, f, fs0, euler::JminBc::Wall);
        double worst = 0.0; // largest asymmetry, relative to the largest |R|
        for (euler::Order o : orders) {
            euler::compute_residual(m, f, o);
            double scale = 0.0;
            double diff = 0.0;
            for (int j = 0; j < m.njc; ++j) {
                for (int i = 0; i < m.nic; ++i) {
                    const int a = m.cell(i, j);
                    const int b = m.cell(m.nic - 1 - i, j);
                    const double sgn[NEQ] = {1.0, 1.0, -1.0, 1.0}; // R2 (y momentum) flips sign
                    for (int k = 0; k < NEQ; ++k) {
                        scale = std::max(scale, std::abs(f.R[k][a]));
                        diff = std::max(diff, std::abs(f.R[k][a] - sgn[k] * f.R[k][b]));
                    }
                }
            }
            worst = std::max(worst, diff / scale);
        }
        check(worst < 1e-12, "E7 symmetric field -> mirror-symmetric residual (1st+2nd): max rel. asymmetry %.2e", worst);
    }

    // Info: farfield regimes taken on a free-stream field.
    for (double M : {0.5, 1.5}) {
        const euler::Freestream fsM = euler::make_freestream(M, 1.25);
        euler::init_freestream(f, fsM);
        int count[4] = {0, 0, 0, 0}; // sub in, sub out, super in, super out
        for (int i = 0; i < m.nic; ++i) {
            const int fc = m.node(i, m.nj - 1);
            const double V = fsM.q[1] * m.njx[fc] + fsM.q[2] * m.njy[fc]; // normal velocity, n outward
            const bool super = std::abs(V) >= 1.0;                       // c_inf = 1
            ++count[(super ? 2 : 0) + (V < 0.0 ? 0 : 1)];
        }
        std::printf("  [INFO] M = %.1f farfield faces: subsonic in %d, subsonic out %d, supersonic in %d, supersonic out %d\n",
                    M, count[0], count[1], count[2], count[3]);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x> [<mesh.x> ...]\n", argv[0]);
        return 1;
    }
    try {
        test_states();
        std::printf("\n=== E. Boundary conditions on the grids\n");
        for (int a = 1; a < argc; ++a) {
            euler::Mesh m = euler::read_plot3d(argv[a]);
            euler::compute_metrics(m);
            euler::Flow f;
            euler::allocate_flow(m, f);
            std::printf("--- %s (%d x %d nodes)\n", argv[a], m.ni, m.nj);
            test_grid(m, f);
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
