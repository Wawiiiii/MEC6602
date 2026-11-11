/*
   Step 2 check (Euler2D slides 21-22).
     1. Small Cartesian mesh (3x3 nodes, unit spacing): print every cell area and
        every face with its left/right cells, normal and length. All known by hand:
        area 1, length 1, i-face normal (1, 0), j-face normal (0, 1).
     2. Airfoil mesh: compute the metrics, write the mesh back with the areas, and check
          - every area > 0 (compute_metrics throws otherwise),
          - sum of the areas = area between the farfield and the airfoil (shoelace),
          - |n| = 1 on every face,
          - constant flux on every face -> residual = 0 in every cell
            (the integral of the normal over a closed surface is 0, slide 22).
   Usage: test_metrics <mesh.x> <out.dat>
*/

#include "mesh.h"
#include "output.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <vector>

// Build a Cartesian mesh by hand (not an O-grid, so it cannot come from read_plot3d).
static euler::Mesh cartesian_mesh(int ni, int nj)
{
    using euler::NG;
    euler::Mesh m;
    m.ni = ni;
    m.nj = nj;
    m.nic = ni - 1;
    m.njc = nj - 1;
    m.ni_tot = ni + 2 * NG;
    m.nj_tot = nj + 2 * NG;
    m.nic_tot = m.nic + 2 * NG;
    m.njc_tot = m.njc + 2 * NG;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    m.x.assign(static_cast<size_t>(m.ni_tot) * m.nj_tot, nan);
    m.y.assign(m.x.size(), nan);
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            m.x[m.node(i, j)] = i; // x = i
            m.y[m.node(i, j)] = j; // y = j
        }
    }
    return m;
}

/*
   Residual of a constant flux vector (fx, fy) through every face, computed like the
   real solver will do it (Euler2D_basic_algo slide 6): one flux per face, added to the
   left cell and subtracted from the right cell (what leaves one enters the other).
   Faces on the boundary also feed a ghost cell; those slots are simply ignored.
*/
static std::vector<double> constant_flux_residual(const euler::Mesh& m, double fx, double fy)
{
    std::vector<double> R(m.vol.size(), 0.0);

    for (int j = 0; j < m.njc; ++j) { // i-faces
        for (int i = 0; i < m.ni; ++i) {
            const int f = m.node(i, j);
            const double flux = (fx * m.nix[f] + fy * m.niy[f]) * m.si[f]; // (F . n) Delta S
            R[m.cell(i - 1, j)] += flux; // left cell
            R[m.cell(i, j)] -= flux;     // right cell
        }
    }
    for (int j = 0; j < m.nj; ++j) { // j-faces
        for (int i = 0; i < m.nic; ++i) {
            const int f = m.node(i, j);
            const double flux = (fx * m.njx[f] + fy * m.njy[f]) * m.sj[f];
            R[m.cell(i, j - 1)] += flux; // left cell
            R[m.cell(i, j)] -= flux;     // right cell
        }
    }
    return R;
}

// Largest |R| over the physical cells, divided by the cell perimeter (so ~1e-16 is perfect).
static double max_relative_residual(const euler::Mesh& m, const std::vector<double>& R)
{
    double worst = 0.0;
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double perimeter = m.si[m.node(i, j)] + m.si[m.node(i + 1, j)]
                                   + m.sj[m.node(i, j)] + m.sj[m.node(i, j + 1)];
            worst = std::max(worst, std::abs(R[m.cell(i, j)]) / perimeter);
        }
    }
    return worst;
}

// Area enclosed by the closed node line j (shoelace formula), absolute value.
static double enclosed_area(const euler::Mesh& m, int j)
{
    double a = 0.0;
    for (int i = 0; i < m.ni - 1; ++i) { // node ni-1 = node 0 closes the polygon
        const int p = m.node(i, j);
        const int q = m.node(i + 1, j);
        a += m.x[p] * m.y[q] - m.x[q] * m.y[p];
    }
    return std::abs(0.5 * a);
}

static void cartesian_check()
{
    euler::Mesh m = cartesian_mesh(3, 3);
    euler::compute_metrics(m);

    std::printf("=== Cartesian 3x3 nodes (2x2 cells), unit spacing\n");
    std::printf("cells (expect area 1):\n");
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            std::printf("  cell (%d,%d)  area %.3f\n", i, j, m.vol[m.cell(i, j)]);
        }
    }
    std::printf("i-faces (expect n = (1, 0), dS = 1):\n");
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.ni; ++i) {
            const int f = m.node(i, j);
            std::printf("  i-face (%d,%d)  left cell (%2d,%d) -> right cell (%d,%d)  n = (%5.2f, %5.2f)  dS = %.3f\n",
                        i, j, i - 1, j, i, j, m.nix[f], m.niy[f], m.si[f]);
        }
    }
    std::printf("j-faces (expect n = (0, 1), dS = 1):\n");
    for (int j = 0; j < m.nj; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const int f = m.node(i, j);
            std::printf("  j-face (%d,%d)  left cell (%d,%2d) -> right cell (%d,%d)  n = (%5.2f, %5.2f)  dS = %.3f\n",
                        i, j, i, j - 1, i, j, m.njx[f], m.njy[f], m.sj[f]);
        }
    }
    const double rx = max_relative_residual(m, constant_flux_residual(m, 1.0, 0.0));
    const double ry = max_relative_residual(m, constant_flux_residual(m, 0.0, 1.0));
    std::printf("constant flux residual: F = (1,0): %.2e   F = (0,1): %.2e   (expect 0)\n\n", rx, ry);
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <mesh.x> <out.dat>\n", argv[0]);
        return 1;
    }

    try {
        cartesian_check();

        euler::Mesh m = euler::read_plot3d(argv[1]);
        euler::compute_metrics(m);
        std::printf("=== %s : %d x %d nodes\n", argv[1], m.ni, m.nj);

        // Areas: range, and total against the shoelace area of the ring.
        double amin = 1e300, amax = 0.0, asum = 0.0;
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const double a = m.vol[m.cell(i, j)];
                amin = std::min(amin, a);
                amax = std::max(amax, a);
                asum += a;
            }
        }
        const double ring = enclosed_area(m, m.nj - 1) - enclosed_area(m, 0);
        std::printf("cell area      : min %.3e  max %.3e  (all > 0)\n", amin, amax);
        std::printf("sum of areas   : %.12e\n", asum);
        std::printf("farfield-wall  : %.12e   (shoelace, rel. diff %.2e)\n", ring, std::abs(asum - ring) / ring);
        std::printf("airfoil area   : %.6f   (shoelace on the wall line)\n", enclosed_area(m, 0));

        // Unit normals: largest deviation of |n| from 1 on all physical faces.
        double ndev = 0.0;
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.ni; ++i) {
                const int f = m.node(i, j);
                ndev = std::max(ndev, std::abs(std::hypot(m.nix[f], m.niy[f]) - 1.0));
            }
        }
        for (int j = 0; j < m.nj; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int f = m.node(i, j);
                ndev = std::max(ndev, std::abs(std::hypot(m.njx[f], m.njy[f]) - 1.0));
            }
        }
        std::printf("max | |n| - 1 |: %.2e\n", ndev);

        // Slide 22 item 5: constant flux -> zero residual, for both flux directions.
        const double rx = max_relative_residual(m, constant_flux_residual(m, 1.0, 0.0));
        const double ry = max_relative_residual(m, constant_flux_residual(m, 0.0, 1.0));
        std::printf("constant flux residual / perimeter: F = (1,0): %.2e   F = (0,1): %.2e   (expect ~1e-16)\n", rx, ry);

        // Small meshes: wall normals, they must point away from the airfoil into the fluid.
        if (m.ni <= 17) {
            std::printf("\nwall j-faces (j = 0): midpoint and normal\n");
            std::printf("   i        xm        ym        nx        ny        dS\n");
            for (int i = 0; i < m.nic; ++i) {
                const int a = m.node(i, 0);
                const int b = m.node(i + 1, 0);
                std::printf("%4d %9.5f %9.5f %9.5f %9.5f %9.5f\n", i, 0.5 * (m.x[a] + m.x[b]),
                            0.5 * (m.y[a] + m.y[b]), m.njx[a], m.njy[a], m.sj[a]);
            }
        }

        // Slide 22 item 4: write the mesh back with the cell areas.
        euler::write_tecplot_mesh(m, argv[2]);
        std::printf("\nwrote %s\n", argv[2]);
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
