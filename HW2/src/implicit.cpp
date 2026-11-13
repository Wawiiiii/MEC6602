#include "implicit.h"

#include "blockmat.h"
#include "jacobian.h"
#include "linsolve.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace euler {

    // Physical cells of a cell array -> compact vector [p * NEQ + k], p = j * nic + i.
    static void gather(const Mesh& m, const std::array<std::vector<double>, NEQ>& a, std::vector<double>& x)
    {
        x.resize(static_cast<size_t>(NEQ) * m.nic * m.njc);
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);    // cell array index
                const int p = j * m.nic + i;   // compact index
                for (int k = 0; k < NEQ; ++k) {
                    x[static_cast<size_t>(p) * NEQ + k] = a[k][c];
                }
            }
        }
    }

    // a (physical cells) = base (physical cells) + s * x (compact vector).
    static void set_plus(const Mesh& m, const std::array<std::vector<double>, NEQ>& base, double s,
                         const std::vector<double>& x, std::array<std::vector<double>, NEQ>& a)
    {
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const int p = j * m.nic + i;
                for (int k = 0; k < NEQ; ++k) {
                    a[k][c] = base[k][c] + s * x[static_cast<size_t>(p) * NEQ + k];
                }
            }
        }
    }

    /*
       Largest factor omega <= 1 such that W + omega dW changes rho and p by at most max_change
       (relative) in every cell and keeps them positive. Startup safeguard only: near convergence
       the updates are small and omega = 1 (full Newton / implicit step).
    */
    static double update_factor(const Mesh& m, const Flow& f, const std::vector<double>& dW, double max_change)
    {
        double omega = 1.0;
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const size_t p = static_cast<size_t>(j * m.nic + i) * NEQ;
                double w[NEQ];  // current state
                double w1[NEQ]; // updated state
                double q[NEQ];
                double q1[NEQ];
                for (int k = 0; k < NEQ; ++k) {
                    w[k] = f.W[k][c];
                    w1[k] = w[k] + dW[p + k];
                }
                cons_to_prim(w, q);
                cons_to_prim(w1, q1);
                const double drho = std::abs(q1[0] - q[0]); // density change
                const double dp = std::abs(q1[3] - q[3]);   // pressure change (full step)
                if (drho > max_change * q[0]) {
                    omega = std::min(omega, max_change * q[0] / drho);
                }
                if (dp > max_change * q[3] || !(q1[3] > 0.0)) {
                    omega = std::min(omega, max_change * q[3] / std::max(dp, 1e-300));
                }
            }
        }
        return omega;
    }

    RunResult run_implicit(const Mesh& m, Flow& f, const Freestream& fs, const ImplicitSettings& s)
    {
        RunResult res;
        const int nc = m.nic * m.njc;                 // physical cells
        const size_t n = static_cast<size_t>(NEQ) * nc; // unknowns
        const double gamma_row = 1.0 + 1.0 / std::sqrt(2.0); // ROS2 gamma

        Flow fp = f;                 // scratch flow for perturbed states (matrix-free products, ROW stage 2)
        BlockMatrix M;               // first-order implicit operator
        std::vector<double> Dfact;   // factorised diagonal blocks
        std::vector<double> diag(nc); // c Omega / dt of each cell
        std::vector<double> dt;      // local time steps (cell array)
        std::vector<double> Rn(n);   // R(W^n), compact
        std::vector<double> Wn(n);   // W^n, compact (only its norm is used)
        std::vector<double> rhs(n);  // right-hand side
        std::vector<double> dW(n);   // update
        std::vector<double> dWold(n, 0.0); // previous update (BDF2)
        std::vector<double> k1(n);   // ROW stage 1
        std::vector<double> k2(n);   // ROW stage 2
        std::vector<double> Rp(n);   // perturbed residual, compact
        Monitor mon(s.max_iter, s.tol_drop, s.cauchy_window, s.cauchy_tol, s.print_every, s.history_file);
        WALL_EXTRAP_ACTIVE = WALL_EXTRAP_SWITCH >= 1.0; // extrapolated wall ghost from the start or later (bc.h)

        for (int it = 1; it <= s.max_iter; ++it) {
            // Residual, forces, convergence check of W^n.
            apply_bc(m, f, fs, s.jmin);
            compute_residual(m, f, s.order);
            const double norm = density_residual_norm(m, f);
            if (mon.record(it, norm, compute_forces(m, f, fs, s.ref), res)) {
                if (!WALL_EXTRAP_ACTIVE && WALL_PRESSURE_EXTRAPOLATE && !res.diverged) {
                    std::printf("  (stopped before the extrapolated wall pressure was switched on)\n");
                }
                break;
            }
            if (WALL_PRESSURE_EXTRAPOLATE && !WALL_EXTRAP_ACTIVE && norm <= WALL_EXTRAP_SWITCH * res.res0) {
                WALL_EXTRAP_ACTIVE = true; // from the next residual on (the CFL ramp follows the residual jump)
                if (s.print_every > 0) {
                    std::printf("           extrapolated wall pressure switched on\n");
                }
            }

            // CFL ramp (slide 40) and local time step.
            const double cfl = std::min(s.cfl_max, std::max(s.cfl0, s.cfl0 * res.res0 / norm));
            local_time_step(m, f, cfl, dt);

            // Time-term coefficient c of c Omega / dt.
            double c = 1.0; // BDF1
            if (s.scheme == TimeScheme::BDF2 && it > 1) {
                c = 1.5;
            } else if (s.scheme == TimeScheme::ROW) {
                c = 1.0 / gamma_row;
            }
            for (int j = 0; j < m.njc; ++j) {
                for (int i = 0; i < m.nic; ++i) {
                    const int cc = m.cell(i, j);
                    diag[j * m.nic + i] = c * m.vol[cc] / dt[cc];
                }
            }

            // First-order Jacobian + time term (slides 26-29), diagonal factorisation.
            assemble_first_order_jacobian(m, f, fs, s.jmin, diag, M);
            factor_diagonal(M, Dfact);

            gather(m, f.R, Rn);
            gather(m, f.W, Wn);
            const double normW = norm2(Wn);

            /*
            Exact operator for GMRES: A v = c Omega/dt v + J v, with the Jacobian-free product
            J v ~ (R(W^n + eps v) - R(W^n)) / eps, eps = sqrt(eps_machine) (1 + ||W||) / ||v|| (slide 40).
            */
            auto apply_A = [&](const std::vector<double>& v, std::vector<double>& y) {
                y.resize(n);
                const double nv = norm2(v);
                if (nv == 0.0) {
                    std::fill(y.begin(), y.end(), 0.0);
                    return;
                }
                const double eps = 1.4901161193847656e-08 * (1.0 + normW) / nv; // sqrt(2^-52) (1 + ||W||) / ||v||
                for (int k = 0; k < NEQ; ++k) {
                    fp.W[k] = f.W[k];
                }
                set_plus(m, f.W, eps, v, fp.W);
                apply_bc(m, fp, fs, s.jmin);
                compute_residual(m, fp, s.order);
                gather(m, fp.R, Rp);
#pragma omp parallel for schedule(static)
                for (int p = 0; p < nc; ++p) {
                    for (int k = 0; k < NEQ; ++k) {
                        const size_t a = static_cast<size_t>(p) * NEQ + k;
                        y[a] = diag[p] * v[a] + (Rp[a] - Rn[a]) / eps;
                    }
                }
            };
            // Preconditioner: SGS sweeps on the first-order operator, from zero.
            auto apply_Minv = [&](const std::vector<double>& r, std::vector<double>& z) {
                z.assign(n, 0.0);
                block_sgs(M, Dfact, r, z, s.sgs_sweeps, s.sgs_omega);
            };
            int lin_its = 0;       // linear iterations of this nonlinear iteration
            double lin_res = 0.0;  // last linear relative residual
            auto solve = [&](const std::vector<double>& b, std::vector<double>& x) {
                x.assign(n, 0.0);
                if (s.linsolve == LinearSolver::GMRES) {
                    const GmresResult g = gmres(apply_A, apply_Minv, b, x, s.gmres_restart, s.gmres_max_iter, s.gmres_rtol);
                    lin_its += g.iterations;
                    lin_res = g.rel_residual;
                } else {
                    block_sgs(M, Dfact, b, x, s.sgs_sweeps, s.sgs_omega);
                    lin_its += s.sgs_sweeps;
                }
            };

            // Right-hand side(s) and update.
            if (s.scheme == TimeScheme::ROW) {
                for (size_t a = 0; a < n; ++a) {
                    rhs[a] = -Rn[a] / gamma_row;
                }
                solve(rhs, k1);
                // R(W^n + k1) for the second stage.
                for (int k = 0; k < NEQ; ++k) {
                    fp.W[k] = f.W[k];
                }
                set_plus(m, f.W, 1.0, k1, fp.W);
                apply_bc(m, fp, fs, s.jmin);
                compute_residual(m, fp, s.order);
                gather(m, fp.R, Rp);
                for (int p = 0; p < nc; ++p) {
                    for (int k = 0; k < NEQ; ++k) {
                        const size_t a = static_cast<size_t>(p) * NEQ + k;
                        rhs[a] = -Rp[a] / gamma_row - 2.0 * diag[p] * k1[a]; // 2 Omega/(gamma dt) = 2 diag
                    }
                }
                solve(rhs, k2);
                for (size_t a = 0; a < n; ++a) {
                    dW[a] = 1.5 * k1[a] + 0.5 * k2[a];
                }
            } else {
                for (int p = 0; p < nc; ++p) {
                    for (int k = 0; k < NEQ; ++k) {
                        const size_t a = static_cast<size_t>(p) * NEQ + k;
                        rhs[a] = -Rn[a];
                        if (s.scheme == TimeScheme::BDF2 && it > 1) {
                            rhs[a] += (diag[p] / 3.0) * dWold[a]; // Omega/(2 dt) = (3/2 Omega/dt) / 3
                        }
                    }
                }
                solve(rhs, dW);
            }

            // Positivity / size limiter, then update W.
            const double omega = update_factor(m, f, dW, s.max_change);
            for (size_t a = 0; a < n; ++a) {
                dW[a] *= omega;
            }
            set_plus(m, f.W, 1.0, dW, f.W);
            dWold = dW;

            if (s.print_every > 0 && (it == 1 || it % s.print_every == 0)) {
                std::printf("           CFL %.3e  linear its %d  lin. rel. res %.2e  update factor %.3f\n", cfl, lin_its,
                            lin_res, omega);
            }
        }
        apply_bc(m, f, fs, s.jmin);
        return res;
    }

}
