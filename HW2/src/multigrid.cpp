#include "multigrid.h"

#include "blockmat.h"
#include "jacobian.h"
#include "linsolve.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace euler {

    // ===================================================================== grids and transfer operators

    bool can_coarsen(const Mesh& m)
    {
        return m.nic % 2 == 0 && m.njc % 2 == 0 && m.nic / 2 >= 8 && m.njc / 2 >= 4;
    }

    Mesh coarsen_mesh(const Mesh& m)
    {
        if (!can_coarsen(m)) {
            throw std::runtime_error("coarsen_mesh: grid cannot be coarsened (odd or too few cells)");
        }
        Mesh c;
        c.ni = (m.ni - 1) / 2 + 1;  // coarse nodes in i
        c.nj = (m.nj - 1) / 2 + 1;  // coarse nodes in j
        c.nic = c.ni - 1;
        c.njc = c.nj - 1;
        c.ni_tot = c.ni + 2 * NG;
        c.nj_tot = c.nj + 2 * NG;
        c.nic_tot = c.nic + 2 * NG;
        c.njc_tot = c.njc + 2 * NG;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        c.x.assign(static_cast<size_t>(c.ni_tot) * c.nj_tot, nan);
        c.y.assign(c.x.size(), nan);
        for (int J = 0; J < c.nj; ++J) {
            for (int I = 0; I < c.ni; ++I) {
                c.x[c.node(I, J)] = m.x[m.node(2 * I, 2 * J)]; // every other node (slide 41)
                c.y[c.node(I, J)] = m.y[m.node(2 * I, 2 * J)];
            }
        }
        compute_metrics(c);
        return c;
    }

    void restrict_volume_average(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine,
                                 std::vector<double>& coarse)
    {
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                double sv = 0.0;  // sum of the fine areas
                double svw = 0.0; // sum of area * value
                for (int b = 0; b < 2; ++b) {
                    for (int a = 0; a < 2; ++a) {
                        const int cf = mf.cell(2 * I + a, 2 * J + b);
                        sv += mf.vol[cf];
                        svw += mf.vol[cf] * fine[cf];
                    }
                }
                coarse[mc.cell(I, J)] = svw / sv;
            }
        }
    }

    void restrict_sum(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine, std::vector<double>& coarse)
    {
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                double s = 0.0; // sum of the 4 fine values
                for (int b = 0; b < 2; ++b) {
                    for (int a = 0; a < 2; ++a) {
                        s += fine[mf.cell(2 * I + a, 2 * J + b)];
                    }
                }
                coarse[mc.cell(I, J)] = s;
            }
        }
    }

    void prolong_bilinear(const Mesh& mf, const Mesh& mc, const std::vector<double>& coarse, std::vector<double>& fine)
    {
        // Coarse value with the boundary rules: periodic in I across the cut, J clamped at wall / farfield.
        auto cval = [&](int I, int J) {
            I = (I + mc.nic) % mc.nic;
            J = std::min(std::max(J, 0), mc.njc - 1);
            return coarse[mc.cell(I, J)];
        };
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 0; i < mf.nic; ++i) {
                const int I = i / 2;                       // coarse cell containing fine cell (i, j)
                const int J = j / 2;
                const int In = (i % 2 == 0) ? I - 1 : I + 1; // nearest coarse neighbour in i
                const int Jn = (j % 2 == 0) ? J - 1 : J + 1; // nearest coarse neighbour in j
                // Slide 50: 9/16 own cell, 3/16 the two side neighbours, 1/16 the corner neighbour.
                fine[mf.cell(i, j)] = (9.0 * cval(I, J) + 3.0 * cval(In, J) + 3.0 * cval(I, Jn) + cval(In, Jn)) / 16.0;
            }
        }
    }

    std::vector<Mesh> coarse_hierarchy(const Mesh& m, int n_levels)
    {
        std::vector<Mesh> h;
        Mesh cur = m;
        for (int l = 1; l < n_levels; ++l) {
            cur = coarsen_mesh(cur);
            h.push_back(cur);
        }
        std::reverse(h.begin(), h.end()); // coarsest first
        return h;
    }

    void inject_piecewise_constant(const Mesh& mf, const Mesh& mc, const Flow& coarse, Flow& fine)
    {
        for (int k = 0; k < NEQ; ++k) {
            for (int j = 0; j < mf.njc; ++j) {
                for (int i = 0; i < mf.nic; ++i) {
                    fine.W[k][mf.cell(i, j)] = coarse.W[k][mc.cell(i / 2, j / 2)];
                }
            }
        }
    }

    int multigrid_levels(const Mesh& m, const MultigridSettings& s)
    {
        int n = 1;    // levels
        int nic = m.nic;
        int njc = m.njc;
        while ((s.levels <= 0 || n < s.levels) && nic % 2 == 0 && njc % 2 == 0 && nic / 2 >= 8 && njc / 2 >= 4) {
            nic /= 2;
            njc /= 2;
            ++n;
        }
        return n;
    }

    // ===================================================================== the cycle

    namespace {

        // One grid level of the hierarchy.
        struct Level {
            Mesh m;                                   // grid of the level
            Flow f;                                   // solution, residual (cell arrays)
            std::array<std::vector<double>, NEQ> Q;   // forcing function (0 on the finest level)
            std::array<std::vector<double>, NEQ> W0;  // solution right after the restriction
            std::array<std::vector<double>, NEQ> Ws;  // W(0) of the multistage stages
            std::vector<double> dt;                   // local time steps
            Order order = Order::First;               // spatial order of the level
            RkScheme rk;                              // multistage scheme of the level
            BlockMatrix M;                            // first-order Jacobian (SGS smoother)
            std::vector<double> Dinv;                 // its inverted diagonal blocks
        };

        // R + Q on the level (ghosts refreshed first).
        void residual_with_forcing(Level& L, const Freestream& fs, JminBc jmin)
        {
            apply_bc(L.m, L.f, fs, jmin);
            compute_residual(L.m, L.f, L.order);
            for (int k = 0; k < NEQ; ++k) {
                for (int j = 0; j < L.m.njc; ++j) {
                    for (int i = 0; i < L.m.nic; ++i) {
                        const int c = L.m.cell(i, j);
                        L.f.R[k][c] += L.Q[k][c];
                    }
                }
            }
        }

        /*
           Largest omega <= 1 such that W + omega dW keeps rho and p within a relative change of
           max_change in every cell (dW: cell arrays). Same safeguard as the implicit solver.
        */
        double change_factor(const Level& L, const std::array<std::vector<double>, NEQ>& dW, double max_change)
        {
            double omega = 1.0;
            for (int j = 0; j < L.m.njc; ++j) {
                for (int i = 0; i < L.m.nic; ++i) {
                    const int c = L.m.cell(i, j);
                    double w[NEQ];
                    double w1[NEQ];
                    double q[NEQ];
                    double q1[NEQ];
                    for (int k = 0; k < NEQ; ++k) {
                        w[k] = L.f.W[k][c];
                        w1[k] = w[k] + dW[k][c];
                    }
                    cons_to_prim(w, q);
                    cons_to_prim(w1, q1);
                    const double drho = std::abs(q1[0] - q[0]);
                    const double dp = std::abs(q1[3] - q[3]);
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

        // n multistage iterations on R + Q (slide 47).
        void smooth_rk(Level& L, const Freestream& fs, const MultigridSettings& s, int n)
        {
            for (int it = 0; it < n; ++it) {
                residual_with_forcing(L, fs, s.jmin);
                local_time_step(L.m, L.f, L.rk.cfl, L.dt);
                for (int k = 0; k < NEQ; ++k) {
                    L.Ws[k] = L.f.W[k];
                }
                for (int st = 0; st < L.rk.stages; ++st) {
                    if (st > 0) {
                        residual_with_forcing(L, fs, s.jmin);
                    }
                    const double a = L.rk.alpha[st]; // alpha_k
                    for (int j = 0; j < L.m.njc; ++j) {
                        for (int i = 0; i < L.m.nic; ++i) {
                            const int c = L.m.cell(i, j);
                            const double fac = a * L.dt[c] / L.m.vol[c];
                            for (int k = 0; k < NEQ; ++k) {
                                L.f.W[k][c] = L.Ws[k][c] - fac * L.f.R[k][c];
                            }
                        }
                    }
                }
            }
        }

        // n implicit Euler steps (Omega/dt + J1) dW = -(R + Q), solved by block SGS (slides 25-35).
        void smooth_sgs(Level& L, const Freestream& fs, const MultigridSettings& s, int n)
        {
            const int nc = L.m.nic * L.m.njc;        // physical cells
            const size_t nu = static_cast<size_t>(NEQ) * nc; // unknowns
            std::vector<double> diag(nc);           // Omega / dt
            std::vector<double> rhs(nu);            // -(R + Q), compact
            std::vector<double> x(nu);              // update, compact
            std::array<std::vector<double>, NEQ> dW; // update, cell arrays
            for (int k = 0; k < NEQ; ++k) {
                dW[k].assign(L.f.W[k].size(), 0.0);
            }
            for (int it = 0; it < n; ++it) {
                residual_with_forcing(L, fs, s.jmin); // also fills f.Q for the Jacobian
                local_time_step(L.m, L.f, s.sgs_cfl, L.dt);
                for (int j = 0; j < L.m.njc; ++j) {
                    for (int i = 0; i < L.m.nic; ++i) {
                        const int c = L.m.cell(i, j);
                        const int p = j * L.m.nic + i;
                        diag[p] = L.m.vol[c] / L.dt[c];
                        for (int k = 0; k < NEQ; ++k) {
                            rhs[static_cast<size_t>(p) * NEQ + k] = -L.f.R[k][c];
                        }
                    }
                }
                assemble_first_order_jacobian(L.m, L.f, fs, s.jmin, diag, L.M);
                factor_diagonal(L.M, L.Dinv);
                std::fill(x.begin(), x.end(), 0.0);
                block_sgs(L.M, L.Dinv, rhs, x, s.sgs_sweeps, 1.0);
                for (int j = 0; j < L.m.njc; ++j) {
                    for (int i = 0; i < L.m.nic; ++i) {
                        const int c = L.m.cell(i, j);
                        const int p = j * L.m.nic + i;
                        for (int k = 0; k < NEQ; ++k) {
                            dW[k][c] = x[static_cast<size_t>(p) * NEQ + k];
                        }
                    }
                }
                const double omega = change_factor(L, dW, s.max_change);
                for (int j = 0; j < L.m.njc; ++j) {
                    for (int i = 0; i < L.m.nic; ++i) {
                        const int c = L.m.cell(i, j);
                        for (int k = 0; k < NEQ; ++k) {
                            L.f.W[k][c] += omega * dW[k][c];
                        }
                    }
                }
            }
        }

        void smooth(Level& L, const Freestream& fs, const MultigridSettings& s, int n)
        {
            if (n <= 0) {
                return;
            }
            if (s.smoother == MgSmoother::RK) {
                smooth_rk(L, fs, s, n);
            } else {
                smooth_sgs(L, fs, s, n);
            }
        }

        // Restriction from level l to level l+1: solution, then forcing function (slides 45-46, 51).
        void go_down(std::vector<Level>& L, int l, const Freestream& fs, const MultigridSettings& s)
        {
            Level& F = L[l];
            Level& C = L[l + 1];
            residual_with_forcing(F, fs, s.jmin); // F.f.R = R_l + Q_l
            for (int k = 0; k < NEQ; ++k) {
                restrict_volume_average(F.m, C.m, F.f.W[k], C.f.W[k]);
                C.W0[k] = C.f.W[k];
                restrict_sum(F.m, C.m, F.f.R[k], C.Q[k]); // I (R_l + Q_l), forcing completed below
            }
            apply_bc(C.m, C.f, fs, s.jmin);
            compute_residual(C.m, C.f, C.order); // R_{l+1}(W_{l+1}^(0))
            for (int k = 0; k < NEQ; ++k) {
                for (int j = 0; j < C.m.njc; ++j) {
                    for (int i = 0; i < C.m.nic; ++i) {
                        const int c = C.m.cell(i, j);
                        C.Q[k][c] -= C.f.R[k][c];
                    }
                }
            }
        }

        // Prolongation of the correction from level l+1 to level l (slides 48, 50).
        void come_up(std::vector<Level>& L, int l, const MultigridSettings& s)
        {
            Level& F = L[l];
            Level& C = L[l + 1];
            std::array<std::vector<double>, NEQ> dW; // prolonged correction, fine cell arrays
            std::vector<double> dc(C.f.W[0].size()); // coarse correction of one variable
            for (int k = 0; k < NEQ; ++k) {
                for (int j = 0; j < C.m.njc; ++j) {
                    for (int i = 0; i < C.m.nic; ++i) {
                        const int c = C.m.cell(i, j);
                        dc[c] = C.f.W[k][c] - C.W0[k][c];
                    }
                }
                dW[k].assign(F.f.W[k].size(), 0.0);
                prolong_bilinear(F.m, C.m, dc, dW[k]);
            }
            const double omega = change_factor(F, dW, s.max_change);
            for (int k = 0; k < NEQ; ++k) {
                for (int j = 0; j < F.m.njc; ++j) {
                    for (int i = 0; i < F.m.nic; ++i) {
                        const int c = F.m.cell(i, j);
                        F.f.W[k][c] += omega * dW[k][c];
                    }
                }
            }
        }

        // One cycle starting on level l, with level last as the coarsest (gamma = 1 V, 2 W).
        void cycle(std::vector<Level>& L, int l, int last, const Freestream& fs, const MultigridSettings& s)
        {
            if (l == last) {
                smooth(L[l], fs, s, s.n_coarse);
                return;
            }
            smooth(L[l], fs, s, s.n_pre);
            go_down(L, l, fs, s);
            for (int g = 0; g < s.gamma; ++g) {
                cycle(L, l + 1, last, fs, s);
            }
            come_up(L, l, s);
            smooth(L[l], fs, s, s.n_post);
        }

    }

    RunResult run_multigrid(const Mesh& m, Flow& f, const Freestream& fs, const MultigridSettings& s)
    {
        const int nlev = multigrid_levels(m, s); // number of levels
        std::vector<Level> L(nlev);
        L[0].m = m;
        L[0].f = f;
        for (int l = 1; l < nlev; ++l) {
            L[l].m = coarsen_mesh(L[l - 1].m);
            allocate_flow(L[l].m, L[l].f);
            init_freestream(L[l].f, fs);
        }
        for (int l = 0; l < nlev; ++l) {
            L[l].order = l == 0 ? s.fine_order : s.coarse_order;
            L[l].rk = l == 0 ? s.fine_rk : s.coarse_rk;
            for (int k = 0; k < NEQ; ++k) {
                L[l].Q[k].assign(L[l].f.W[k].size(), 0.0);
                L[l].W0[k].assign(L[l].f.W[k].size(), 0.0);
                L[l].Ws[k].assign(L[l].f.W[k].size(), 0.0);
            }
        }
        if (s.print_every > 0) {
            std::printf("  multigrid: %d levels (finest %d x %d, coarsest %d x %d nodes), %s-cycle, smoother %s\n", nlev,
                        L[0].m.ni, L[0].m.nj, L[nlev - 1].m.ni, L[nlev - 1].m.nj, s.gamma == 1 ? "V" : "W",
                        s.smoother == MgSmoother::RK ? "RK" : "SGS");
        }
        WALL_EXTRAP_ACTIVE = true;

        // Full multigrid start (slide 42): cycles on each coarse level, coarsest first, then inject.
        if (s.fmg_cycles > 0 && nlev > 1) {
            for (int l = nlev - 1; l >= 1; --l) {
                for (int c = 0; c < s.fmg_cycles; ++c) {
                    cycle(L, l, nlev - 1, fs, s); // level l is the finest of this sub-cycle (Q_l = 0)
                }
                inject_piecewise_constant(L[l - 1].m, L[l].m, L[l].f, L[l - 1].f);
                if (s.print_every > 0) {
                    std::printf("  FMG: %d cycles on %d x %d, injected into %d x %d\n", s.fmg_cycles, L[l].m.ni, L[l].m.nj,
                                L[l - 1].m.ni, L[l - 1].m.nj);
                }
            }
        }

        RunResult res;
        Monitor mon(s.max_iter, s.tol_drop, s.cauchy_window, s.cauchy_tol, s.print_every, s.history_file);
        for (int it = 1; it <= s.max_iter; ++it) {
            // Residual of the finest level (Q_0 = 0) and forces, then one cycle.
            apply_bc(L[0].m, L[0].f, fs, s.jmin);
            compute_residual(L[0].m, L[0].f, L[0].order);
            const double norm = density_residual_norm(L[0].m, L[0].f);
            if (mon.record(it, norm, compute_forces(L[0].m, L[0].f, fs, s.ref), res)) {
                break;
            }
            cycle(L, 0, nlev - 1, fs, s);
        }
        f = L[0].f;
        apply_bc(m, f, fs, s.jmin);
        return res;
    }

}
