#include "solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace euler {

    // Runge-Kutta coefficients and CFL for the spatial order and number of stages.
    RkScheme rk_scheme(Order spatial, int stages)
    {
        RkScheme s;
        if (stages == 1) {
            s.stages = 1;
            s.alpha[0] = 1.0;
            s.cfl = 0.9; // explicit Euler: CFL lower than 1
            return s;
        }

        // Coefficient table, upwind schemes. Rows: 3, 4, 5 stages; columns: CFL, then alpha_1 .. alpha_m.
        static const double first[3][6] = {
            {1.5, 0.1481, 0.4000, 1.0000, 0.0, 0.0},     // 1st-order scheme, 3 stages
            {2.0, 0.0833, 0.2069, 0.4265, 1.0000, 0.0},  // 1st-order scheme, 4 stages
            {2.5, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000} // 1st-order scheme, 5 stages
        };
        static const double second[3][6] = {
            {0.69, 0.1918, 0.4929, 1.0000, 0.0, 0.0},     // 2nd-order scheme, 3 stages
            {0.92, 0.1084, 0.2602, 0.5052, 1.0000, 0.0},  // 2nd-order scheme, 4 stages
            {1.15, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000} // 2nd-order scheme, 5 stages
        };
        if (stages < 3 || stages > 5) {
            throw std::runtime_error("rk_scheme: stages must be 1, 3, 4 or 5");
        }
        const double* row = (spatial == Order::First ? first : second)[stages - 3]; // table row
        s.stages = stages;
        s.cfl = row[0];
        for (int k = 0; k < stages; ++k) {
            s.alpha[k] = row[1 + k];
        }
        return s;
    }

    // Local time step dt = CFL Omega / (Lambda_I + Lambda_J) of every physical cell.
    void local_time_step(const Mesh& m, const Flow& f, double cfl, std::vector<double>& dt)
    {
        dt.resize(m.vol.size());
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const double rho = f.Q[0][c];                    // density
                const double u = f.Q[1][c];                      // x velocity
                const double v = f.Q[2][c];                      // y velocity
                const double a = std::sqrt(GAMMA * f.Q[3][c] / rho); // speed of sound

                // I direction: i-faces at node(i, j) and node(i+1, j).
                const int fa = m.node(i, j);
                const int fb = m.node(i + 1, j);
                const double nIx = 0.5 * (m.nix[fa] + m.nix[fb]); // averaged i normal, x
                const double nIy = 0.5 * (m.niy[fa] + m.niy[fb]); // averaged i normal, y
                const double sI = 0.5 * (m.si[fa] + m.si[fb]);    // averaged i face length
                const double lamI = (std::abs(u * nIx + v * nIy) + a) * sI; // Lambda_I

                // J direction: j-faces at node(i, j) and node(i, j+1).
                const int fc = m.node(i, j);
                const int fd = m.node(i, j + 1);
                const double nJx = 0.5 * (m.njx[fc] + m.njx[fd]); // averaged j normal, x
                const double nJy = 0.5 * (m.njy[fc] + m.njy[fd]); // averaged j normal, y
                const double sJ = 0.5 * (m.sj[fc] + m.sj[fd]);    // averaged j face length
                const double lamJ = (std::abs(u * nJx + v * nJy) + a) * sJ; // Lambda_J

                dt[c] = cfl * m.vol[c] / (lamI + lamJ);
            }
        }
    }

    /*
       RMS of R_rho / Omega over the physical cells.
       Rows are summed in parallel, then added in order, so the result does not depend on the
       number of threads.
    */
    double density_residual_norm(const Mesh& m, const Flow& f)
    {
        std::vector<double> row(m.njc, 0.0); // sum of (R_rho / Omega)^2 over each j row
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            double s = 0.0;
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const double r = f.R[0][c] / m.vol[c];
                s += r * r;
            }
            row[j] = s;
        }
        double sum = 0.0; // sum of (R_rho / Omega)^2
        for (int j = 0; j < m.njc; ++j) {
            sum += row[j];
        }
        return std::sqrt(sum / (static_cast<double>(m.nic) * m.njc));
    }

    Monitor::Monitor(int max_iter, double tol_drop, int cauchy_window, double cauchy_tol, int print_every,
                     const std::string& history_file)
        : max_iter_(max_iter), tol_drop_(tol_drop), window_(cauchy_window), cauchy_tol_(cauchy_tol),
          print_every_(print_every), t0_(std::chrono::steady_clock::now())
    {
        if (!history_file.empty()) {
            hist_ = std::fopen(history_file.c_str(), "w");
            if (!hist_) {
                throw std::runtime_error("cannot open history file: " + history_file);
            }
            std::fprintf(hist_, "# iteration  ||R_rho||  ||R_rho||/||R_rho||_0  CL  CD  CM  seconds\n");
        }
    }

    Monitor::~Monitor()
    {
        if (hist_) {
            std::fclose(hist_);
        }
    }

    double Monitor::seconds() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    }

    void Monitor::set_extra_columns(const std::string& header)
    {
        extra_header_ = header;
    }

    void Monitor::end_row()
    {
        if (row_open_) {
            std::printf("\n");
            std::fflush(stdout); // show the row at once, also when the output is redirected to a file
            row_open_ = false;
        }
    }

    // Store and print iteration it, then apply the stopping tests.
    bool Monitor::record(int it, double norm, const Forces& F, RunResult& r)
    {
        if (it == 1) {
            r.res0 = norm;
        }
        r.res = norm;
        r.iterations = it;
        r.forces = F;
        r.history.push_back(norm);
        r.cl.push_back(F.cl);
        r.cd.push_back(F.cd);
        r.cm.push_back(F.cm);
        r.seconds = seconds();

        if (hist_) {
            std::fprintf(hist_, "%d %.10e %.10e %.12e %.12e %.12e %.4f\n", it, norm, norm / r.res0, F.cl, F.cd, F.cm,
                         r.seconds);
        }
        // Terminal: one row per printed iteration, under a header printed once.
        if (print_every_ > 0 && (it == 1 || it % print_every_ == 0)) {
            if (!header_printed_) {
                std::printf("\n%7s %11s %10s %14s %14s %14s %9s%s\n", "iter", "||R_rho||", "drop", "CL", "CD", "CM",
                            "time_s", extra_header_.c_str());
                header_printed_ = true;
            }
            std::printf("%7d %11.4e %10.3e %14.6e %14.6e %14.6e %9.2f", it, norm, norm / r.res0, F.cl, F.cd, F.cm,
                        r.seconds);
            if (extra_header_.empty()) {
                std::printf("\n");
                std::fflush(stdout); // see end_row()
            } else {
                row_open_ = true; // the solver adds its columns, then calls end_row()
            }
        }

        // Stopping tests: divergence, residual drop, Cauchy criterion.
        if (!std::isfinite(norm)) {
            r.diverged = true;
            end_row();
            return true;
        }
        if (norm <= tol_drop_ * r.res0 || norm == 0.0) {
            r.converged = true;
            end_row();
            return true;
        }
        // Cauchy criterion: max - min of CL, CD, CM over the last window_ iterations.
        if (window_ > 0 && static_cast<int>(r.cl.size()) >= window_) {
            bool ok = true;
            for (const std::vector<double>* h : {&r.cl, &r.cd, &r.cm}) {
                double lo = 1e300;  // smallest value in the window
                double hi = -1e300; // largest value in the window
                for (size_t k = h->size() - window_; k < h->size(); ++k) {
                    lo = std::min(lo, (*h)[k]);
                    hi = std::max(hi, (*h)[k]);
                }
                ok = ok && (hi - lo < cauchy_tol_);
            }
            if (ok) {
                r.cauchy = true;
                end_row();
                return true;
            }
        }
        // The iteration limit is the solvers' loop bound: iteration max_iter still does its update.
        return false;
    }

    // Explicit multi-stage Runge-Kutta run with local time steps.
    RunResult run_explicit(const Mesh& m, Flow& f, const Freestream& fs, const ExplicitSettings& s)
    {
        RunResult res;
        std::vector<double> dt;                  // local time step of each cell
        std::array<std::vector<double>, NEQ> W0; // W(0) = W^n, kept over the stages
        for (int k = 0; k < NEQ; ++k) {
            W0[k].resize(f.W[k].size());
        }
        Monitor mon(s.max_iter, s.tol_drop, s.cauchy_window, s.cauchy_tol, s.print_every, s.history_file);
        WALL_EXTRAP_ACTIVE = WALL_EXTRAP_SWITCH >= 1.0; // extrapolated wall pressure from the start or later (bc.h)

        for (int it = 1; it <= s.max_iter; ++it) {
            // 1. Residual, time step and forces of W^n (ghosts up to date for the wall pressure).
            apply_bc(m, f, fs, s.jmin);
            compute_residual(m, f, s.order);
            local_time_step(m, f, s.rk.cfl, dt);
            const double norm = density_residual_norm(m, f);
            if (mon.record(it, norm, compute_forces(m, f, fs, s.ref), res)) {
                break;
            }
            // Switch to the extrapolated wall pressure once the residual has dropped enough.
            if (WALL_PRESSURE_EXTRAPOLATE && !WALL_EXTRAP_ACTIVE && norm <= WALL_EXTRAP_SWITCH * res.res0) {
                WALL_EXTRAP_ACTIVE = true; // from the next stage on
            }

            // 2. Stages: W(k) = W(0) - alpha_k dt/Omega R(W(k-1)), W(0) = W^n.
            for (int k = 0; k < NEQ; ++k) {
                W0[k] = f.W[k];
            }
            for (int st = 0; st < s.rk.stages; ++st) {
                if (st > 0) {
                    apply_bc(m, f, fs, s.jmin);
                    compute_residual(m, f, s.order); // R(W(k-1))
                }
                const double a = s.rk.alpha[st]; // alpha_k
#pragma omp parallel for schedule(static)
                for (int j = 0; j < m.njc; ++j) {
                    for (int i = 0; i < m.nic; ++i) {
                        const int c = m.cell(i, j);
                        const double fac = a * dt[c] / m.vol[c]; // alpha_k dt / Omega
                        for (int k = 0; k < NEQ; ++k) {
                            f.W[k][c] = W0[k][c] - fac * f.R[k][c];
                        }
                    }
                }
            }
        }
        // Ghosts consistent with the final W (for outputs).
        apply_bc(m, f, fs, s.jmin);
        return res;
    }

}
