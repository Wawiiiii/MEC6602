#include "flow.h"

#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace euler {

    void allocate_flow(const Mesh& m, Flow& f)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const size_t ncell = static_cast<size_t>(m.nic_tot) * m.njc_tot; // stored cells, ghosts included
        const size_t nnode = static_cast<size_t>(m.ni_tot) * m.nj_tot;   // stored nodes = face slots

        for (int k = 0; k < NEQ; ++k) {
            f.W[k].assign(ncell, nan);
            f.Q[k].assign(ncell, nan);
            f.R[k].assign(ncell, nan);
            f.Fi[k].assign(nnode, 0.0);
            f.Fj[k].assign(nnode, 0.0);
        }
    }

    void compute_primitives(Flow& f)
    {
        const int ncell = static_cast<int>(f.W[0].size()); // stored cells, ghosts included

        // Every cell is independent: OpenMP parallel loop (slides 68-70).
#pragma omp parallel for schedule(static)
        for (int c = 0; c < ncell; ++c) {
            const double w[NEQ] = {f.W[0][c], f.W[1][c], f.W[2][c], f.W[3][c]}; // conservative state of cell c
            double q[NEQ];                                                       // primitive state of cell c
            cons_to_prim(w, q);
            for (int k = 0; k < NEQ; ++k) {
                f.Q[k][c] = q[k];
            }
        }
    }

    void set_threads(int n)
    {
#ifdef _OPENMP
        if (n > 0) {
            omp_set_num_threads(n);
        }
#else
        (void)n;
#endif
    }

    int get_threads()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

}
