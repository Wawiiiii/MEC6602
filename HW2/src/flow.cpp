#include "flow.h"

#include <limits>

namespace euler {

    void allocate_flow(const Mesh& m, Flow& f)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const size_t ncell = static_cast<size_t>(m.nic_tot) * m.njc_tot; // stored cells, ghosts included

        for (int k = 0; k < NEQ; ++k) {
            f.W[k].assign(ncell, nan);
            f.Q[k].assign(ncell, nan);
            f.R[k].assign(ncell, nan);
        }
    }

    void compute_primitives(Flow& f)
    {
        const size_t ncell = f.W[0].size(); // stored cells, ghosts included

        for (size_t c = 0; c < ncell; ++c) {
            const double w[NEQ] = {f.W[0][c], f.W[1][c], f.W[2][c], f.W[3][c]}; // conservative state of cell c
            double q[NEQ];                                                       // primitive state of cell c
            cons_to_prim(w, q);
            for (int k = 0; k < NEQ; ++k) {
                f.Q[k][c] = q[k];
            }
        }
    }

}
