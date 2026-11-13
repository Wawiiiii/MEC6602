#pragma once

#include "flow.h"

#include <vector>

namespace euler {

    /*
       Block-structured sparse matrix of the implicit operator (Roe_scheme_accel slides 28-29):
       4x4 blocks (NEQ = 4), one block row per PHYSICAL cell, five blocks per row:
         D  : diagonal block            (cell (i, j) itself)
         Wb : west  neighbour (i-1, j)  (i = 0     -> i = nic-1, across the TE cut)
         Eb : east  neighbour (i+1, j)  (i = nic-1 -> i = 0,     across the TE cut)
         Sb : south neighbour (i, j-1)  (zero at j = 0: wall, boundary lumped on D, slide 30)
         Nb : north neighbour (i, j+1)  (zero at j = njc-1: farfield, lumped on D)

       Compact numbering of the physical cells (no ghosts): p = j * nic + i.
       A vector of unknowns has size NEQ * nic * njc, entry [p * NEQ + k], k = equation
       (conservative variables rho, rho u, rho v, rho E).
       A block is stored row-major: block(p)[r * NEQ + c] = d R_p,r / d W_q,c  (r = row, c = column).
       Each array has size nic * njc * NEQ * NEQ; block of cell p starts at p * NEQ * NEQ.

       Row p of y = M x:
         y_p = D_p x_p + Wb_p x_west(p) + Eb_p x_east(p) + Sb_p x_south(p) + Nb_p x_north(p)
       (the south / north terms are absent at j = 0 / j = njc-1).
    */
    struct BlockMatrix {
        int nic = 0;            // physical cells in i
        int njc = 0;            // physical cells in j
        std::vector<double> D;  // diagonal blocks
        std::vector<double> Wb; // west blocks  (i-1)
        std::vector<double> Eb; // east blocks  (i+1)
        std::vector<double> Sb; // south blocks (j-1)
        std::vector<double> Nb; // north blocks (j+1)

        // Allocate for nic x njc cells, all blocks set to zero.
        void resize(int nic_, int njc_)
        {
            nic = nic_;
            njc = njc_;
            const size_t n = static_cast<size_t>(nic) * njc * NEQ * NEQ;
            D.assign(n, 0.0);
            Wb.assign(n, 0.0);
            Eb.assign(n, 0.0);
            Sb.assign(n, 0.0);
            Nb.assign(n, 0.0);
        }

        int ncell() const { return nic * njc; }                         // number of block rows
        int index(int i, int j) const { return j * nic + i; }          // compact cell number p
        int west(int i) const { return i == 0 ? nic - 1 : i - 1; }     // i of the west neighbour (periodic)
        int east(int i) const { return i == nic - 1 ? 0 : i + 1; }     // i of the east neighbour (periodic)
    };

}
