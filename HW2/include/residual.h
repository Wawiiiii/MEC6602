#pragma once

#include "flow.h"
#include "mesh.h"

namespace euler {

    // Spatial order of the convective flux: cell values (first) or MUSCL values (second), slide 12 step 1.
    enum class Order {
        First,
        Second
    };

    /*
       Roe flux times the face length, (Fc dS), on one face (slides 12-16):
         i-face (i, j): between cells (i-1, j) and (i, j), i = 0..ni-1, j = 0..njc-1,
                        MUSCL along i with cells i-2 .. i+1
         j-face (i, j): between cells (i, j-1) and (i, j), i = 0..nic-1, j = 0..nj-1,
                        MUSCL along j with cells j-2 .. j+1
       The normal goes from the left cell to the right cell (mesh.h).
       f.Q must be up to date (compute_primitives), ghost cells included.
       Physical boundary faces (j-faces j = 0 wall and j = nj-1 farfield) do not use Roe:
       they use the flux of the averaged state Fc(1/2 (W_ghost + W_interior)) (basic_algo
       slide 7, method #2). The i-cut faces are connectivity and use Roe like interior faces.
    */
    void i_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double F[NEQ]);
    void j_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double F[NEQ]);

    /*
       Residual of every physical cell (Euler2D slide 9, basic_algo slide 6):
         R = sum over the 4 faces of (Fc dS), normal pointing out of the cell.
       One flux per face: added to the left cell, subtracted from the right cell.
       Updates f.Q from f.W first. Ghost slots of f.R get partial sums and must be ignored.
    */
    void compute_residual(const Mesh& m, Flow& f, Order order);

}
