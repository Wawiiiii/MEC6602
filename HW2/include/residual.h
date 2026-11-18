#pragma once

#include "flow.h"
#include "mesh.h"

namespace euler {

    // on a start en ordre 1 pour verif de Roe mais le solveur est en 2 par defaut
    enum class Order {
        First,
        Second
    };

    void i_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double flux[NEQ]);
    void j_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double flux[NEQ]);

    void compute_residual(const Mesh& m, Flow& f, Order order);

}
