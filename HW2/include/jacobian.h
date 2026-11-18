#pragma once

#include "bc.h"
#include "blockmat.h"
#include "flow.h"
#include "mesh.h"
#include "roe.h"

#include <vector>

namespace euler {

    void flux_jacobian(const double q[NEQ], double nx, double ny, double A[NEQ * NEQ]);


    void roe_abs_matrix(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double absA[NEQ * NEQ],
                        double harten_fraction = HARTEN_FRACTION);

    void roe_face_jacobians(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double dFdL[NEQ * NEQ],
                            double dFdR[NEQ * NEQ]);

    void assemble_first_order_jacobian(const Mesh& m, const Flow& f, const Freestream& fs, JminBc jmin,
                                       const std::vector<double>& diag_add, BlockMatrix& M);

}
