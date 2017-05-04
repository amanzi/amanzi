/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin method.
*/

#ifndef AMANZI_WHETSTONE_DG_HH_
#define AMANZI_WHETSTONE_DG_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "WhetStone_typedefs.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

// Gauss quadrature on interval (0,1)
const double q1d_weights[4][4] = {
    1.0, 0.0, 0.0, 0.0,
    0.5, 0.5, 0.0, 0.0,
    0.277777777777778, 0.444444444444444, 0.277777777777778, 0.0,
    0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
};
const double q1d_points[4][4] = {
    0.5, 0.0, 0.0, 0.0,
    0.211324865405187, 0.788675134594813, 0.0, 0.0,
    0.112701665379258, 0.5, 0.887298334620742, 0.0,
    0.0694318442029737, 0.330009478207572, 0.669990521792428, 0.930568155797026
};

class DG { 
 public:
  DG() : mesh_(Teuchos::null), d_(0) {};
  DG(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), d_(mesh_->space_dimension()) {};
  ~DG() {};

  int TaylorMassMatrix(int c, int order, DenseMatrix& M);

 private:
  void IntegrateMonomialsCell_(int c, int k, double* monomials);
  void IntegrateMonomialsFace_(int f, int k, double factor, double* monomials);
  void IntegrateMonomialsEdge_(const AmanziGeometry::Point& x1,
                               const AmanziGeometry::Point& x2,
                               int k, double factor, double* monomials);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

