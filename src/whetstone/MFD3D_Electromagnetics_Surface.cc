/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The mimetic finite difference method.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "MFD3D_Electromagnetics.hh"
#include "MFD3D_Diffusion.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for the mass matrix in electromagnetics for
* face f of a 3D cell.
****************************************************************** */
int
MFD3D_Electromagnetics::L2consistencyBoundary(int f,
                                              const Tensor<>& T,
                                              DenseMatrix<>& N,
                                              DenseMatrix<>& Mc)
{
  auto [edges, dirs] = mesh_->getFaceEdgesAndDirections(f);
  int nedges = edges.size();

  N.reshape(nedges, d_ - 1);
  Mc.reshape(nedges, nedges);

  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);
  const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
  const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
  double area = mesh_->getFaceArea(f);

  // calculate rotation matrix
  Tensor P(d_, 2);

  v3 = normal / area;
  for (int i = 0; i < d_; i++) {
    for (int j = 0; j < d_; j++) { P(i, j) = v3[i] * v3[j]; }
  }
  P(0, 1) -= v3[2];
  P(0, 2) += v3[1];

  P(1, 0) += v3[2];
  P(1, 2) -= v3[0];

  P(2, 0) -= v3[1];
  P(2, 1) += v3[0];

  // define rotated tensor
  Tensor PTP(d_, 2), Pt(P), Tinv(T);
  Tinv.Inverse();
  Pt.Transpose();
  PTP.assign(Pt * Tinv * P);

  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);
    double a1 = mesh_->getEdgeLength(e);
    v2 = PTP * (xe - xf);

    for (int j = i; j < nedges; j++) {
      e = edges[j];
      const AmanziGeometry::Point& ye = mesh_->getEdgeCentroid(e);
      double a2 = mesh_->getEdgeLength(e);

      v1 = ye - xf;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / area;
    }
  }

  // Rows of matrix N are normal vectors in the plane of face f.
  v1 = mesh_->getEdgeVector(edges[0]) / mesh_->getEdgeLength(edges[0]);
  v2 = v3 ^ v1;
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->getEdgeVector(e);
    double len = mesh_->getEdgeLength(e);
    N(i, 0) = -(tau * v2) * dirs[i] / len;
    N(i, 1) = (tau * v1) * dirs[i] / len;
  }

  return 0;
}


/* ******************************************************************
* Matrix matrix for edge-based discretization.
****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixBoundary(int f, const Tensor<>& T, DenseMatrix<>& M)
{
  DenseMatrix N;

  int ok = L2consistencyBoundary(f, T, N, M);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
