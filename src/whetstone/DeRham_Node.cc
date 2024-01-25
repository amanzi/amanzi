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

  Derham complex: mimetic inner products for nodal DOFs.
*/

#include "Mesh.hh"

#include "DeRham_Node.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
* Non-symmetric tensor is not yet used.
****************************************************************** */
int
DeRham_Node::L2consistency(int c, const Tensor<>& T, DenseMatrix<>& N, DenseMatrix<>& Mc, bool symmetry)
{

  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  N.reshape(nnodes, 1);
  Mc.reshape(nnodes, nnodes);

  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  double volume = mesh_->getCellVolume(c);
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

  // to calculate matrix R, we use temporary matrix N
  N.putScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    double tmp = (xf - xc) * normal;

    if (d_ == 2) {
      int m = (n + 1) % nnodes;
      N(n, 0) += tmp / 4;
      N(m, 0) += tmp / 4;
    }
  }

  // calculate upper part of R T R^T / volume
  for (int i = 0; i < nnodes; i++) {
    double a = N(i, 0) * T(0, 0) / volume;
    for (int j = i; j < nnodes; j++) { Mc(i, j) = a * N(j, 0); }
  }

  // populate matrix N
  for (int i = 0; i < nnodes; i++) { N(i, 0) = 1.0; }

  return 0;
}


/* ******************************************************************
* Mass matrix: adding stability matrix to the consistency matrix.
****************************************************************** */
int
DeRham_Node::MassMatrix(int c, const Tensor<>& T, DenseMatrix<>& M)
{
  DenseMatrix<> N;

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  // StabilityScalar_(N, M);
  StabilityOptimized_(T, N, M);
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
