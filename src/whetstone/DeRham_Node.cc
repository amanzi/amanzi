/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products for nodal DOFs.
*/

#include "MeshLight.hh"

#include "DeRham_Node.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
* Non-symmetric tensor is not yet used.
****************************************************************** */
int
DeRham_Node::L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List nodes, face_nodes;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  N.Reshape(nnodes, 1);
  Mc.Reshape(nnodes, nnodes);

  const auto& faces = mesh_->cell_get_faces(c);
  int nfaces = faces.size();

  double volume = mesh_->cell_volume(c);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // to calculate matrix R, we use temporary matrix N
  N.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

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
DeRham_Node::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  // StabilityScalar_(N, M);
  StabilityOptimized_(T, N, M);
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
