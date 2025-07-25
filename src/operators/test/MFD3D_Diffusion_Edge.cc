/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, version 2.1
  Release name: naka-to.

  Mimetic discretization of elliptic operator using edge-based
  degrees of freedom shows flexibility of the discretization framework.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion_Edge.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

RegisteredFactory<MFD3D_Diffusion_Edge> MFD3D_Diffusion_Edge::reg_("diffusion edge");

/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction.
* Only the upper triangular part of Ac is calculated.
* The degrees of freedom are at nodes.
****************************************************************** */
int
MFD3D_Diffusion_Edge::H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  N.Reshape(nedges, d_ + 1);
  Ac.Reshape(nedges, nedges);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  // calculate matrix R (we re-use matrix N)
  if (d_ == 3) N.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    if (d_ == 2) {
      for (int k = 0; k < d_; k++) N(n, k) = normal[k] * dirs[n];
    } else {
      const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
      double area = mesh_->getFaceArea(f);

      auto [fedges, edirs] = mesh_->getFaceEdgesAndDirections(f);
      auto map = mesh_->getFaceCellEdgeMap(f, c);
      int nfedges = fedges.size();

      int e0 = fedges[0];
      const AmanziGeometry::Point& xe0 = mesh_->getEdgeCentroid(e0);

      for (int k = 0; k < d_; ++k) N(map[0], k) += normal[k] * dirs[n];

      for (int m = 0; m < nfedges; ++m) {
        int e = fedges[m];
        const AmanziGeometry::Point& tau = mesh_->getEdgeVector(e);

        double tmp = ((tau ^ normal) * (xf - xe0)) * dirs[n] * edirs[m] / area;

        for (int k = 0; k < d_; ++k) {
          N(map[m], k) += normal[k] * tmp / area;
        }
      }
    }
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);

  for (int n = 0; n < nedges; ++n) {
    for (int k = 0; k < d_; k++) v1[k] = N(n, k);
    v2 = K * v1;

    for (int m = n; m < nedges; m++) {
      for (int k = 0; k < d_; k++) v1[k] = N(m, k);
      Ac(n, m) = (v1 * v2) / volume;
    }
  }

  // calculate N
  for (int n = 0; n < nedges; n++) {
    int e = edges[n];
    const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);
    for (int k = 0; k < d_; k++) N(n, k) = xe[k] - xc[k];
    N(n, d_) = 1.0;
  }

  // Internal verification
  // DenseMatrix NtR(d_ + 1, d_ + 1);
  // NtR.Multiply(N, R, true);
  // std::cout << NtR << std::endl;

  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int
MFD3D_Diffusion_Edge::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
