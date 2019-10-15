/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion_Edge.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

RegisteredFactory<MFD3D_Diffusion_Edge>
  MFD3D_Diffusion_Edge::factory_("diffusion edge");

/* ******************************************************************
 * Consistency condition for stiffness matrix in heat conduction.
 * Only the upper triangular part of Ac is calculated.
 * The degrees of freedom are at nodes.
 ****************************************************************** */
int
MFD3D_Diffusion_Edge::H1consistency(int c, const Tensor& K, DenseMatrix& N,
                                    DenseMatrix& Ac)
{
  Entity_ID_List faces, edges, fedges;
  std::vector<int> dirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
  int nfaces = faces.size();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  N.Reshape(nedges, d_ + 1);
  Ac.Reshape(nedges, nedges);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  // calculate matrix R (we re-use matrix N)
  if (d_ == 3) N.putScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    if (d_ == 2) {
      for (int k = 0; k < d_; k++) N(n, k) = normal[k] * dirs[n];
    } else {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);

      mesh_->face_get_edges_and_dirs(f, &fedges, &edirs);
      mesh_->face_to_cell_edge_map(f, c, &map);
      int nfedges = fedges.size();

      int e0 = fedges[0];
      const AmanziGeometry::Point& xe0 = mesh_->edge_centroid(e0);

      for (int k = 0; k < d_; ++k) N(map[0], k) += normal[k] * dirs[n];

      for (int m = 0; m < nfedges; ++m) {
        int e = fedges[m];
        const AmanziGeometry::Point& tau = mesh_->edge_vector(e);

        double tmp = ((tau ^ normal) * (xf - xe0)) * dirs[n] * edirs[m] / area;

        for (int k = 0; k < d_; ++k) { N(map[m], k) += normal[k] * tmp / area; }
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
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    for (int k = 0; k < d_; k++) N(n, k) = xe[k] - xc[k];
    N(n, d_) = 1.0;
  }

  // Internal verification
  // DenseMatrix NtR(d_ + 1, d_ + 1);
  // NtR.elementWiseMultiply(N, R, true);
  // std::cout << NtR << std::endl;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
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
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

} // namespace WhetStone
} // namespace Amanzi
