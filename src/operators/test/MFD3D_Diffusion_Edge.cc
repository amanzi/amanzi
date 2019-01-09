/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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

namespace Amanzi {
namespace WhetStone {

RegisteredFactory<MFD3D_Diffusion_Edge> MFD3D_Diffusion_Edge::factory_("diffusion edge");

/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
* The degrees of freedom are at nodes.
****************************************************************** */
int MFD3D_Diffusion_Edge::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List faces, edges;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  N.Reshape(nedges, d_ + 1);
  Ac.Reshape(nedges, nedges);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

   AmanziGeometry::Point v1(d_), v2(d_);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    double a1 = mesh_->face_area(f);
    v2 = K * (xf - xc);

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& yf = mesh_->face_centroid(f);
      double a2 = mesh_->face_area(f);

      v1 = yf - xc;
      Ac(i, j) = (v1 * v2) * (a1 * a2) / volume;
    }
  }

  // calculate N
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    for (int k = 0; k < d_; k++) N(n, k) = xf[k] - xc[k];
    N(n, d_) = 1.0;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Diffusion_Edge::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



