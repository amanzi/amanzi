/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
int MFD3D_Electromagnetics::L2consistencyBoundary(
    int f, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc)
{
  Entity_ID_List edges;
  std::vector<int> dirs;

  mesh_->face_get_edges_and_dirs(f, &edges, &dirs);
  int nedges = edges.size();

  N.Reshape(nedges, d_ - 1);
  Mc.Reshape(nedges, nedges);

  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
  double area = mesh_->face_area(f);

  // calculate rotation matrix
  Tensor P(d_, 2); 

  v3 = normal / area;
  for (int i = 0; i < d_; i++) {
    for (int j = 0; j < d_; j++) { 
      P(i, j) = v3[i] * v3[j];
    }
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
  PTP = Pt * Tinv * P;

  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    double a1 = mesh_->edge_length(e) * dirs[i];
    v2 = PTP * (xe - xf);

    for (int j = i; j < nedges; j++) {
      e = edges[j];
      const AmanziGeometry::Point& ye = mesh_->edge_centroid(e);
      double a2 = mesh_->edge_length(e) * dirs[j];

      v1 = ye - xf;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / area;
    }
  }

  // Rows of matrix N are edge vectors in a coordinate system of face f.
  v1 = mesh_->edge_vector(edges[0]) / mesh_->edge_length(edges[0]);
  v2 = v3^v1;
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);
    N(i, 0) =-(tau * v2) / len; 
    N(i, 1) = (tau * v1) / len; 
  }

  return 0;
}


/* ******************************************************************
* Matrix matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrixBoundary(int f, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistencyBoundary(f, T, N, M);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi



