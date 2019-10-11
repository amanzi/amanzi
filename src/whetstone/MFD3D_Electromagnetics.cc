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
 * Efficient implementation is possible in 2D. Hence, we fork the code.
 * This is the experimental algorithm.
 ****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency(int c, const Tensor& T, DenseMatrix& N,
                                      DenseMatrix& Ac)
{
  int ok;
  if (d_ == 2) {
    ok = H1consistency2D_(c, T, N, Ac);
  } else {
    ok = H1consistency3D_(c, T, N, Ac);
  }
  return ok;
}


/* ******************************************************************
 * Consistency condition for a stiffness matrix.
 ****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency2D_(int c, const Tensor& T, DenseMatrix& N,
                                         DenseMatrix& Ac)
{
  Kokkos::View<Entity_ID*> faces;
  Kokkos::View<int*> fdirs;

  mesh_->cell_get_faces_and_dirs(c, faces, fdirs);
  int nfaces = faces.extent(0);

  int nd = 3;
  N.Reshape(nfaces, nd);
  Ac.Reshape(nfaces, nfaces);

  // calculate Ac = R (R^T N)^{+} R^T
  double T00 = T(0, 0);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces(i);
    double a1 = mesh_->face_area(f);

    for (int j = i; j < nfaces; ++j) {
      f = faces(j);
      double a2 = mesh_->face_area(f);
      Ac(i, j) = a1 * a2 / (T00 * fdirs(i) * fdirs(j) * volume);
    }
  }

  // Matrix N(:, 1:2) are simply tangents
  for (int i = 0; i < nfaces; i++) {
    int f = faces(i);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    double len = mesh_->face_area(f);

    for (int k = 0; k < d_; ++k) {
      len = -len;
      N(i, k) = normal[1 - k] / len;
    }

    N(i, d_) = (xf - xc) * normal / len;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Consistency condition for a stiffness matrix.
 ****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency3D_(int c, const Tensor& T, DenseMatrix& N,
                                         DenseMatrix& Ac)
{
  Kokkos::View<Entity_ID*> faces, fedges, edges;
  std::vector<int> map;
  Kokkos::View<int*> fdirs, edirs;

  mesh_->cell_get_faces_and_dirs(c, faces, fdirs);
  int nfaces = faces.extent(0);

  mesh_->cell_get_edges(c, edges);
  int nedges = edges.extent(0);

  int nd = 6; // order_ * (order_ + 2) * (order_ + 3) / 2;
  N.Reshape(nedges, nd);
  Ac.Reshape(nedges, nedges);

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces(i);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    AmanziGeometry::Point normal = mesh_->face_normal(f);
    normal /= mesh_->face_area(f);

    mesh_->face_get_edges_and_dirs(f, fedges, &edirs);
    int nfedges = fedges.extent(0);

    mesh_->face_to_cell_edge_map(f, c, &map);

    for (int m = 0; m < nfedges; ++m) {
      int e = fedges(m);
      const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
      double len = mesh_->edge_length(e);

      len *= 2 * fdirs(i) * edirs(m);
      v2 = xe - xf;

      for (int k = 0; k < d_; ++k) N(map[m], k) += len * v2[k];
    }
  }

  // calculate Ac = R (R^T N)^{+} R^T
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Ac(i, j) = (v2 * v3) / (4 * volume);
    }
  }

  // Matrix N(:, 1:3) are simply tangents
  for (int i = 0; i < nedges; i++) {
    int e = edges(i);
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    AmanziGeometry::Point tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);

    tau /= len;
    v1 = xe - xc;
    v3 = tau ^ v1;

    for (int k = 0; k < d_; ++k) N(i, k) = tau[k];
    for (int k = 0; k < d_; ++k) N(i, d_ + k) = v3[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Mass matrix optimized for sparsity.
 ****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixOptimized(int c, const Tensor& T,
                                            DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  ok = StabilityOptimized_(T, N, M);
  return ok;
}


/* ******************************************************************
 * Inverse of matrix matrix for edge-based discretization optimized
 * for starcity.
 ****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixInverseOptimized(int c, const Tensor& T,
                                                   DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, T, R, W, true);
  if (ok) return ok;

  ok = StabilityOptimized_(T, R, W);
  return ok;
}


/* ******************************************************************
 * A simple mass matrix for testing.
 ****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixDiagonal(int c, const Tensor& T,
                                           DenseMatrix& M)
{
  double volume = mesh_->cell_volume(c, false);

  Kokkos::View<Entity_ID*> edges;
  mesh_->cell_get_edges(c, edges);
  int nedges = edges.extent(0);

  M.PutScalar(0.0);
  for (int n = 0; n < nedges; n++) {
    int e = edges(n);
    M(n, n) = d_ * volume / (nedges * T(0, 0));
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix: the standard algorithm.
 ****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix M, C;

  int ok = StiffnessMatrix(c, T, A, M, C);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix: the standard algorithm. Curls in 2D and 3D are
 * defined using exterior face normals.
 ****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A,
                                        DenseMatrix& M, DenseMatrix& C)
{
  MFD3D_Diffusion diffusion(mesh_);
  int ok = diffusion.MassMatrixScaledArea(c, T, M);

  // populate curl matrix
  CurlMatrix(c, C);
  int nfaces = C.NumRows();
  int nedges = C.NumCols();

  A.Reshape(nedges, nedges);
  DenseMatrix MC(nfaces, nedges);

  MC.Multiply(M, C, false);
  A.Multiply(C, MC, true);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix: the algorithm based on new
 ****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrixGeneralized(int c, const Tensor& T,
                                                   DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  AddGradientToProjector_(c, N);

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix: the algorithm based on new
 ****************************************************************** */
void
MFD3D_Electromagnetics::AddGradientToProjector_(int c, DenseMatrix& N)
{
  Kokkos::View<Entity_ID*> edges, nodes;

  mesh_->cell_get_edges(c, edges);
  int nedges = edges.extent(0);

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  // reserve map: gid -> lid
  std::map<int, int> lid;
  for (int n = 0; n < nnodes; ++n) lid[nodes(n)] = n;

  // populate discrete gradient
  int v1, v2;
  DenseMatrix G(nedges, nnodes);
  G.PutScalar(0.0);

  for (int m = 0; m < nedges; ++m) {
    int e = edges(m);
    double length = mesh_->edge_length(e);
    mesh_->edge_get_nodes(e, &v1, &v2);

    G(m, lid[v1]) += 1.0 / length;
    G(m, lid[v2]) -= 1.0 / length;
  }

  // create matrix [N G] with possibly linearly dependent vectors
  int ncols = N.NumCols();
  DenseMatrix NG(nedges, ncols + nnodes - 1);
  for (int n = 0; n < ncols; ++n) {
    for (int m = 0; m < nedges; ++m) NG(m, n) = N(m, n);
  }
  for (int n = 0; n < nnodes - 1; ++n) {
    for (int m = 0; m < nedges; ++m) NG(m, ncols + n) = G(m, n);
  }

  // identify linearly independent vectors by using the
  // Gramm-Schmidt orthogonalization process
  int ngcols = ncols + nnodes - 1;
  double scale = G.Norm2();
  double tol = 1e-24 * scale * scale;

  std::vector<int> map;
  for (int n = 0; n < ngcols; ++n) {
    double l22 = 0.0;
    for (int m = 0; m < nedges; m++) l22 += NG(m, n) * NG(m, n);

    // skip column of matrix G
    if (n >= ncols && l22 < tol) continue;

    map.push_back(n);
    l22 = 1.0 / sqrt(l22);
    for (int m = 0; m < nedges; ++m) NG(m, n) *= l22;

    for (int k = n + 1; k < ngcols; ++k) {
      double s = 0.0;
      for (int m = 0; m < nedges; ++m) s += NG(m, n) * NG(m, k);
      for (int m = 0; m < nedges; ++m) NG(m, k) -= s * NG(m, n);
    }
  }

  // create new matrix N
  int nmap = map.size();
  N.Reshape(nedges, nmap);
  for (int n = 0; n < nmap; ++n) {
    for (int m = 0; m < nedges; ++m) N(m, n) = NG(m, map[n]);
  }
}


/* ******************************************************************
 * Curl matrix acts onto the space of total fluxes; hence, there is
 * no face area scaling below.
 ****************************************************************** */
void
MFD3D_Electromagnetics::CurlMatrix(int c, DenseMatrix& C)
{
  std::vector<int> map;
  Kokkos::View<Entity_ID*> faces, fedges, edges, nodes;
  Kokkos::View<int*> fdirs, edirs;

  mesh_->cell_get_edges(c, edges);
  int nedges = edges.extent(0);

  mesh_->cell_get_faces_and_dirs(c, faces, fdirs);
  int nfaces = faces.extent(0);

  C.Reshape(nfaces, nedges);
  C.PutScalar(0.0);

  if (d_ == 2) {
    mesh_->cell_get_nodes(c, nodes);

    for (int i = 0; i < nfaces; ++i) {
      int j = (i + 1) % nfaces;
      C(i, i) = -1.0;
      C(i, j) = 1.0;
    }
  } else {
    for (int i = 0; i < nfaces; ++i) {
      int f = faces(i);

      mesh_->face_to_cell_edge_map(f, c, &map);
      mesh_->face_get_edges_and_dirs(f, fedges, &edirs);
      int nfedges = fedges.extent(0);

      for (int j = 0; j < nfedges; ++j) {
        int e = fedges(j);
        double len = mesh_->edge_length(e);
        C(i, map[j]) = len * edirs(j) * fdirs(i);
      }
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
