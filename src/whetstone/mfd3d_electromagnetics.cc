/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Release name: ara-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "tensor.hh"
#include "mfd3d_electromagnetics.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistency(int c, const Tensor& T,
                                          DenseMatrix& N, DenseMatrix& Mc)
{
  int ok, d = mesh_->space_dimension();
  if (d == 2) {
    ok = L2consistency2D_(c, T, N, Mc);
  } else {
    ok = L2consistency3D_(c, T, N, Mc);
  }

  return ok;
}


/* ******************************************************************
* Consistency condition for the mass matrix in electromagnetics.
* Only the upper triangular part of Mc = R (R^T N)^{-1} R^T is 
* calculated. Here R^T N = |c| T.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistency2D_(int c, const Tensor& T,
                                             DenseMatrix& N, DenseMatrix& Mc)
{
  int n1, n2, d(2);
  Entity_ID_List edges;
  std::vector<int> edirs;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  AmanziGeometry::Point v1(d), v2(d), v3(d), tau(d), p1(d), p2(d);

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < nedges; ++i) {
    int e = edges[i];
    double len = mesh_->edge_length(e);

    mesh_->edge_get_nodes(e, &n1, &n2);
    mesh_->node_get_coordinates(n1, &p1);
    mesh_->node_get_coordinates(n2, &p2);

    v1 = ((p1 + p2) / 2) - xc;

    for (int k = 0; k < d; ++k) {
      len = -len;
      N(i, k) = len * v1[1 - k];
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T 
  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d; ++k) v1[k] = N(i, k);
    v2 = Tinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d; ++k) v3[k] = N(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

  // Rows of matrix N are oriented tangent vectors. Since N goes to 
  // the Gramm-Schmidt orthogonalizetion procedure, we can skip scaling
  // tensorial factor T.
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);
    for (int k = 0; k < d; ++k) N(i, k) = tau[k] / len;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for the mass matrix in electromagnetics.
* Only the upper triangular part of Mc = R (R^T N)^{-1} R^T is 
* calculated. Here R^T N = |c| T.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistency3D_(int c, const Tensor& T,
                                             DenseMatrix& N, DenseMatrix& Mc)
{
  int n1, n2, d(3);
  Entity_ID_List edges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  AmanziGeometry::Point v1(d), v2(d), v3(d), tau(d), p1(d), p2(d);
  AmanziGeometry::Point vv[3];

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);

    v1 = xc - xf; 
 
    double a1 = normal * v1;
    for (int k = 0; k < d; ++k) {
      vv[k].set(normal[k] * v1);
      vv[k][k] -= a1;
    }

    mesh_->face_get_edges_and_dirs(f, &edges, &edirs);
    int nedges = edges.size();

    mesh_->face_to_cell_edge_map(f, c, &map);

    for (int m = 0; m < nedges; ++m) {
      int e = edges[m];
      mesh_->edge_get_nodes(e, &n1, &n2);
      mesh_->node_get_coordinates(n1, &p1);
      mesh_->node_get_coordinates(n2, &p2);
 
      v3 = ((p1 + p2) / 2) - xf;

      double len = mesh_->edge_length(e);
      len /= 2.0 * area * area * fdirs[i] * edirs[m];

      for (int k = 0; k < d; ++k) {
        N(map[m], k) -= len * ((vv[k]^v3) * normal);
      }
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T 
  Tensor Tinv(T);
  Tinv.Inverse();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d; ++k) v1[k] = N(i, k);
    v2 = Tinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d; ++k) v3[k] = N(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

  // Rows of matrix N are simply tangents. Since N goes to the
  // Gramm-Schmidt orthogonalizetion procedure, we can skip scaling
  // tensorial factor T.
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);
    for (int k = 0; k < d; ++k) N(i, k) = tau[k] / len;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of the mass matrix in 
* electromagnetics. Only the upper triangular part of matrix 
* Wc = N (N^T R)^{-1} N^T is calculated. Here N^T R = |c| T.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistencyInverse(
    int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List edges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  int n1, n2, d = mesh_->space_dimension();
  ASSERT(d == 3);
  AmanziGeometry::Point v1(d), v2(d), v3(d), tau(d), p1(d), p2(d);
  AmanziGeometry::Point vv[3];

  // Since the matrix N is the matrix of scaled tangent vextors,
  // N = N0 T, we use the tensor T.
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  for (int i = 0; i < nedges; i++) {
    int e1 = edges[i];
    const AmanziGeometry::Point& v1 = mesh_->edge_vector(e1);
    double a1 = mesh_->edge_length(e1);
    v3 = T * v1;

    for (int j = i; j < nedges; j++) {
      int e2 = edges[j];
      const AmanziGeometry::Point& v2 = mesh_->edge_vector(e2);
      double a2 = mesh_->edge_length(e2);
      Wc(i, j) = (v2 * v3) / (a1 * a2 * volume);
    }
  }

  // Calculate matrix R
  R.PutScalar(0.0);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);

    v1 = xc - xf; 
 
    double a1 = normal * v1;
    for (int k = 0; k < d; ++k) {
      vv[k].set(normal[k] * v1);
      vv[k][k] -= a1;
    }

    mesh_->face_get_edges_and_dirs(f, &edges, &edirs);
    int nedges = edges.size();

    mesh_->face_to_cell_edge_map(f, c, &map);

    for (int m = 0; m < nedges; ++m) {
      int e = edges[m];
      mesh_->edge_get_nodes(e, &n1, &n2);
      mesh_->node_get_coordinates(n1, &p1);
      mesh_->node_get_coordinates(n2, &p2);
 
      v3 = ((p1 + p2) / 2) - xf;

      double len = mesh_->edge_length(e);
      len /= 2.0 * area * area * fdirs[i] * edirs[m];

      for (int k = 0; k < d; ++k) {
        R(map[m], k) -= len * ((vv[k]^v3) * normal);
      }
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Matrix matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::H1consistency(int c, const Tensor& T,
                                          DenseMatrix& N, DenseMatrix& Ac)
{
  int n1, n2, d = mesh_->space_dimension();
  ASSERT(d == 3);

  Entity_ID_List edges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  AmanziGeometry::Point v1(d), v2(d), v3(d), p1(d), p2(d);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);

    v1 = xc - xf; 
 
    mesh_->face_get_edges_and_dirs(f, &edges, &edirs);
    int nedges = edges.size();

    mesh_->face_to_cell_edge_map(f, c, &map);

    for (int m = 0; m < nedges; ++m) {
      int e = edges[m];

      double len = mesh_->edge_length(e);
      len *= 2 * fdirs[i] * edirs[m];

      for (int k = 0; k < d; ++k) {
        N(map[m], k) += len * v1[k];
      }
    }
  }
  
  // calculate Ac = N (R^T N)^{+} N^T
  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d; ++k) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d; ++k) v3[k] = N(j, k);
      Ac(i, j) = (v2 * v3) / (4 * volume);
    }
  }

  // Matrix N(:, 1:3) are simply tangents
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);

    for (int k = 0; k < d; ++k) N(i, k) = tau[k] / len;

    mesh_->edge_get_nodes(e, &n1, &n2);
    mesh_->node_get_coordinates(n1, &p1);
    mesh_->node_get_coordinates(n2, &p2);
 
    v1 = ((p1 + p2) / 2) - xc;
    v2 = v1^tau;

    for (int k = 0; k < d; ++k) {
      N(i, k + d) = v2[k] / len;
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Matrix matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = M.NumRows();

  DenseMatrix N(nrows, d);
  DenseMatrix Mc(nrows, nrows);

  int ok = L2consistency(c, T, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix optimized for sparsity.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = M.NumRows();

  DenseMatrix N(nrows, d);
  DenseMatrix Mc(nrows, nrows);

  int ok = L2consistency(c, T, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(T, N, Mc, M);
  return ok;
}


/* ******************************************************************
* Inverse of the mass matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nrows = W.NumRows();

  DenseMatrix R(nrows, d);
  DenseMatrix Wc(nrows, nrows);

  int ok = L2consistencyInverse(c, T, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse of matrix matrix for edge-based discretization optimized
* for starcity.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrixInverseOptimized(
    int c, const Tensor& T, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nrows = W.NumRows();

  DenseMatrix R(nrows, d);
  DenseMatrix Wc(nrows, nrows);

  int ok = L2consistencyInverse(c, T, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(T, R, Wc, W);
  return ok;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Electromagnetics::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nedges = A.NumRows();

  int nd = d * (d + 1) / 2;
  DenseMatrix N(nedges, nd);
  DenseMatrix Ac(nedges, nedges);

  int ok = H1consistency(c, T, N, Ac);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, Ac, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix optimized for sparsity.
****************************************************************** */
int MFD3D_Electromagnetics::StiffnessMatrixOptimized(
    int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nedges = A.NumRows();

  int nd = d * (d + 1) / 2;
  DenseMatrix N(nedges, nd);
  DenseMatrix Ac(nedges, nedges);

  int ok = H1consistency(c, T, N, Ac);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(T, N, Ac, A);
  return ok;
}

}  // namespace WhetStone
}  // namespace Amanzi



