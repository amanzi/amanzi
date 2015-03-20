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
* Consistency condition for the mass matrix in electromagnetics.
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistency(int c, const Tensor& T,
                                          DenseMatrix& N, DenseMatrix& Mc)
{
  Entity_ID_List edges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  int n1, n2, d = mesh_->space_dimension();
  ASSERT(d == 3);
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

  // Matrix N are simply tangents
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);
    for (int k = 0; k < d; ++k) N(i, k) = tau[k] / len;
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


}  // namespace WhetStone
}  // namespace Amanzi



