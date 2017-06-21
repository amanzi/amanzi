/*
  WhetStone, version 2.1
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
int MFD3D_Electromagnetics::H1consistency(
    int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  int ok, d = mesh_->space_dimension();
  if (d == 2) {
    ok = H1consistency2DExperimental_(c, T, N, Ac);
  } else {
    ok = H1consistency3DExperimental_(c, T, N, Ac);
  }
  return ok;
}


/* ******************************************************************
* Stiffness matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::H1consistency2DExperimental_(
    int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  int d(2);
  Entity_ID_List faces;
  std::vector<int> fdirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  // calculate Ac = R (R^T N)^{+} R^T
  double T00 = T(0, 0);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    double a1 = mesh_->face_area(f);

    for (int j = i; j < nfaces; ++j) {
      f = faces[j];
      double a2 = mesh_->face_area(f);
      Ac(i, j) = a1 * a2 / (T00 * fdirs[i] * fdirs[j] * volume);
    }
  }

  // Matrix N(:, 1:2) are simply tangents
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    double len = mesh_->face_area(f);

    for (int k = 0; k < d; ++k) {
      len = -len;
      N(i, k) = normal[1 - k] / len;
    }

    N(i, d) = (xf - xc) * normal / len; 
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::H1consistency3DExperimental_(
    int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  int d(3);
  Entity_ID_List edges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  AmanziGeometry::Point v1(d), v2(d), v3(d);

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
  
  // calculate Ac = R (R^T N)^{+} R^T
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
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);

    for (int k = 0; k < d; ++k) N(i, k) = tau[k] / len;

    v1 = xe - xc;
    v2 = v1^tau;

    for (int k = 0; k < d; ++k) {
      N(i, k + d) = v2[k] / len;
    }
  }

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

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  ok = StabilityOptimized_(T, N, M);
  return ok;
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

  int ok = L2consistencyInverse(c, T, R, W, true);
  if (ok) return ok;

  ok = StabilityOptimized_(T, R, W);
  return ok;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Electromagnetics::StiffnessMatrix(
    int c, const Tensor& T, DenseMatrix& A)
{
  Entity_ID_List faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();
  int nedges = A.NumRows();

  DenseMatrix M(nfaces, nfaces), C(nfaces, nedges);

  int ok = StiffnessMatrix(c, T, A, M, C);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. Curl in 2D and 3D is 
* defined using the exterior face normal.
****************************************************************** */
int MFD3D_Electromagnetics::StiffnessMatrix(
    int c, const Tensor& T, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C)
{
  Entity_ID_List faces, nodes, fnodes, edges;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();
  int nedges = A.NumRows();

  DenseMatrix MC(nfaces, nedges);

  MFD3D_Diffusion diffusion(mesh_);
  int ok = diffusion.MassMatrixScaledArea(c, T, M);

  // populate curl matrix. Curl acts to the space of total
  // fluxes; hence, there is no face area scaling.
  C.PutScalar(0.0);

  int d = mesh_->space_dimension();
  if (d == 2) {
    mesh_->cell_get_nodes(c, &nodes);

    for (int i = 0; i < nfaces; ++i) {
      int j = (i + 1) % nfaces;
      C(i, i) = -1.0;
      C(i, j) = 1.0;
    }
  } else {
    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];

      mesh_->face_to_cell_edge_map(f, c, &map);
      mesh_->face_get_edges_and_dirs(f, &edges, &edirs);
      int medges = edges.size();

      for (int j = 0; j < medges; ++j) {
        int e = edges[j]; 
        double len = mesh_->edge_length(e);
        C(i, map[j]) = len * edirs[j] * fdirs[i];
      }
    }
  }

  MC.Multiply(M, C, false);
  A.Multiply(C, MC, true); 

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the experimental algorithm.
****************************************************************** */
int MFD3D_Electromagnetics::StiffnessMatrixExperimental(
    int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nedges = A.NumRows();

  int nd = d * (d + 1) / 2;
  DenseMatrix N(nedges, nd);

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



