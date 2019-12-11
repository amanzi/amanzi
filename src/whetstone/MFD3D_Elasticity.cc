/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method for elasticity.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_Elasticity.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Consistency condition for mass matrix in mechanics.
 * Only the upper triangular part of Ac is calculated.
 * Requires mesh_get_edges to complete the implementation.
 ****************************************************************** */
int
MFD3D_Elasticity::L2consistency(int c, const Tensor& T, DenseMatrix& N,
                                DenseMatrix& Mc, bool symmetry)
{
  Kokkos::View<Entity_ID*> faces;

  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  double volume = mesh_->cell_volume(c, false);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nfaces; i++) {
    int f = faces(i);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = fm - cm;
    double a = normal * v1;
    for (int k = 0; k < d_; k++) v1[k] = a - normal[k] * v1[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Consistency condition for stiffness matrix in mechanics.
 * Only the upper triangular part of Ac is calculated.
 ****************************************************************** */
int
MFD3D_Elasticity::H1consistency(int c, const Tensor& T, DenseMatrix& N,
                                DenseMatrix& Ac)
{
  Kokkos::View<Entity_ID*> faces, nodes;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.extent(0);

  int nrows = d_ * nnodes;
  int md = d_ * (d_ + 1);
  N.Reshape(nrows, md);
  Ac.Reshape(nrows, nrows);

  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  // convolution of tensors
  std::vector<Tensor> TE;

  for (int k = 0; k < d_; k++) {
    Tensor E(d_, 2);
    E(k, k) = 1.0;
    TE.push_back(T * E);
  }

  for (int k = 0; k < d_; k++) {
    for (int l = k + 1; l < d_; l++) {
      Tensor E(d_, 2);
      E(k, l) = E(l, k) = 1.0;
      TE.push_back(T * E);
    }
  }
  int nd = TE.size();

  // to calculate matrix R, we use temporary matrix N
  DenseMatrix R(nrows, nd, N.Values(), WHETSTONE_DATA_ACCESS_VIEW);
  R.PutScalar(0.0);

  for (int i = 0; i < nfaces; i++) {
    int f = faces(i);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    Kokkos::View<Entity_ID*> face_nodes;
    mesh_->face_get_nodes(f, face_nodes);
    int num_face_nodes = face_nodes.extent(0);

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes(j);
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs(i);
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes(jnext);
        int vprev = face_nodes(jprev);

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs(i) * norm(v3) / (4 * area);
      }

      int pos = 0;
      for (pos = 0; pos < nodes.extent(0); ++pos) {
        if (nodes(pos) == v) { break; }
      }
      for (int k = 0; k < nd; k++) {
        v1 = TE[k] * normal;
        for (int l = 0; l < d_; l++) R(l * nnodes + pos, k) += v1[l] * u;
      }
    }
  }

  // calculate R inv(T) R^T / volume
  Tensor Tinv(T);
  Tinv.Inverse();

  double volume = mesh_->cell_volume(c, false);
  Tinv *= 1.0 / volume;

  DenseMatrix RT(nrows, nd);
  if (Tinv.rank() == 1) {
    double* data_N = R.Values(); // We stress that memory belongs to matrix N.
    double* data_RT = RT.Values();
    double s = Tinv(0, 0);
    for (int i = 0; i < nrows * nd; i++) data_RT[i] = data_N[i] * s;
  } else if (Tinv.rank() == 4) {
    DenseMatrix Ttmp(nd, nd, Tinv.data_ptr(), WHETSTONE_DATA_ACCESS_VIEW);
    MatrixMatrixProduct_(R, Ttmp, false, RT);
  }
  DenseMatrix AcAc(nrows, nrows);
  MatrixMatrixProduct_(RT, R, true, AcAc);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) Ac(i, j) = AcAc(i, j);
  }

  // calculate matrix N
  N.PutScalar(0.0);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < nnodes; i++) {
    int v = nodes(i);
    mesh_->node_get_coordinates(v, &p);
    v1 = p - cm;

    int md = 0;
    for (int k = 0; k < d_; k++) {
      N(k * nnodes + i, md) = v1[k];
      md++;
    }
    for (int k = 0; k < d_; k++) {
      for (int l = k + 1; l < d_; l++) {
        N(k * nnodes + i, md) = v1[l];
        N(l * nnodes + i, md) = v1[k];
        md++;
      }
    }
    for (int k = 0; k < d_; k++) { // additional columns correspod to kernel
      N(k * nnodes + i, md) = 1.0;
      md++;
    }
    for (int k = 0; k < d_; k++) {
      for (int l = k + 1; l < d_; l++) {
        N(k * nnodes + i, md) = v1[l];
        N(l * nnodes + i, md) = -v1[k];
        md++;
      }
    }
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Lame stiffness matrix: a wrapper for other low-level routines
 ****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Lame stiffness matrix: a wrapper for other low-level routines
 ****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrixOptimized(int c, const Tensor& T,
                                           DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  StabilityOptimized_(T, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Lame stiffness matrix: a wrapper for other low-level routines
 * For education purpose only: ther are no M-matrices in elasticity.
 ****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE;
  ok = StabilityMMatrix_(c, N, A, objective);
  return ok;
}


/* ******************************************************************
 * Classical matrix-matrix product
 ****************************************************************** */
void
MFD3D_Elasticity::MatrixMatrixProduct_(const DenseMatrix& A,
                                       const DenseMatrix& B, bool transposeB,
                                       DenseMatrix& AB)
{
  int nrows = A.NumRows();
  int ncols = A.NumCols();

  int mrows = B.NumRows();
  int mcols = B.NumCols();

  if (transposeB) {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < mrows; j++) {
        double s = 0.0;
        for (int k = 0; k < ncols; k++) s += A(i, k) * B(j, k);
        AB(i, j) = s;
      }
    }
  } else {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < mcols; j++) {
        double s = 0.0;
        for (int k = 0; k < ncols; k++) s += A(i, k) * B(k, j);
        AB(i, j) = s;
      }
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
