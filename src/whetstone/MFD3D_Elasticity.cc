/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method for elasticity.
*/

#include <cmath>
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
int MFD3D_Elasticity::L2consistency(int c, const Tensor& T,
                                    DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List faces;

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(c);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = fm - cm;
    double a = normal * v1;
    for (int k = 0; k < d; k++) v1[k] = a - normal[k] * v1[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in mechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Elasticity::L2consistencyInverse(
    int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) {
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in mechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Elasticity::H1consistency(int c, const Tensor& T,
                                    DenseMatrix& N, DenseMatrix& Ac)
{
  int nrows = N.NumRows();

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int num_nodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int d = mesh_->space_dimension();
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  // convolution of tensors
  std::vector<Tensor> TE;

  for (int k = 0; k < d; k++) {
    Tensor E(d, 2);
    E(k, k) = 1.0;
    TE.push_back(T * E);
  }

  for (int k = 0; k < d; k++) {
    for (int l = k + 1; l < d; l++) {
      Tensor E(d, 2);
      E(k, l) = E(l, k) = 1.0;
      TE.push_back(T * E);
    }
  }
  int nd = TE.size();

  // to calculate matrix R, we use temporary matrix N
  DenseMatrix R(nrows, nd, N.Values(), WHETSTONE_DATA_ACCESS_VIEW);
  R.PutScalar(0.0);

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d == 2) {
        u = 0.5 * dirs[i]; 
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1^v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = FindPosition_(v, nodes);
      for (int k = 0; k < nd; k++) {
        v1 = TE[k] * normal;
        for (int l = 0; l < d; l++) R(l * num_nodes + pos, k) += v1[l] * u;
      }
    }
  }

  // calculate R inv(T) R^T / volume
  Tensor Tinv(T);
  Tinv.Inverse();

  double volume = mesh_->cell_volume(c);
  Tinv *= 1.0 / volume;

  DenseMatrix RT(nrows, nd);
  if (Tinv.rank() == 1) {
    double* data_N = R.Values();  // We stress that memory belongs to matrix N.
    double* data_RT = RT.Values();
    double s = Tinv(0, 0);
    for (int i = 0; i < nrows * nd; i++) data_RT[i] = data_N[i] * s;
  } else if (Tinv.rank() == 4) {
    DenseMatrix Ttmp(nd, nd, Tinv.data(), WHETSTONE_DATA_ACCESS_VIEW);
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

  for (int i = 0; i < num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    v1 = p - cm;

    int md = 0;
    for (int k = 0; k < d; k++) {
      N(k * num_nodes + i, md) = v1[k];
      md++;
    }
    for (int k = 0; k < d; k++) {
      for (int l = k + 1; l < d; l++) {
        N(k * num_nodes + i, md) = v1[l];
        N(l * num_nodes + i, md) = v1[k];
        md++;
      }
    }
    for (int k = 0; k < d; k++) {  // additional columns correspod to kernel
      N(k * num_nodes + i, md) = 1.0;
      md++;
    }
    for (int k = 0; k < d; k++) {
      for (int l = k + 1; l < d; l++) {
        N(k * num_nodes + i, md) =  v1[l];
        N(l * num_nodes + i, md) = -v1[k];
        md++;
      }
    }
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nd = d * (d + 1);
  int nrows = A.NumRows();

  DenseMatrix N(nrows, nd);

  int ok = H1consistency(c, T, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrixOptimized(int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nd = d * (d + 1);
  int nrows = A.NumRows();

  DenseMatrix N(nrows, nd);

  int ok = H1consistency(c, T, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityOptimized(T, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
* For education purpose only: ther are no M-matrices in elasticity.
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nd = d * (d + 1);
  int nrows = A.NumRows();

  DenseMatrix N(nrows, nd);

  int ok = H1consistency(c, T, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE;
  ok = StabilityMMatrix_(c, N, A, objective);

  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Classical matrix-matrix product
****************************************************************** */
void MFD3D_Elasticity::MatrixMatrixProduct_(
    const DenseMatrix& A, const DenseMatrix& B, bool transposeB,
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

}  // namespace WhetStone
}  // namespace Amanzi



