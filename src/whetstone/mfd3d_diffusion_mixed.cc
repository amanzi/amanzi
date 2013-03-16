/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: ara-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>
#include <vector>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::DarcyMass(int cell, const Tensor& permeability,
                     Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();
  int nfaces = M.numRows();

  Teuchos::SerialDenseMatrix<int, double> N(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Mc(nfaces, nfaces);

  Tensor permeability_inv(permeability);
  permeability_inv.inverse();

  int ok = L2consistency(cell, permeability_inv, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::DarcyMassInverse(int cell, const Tensor& permeability,
                            Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.numRows();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D::DarcyMassInverseHex(int cell, const Tensor& permeability,
                               Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.numRows();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  int flag = StabilityMonotoneHex(cell, permeability, Wc, W);
  if (flag) StabilityScalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* This is a debug version of the above routine for a scalar tensor
* and an orthogonal brick element.
****************************************************************** */
int MFD3D::DarcyMassInverseDiagonal(int cell, const Tensor& permeability,
                                    Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
  int nfaces = faces.size();

  W.putScalar(0.0);
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    W(n, n) = nfaces * permeability(0, 0) * area * area / (d * volume);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Second-generation MFD method as inlemented in RC1.
****************************************************************** */
int MFD3D::DarcyMassInverseSO(int cell, const Tensor& permeability,
                              Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);

  AmanziMesh::Entity_ID_List nodes, corner_faces;
  mesh_->cell_get_nodes(cell, &nodes);
  int nnodes = nodes.size();

  Tensor K(permeability);
  K.inverse();

  // collect all corner matrices
  std::vector<Tensor> Mv;
  std::vector<double> cwgt;

  Tensor N(d, 2), NK(d, 2), Mv_tmp(d, 2);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, cell, AmanziMesh::USED, &corner_faces);
    int nfaces = corner_faces.size();
    if (nfaces < d) {
      Errors::Message msg;
      msg << "WhetStone MFD3D: number of faces forming a corner is small.";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i < d; i++) {
      int f = corner_faces[i];
      N.add_column(i, mesh_->face_normal(f));
    }
    double cwgt_tmp = fabs(N.determinant());

    N.inverse();
    NK = N * K;

    N.transpose();
    Mv_tmp = NK * N;
    Mv.push_back(Mv_tmp);

    for (int i = 0; i < d; i++) {
      int f = corner_faces[i];
      cwgt_tmp /= mesh_->face_area(f);
    }
    cwgt.push_back(cwgt_tmp);
  }

  // rescale corner weights
  double factor = 0.0;
  for (int n = 0; n < nnodes; n++) factor += cwgt[n];
  factor = mesh_->cell_volume(cell) / factor;

  for (int n = 0; n < nnodes; n++) cwgt[n] *= factor;

  // assemble corner matrices
  W.putScalar(0.0);
  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, cell, AmanziMesh::USED, &corner_faces);

    Tensor& Mv_tmp = Mv[n];
    for (int i = 0; i < d; i++) {
      int k = FindPosition_(corner_faces[i], faces);
      for (int j = i; j < d; j++) {
        int l = FindPosition_(corner_faces[j], faces);
        W(k, l) += Mv_tmp(i, j) * cwgt[n] * fdirs[k] * fdirs[l];
        W(l, k) = W(k, l);
      }
    }
  }
 
  // invert matrix W
  Teuchos::LAPACK<int, double> lapack;
  int info, size = W.numRows();

  int ipiv[size];
  double work[size];

  lapack.GETRF(size, size, W.values(), size, ipiv, &info);
  lapack.GETRI(size, W.values(), size, ipiv, work, size, &info);
  if (info != 0) {
    Errors::Message msg;
    msg << "WhetStone MFD3D: support operator generated bad elemental mass matrix.";
    Exceptions::amanzi_throw(msg);
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy inverse mass matrix via optimziation.
****************************************************************** */
int MFD3D::DarcyMassInverseOptimized(int cell, const Tensor& permeability,
                                     Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.numRows();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(permeability, R, Wc, W);
  return ok;
}


/* ******************************************************************
* Darcy inverse mass matrix via optimization, experimental.
****************************************************************** */
int MFD3D::DarcyMassInverseOptimizedScaled(int cell, const Tensor& permeability,
                                           Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.numRows();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverseScaled(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  ok = StabilityOptimized(permeability, R, Wc, W);
  RescaleDarcyMassInverse_(cell, W);

  return ok;
}


/* ******************************************************************
* Consistency condition for inner product in space of Darcy fluxes. 
* Only upper triangular part of Mc is calculated.
* Darcy flux is scaled by area!
****************************************************************** */
int MFD3D::L2consistency(int cell, const Tensor& T,
                         Teuchos::SerialDenseMatrix<int, double>& N,
                         Teuchos::SerialDenseMatrix<int, double>& Mc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);

  int num_faces = faces.size();
  if (num_faces != N.numRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  Tensor Tinv(T);
  Tinv.inverse();

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    v2 = Tinv * (fm - cm);

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
      v1 = fm - cm;
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    for (int k = 0; k < d; k++) N(i, k) = normal[k];
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of Darcy 
* fluxes. Only the upper triangular part of Wc is calculated.
* Darcy flux is scaled by area!
****************************************************************** */
int MFD3D::L2consistencyInverse(int cell, const Tensor& T,
                                Teuchos::SerialDenseMatrix<int, double>& R,
                                Teuchos::SerialDenseMatrix<int, double>& Wc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.numRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = T * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume);
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d; k++) R(i, k) = fm[k] - cm[k];
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of Darcy 
* fluxes. Only the upper triangular part of Wc is calculated.
****************************************************************** */
int MFD3D::L2consistencyInverseScaled(int cell, const Tensor& T,
                                      Teuchos::SerialDenseMatrix<int, double>& R,
                                      Teuchos::SerialDenseMatrix<int, double>& Wc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.numRows()) return num_faces;  // matrix was not reshaped

  // calculate areas of possibly curved faces
  std::vector<double> areas(num_faces, 0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // populate matrix W_0
  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = T * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d; k++) R(i, k) = (fm[k] - cm[k]) * mesh_->face_area(f);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Rescale matrix to area-weighted fluxes.
****************************************************************** */
void MFD3D::RescaleDarcyMassInverse_(int cell, Teuchos::SerialDenseMatrix<int, double>& W)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
  int num_faces = faces.size();

  // calculate areas of possibly curved faces
  std::vector<double> areas(num_faces, 0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // back to area-weighted fluxes
  for (int i = 0; i < num_faces; i++) {
    for (int j = 0; j < num_faces; j++) W(i, j) *= areas[i] * areas[j];
  }
}


/* ******************************************************************
* A simple monotone stability term for a 2D or 3D brick element. 
****************************************************************** */
int MFD3D::StabilityMonotoneHex(int cell, const Tensor& T,
                                Teuchos::SerialDenseMatrix<int, double>& Mc,
                                Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();
  int nrows = 2 * d;

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j) = Mc(i, j);
  }

  // create groups of quasi-parallel faces
  int map[nrows];
  for (int i = 0; i < nrows; i++) map[i] = i;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int i1, i2, k, l;
  double s1, s2, area1, area2;
  for (int i = 0; i < d-1; i++) {
    int i1 = 2*i;
    int i2 = i1 + 1;

    int f = faces[i1];
    const AmanziGeometry::Point& normal1 = mesh_->face_normal(f);
    area1 = mesh_->face_area(f);

    s1 = 1.0;
    for (int j = i2; j < nrows; j++) {
      int f = faces[j];
      const AmanziGeometry::Point& normal2 = mesh_->face_normal(f);
      area2 = mesh_->face_area(f);

      s2 = (normal1 * normal2) * (dirs[i] * dirs[j]) / (area1 * area2);
      if (s2 < s1) {  // swap map values in positions i2 and j
        k = map[i2];
        map[i2] = j;
        map[j] = k;
        s1 = s2;
      }
    }
    // if (s1 >= 0.0) return -1;  // hex is too disturted
  }

  // define transformed tensor
  Tensor T1(d, 2);
  for (int i = 0; i < d; i++) {
    k = map[2*i];
    int f = faces[k];
    const AmanziGeometry::Point& normal1 = mesh_->face_normal(f);
    area1 = mesh_->face_area(f);

    for (int j = i; j < d; j++) {
      l = map[2*j];
      f = faces[l];
      const AmanziGeometry::Point& normal2 = mesh_->face_normal(f);
      area2 = mesh_->face_area(f);

      s1 = (T * normal1) * normal2 * (dirs[k] * dirs[l]) / (area1 * area2);
      if (i-j)
        T1(i, j) = T1(j, i) = -fabs(s1);
      else
        T1(i, i) = s1;
    }
  }

  // add stability term D_ik T1_kl D_il
  double volume = mesh_->cell_volume(cell);
  for (int i = 0; i < nrows; i++) {
    i1 = i / 2;
    k = map[i];
    area1 = mesh_->face_area(faces[k]);

    for (int j = i; j < nrows; j++) {
      i2 = j / 2;
      l = map[j];
      area2 = mesh_->face_area(faces[l]);
      M(l, k) = M(k, l) += T1(i1, i2) * area1 * area2 / volume;  // Fix (lipnikov@lanl.gov)
    }
  }
  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi



