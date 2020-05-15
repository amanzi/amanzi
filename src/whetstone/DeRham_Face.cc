/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products on faces.
*/

#include "Mesh.hh"

#include "DeRham_Face.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Only upper triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
****************************************************************** */
int DeRham_Face::L2consistency(
    int c, const Tensor<>& K, DenseMatrix<>& N, DenseMatrix<>& Mc, bool symmetry)
{
  AmanziMesh::Entity_ID_View faces;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.size();

  N.reshape(nfaces, d_);
  Mc.reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  Tensor<> Kinv(K);
  Kinv.Inverse();
  Kinv.Transpose();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double a1 = mesh_->face_area(f);
    v2 = Kinv * (fm - cm);

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
      double a2 = mesh_->face_area(f);
      v1 = fm - cm;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / volume;
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);
    for (int k = 0; k < d_; k++) N(i, k) = normal[k] * dirs[i] / area;
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix: adding stability matrix to the consistency matrix.
****************************************************************** */
int DeRham_Face::MassMatrix(int c, const Tensor<>& K, DenseMatrix<>& M)
{
  DenseMatrix<> N;

  Tensor<> Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix.
* Only the upper triangular part of Wc is calculated.
****************************************************************** */
int DeRham_Face::L2consistencyInverse(
    int c, const Tensor<>& K, DenseMatrix<>& R, DenseMatrix<>& Wc, bool symmetry)
{
  AmanziMesh::Entity_ID_View faces;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.size();

  R.reshape(nfaces, d_);
  Wc.reshape(nfaces, nfaces);

  // calculate areas of possibly curved faces
  std::vector<double> areas(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // populate matrix W_0
  AmanziGeometry::Point v1(d_);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = K * normal;

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d_; k++) R(i, k) = (fm[k] - cm[k]) * areas[i];
  }

  /* Internal verification 
  DenseMatrix NtR(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      NtR(i, j) = 0.0;
      for (int k = 0; k < nfaces; k++) {
        const AmanziGeometry::Point& v1 = mesh_->face_normal(faces[k]);
        NtR(i, j) += v1[i] * R(k, j) / areas[k] * dirs[k];
      }
    }
  }
  */
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix: adding stability to the consistency
****************************************************************** */
int DeRham_Face::MassMatrixInverse(int c, const Tensor<>& K, DenseMatrix<>& W)
{
  DenseMatrix<> R;

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return ok;
 
  StabilityScalar_(R, W);

  return ok;
}

}  // namespace WhetStone
}  // namespace Amanzi

