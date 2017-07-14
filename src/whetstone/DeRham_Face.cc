/*
  WhetStone, version 2.1
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
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Only upper triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
****************************************************************** */
int DeRham_Face::L2consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  Tensor Kinv(K);
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
int DeRham_Face::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
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
    int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_);
  Wc.Reshape(nfaces, nfaces);

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
int DeRham_Face::MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return ok;
 
  StabilityScalar_(R, W);

  return ok;
}


/* ******************************************************************
* Consistency condition for inner product on a generized polyhedron.
****************************************************************** */
int DeRham_Face::L2consistencyGeneralized(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List faces, nodes;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int nfaces = faces.size();
  int nx(d_ * nfaces);

  N.Reshape(nx, d_);
  Mc.Reshape(nx, nx);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  AmanziGeometry::Point v1(d_), v2(d_);
  std::vector<AmanziGeometry::Point> vv(3), xm(3);

  // populate matrices R and N
  DenseMatrix R(N);
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);  
    double area_div = norm(normal);

    CurvedFaceGeometry_(f, dirs[i], vv, xm);

    for (int k = 0; k < d_; ++k) {
      R(d_ * i, k) = area * xm[0][k] - area_div * xc[k];
      R(d_ * i + 1, k) = area * xm[1][k];
      R(d_ * i + 2, k) = area * xm[2][k];

      for (int l = 0; l < d_; ++l) {
        N(d_ * i + l, k) = vv[l][k];
      }
    }
  }

  // upper triangular part of the consistency term
  Tensor Kinv(K);
  Kinv.Inverse();

  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < d_; ++k) v1[k] = R(i, k);
    v2 = Kinv * v1;

    for (int j = i; j < nx; ++j) {
      for (int k = 0; k < d_; ++k) v1[k] = R(j, k);
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix for genelized polyhedron
****************************************************************** */
int DeRham_Face::MassMatrixGeneralized(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistencyGeneralized(c, Kinv, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of inner product on a generized 
* polyhedron.
****************************************************************** */
int DeRham_Face::L2consistencyInverseGeneralized(
    int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry)
{
  Entity_ID_List faces, nodes;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int nfaces = faces.size();
  int nx(d_ * nfaces);

  R.Reshape(nx, d_);
  Wc.Reshape(nx, nx);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  AmanziGeometry::Point v1(d_), v2(d_);
  std::vector<AmanziGeometry::Point> vv(3), xm(3);

  // populate matrices R and N
  DenseMatrix N(R);
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);  
    double area_div = norm(normal);

    CurvedFaceGeometry_(f, dirs[i], vv, xm);

    for (int k = 0; k < d_; ++k) {
      R(d_ * i, k) = area * xm[0][k] - area_div * xc[k];
      R(d_ * i + 1, k) = area * xm[1][k];
      R(d_ * i + 2, k) = area * xm[2][k];

      for (int l = 0; l < d_; ++l) {
        N(d_ * i + l, k) = vv[l][k];
      }
    }
  }

  // upper triangular part of the consistency term
  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nx; ++j) {
      for (int k = 0; k < d_; ++k) v1[k] = N(j, k);
      Wc(i, j) = (v1 * v2) / volume;
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix for generalized polyhedron
****************************************************************** */
int DeRham_Face::MassMatrixInverseGeneralized(
    int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseGeneralized(c, K, R, W, true);
  if (ok) return ok;

  StabilityScalar_(R, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi


