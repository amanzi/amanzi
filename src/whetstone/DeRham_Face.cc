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
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  const auto faces_dirs = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces_dirs.first.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  Tensor Kinv(K);
  Kinv.Inverse();
  Kinv.Transpose();

  for (int i = 0; i < nfaces; i++) {
    int f = faces_dirs.first[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    double a1 = mesh_->getFaceArea(f);
    v2 = Kinv * (fm - cm);

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < nfaces; j++) {
      f = faces_dirs.first[j];
      const AmanziGeometry::Point& fm2 = mesh_->getFaceCentroid(f);
      double a2 = mesh_->getFaceArea(f);
      v1 = fm2 - cm;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / volume;
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces_dirs.first[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);
    for (int k = 0; k < d_; k++) N(i, k) = normal[k] * faces_dirs.second[i] / area;
  }
  return 0;
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
  return 0;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix.
* Only the upper triangular part of Wc is calculated.
****************************************************************** */
int DeRham_Face::L2consistencyInverse(
    int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry)
{
  const auto faces_dirs = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces_dirs.first.size();

  R.Reshape(nfaces, d_);
  Wc.Reshape(nfaces, nfaces);

  // calculate areas of possibly curved faces
  std::vector<double> areas(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces_dirs.first[i];
    areas[i] = norm(mesh_->getFaceNormal(f));
  }

  // populate matrix W_0
  AmanziGeometry::Point v1(d_);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces_dirs.first[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    v1 = K * normal;

    for (int j = i; j < nfaces; j++) {
      f = faces_dirs.first[j];
      const AmanziGeometry::Point& v2 = mesh_->getFaceNormal(f);
      Wc(i, j) = (v1 * v2) / (faces_dirs.second[i] * faces_dirs.second[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces_dirs.first[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    for (int k = 0; k < d_; k++) R(i, k) = (fm[k] - cm[k]) * areas[i];
  }

  /* Internal verification 
  DenseMatrix NtR(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      NtR(i, j) = 0.0;
      for (int k = 0; k < nfaces; k++) {
        const AmanziGeometry::Point& v1 = mesh_->getFaceNormal(faces_dirs.first[k]);
        NtR(i, j) += v1[i] * R(k, j) / areas[k] * faces_dirs.second[k];
      }
    }
  }
  */
  return 0;
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

}  // namespace WhetStone
}  // namespace Amanzi


