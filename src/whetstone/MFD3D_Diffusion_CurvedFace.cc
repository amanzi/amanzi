/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The mimetic finite difference method for diffusion.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "MeshLight.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion_CurvedFace.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Only upper triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
****************************************************************** */
int
MFD3D_Diffusion_CurvedFace::L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc)
{
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // populate areas
  DenseVector area(nfaces);
  for (int i = 0; i < nfaces; ++i) area(i) = norm(mesh_->face_normal(faces[i]));

  Tensor Kinv(K);
  Kinv.Inverse();

  for (int i = 0; i < nfaces; i++) {
    double a1 = area(i);
    v2 = Kinv * ((*bf_)[faces[i]] - cm);

    for (int j = i; j < nfaces; j++) {
      double a2 = area(j);
      v1 = (*bf_)[faces[j]] - cm;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / volume;
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    for (int k = 0; k < d_; k++) N(i, k) = normal[k] * dirs[i] / area(i);
  }
  return 0;
}


/* ******************************************************************
* Mass matrix: adding stability matrix to the consistency matrix.
*e***************************************************************** */
int
MFD3D_Diffusion_CurvedFace::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Inverse mass matrix in flux space via optimization, experimental.
****************************************************************** */
int
MFD3D_Diffusion_CurvedFace::MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
{
  int ok = MassMatrix(c, K, W);
  W.Inverse();
  return ok;
}


/* ******************************************************************
* Stiffness matrix is calculated by a hybridization algorithm.
****************************************************************** */
int
MFD3D_Diffusion_CurvedFace::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix M;
  MassMatrixInverse(c, K, M);

  const auto& faces = mesh_->cell_get_faces(c);
  int nfaces = faces.size();

  // populate areas
  DenseVector area(nfaces);
  for (int i = 0; i < nfaces; ++i) area(i) = norm(mesh_->face_normal(faces[i]));

  // populate stiffness matrix
  A.Reshape(nfaces + 1, nfaces + 1);

  double cntr(0.0);
  for (int i = 0; i < nfaces; ++i) {
    for (int j = 0; j < nfaces; ++j) { A(i, j) = M(i, j) * area(i) * area(j); }

    double add(0.0);
    for (int j = 0; j < nfaces; ++j) { add -= M(i, j) * area(j); }
    A(nfaces, i) = A(i, nfaces) = add * area(i);

    cntr -= add * area(i);
  }
  A(nfaces, nfaces) = cntr;

  return 0;
}


/* *****************************************************************
* Low-order L2 projector.
* NOTE: we abuse the interface and return a linear polynomial.
***************************************************************** */
void
MFD3D_Diffusion_CurvedFace::L2Cell(int c,
                                   const std::vector<Polynomial>& ve,
                                   const std::vector<Polynomial>& vf,
                                   const Polynomial* moments,
                                   Polynomial& vc)
{
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  vc.Reshape(d_, 1, true);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);

    for (int k = 0; k < d_; k++) {
      double Rik = fm[k] - cm[k];
      vc(k + 1) += Rik * vf[i](0) * dirs[i];
    }
  }

  vc *= -1.0 / mesh_->cell_volume(c);
}

} // namespace WhetStone
} // namespace Amanzi
