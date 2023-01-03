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

  The mimetic finite difference method on 2D surfaces in 3D.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of
* fluxes for a non-flat surface. Only the upper triangular part of
* Wc is calculated. Darcy flux is scaled by the area!
* WARNING: routine works for scalar T only.
****************************************************************** */
int
MFD3D_Diffusion::L2consistencyInverseSurface(int c,
                                             const Tensor& T,
                                             DenseMatrix& R,
                                             DenseMatrix& Wc)
{
  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_ - 1);
  Wc.Reshape(nfaces, nfaces);

  double volume = mesh_->getCellVolume(c);

  // calculate cell normal
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  const AmanziGeometry::Point& xf1 = mesh_->getFaceCentroid(faces[0]);
  const AmanziGeometry::Point& xf2 = mesh_->getFaceCentroid(faces[1]);
  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  v1 = (xf1 - xc) ^ (xf2 - xc);
  v1 /= norm(v1);

  // calculate projector
  Tensor P(d_, 2);
  for (int i = 0; i < d_; i++) {
    P(i, i) = 1.0;
    for (int j = 0; j < d_; j++) { P(i, j) -= v1[i] * v1[j]; }
  }

  // cell-based coordinate system
  v2 = xf1 - xc;
  v2 /= norm(v2);
  v3 = v1 ^ v2;

  // define new tensor
  Tensor PTP(d_, 2);
  PTP = P * T * P;

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_getFaceNormal(f, c);

    v1 = PTP * normal;

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v4 = mesh_getFaceNormal(f, c);
      Wc(i, j) = (v1 * v4) / volume;
    }
  }

  // calculate matrix R
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);

    R(i, 0) = v2 * (fm - cm);
    R(i, 1) = v3 * (fm - cm);
  }

  return 0;
}


/* ******************************************************************
* Darcy inverse mass matrix for surface: the standard algorithm
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseSurface(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseSurface(c, K, R, W);
  if (ok) return ok;

  StabilityOptimized_(K, R, W);
  return 0;
}


/* ******************************************************************
* Mass matrix for a polyhedral element via simplex method.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseSurfaceMMatrix(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  // use boolean flag to populate the whole matrix
  int ok = L2consistencyInverseSurface(c, K, R, W);
  if (ok) return ok;

  // scaling of matrix W for numerical stability
  double s = W.Trace() / W.NumRows();
  W /= s;

  ok = StabilityMMatrix_(c, R, W);
  if (ok) return ok;

  W *= s;
  simplex_functional_ *= s;
  return 0;
}


/* ******************************************************************
* Exterior normal to 2D face in 3D space.
****************************************************************** */
AmanziGeometry::Point
MFD3D_Diffusion::mesh_getFaceNormal(int f, int c)
{
  AmanziGeometry::Point v0(d_), v1(d_);
  Entity_ID_List nodes;

  nodes = mesh_->getFaceNodes(f);
  v0 = mesh_->getNodeCoordinate(nodes[0]);
  v1 = mesh_->getNodeCoordinate(nodes[1]);

  AmanziGeometry::Point tau(v1 - v0);
  AmanziGeometry::Point normal = v0 - mesh_->getCellCentroid(c);

  // orthogonalize and rescale normal
  double len = norm(tau);
  double s = (normal * tau) / len / len;
  normal -= s * tau;
  normal *= len / norm(normal);

  return normal;
}


/* ******************************************************************
* The conventional FV scheme for a general mesh.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseSurfaceTPFA(int c, const Tensor& K, DenseMatrix& W)
{
  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  W.Reshape(nfaces, nfaces);
  W.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  AmanziGeometry::Point a(d_);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f, c);

    a = xf - xc;
    double s = mesh_->getFaceArea(f) / norm(a);
    double Knn = ((K * a) * normal) * s;
    double dxn = a * normal;
    W(n, n) = Knn / fabs(dxn);
  }

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
