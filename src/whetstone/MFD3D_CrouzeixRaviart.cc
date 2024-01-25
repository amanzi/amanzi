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

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#include <cmath>
#include <tuple>
#include <vector>

// Amanzi
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SingleFaceMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_CrouzeixRaviart::MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Consistency condition for stiffness matrix.
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int
MFD3D_CrouzeixRaviart::H1consistency(int c, const Tensor<>& K, DenseMatrix<>& N, DenseMatrix<>& Ac)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  N.reshape(nfaces, d_ + 1);
  R_.reshape(nfaces, d_ + 1);
  Ac.reshape(nfaces, nfaces);

  // calculate matrix R
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    for (int k = 0; k < d_; k++) R_(n, k) = normal[k] * dirs[n];
    R_(n, d_) = 0.0;
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->getCellVolume(c);
  for (int n = 0; n < nfaces; ++n) {
    for (int k = 0; k < d_; k++) v1[k] = R_(n, k);
    v2 = K * v1;

    for (int m = n; m < nfaces; m++) {
      for (int k = 0; k < d_; k++) v1[k] = R_(m, k);
      Ac(n, m) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    for (int k = 0; k < d_; k++) N(n, k) = xf[k] - xc[k];
    N(n, d_) = 1.0; // additional column is added to the consistency condition
  }

  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int
MFD3D_CrouzeixRaviart::StiffnessMatrix(int c, const Tensor<>& K, DenseMatrix<>& A)
{
  DenseMatrix<> N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  MFD3D::StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Projector on a manifold
****************************************************************** */
void
MFD3D_CrouzeixRaviart::H1Face(int f,
        const std::vector<Polynomial<>>& ve,
        const Polynomial<>* moments,
        Polynomial<>& vf)
{
  const auto& origin = mesh_->getFaceCentroid(f);
  const auto& normal = mesh_->getFaceNormal(f);
  AmanziGeometry::SurfaceCoordinateSystem coordsys(origin, normal);

  Teuchos::RCP<AmanziMesh::SingleFaceMesh> surf_mesh =
    Teuchos::rcp(new AmanziMesh::SingleFaceMesh(mesh_, f, coordsys));
  ProjectorCell_(surf_mesh, f, ve, ve, vf);
}


/* ******************************************************************
* Energy projector on the space of linear polynomials in cell c.
****************************************************************** */
void
MFD3D_CrouzeixRaviart::ProjectorCell_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
        int c,
        const std::vector<Polynomial<>>& ve,
        const std::vector<Polynomial<>>& vf,
        Polynomial<>& uc)
{
  const auto& [faces,dirs] = mymesh->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mymesh->getCellCentroid(c);
  double vol = mymesh->getCellVolume(c);

  // create zero vector polynomial
  uc.reshape(d_, 1, true);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mymesh->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mymesh->getFaceNormal(f);

    double tmp = vf[n].Value(xf) * dirs[n] / vol;

    for (int j = 0; j < d_; ++j) { uc(1, j) += tmp * normal[j]; }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_);
  for (int j = 0; j < d_; ++j) { grad[j] = uc(1, j); }

  double a1(0.0), a2(0.0), tmp;
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mymesh->getFaceCentroid(f);
    double area = mymesh->getFaceArea(f);

    tmp = vf[n].Value(xf) - grad * (xf - xc);
    a1 += tmp * area;
    a2 += area;
  }

  uc(0) = a1 / a2;

  // set the correct origin
  uc.set_origin(xc);
}

} // namespace WhetStone
} // namespace Amanzi
