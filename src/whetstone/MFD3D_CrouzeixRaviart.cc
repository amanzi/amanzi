/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#include <cmath>
#include <tuple>
#include <vector>

// Amanzi
#include "MeshLight.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SurfaceMeshLight.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_CrouzeixRaviart::MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                                             const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
  : BilinearForm(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_ + 1);
  R_.Reshape(nfaces, d_ + 1);
  Ac.Reshape(nfaces, nfaces);

  // calculate matrix R
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int k = 0; k < d_; k++) R_(n, k) = normal[k] * dirs[n];
    R_(n, d_) = 0.0;
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->cell_volume(c);
  for (int n = 0; n < nfaces; ++n) {
    for (int k = 0; k < d_; k++) v1[k] = R_(n, k);
    v2 = K * v1;

    for (int m = n; m < nfaces; m++) {
      for (int k = 0; k < d_; k++) v1[k] = R_(m, k);
      Ac(n, m) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    for (int k = 0; k < d_; k++) N(n, k) = xf[k] - xc[k];
    N(n, d_) = 1.0;  // additional column is added to the consistency condition
  }

  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. 
****************************************************************** */
int MFD3D_CrouzeixRaviart::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Projector on a manifold
****************************************************************** */
void MFD3D_CrouzeixRaviart::H1Face(int f, const std::vector<Polynomial>& ve,
                                   const Polynomial* moments, Polynomial& vf)
{
  const auto& origin = mesh_->face_centroid(f);
  const auto& normal = mesh_->face_normal(f);
  SurfaceCoordinateSystem coordsys(origin, normal);

  Teuchos::RCP<const SurfaceMeshLight> surf_mesh = Teuchos::rcp(new SurfaceMeshLight(mesh_, f, coordsys));
  ProjectorCell_<SurfaceMeshLight>(surf_mesh, f, ve, ve, vf);
}

}  // namespace WhetStone
}  // namespace Amanzi

