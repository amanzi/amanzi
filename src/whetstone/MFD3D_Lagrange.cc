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

  Lagrange-type element: degrees of freedom are nodal values.
*/

#include <cmath>
#include <tuple>
#include <vector>

// Amanzi
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_Lagrange::MFD3D_Lagrange(const Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* High-order consistency condition for the stiffness matrix.
****************************************************************** */
int
MFD3D_Lagrange::H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  N.Reshape(nnodes, d_ + 1);
  Ac.Reshape(nnodes, nnodes);

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

  double volume = mesh_->getCellVolume(c);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  // calculate matrix R. We temporaryly reuse matrix N
  N.PutScalar(0.0);

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    double area = mesh_->getFaceArea(f);

    auto face_nodes = mesh_->getFaceNodes(f);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs[i];
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        p = mesh_->getNodeCoordinate(v);
        pnext = mesh_->getNodeCoordinate(vnext);
        pprev = mesh_->getNodeCoordinate(vprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      for (int k = 0; k < d_; k++) N(pos, k) += normal[k] * u;
    }
  }

  // calculate upper part of R K R^T / volume
  for (int i = 0; i < nnodes; i++) {
    for (int k = 0; k < d_; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nnodes; j++) {
      for (int k = 0; k < d_; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh_->getNodeCoordinate(v);
    for (int k = 0; k < d_; k++) N(i, k) = p[k] - cm[k];
    N(i, d_) = 1.0; // additional column is added to the consistency condition
  }

  return 0;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int
MFD3D_Lagrange::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* *****************************************************************
* N and R are used. Here we use simplified versions.
***************************************************************** */
void
MFD3D_Lagrange::ProjectorCell_(int c,
                               const std::vector<Polynomial>& ve,
                               const std::vector<Polynomial>& vf,
                               Polynomial& uc)
{
  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int num_faces = faces.size();

  // populate matrix R (should be a separate routine lipnikov@lanl.gv)
  DenseMatrix R(nnodes, d_);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  R.PutScalar(0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    double area = mesh_->getFaceArea(f);

    auto face_nodes = mesh_->getFaceNodes(f);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs[i];
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        p = mesh_->getNodeCoordinate(v);
        pnext = mesh_->getNodeCoordinate(vnext);
        pprev = mesh_->getNodeCoordinate(vprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      for (int k = 0; k < d_; k++) R(pos, k) += normal[k] * u;
    }
  }

  uc.Reshape(d_, 1, true);
  for (int i = 0; i < nnodes; i++) {
    for (int k = 0; k < d_; k++) { uc(k + 1) += R(i, k) * vf[i](0); }
  }
  uc *= 1.0 / mesh_->getCellVolume(c);
}

} // namespace WhetStone
} // namespace Amanzi
