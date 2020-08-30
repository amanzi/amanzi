/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
  : BilinearForm(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
****************************************************************** */
int MFD3D_Lagrange::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  N.Reshape(nnodes, d_ + 1);
  Ac.Reshape(nnodes, nnodes);

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  double volume = mesh_->cell_volume(c);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  // calculate matrix R. We temporaryly reuse matrix N
  N.PutScalar(0.0);

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

      if (d_ == 2) {
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

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    for (int k = 0; k < d_; k++) N(i, k) = p[k] - cm[k];
    N(i, d_) = 1.0;  // additional column is added to the consistency condition
  }

  return 0;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_Lagrange::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
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
void MFD3D_Lagrange::ProjectorCell_(
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf, Polynomial& uc)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  // populate matrix R (should be a separate routine lipnikov@lanl.gv)
  DenseMatrix R(nnodes, d_);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  R.PutScalar(0.0);
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

      if (d_ == 2) {
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

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      for (int k = 0; k < d_; k++) R(pos, k) += normal[k] * u;
    }
  }

  uc.Reshape(d_, 1, true);
  for (int i = 0; i < nnodes; i++) {
    for (int k = 0; k < d_; k++) {
      uc(k + 1) += R(i, k) * vf[i](0);
    }
  }
  uc *= 1.0 / mesh_->cell_volume(c);
}

}  // namespace WhetStone
}  // namespace Amanzi

