/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two states 
  of a deformable mesh: finite element implementation.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps_FEM.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculation of Jacobian for linearized map xi + t (F(xi) - xi).
****************************************************************** */
void MeshMaps_FEM::JacobianCellValue(
    int c, double t, const AmanziGeometry::Point& xref, Tensor& J) const
{
  J = JacobianValueInternal_(mesh1_, c, xref);
 
  // convolution of two maps.
  Tensor J0 = JacobianValueInternal_(mesh0_, c, xref);
  J0.Inverse();
  J = J * J0;

  J *= t;
  J += 1.0 - t;
}


/* ******************************************************************
* Calculate determinant of the Jacobian at time t.
* NOTE: We assume that cell c is a parallepiped.
****************************************************************** */
void MeshMaps_FEM::JacobianDet(int c, double t, Polynomial& v) const
{
  Entity_ID_List nodes;

  mesh1_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  ASSERT(nnodes == 4);

  AmanziGeometry::Point p1(d_), p2(d_), p3(d_), p4(d_), j0(d_), j1(d_);
  mesh1_->node_get_coordinates(nodes[0], &p1);
  mesh1_->node_get_coordinates(nodes[1], &p2);
  mesh1_->node_get_coordinates(nodes[2], &p3);
  mesh1_->node_get_coordinates(nodes[3], &p4);

  AmanziGeometry::Point a0(p2 - p1), a1(p4 - p1);
  AmanziGeometry::Point b0(p3 - p4 - a0), b1(p3 - p2 - a1);

  // By our assumption, the Jacobian is constant
  AmanziGeometry::Point xref(0.0, 0.0);
  Tensor J0 = JacobianValueInternal_(mesh0_, c, xref);
  J0.Inverse();

  a0 *= t * J0(0, 0);
  a1 *= t * J0(1, 1);
  b0 *= t * J0(0, 0);
  b1 *= t * J0(1, 1);

  a0[0] += 1.0 - t;
  a1[1] += 1.0 - t;

  v.Reshape(2, 2);
  v.monomials(0).coefs()[0] = a0[0] * a1[1] - a0[1] * a1[0];

  v.monomials(1).coefs()[0] = a0[0] * b1[1] - a0[1] * b1[0];
  v.monomials(1).coefs()[1] = b0[0] * a1[1] - b0[1] * a1[0];

  v.monomials(2).coefs()[1] = b0[0] * b1[1] - b0[1] * b1[0];
}


/* ******************************************************************
* Calculate mesh velocity on face f.
****************************************************************** */
void MeshMaps_FEM::VelocityFace(int c, int f, std::vector<Polynomial>& v) const
{
  AmanziMesh::Entity_ID_List nodes;
  AmanziGeometry::Point x0, x1;

  const AmanziGeometry::Point& xf0 = mesh0_->face_centroid(f);
  const AmanziGeometry::Point& xf1 = mesh1_->face_centroid(f);

  // velocity order 0
  for (int i = 0; i < d_; ++i) {
    v[i].Reset();
    v[i].monomials(0).coefs()[0] = xf1[i] - xf0[i];
  }

  // velocity order 1 (2D algorithm)
  mesh0_->face_get_nodes(f, &nodes);
  mesh0_->node_get_coordinates(nodes[0], &x0);
  mesh1_->node_get_coordinates(nodes[0], &x1);

  x0 -= xf0;
  x1 -= xf1;

  WhetStone::Tensor A(2, 2);
  AmanziGeometry::Point b(2);

  A(0, 0) = x0[0];
  A(0, 1) = A(1, 0) = x0[1];
  A(1, 1) = -x0[0];

  A.Inverse();
  b = A * (x1 - x0);

  v[0].monomials(1).coefs() = { b[0], b[1]};
  v[1].monomials(1).coefs() = {-b[1], b[0]};

  // we change to the global coordinate system
  v[0].monomials(0).coefs()[0] -= b * xf0;
  v[1].monomials(0).coefs()[0] -= (b^xf0)[0];
}


/* ******************************************************************
* Calculate Jacobian at point x of face f.
****************************************************************** */
void MeshMaps_FEM::JacobianFaceValue(
    int c, int f, const std::vector<Polynomial>& v,
    const AmanziGeometry::Point& x, Tensor& J) const
{
  Entity_ID_List faces;
  mesh0_->cell_get_faces(c, &faces);
  
  AmanziGeometry::Point xref;
  for (int n = 0; n < faces.size(); n++) {
    if (faces[n] == f) {
      if (n == 0) { 
        xref.set(0.5, 0.0);
      } else if (n == 1) {
        xref.set(1.0, 0.5);
      } else if (n == 2) { 
        xref.set(0.5, 1.0);
      } else if (n == 3) {
        xref.set(0.0, 0.5);
      }
      break;
    }
  }

  JacobianCellValue(c, 1.0, xref, J);
}


/* ******************************************************************
* Supporting routine for calculating Jacobian.
****************************************************************** */
Tensor MeshMaps_FEM::JacobianValueInternal_(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
    int c, const AmanziGeometry::Point& xref) const
{
  Entity_ID_List nodes;

  mesh->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  ASSERT(nnodes == 4);

  AmanziGeometry::Point p1(d_), p2(d_), p3(d_), p4(d_), j0(d_), j1(d_);
  mesh->node_get_coordinates(nodes[0], &p1);
  mesh->node_get_coordinates(nodes[1], &p2);
  mesh->node_get_coordinates(nodes[2], &p3);
  mesh->node_get_coordinates(nodes[3], &p4);

  j0 = (1.0 - xref[1]) * (p2 - p1) + xref[1] * (p3 - p4);
  j1 = (1.0 - xref[0]) * (p4 - p1) + xref[0] * (p3 - p2);

  Tensor jac(d_, 2);
  jac.SetColumn(0, j0);
  jac.SetColumn(1, j1);

  return jac;
}


/* ******************************************************************
* Bilinear map (2D algorithm)
****************************************************************** */
AmanziGeometry::Point MeshMaps_FEM::Map_(int c, const AmanziGeometry::Point& xref) const
{
  Entity_ID_List nodes;

  mesh1_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  ASSERT(nnodes == 4);

  AmanziGeometry::Point p1(d_), p2(d_), p3(d_), p4(d_), f(d_);
  mesh1_->node_get_coordinates(nodes[0], &p1);
  mesh1_->node_get_coordinates(nodes[1], &p2);
  mesh1_->node_get_coordinates(nodes[2], &p3);
  mesh1_->node_get_coordinates(nodes[3], &p4);

  double x(xref[0]), y(xref[1]);
  f = (1.0 - x) * (1.0 - y) * p1 + x * (1.0 - y) * p2
    + x * y * p3 + (1.0 - x) * y * p4;

  return f;
}


}  // namespace WhetStone
}  // namespace Amanzi

