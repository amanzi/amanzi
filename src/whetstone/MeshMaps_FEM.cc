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
* Calculate mesh velocity in cell c. 
****************************************************************** */
void MeshMaps_FEM::VelocityCell(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  Entity_ID_List nodes;

  mesh1_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  ASSERT(nnodes == 4);

  AmanziGeometry::Point p1(d_), p2(d_), p3(d_), p4(d_);
  mesh1_->node_get_coordinates(nodes[0], &p1);
  mesh1_->node_get_coordinates(nodes[1], &p2);
  mesh1_->node_get_coordinates(nodes[2], &p3);
  mesh1_->node_get_coordinates(nodes[3], &p4);

  AmanziGeometry::Point p21(p2 - p1), p41(p4 - p1), p32(p3 - p2);
  AmanziGeometry::Point pp(p32 - p41);

  // By our assumption, the Jacobian is constant and diagonal
  AmanziGeometry::Point xref(0.0, 0.0);
  Tensor J0 = JacobianValueInternal_(mesh0_, c, xref);
  J0.Inverse();

  //calculate map in coordinate system centered at point 1
  vc.resize(d_);
  for (int i = 0; i < d_; ++i) {
    vc[i].Reshape(d_, 2);
    vc[i](0, 0) = p1[i];
    vc[i](1, 0) = p21[i] * J0(0, 0);
    vc[i](1, 1) = p41[i] * J0(1, 1);
    vc[i](2, 1) = pp[i] * J0(0, 0) * J0(1, 1);
  }

  // rebase polynomial to global coordinate system
  mesh0_->node_get_coordinates(nodes[0], &p1);
  AmanziGeometry::Point zero(0.0, 0.0);
  for (int i = 0; i < d_; ++i) {
    vc[i].set_origin(p1);
    vc[i].ChangeOrigin(zero);
  }

  // calculate velocity u(X) = F(X) - X
  for (int i = 0; i < d_; ++i) {
    vc[i](1, i) -= 1.0;
  }
}


/* ******************************************************************
* Transformation of normal is defined completely by face data.
****************************************************************** */
void MeshMaps_FEM::NansonFormula(
    int f, double t, const VectorPolynomial& vc, VectorPolynomial& cn) const
{
  Entity_ID_List cells;
  mesh0_->face_get_cells(f, AmanziMesh::USED, &cells);

  WhetStone::MatrixPolynomial C;
  Cofactors(cells[0], t, vc, C);

  AmanziGeometry::Point normal = mesh0_->face_normal(f);
  cn[0].Multiply(C, normal, cn, false);
}


/* ******************************************************************
* Calculation of matrix of cofactors
****************************************************************** */
void MeshMaps_FEM::Cofactors(
    int c, double t, const VectorPolynomial& vc, MatrixPolynomial& C) const
{
  Entity_ID_List nodes;

  mesh1_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  ASSERT(nnodes == 4);

  AmanziGeometry::Point p1(d_), p2(d_), p3(d_), p4(d_);
  mesh1_->node_get_coordinates(nodes[0], &p1);
  mesh1_->node_get_coordinates(nodes[1], &p2);
  mesh1_->node_get_coordinates(nodes[2], &p3);
  mesh1_->node_get_coordinates(nodes[3], &p4);

  AmanziGeometry::Point p21(p2 - p1), p41(p4 - p1), p32(p3 - p2), p34(p3 - p4);
  
  // By our assumption, the Jacobian is constant and diagonal
  AmanziGeometry::Point xref(0.0, 0.0);
  Tensor J0 = JacobianValueInternal_(mesh0_, c, xref);
  J0.Inverse();

  // allocate memory for matrix
  C.resize(d_);
  for (int i = 0; i < d_; ++i) {
    C[i].resize(2);
    for (int j = 0; j < d_; ++j) {
      C[i][j].Reshape(d_, 1);
    }
  } 

  // constant part
  C[0][0](0, 0) = p41[1];
  C[1][0](0, 0) = -p41[0];

  C[0][1](0, 0) = -p21[1];
  C[1][1](0, 0) = p21[0];

  // linear part
  C[0][0](1, 0) = (p32[1] - p41[1]) * J0(1, 1);
  C[1][0](1, 0) =-(p32[0] - p41[0]) * J0(0, 0);

  C[0][1](1, 1) =-(p34[1] - p21[1]) * J0(0, 0);
  C[1][1](1, 1) = (p34[0] - p21[0]) * J0(1, 1);

  // rebase polynomials to global coordinate system
  mesh0_->node_get_coordinates(nodes[0], &p1);
  AmanziGeometry::Point zero(0.0, 0.0);
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      C[i][j].set_origin(p1);
      C[i][j].ChangeOrigin(zero);
    }
  }

  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      C[i][j] *= t * J0(j, j);
    }
    C[i][i](0, 0) += 1.0 - t;
  }
}


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
* NOTE: We assume that cell c is a rectangle on mesh 0.
****************************************************************** */
void MeshMaps_FEM::JacobianDet(
    int c, double t, const std::vector<VectorPolynomial>& vf, Polynomial& vc) const
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

  // By our assumption, the Jacobian is constant and diagonal
  AmanziGeometry::Point xref(0.0, 0.0);
  Tensor J0 = JacobianValueInternal_(mesh0_, c, xref);
  J0.Inverse();

  mesh0_->node_get_coordinates(nodes[0], &p1);
  a0 -= b0 * (p1[1] * J0(1, 1));
  a1 -= b1 * (p1[0] * J0(0, 0));
  b0 *= J0(1, 1);
  b1 *= J0(0, 0);

  a0 *= t * J0(0, 0);
  a1 *= t * J0(1, 1);
  b0 *= t * J0(0, 0);
  b1 *= t * J0(1, 1);

  a0[0] += 1.0 - t;
  a1[1] += 1.0 - t;

  vc.Reshape(2, 2);
  vc(0, 0) = a0[0] * a1[1] - a0[1] * a1[0];

  vc(1, 0) = a0[0] * b1[1] - a0[1] * b1[0];
  vc(1, 1) = b0[0] * a1[1] - b0[1] * a1[0];

  vc(2, 1) = b0[0] * b1[1] - b0[1] * b1[0];
}


/* ******************************************************************
* Calculate Jacobian at point x of face f.
****************************************************************** */
void MeshMaps_FEM::JacobianFaceValue(
    int f, const VectorPolynomial& v,
    const AmanziGeometry::Point& x, Tensor& J) const
{
  Entity_ID_List faces, cells;
  mesh0_->face_get_cells(f, AmanziMesh::USED, &cells);

  int c(cells[0]);
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

