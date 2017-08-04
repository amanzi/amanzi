/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin methods.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity on face f.
****************************************************************** */
void MeshMaps::FaceVelocity(int c, int f, std::vector<Polynomial>& v) const
{
  ASSERT(mesh0_ != Teuchos::null);
  VEM_FaceVelocity_(c, f, v);
}


/* ******************************************************************
* Calculate mesh velocity on face f: VEM implemenetation
****************************************************************** */
void MeshMaps::VEM_FaceVelocity_(int c, int f, std::vector<Polynomial>& v) const
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

  v[0].monomials(0).coefs()[0] -= b * xf0;
  v[1].monomials(0).coefs()[0] -= (b^xf0)[0];
}


/* ******************************************************************
* Calculate Jacobian at point x of face f 
****************************************************************** */
Tensor MeshMaps::FaceJacobian(int c, int f, const std::vector<Polynomial>& v,
                              const AmanziGeometry::Point& x) const
{
  if (method_ == WHETSTONE_METHOD_FEM) {
    return FEM_FaceJacobian_(c, f, v, x);
  } else if (method_ == WHETSTONE_METHOD_VEM) {
    return VEM_FaceJacobian_(c, f, v, x);
  }
}


/* ******************************************************************
* Calculate Jacobian at point x of face f: VEM implementation
****************************************************************** */
Tensor MeshMaps::VEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                                   const AmanziGeometry::Point& x) const
{
  // FIXME x is not used
  Tensor jac(d_, 2);

  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      jac(i, j) = v[i].monomials(1).coefs()[j];
    }
  }
  jac += 1.0;

  return jac;
}


/* ******************************************************************
* Calculate Jacobian at point x of face f: FEM implementation
****************************************************************** */
Tensor MeshMaps::FEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                                   const AmanziGeometry::Point& x) const
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

  return FEM_Jacobian(c, xref);
}


/* ******************************************************************
* Polynomial approximation of map x2 = F(x1).
* We assume that vectors of vertices have a proper length.
****************************************************************** */
int MeshMaps::LeastSquareFit(int order,
                             const std::vector<AmanziGeometry::Point>& x1, 
                             const std::vector<AmanziGeometry::Point>& x2,
                             std::vector<AmanziGeometry::Point>& u) const
{
  int nk = (order + 1) * (order + 2) / 2;
  int nx = x1.size();

  // evaluate basis functions at given points
  int i1(0);
  DenseMatrix psi(nk, nx);

  for (int k = 0; k <= order; ++k) {
    for (int i = 0; i < k + 1; ++i) {
      for (int n = 0; n < nx; ++n) {
        psi(i1, n) = std::pow(x1[n][0], k - i) * std::pow(x1[n][1], i);
      }
      i1++;
    }
  }
      
  // form linear system
  DenseMatrix A(nk, nk);
  DenseVector bx(nk), by(nk), ux(nk), uy(nk);

  for (int i = 0; i < nk; ++i) {
    for (int j = i; j < nk; ++j) {
      double tmp(0.0);
      for (int n = 0; n < nx; ++n) {
        tmp += psi(i, n) * psi(j, n);
      }
      A(i, j) = A(j, i) = tmp;
    }

    bx(i) = 0.0;
    by(i) = 0.0;
    for (int n = 0; n < nx; ++n) {
      bx(i) += x2[n][0] * psi(i, n);
      by(i) += x2[n][1] * psi(i, n);
    }
  }

  // solver linear systems
  A.Inverse();
  A.Multiply(bx, ux, false);
  A.Multiply(by, uy, false);

  u.clear();
  for (int i = 0; i < nk; ++i) {
    u.push_back(AmanziGeometry::Point(ux(i), uy(i)));
  }
}


/* ******************************************************************
* Support of finite element meshes: bilinear map (2D algorithm)
****************************************************************** */
AmanziGeometry::Point MeshMaps::FEM_Map(int c, const AmanziGeometry::Point& xref) const
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


/* ******************************************************************
* Support of finite element meshes: Jacobian
****************************************************************** */
Tensor MeshMaps::FEM_Jacobian(int c, const AmanziGeometry::Point& xref) const
{
  Tensor jac = FEM_JacobianInternal_(mesh1_, c, xref);
 
  // Jacobian for convolution of two maps.
  if (mesh0_ != Teuchos::null) {
    Tensor jac0 = FEM_JacobianInternal_(mesh0_, c, xref);
    jac0.Inverse();
    jac = jac * jac0;
  }

  return jac;
}


Tensor MeshMaps::FEM_JacobianInternal_(
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

}  // namespace WhetStone
}  // namespace Amanzi

