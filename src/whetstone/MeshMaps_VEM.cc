/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two states 
  of a deformable mesh: virtual element implementation.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps_VEM.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculation of Jacobian.
****************************************************************** */
void MeshMaps_VEM::JacobianCellValue(
    int c, double t, const AmanziGeometry::Point& xref, Tensor& J) const
{
}


/* ******************************************************************
* Calculate determinant of a Jacobian.
****************************************************************** */
void MeshMaps_VEM::JacobianDet(int c, double t, Polynomial& v) const
{
  v.Reshape(d_, 0);
  v.monomials(0).coefs()[0] = mesh1_->cell_volume(c);
}


/* ******************************************************************
* Calculate mesh velocity on face f.
****************************************************************** */
void MeshMaps_VEM::VelocityFace(int c, int f, std::vector<Polynomial>& v) const
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
* Calculate Jacobian at point x of face f 
****************************************************************** */
void MeshMaps_VEM::JacobianFaceValue(
    int c, int f, const std::vector<Polynomial>& v,
    const AmanziGeometry::Point& x, Tensor& J) const
{
  // FIXME x is not used
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      J(i, j) = v[i].monomials(1).coefs()[j];
    }
  }
  J += 1.0;
}

}  // namespace WhetStone
}  // namespace Amanzi

