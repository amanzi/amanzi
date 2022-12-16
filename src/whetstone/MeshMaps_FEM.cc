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
void
MeshMaps_FEM::VelocityCell(int c,
                           const std::vector<VectorPolynomial>& ve,
                           const std::vector<VectorPolynomial>& vf,
                           VectorPolynomial& vc) const
{
  Entity_ID_List nodes;

  mesh1_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  AMANZI_ASSERT(nnodes == 4);

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

  // calculate map in coordinate system centered at point 1
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
  for (int i = 0; i < d_; ++i) { vc[i](1, i) -= 1.0; }
}


/* ******************************************************************
* Supporting routine for calculating Jacobian.
****************************************************************** */
Tensor
MeshMaps_FEM::JacobianValueInternal_(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                     int c,
                                     const AmanziGeometry::Point& xref) const
{
  Entity_ID_List nodes;

  mesh->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  AMANZI_ASSERT(nnodes == 4);

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

} // namespace WhetStone
} // namespace Amanzi
