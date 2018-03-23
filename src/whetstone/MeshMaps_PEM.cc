/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two states 
  of a deformable mesh: polygonal finite element implementation.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps_PEM.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate piecewise polynomial mesh velocity in cell c. 
****************************************************************** */
void MeshMaps_PEM::VelocityCell(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  ASSERT(d_ == 2);

  Entity_ID_List faces, nodes;
  AmanziGeometry::Point p(d_), q(d_);
  Tensor T0(d_, d_), T1(d_, d_);

  AmanziGeometry::Point x0 = cell_geometric_center(0, c);
  AmanziGeometry::Point x1 = cell_geometric_center(1, c);
  // const AmanziGeometry::Point x0 = mesh0_->cell_centroid(c);
  // const AmanziGeometry::Point x1 = mesh1_->cell_centroid(c);

  mesh0_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  int m(0);
  vc.resize(nfaces * d_);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    int m = n * d_;

    // calculate 2x2 map F(X) = q + J * X
    mesh0_->face_get_nodes(f, &nodes);
    for (int i = 0; i < 2; ++i) {
      mesh0_->node_get_coordinates(nodes[i], &p);
      mesh1_->node_get_coordinates(nodes[i], &q);
    
      T0.SetColumn(i, p - x0);
      T1.SetColumn(i, q - x1);
    }

    T0.Inverse();
    Tensor J = T1 * T0;
    q = x1 - J * x0;

    // calculate velocity u(X) = F(X) - X
    J += -1.0;
    for (int i = 0; i < d_; ++i) {
      vc[m].Reshape(d_, 1);

      vc[m](0, 0) = q[i];
      for (int j = 0; j < d_; ++j) {
        vc[m](1, j) = J(i, j);
      }
      m++;
    }
  }
}


/* ******************************************************************
* Transformation of normal is defined completely by face data.
****************************************************************** */
void MeshMaps_PEM::NansonFormula(
    int f, double t, const VectorPolynomial& vf, VectorPolynomial& cn) const
{
  AmanziGeometry::Point p(d_);
  WhetStone::Tensor J(d_, 2);

  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      J(i, j) = t * vf[i](1, j);
    }
  }
  J += 1.0;

  Tensor C = J.Cofactors();
  p = C * mesh0_->face_normal(f);

  cn.resize(d_);
  for (int i = 0; i < d_; ++i) {
    cn[i].Reshape(d_, 0);
    cn[i](0, 0) = p[i];
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

