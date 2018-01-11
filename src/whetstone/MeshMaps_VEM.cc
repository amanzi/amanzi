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

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity in cell c: new algorithm
****************************************************************** */
void MeshMaps_VEM::VelocityCell(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  // LeastSquareProjector_Cell_(order_, c, vf, vc);
  if (order_ < 2) {
    projector.HarmonicCell_CR1(c, vf, vc);
  } else {
    projector.HarmonicCell_CRk(c, order_, vf, vc);
  }
}


/* ******************************************************************
* Calculate mesh velocity on face f.
****************************************************************** */
void MeshMaps_VEM::VelocityFace(int f, VectorPolynomial& vf) const
{
  if (d_ == 2) {
    MeshMaps::VelocityFace(f, vf);
  } else {
    AmanziMesh::Entity_ID_List edges;
    std::vector<int> dirs;

    mesh0_->face_get_edges_and_dirs(f, &edges, &dirs);
    int nedges = edges.size();

    VectorPolynomial v;
    std::vector<VectorPolynomial> ve;

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      VelocityEdge_(e, v);
      ve.push_back(v);
    }

    AmanziGeometry::Point p0(mesh1_->face_centroid(f) - mesh0_->face_centroid(f));
    projector.HarmonicFace_CR1(f, p0, ve, vf);
  }
}


/* ******************************************************************
* Calculate mesh velocity on 2D or 3D edge e.
****************************************************************** */
void MeshMaps_VEM::VelocityEdge_(int e, VectorPolynomial& ve) const
{
  const AmanziGeometry::Point& xe0 = mesh0_->edge_centroid(e);
  const AmanziGeometry::Point& xe1 = mesh1_->edge_centroid(e);

  // velocity order 0
  ve.resize(d_);
  for (int i = 0; i < d_; ++i) {
    ve[i].Reshape(d_, 1);
  }

  // velocity order 1
  int n0, n1;
  AmanziGeometry::Point x0, x1;

  mesh0_->edge_get_nodes(e, &n0, &n1);
  mesh0_->node_get_coordinates(n0, &x0);
  mesh1_->node_get_coordinates(n0, &x1);

  x0 -= xe0;
  x1 -= xe1;

  // operator F(\xi) = x_c + R (\xi - \xi_c) where R = x1 * x0^T
  x0 /= L22(x0);
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      ve[i](1, j) = x1[i] * x0[j];
    }
    ve[i](0, 0) = xe1[i] - x1[i] * (x0 * xe0);
    ve[i](1, i) -= 1.0;
  }
}


/* ******************************************************************
* Transformation of normal is defined completely by face data.
* NOTE: Limited to P1 elements. FIXME
****************************************************************** */
void MeshMaps_VEM::NansonFormula(
    int f, double t, const VectorPolynomial& v, VectorPolynomial& cn) const
{
  AmanziGeometry::Point p(d_);
  WhetStone::Tensor J(d_, 2);

  JacobianFaceValue_(f, v, p, J);
  J *= t;
  J += 1.0 - t;

  Tensor C = J.Cofactors();
  p = C * mesh0_->face_normal(f);

  cn.resize(d_);
  for (int i = 0; i < d_; ++i) {
    cn[i].Reshape(d_, 0);
    cn[i](0, 0) = p[i];
  }
}


/* ******************************************************************
* Calculate Jacobian at point x of face f 
* NOTE: limited to linear velocity FIXME
****************************************************************** */
void MeshMaps_VEM::JacobianFaceValue_(
    int f, const VectorPolynomial& v,
    const AmanziGeometry::Point& x, Tensor& J) const
{
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      J(i, j) = v[i](1, j);
    }
  }
  J += 1.0;
}


/* ******************************************************************
* Calculate mesh velocity in cell c: old algorithm
****************************************************************** */
void MeshMaps_VEM::LeastSquareProjector_Cell_(
    int order, int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  vc.resize(d_);
  for (int i = 0; i < d_; ++i) vc[i].Reshape(d_, 1);
  
  Entity_ID_List nodes;
  mesh0_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  AmanziGeometry::Point px;
  std::vector<AmanziGeometry::Point> x1, x2, u;
  for (int i = 0; i < nnodes; ++i) {
    mesh0_->node_get_coordinates(nodes[i], &px);
    x1.push_back(px);
  }
  for (int i = 0; i < nnodes; ++i) {
    mesh1_->node_get_coordinates(nodes[i], &px);
    x2.push_back(px);
  }

  // calculate velocity u(X) = F(X) - X
  LeastSquareFit(order, x1, x2, vc);

  for (int i = 0; i < d_; ++i) {
    vc[i](1, i) -= 1.0;
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

