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
#include "MFD3DFactory.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity in cell c: new algorithm.
* NOTE: second mesh is not used, so does it belong here?
****************************************************************** */
void MeshMaps_VEM::VelocityCell(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  VectorPolynomial moments(d_, vf.size());

  WhetStone::MFD3DFactory factory;
  auto mfd = factory.CreateMFD3D(mesh0_, method_, order_);

  if (projector_ == "H1") {
    mfd->H1Cell(c, vf, moments, vc);
  }
  else if (projector_ == "L2") {
    mfd->L2Cell(c, vf, moments, vc);
  }
  else if (projector_ == "least square") {
    LeastSquareProjector_Cell_(order_, c, vf, vc);
  } 
  else {
    AMANZI_ASSERT(false);
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

    MFD3D_CrouzeixRaviart mfd(mesh0_);
    mfd.set_order(order_);

    AmanziGeometry::Point p0(mesh1_->face_centroid(f) - mesh0_->face_centroid(f));
    mfd.H1Face(f, p0, ve, vf);
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
    int f, double t, const VectorPolynomial& vf, VectorPolynomial& cn) const
{
  AmanziGeometry::Point p(d_);
  WhetStone::Tensor J(d_, 2);

  JacobianFaceValue_(f, vf, p, J);
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
* Calculate mesh velocity in cell c: new algorithm
****************************************************************** */
void MeshMaps_VEM::JacobianCell(
    int c, const std::vector<VectorPolynomial>& vf, MatrixPolynomial& J) const
{
  auto moments = std::make_shared<DenseVector>();
  MFD3D_CrouzeixRaviart mfd(mesh0_);

  mfd.set_order(order_ + 1);
  mfd.L2GradientCell(c, vf, moments, J);
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
  AMANZI_ASSERT(order == 1 || order == 2);

  vc.resize(d_);
  for (int i = 0; i < d_; ++i) vc[i].Reshape(d_, order);
  
  AmanziGeometry::Point px1, px2;
  std::vector<AmanziGeometry::Point> x1, x2;

  Entity_ID_List nodes, faces;
  mesh0_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  for (int n = 0; n < nnodes; ++n) {
    mesh0_->node_get_coordinates(nodes[n], &px1);
    x1.push_back(px1);

    mesh1_->node_get_coordinates(nodes[n], &px2);
    x2.push_back(px2 - px1);
  }

  // FIXME
  if (order > 1) {
    mesh0_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      const auto& xf = mesh0_->face_centroid(faces[n]);
      x1.push_back(xf);

      for (int i = 0; i < d_; ++i)  {
        px2[i] = vf[n][i].Value(xf); 
      }
      x2.push_back(px2);
    }
  }

  // calculate velocity u(X) = F(X) - X
  LeastSquareFit(order, x1, x2, vc);
}


/* ******************************************************************
* Calculate mesh velocity in cell c: old algorithm
****************************************************************** */
void MeshMaps_VEM::ParseInputParameters_(const Teuchos::ParameterList& plist)
{
  method_ = plist.get<std::string>("method");
  order_ = plist.get<int>("method order");
  projector_ = plist.get<std::string>("projector");
}

}  // namespace WhetStone
}  // namespace Amanzi

