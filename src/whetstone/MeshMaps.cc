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

  Base class for maps between mesh objects located on different
  meshes, e.g. two states of a deformable mesh.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity on 2D or 3D edge e.
****************************************************************** */
void
MeshMaps::VelocityEdge(int e, VectorPolynomial& v) const
{
  AMANZI_ASSERT(d_ == 3);

  auto points0 = mesh0_->getEdgeHOCoordinates(e);
  auto points1 = mesh1_->getEdgeHOCoordinates(e);
  AMANZI_ASSERT(points0.size() == points1.size());

  // local coordinate system (-0.5, 0.5)
  std::vector<AmanziGeometry::Point> tau(1, mesh0_->getEdgeVector(e));

  // polynomial is converted from local to global coordinate system
  int n0, n1;
  AmanziGeometry::Point x0, x1, y0, y1, ye;

  const AmanziGeometry::Point& xe = mesh0_->getEdgeCentroid(e);
  ye = mesh1_->getEdgeCentroid(e);

  mesh0_->getEdgeNodes(e, &n0, &n1);
  x0 = mesh0_->getNodeCoordinate(n0);
  x1 = mesh0_->getNodeCoordinate(n1);

  y0 = mesh1_->getNodeCoordinate(n0);
  y1 = mesh1_->getNodeCoordinate(n1);

  // velocity at points defining the polynomial
  y0 -= x0;
  y1 -= x1;
  ye -= xe;

  for (int i = 0; i < points1.size(); ++i) { points1[i] -= points0[i]; }

  // velocity is transformed from local to global coordinate systems
  int order = points1.size() + 1;

  v.resize(d_);
  for (int i = 0; i < d_; ++i) {
    v[i].Reshape(d_ - 2, order);
    if (order == 1) {
      v[i](0) = ye[i];
      v[i](1) = y1[i] - y0[i];
    } else {
      v[i](0) = points1[0][i];
      v[i](1) = y1[i] - y0[i];
      v[i](2) = 2 * y0[i] + 2 * y1[i] - 4 * points1[0][i];
    }

    v[i].InverseChangeCoordinates(xe, tau);
  }
}


/* ******************************************************************
* Calculate mesh velocity on 2D face f.
****************************************************************** */
void
MeshMaps::VelocityFace(int f, VectorPolynomial& v) const
{
  AMANZI_ASSERT(d_ == 2);

  auto points0 = mesh0_->getFaceHOCoordinates(f);
  auto points1 = mesh1_->getFaceHOCoordinates(f);
  AMANZI_ASSERT(points0.size() == points1.size());

  // local coordinate system
  const AmanziGeometry::Point& xf = mesh0_->getFaceCentroid(f);
  const AmanziGeometry::Point& normal = mesh0_->getFaceNormal(f);

  SurfaceCoordinateSystem coordsys(xf, normal);
  const auto& tau = *coordsys.tau();

  // polynomial is converted from local to global coordinate system
  AmanziMesh::Entity_ID_List nodes;
  AmanziGeometry::Point x0, x1, y0, y1;

  AmanziGeometry::Point yf = mesh1_->getFaceCentroid(f);

  nodes = mesh0_->getFaceNodes(f);
  x0 = mesh0_->getNodeCoordinate(nodes[0]);
  x1 = mesh0_->getNodeCoordinate(nodes[1]);

  y0 = mesh1_->getNodeCoordinate(nodes[0]);
  y1 = mesh1_->getNodeCoordinate(nodes[1]);

  // velocity at points defining the polynomial
  y0 -= x0;
  y1 -= x1;
  yf -= xf;

  for (int i = 0; i < points1.size(); ++i) { points1[i] -= points0[i]; }

  // velocity is transformed from local to global coordinate systems
  int order = points1.size() + 1;

  v.resize(d_);
  for (int i = 0; i < d_; ++i) {
    v[i].Reshape(d_ - 1, order);
    if (order == 1) {
      v[i](0) = yf[i];
      v[i](1) = y1[i] - y0[i];
    } else {
      v[i](0) = points1[0][i];
      v[i](1) = y1[i] - y0[i];
      v[i](2) = 2 * y0[i] + 2 * y1[i] - 4 * points1[0][i];
    }

    v[i].InverseChangeCoordinates(xf, tau);
  }
}


/* ******************************************************************
* Transformation of normal is defined completely by face data.
****************************************************************** */
void
MeshMaps::NansonFormula(int f,
                        const VectorSpaceTimePolynomial& map,
                        VectorSpaceTimePolynomial& cn) const
{
  const auto& normal = mesh0_->getFaceNormal(f);
  cn.resize(d_);

  auto grad = Gradient(map);

  if (d_ == 2) {
    cn[1] = grad(0, 0) * normal[1] - grad(0, 1) * normal[0];
    cn[0] = grad(1, 1) * normal[0] - grad(1, 0) * normal[1];
  } else {
    for (int i = 0; i < d_; ++i) {
      int j = (i + 1) % d_;
      int k = (j + 1) % d_;
      cn[i] = (grad(j, j) * grad(k, k) - grad(j, k) * grad(k, j)) * normal[i] +
              (grad(j, k) * grad(k, i) - grad(j, i) * grad(k, k)) * normal[j] +
              (grad(j, i) * grad(k, j) - grad(j, j) * grad(k, i)) * normal[k];
    }
  }
}


/* ******************************************************************
* Calculation of Jacobian.
* Multiple velocities are packed in a rectagular matrix.
****************************************************************** */
void
MeshMaps::Jacobian(const VectorPolynomial& vc, MatrixPolynomial& J) const
{
  // allocate memory
  int nvc = vc.size();
  J.Reshape(d_, nvc, d_, 0, false);

  // copy velocity gradients to Jacobian
  for (int i = 0; i < nvc; ++i) {
    auto tmp = Gradient(vc[i]);
    for (int j = 0; j < d_; ++j) { J(i, j) = tmp[j]; }
  }
}


/* ******************************************************************
* Polynomial approximation v of map x2 = F(x1).
* We assume that vectors of vertices have a proper length.
****************************************************************** */
int
MeshMaps::LeastSquareFit(int order,
                         const std::vector<AmanziGeometry::Point>& x1,
                         const std::vector<AmanziGeometry::Point>& x2,
                         VectorPolynomial& v) const
{
  Polynomial poly(d_, order);

  int nk = poly.size();
  int nx = x1.size();

  // evaluate basis functions at given points
  DenseMatrix psi(nx, nk);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    int i = it.PolynomialPosition();
    const int* idx = it.multi_index();

    for (int n = 0; n < nx; ++n) {
      double val(1.0);
      for (int k = 0; k < d_; ++k) { val *= std::pow(x1[n][k], idx[k]); }
      psi(n, i) = val;
    }
  }

  // form linear system
  DenseMatrix A(nk, nk);

  A.Multiply(psi, psi, true);
  A.Inverse();

  // solver linear systems
  DenseVector b(nk), u(nk);

  v.resize(d_);
  for (int k = 0; k < d_; ++k) {
    v[k].Reshape(d_, order);
    v[k].set_origin(AmanziGeometry::Point(d_));

    for (int i = 0; i < nk; ++i) {
      b(i) = 0.0;
      for (int n = 0; n < nx; ++n) { b(i) += x2[n][k] * psi(n, i); }
    }

    A.Multiply(b, u, false);

    for (auto i = 0; i < nk; ++i) { v[k](i) = u(i); }
  }

  return 0;
}


/* ******************************************************************
* Project polynomial on mesh0 to polynomial space on mesh1.
****************************************************************** */
void
MeshMaps::ProjectPolynomial(int c, Polynomial& poly) const
{
  int order = poly.order();

  WhetStone::Entity_ID_List nodes;

  const auto& faces = mesh0_->getCellFaces(c);
  int nfaces = faces.size();

  AmanziGeometry::Point v0(d_), v1(d_);
  std::vector<Polynomial> vvf;

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh1_->getFaceCentroid(f);
    nodes = mesh0_->getFaceNodes(f);

    v0 = mesh0_->getNodeCoordinate(nodes[0]);
    v1 = mesh0_->getNodeCoordinate(nodes[1]);
    double f0 = poly.Value(v0);
    double f1 = poly.Value(v1);

    WhetStone::Polynomial vf(d_ - 1, 1);
    vf.Reshape(d_ - 1, order);
    if (order == 1) {
      vf(0, 0) = (f0 + f1) / 2;
      vf(1, 0) = f1 - f0;
    } else if (order == 2) {
      double f2 = poly.Value((v0 + v1) / 2);
      vf(0, 0) = f2;
      vf(1, 0) = f1 - f0;
      vf(2, 0) = -4 * f2 + 2 * f0 + 2 * f1;
    } else {
      AMANZI_ASSERT(0);
    }

    v0 = mesh1_->getNodeCoordinate(nodes[0]);
    v1 = mesh1_->getNodeCoordinate(nodes[1]);

    std::vector<AmanziGeometry::Point> tau;
    tau.push_back(v1 - v0);
    vf.InverseChangeCoordinates(xf, tau);
    vvf.push_back(vf);
  }

  Polynomial moments(d_, 0);

  if (order == 2) {
    NumericalIntegration numi(mesh1_);
    double mass = numi.IntegratePolynomialCell(c, poly);

    moments(0) = mass / mesh1_->getCellVolume(c);
  }

  Teuchos::ParameterList plist;
  plist.set<int>("method order", order);

  MFD3D_LagrangeSerendipity mfd(plist, mesh1_);
  mfd.L2Cell(c, vvf, vvf, &moments, poly);
}

} // namespace WhetStone
} // namespace Amanzi
