/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for maps between mesh objects located on different
  meshes, e.g. two states of a deformable mesh.
*/

#include "Point.hh"

#include "CoordinateSystems.hh"
#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Calculate mesh velocity on 2D face f.
 ****************************************************************** */
void
MeshMaps::VelocityFace(int f, VectorPolynomial& v) const
{
  AMANZI_ASSERT(d_ == 2);
  AmanziGeometry::Point_List points0, points1;
  mesh0_->face_get_ho_nodes(f, &points0);
  mesh1_->face_get_ho_nodes(f, &points1);
  AMANZI_ASSERT(points0.size() == points1.size());

  // local coordinate system
  const AmanziGeometry::Point& normal = mesh0_->face_normal(f);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  FaceCoordinateSystem(normal, tau);

  // polynomial is converted from local to global coordinate system
  Kokkos::View<AmanziMesh::Entity_ID*> nodes;
  AmanziGeometry::Point x0, x1, y0, y1;

  const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
  AmanziGeometry::Point yf = mesh1_->face_centroid(f);

  mesh0_->face_get_nodes(f, nodes);
  mesh0_->node_get_coordinates(nodes(0), &x0);
  mesh0_->node_get_coordinates(nodes(1), &x1);

  mesh1_->node_get_coordinates(nodes(0), &y0);
  mesh1_->node_get_coordinates(nodes(1), &y1);

  // velocity at points defining the polynomial
  y0 -= x0;
  y1 -= x1;
  yf -= xf;

  for (int i = 0; i < points1.size(); ++i) { points1[i] -= points0[i]; }

  // velocity is transformed from local to glocal coordinate systems
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
    v[i].ChangeOrigin(AmanziGeometry::Point(d_));
  }
}


/* ******************************************************************
 * Transformation of normal is defined completely by face data.
 ****************************************************************** */
void
MeshMaps::NansonFormula(int f, const VectorPolynomial& map,
                        VectorPolynomial& cn) const
{
  AMANZI_ASSERT(d_ == 2);

  const AmanziGeometry::Point& normal = mesh0_->face_normal(f);

  cn.resize(d_);
  auto grad = Gradient(map[0]);
  cn[1] = grad[0] * normal[1] - grad[1] * normal[0];
  cn[1](0) += normal[1];

  grad = Gradient(map[1]);
  cn[0] = grad[1] * normal[0] - grad[0] * normal[1];
  cn[0](0) += normal[0];
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
 * Calculation of matrix of cofactors.
 * Multiple cofactors are packed in a rectagular matrix.
 ****************************************************************** */
void
MeshMaps::Cofactors(const MatrixPolynomial& J, MatrixPolynomial& C) const
{
  // allocate memory for matrix of cofactors
  int nJ = J.NumRows();
  C.Reshape(d_, nJ, d_, 0, false);

  // calculate cofactors
  int kJ = nJ / d_;
  for (int n = 0; n < kJ; ++n) {
    int m0 = n * d_;
    int m1 = m0 + 1;
    if (d_ == 2) {
      C(m1, 1) = J(m0, 0);
      C(m1, 0) = J(m0, 1);
      C(m1, 0) *= -1.0;

      C(m0, 0) = J(m1, 1);
      C(m0, 1) = J(m1, 0);
      C(m0, 1) *= -1.0;
    } else if (d_ == 3) {
      int m2 = m0 + 2;
      C(m0, 0) = J(m1, 1) * J(m2, 2) - J(m2, 1) * J(m1, 2);
      C(m1, 0) = J(m2, 1) * J(m0, 2) - J(m0, 1) * J(m2, 2);
      C(m2, 0) = J(m0, 1) * J(m1, 2) - J(m1, 1) * J(m0, 2);

      C(m0, 1) = J(m2, 0) * J(m1, 2) - J(m1, 0) * J(m2, 2);
      C(m1, 1) = J(m0, 0) * J(m2, 2) - J(m2, 0) * J(m0, 2);
      C(m2, 1) = J(m1, 0) * J(m0, 2) - J(m0, 0) * J(m1, 2);

      C(m0, 2) = J(m1, 0) * J(m2, 1) - J(m2, 0) * J(m1, 1);
      C(m1, 2) = J(m2, 0) * J(m0, 1) - J(m0, 0) * J(m2, 1);
      C(m2, 2) = J(m0, 0) * J(m1, 1) - J(m1, 0) * J(m0, 1);
    }
  }
}


/* ******************************************************************
 * Calculate detminant at time t.
 * Multiple determinatds are packed in a vector.
 ****************************************************************** */
void
MeshMaps::Determinant(const MatrixPolynomial& J, VectorPolynomial& det) const
{
  int ndet = J.NumRows() / d_;
  det.resize(ndet);

  for (int n = 0; n < ndet; ++n) {
    int m0 = n * d_;
    int m1 = m0 + 1;

    if (d_ == 2) {
      det[n] = J(m0, 0) * J(m1, 1) - J(m0, 1) * J(m1, 0);
    } else if (d_ == 3) {
      int m2 = m0 + 2;
      det[n] = J(m0, 0) * J(m1, 1) * J(m2, 2) + J(m2, 0) * J(m0, 1) * J(m1, 2) +
               J(m1, 0) * J(m2, 1) * J(m0, 2) - J(m2, 0) * J(m1, 1) * J(m0, 2) -
               J(m1, 0) * J(m0, 1) * J(m2, 2) - J(m0, 0) * J(m2, 1) * J(m1, 2);
    }
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

  Kokkos::View<WhetStone::Entity_ID*> faces, nodes;

  mesh0_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  AmanziGeometry::Point v0(d_), v1(d_);
  std::vector<Polynomial> vvf;

  for (int i = 0; i < nfaces; ++i) {
    int f = faces(i);
    const AmanziGeometry::Point& xf = mesh1_->face_centroid(f);
    mesh0_->face_get_nodes(f, nodes);

    mesh0_->node_get_coordinates(nodes(0), &v0);
    mesh0_->node_get_coordinates(nodes(1), &v1);
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

    mesh1_->node_get_coordinates(nodes(0), &v0);
    mesh1_->node_get_coordinates(nodes(1), &v1);

    std::vector<AmanziGeometry::Point> tau;
    tau.push_back(v1 - v0);
    vf.InverseChangeCoordinates(xf, tau);
    vvf.push_back(vf);
  }

  Polynomial moments(d_, 0);

  if (order == 2) {
    NumericalIntegration numi(mesh1_);
    double mass = numi.IntegratePolynomialCell(c, poly);

    moments(0) = mass / mesh1_->cell_volume(c, false);
  }

  Teuchos::ParameterList plist;
  plist.set<int>("method order", order);

  MFD3D_LagrangeSerendipity mfd(plist, mesh1_);
  mfd.L2Cell(c, vvf, &moments, poly);
}

} // namespace WhetStone
} // namespace Amanzi
