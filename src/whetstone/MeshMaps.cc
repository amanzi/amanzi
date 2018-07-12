/*
  WhetStone, version 2.1
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

#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity on face f.
* NOTE: 2D algorithm for linear velocity.
****************************************************************** */
void MeshMaps::VelocityFace(int f, VectorPolynomial& v) const
{
  AmanziMesh::Entity_ID_List nodes;
  AmanziGeometry::Point x0, x1;

  const AmanziGeometry::Point& xf0 = mesh0_->face_centroid(f);
  const AmanziGeometry::Point& xf1 = mesh1_->face_centroid(f);

  mesh0_->face_get_nodes(f, &nodes);
  mesh0_->node_get_coordinates(nodes[0], &x0);
  mesh1_->node_get_coordinates(nodes[0], &x1);

  x0 -= xf0;
  x1 -= xf1;

  v.resize(d_);
  for (int i = 0; i < d_; ++i) {
    v[i].Reshape(d_, 1);
  }

  x0 /= L22(x0);
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      v[i](1, j) = x1[i] * x0[j];
    }
    v[i](0, 0) = xf1[i] - x1[i] * (x0 * xf0);
    v[i](1, i) -= 1.0;
  }
}


/* ******************************************************************
* Calculation of Jacobian.
* Multiple velocities are packed in a rectagular matrix.
****************************************************************** */
void MeshMaps::Jacobian(const VectorPolynomial& vc, MatrixPolynomial& J) const
{
  // allocate memory
  int nvc = vc.size();
  J.resize(nvc);
  for (int i = 0; i < nvc; ++i) {
    J[i].resize(d_);
  }

  // copy velocity gradients to Jacobian
  VectorPolynomial tmp(d_, 0);
  for (int i = 0; i < nvc; ++i) {
    tmp.Gradient(vc[i]);
    for (int j = 0; j < d_; ++j) {
      J[i][j] = tmp[j];
    }
  }
}


/* ******************************************************************
* Calculation of matrix of cofactors.
* Multiple cofactors are packed in a rectagular matrix.
****************************************************************** */
void MeshMaps::Cofactors(
    double t, const MatrixPolynomial& J, MatrixPolynomial& C) const
{
  // allocate memory for matrix of cofactors
  int nJ = J.size();
  C.resize(nJ);
  for (int i = 0; i < nJ; ++i) {
    C[i].resize(d_);
  }

  // calculate cofactors
  int kJ = nJ / d_;
  for (int n = 0; n < kJ; ++n) {
    int m0 = n * d_;
    int m1 = m0 + 1;
    if (d_ == 2) {
      C[m1][1] = J[m0][0];
      C[m1][0] = J[m0][1];
      C[m1][0] *= -1.0;

      C[m0][0] = J[m1][1];
      C[m0][1] = J[m1][0];
      C[m0][1] *= -1.0;
    }
    else if (d_ == 3) {
      int m2 = m0 + 2;
      C[m0][0] = J[m1][1] * J[m2][2] - J[m2][1] * J[m1][2];
      C[m1][0] = J[m2][1] * J[m0][2] - J[m0][1] * J[m2][2];
      C[m2][0] = J[m0][1] * J[m1][2] - J[m1][1] * J[m0][2];

      C[m0][1] = J[m2][0] * J[m1][2] - J[m1][0] * J[m2][2];
      C[m1][1] = J[m0][0] * J[m2][2] - J[m2][0] * J[m0][2];
      C[m2][1] = J[m1][0] * J[m0][2] - J[m0][0] * J[m1][2];

      C[m0][2] = J[m1][0] * J[m2][1] - J[m2][0] * J[m1][1];
      C[m1][2] = J[m2][0] * J[m0][1] - J[m0][0] * J[m2][1];
      C[m2][2] = J[m0][0] * J[m1][1] - J[m1][0] * J[m0][1];
    }
  }

  // add time dependence
  for (int i = 0; i < nJ; ++i) {
    for (int j = 0; j < d_; ++j) {
      C[i][j] *= t;
    }
    C[i][i % d_](0, 0) += 1.0;
  }
}


/* ******************************************************************
* Calculate detminant at time t.
* Multiple determinatds are packed in a vector.
****************************************************************** */
void MeshMaps::Determinant(
   double t, const MatrixPolynomial& J, VectorPolynomial& det) const
{
  int ndet = J.size() / d_;
  det.resize(ndet);

  MatrixPolynomial Jt;
  Jt.resize(d_);
  for (int i = 0; i < d_; ++i) {
    Jt[i].resize(d_);
  }

  for (int n = 0; n < ndet; ++n) {
    int m = n * d_;
    for (int i = 0; i < d_; ++i) {
      for (int j = 0; j < d_; ++j) {
        Jt[i][j] = J[m + i][j] * t;
      }
      Jt[i][i](0, 0) += 1.0;
    }

    if (d_ == 2) {
      det[n] = Jt[0][0] * Jt[1][1] - Jt[0][1] * Jt[1][0];
    }
    else if (d_ == 3) {
      det[n] = Jt[0][0] * Jt[1][1] * Jt[2][2] 
             + Jt[2][0] * Jt[0][1] * Jt[1][2] 
             + Jt[1][0] * Jt[2][1] * Jt[0][2] 
             - Jt[2][0] * Jt[1][1] * Jt[0][2] 
             - Jt[1][0] * Jt[0][1] * Jt[2][2] 
             - Jt[0][0] * Jt[2][1] * Jt[1][2]; 
    }
  }
}


/* ******************************************************************
* Polynomial approximation v of map x2 = F(x1).
* We assume that vectors of vertices have a proper length.
****************************************************************** */
int MeshMaps::LeastSquareFit(int order,
                             const std::vector<AmanziGeometry::Point>& x1, 
                             const std::vector<AmanziGeometry::Point>& x2,
                             VectorPolynomial& v) const
{
  Polynomial poly(d_, order);

  int nk = poly.size();
  int nx = x1.size();

  // evaluate basis functions at given points
  DenseMatrix psi(nk, nx);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int i = it.PolynomialPosition();
    const int* idx = it.multi_index();

    for (int n = 0; n < nx; ++n) {
      double val(1.0);
      for (int k = 0; k < d_; ++k) {
        val *= std::pow(x1[n][k], idx[k]);
      }
      psi(i, n) = val;
    }
  }
      
  // form matrix of linear system
  DenseMatrix A(nk, nk);

  for (int i = 0; i < nk; ++i) {
    for (int j = i; j < nk; ++j) {
      double tmp(0.0);
      for (int n = 0; n < nx; ++n) {
        tmp += psi(i, n) * psi(j, n);
      }
      A(i, j) = A(j, i) = tmp;
    }
  }

  A.Inverse();

  // solver linear systems
  DenseVector b(nk), u(nk);

  v.resize(d_);
  for (int k = 0; k < d_; ++k) { 
    v[k].Reshape(d_, order);
    v[k].set_origin(AmanziGeometry::Point(d_));

    for (int i = 0; i < nk; ++i) {
      b(i) = 0.0;
      for (int n = 0; n < nx; ++n) {
        b(i) += x2[n][k] * psi(i, n);
      }
    }

    A.Multiply(b, u, false);

    for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
      int n = it.MonomialSetOrder();
      int m = it.MonomialSetPosition();
      int i = it.PolynomialPosition();
      v[k](n, m) = u(i);
    }
  }

  return 0;
}


/* ******************************************************************
* Project polynomial on mesh0 to polynomial space on mesh1.
****************************************************************** */
void MeshMaps::ProjectPolynomial(int c, Polynomial& poly) const
{
  int order = poly.order();

  WhetStone::Entity_ID_List faces, nodes;
  std::vector<int> dirs;

  mesh0_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();  

  AmanziGeometry::Point v0(d_), v1(d_);
  std::vector<VectorPolynomial> vvf;

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh1_->face_centroid(f);
    mesh0_->face_get_nodes(f, &nodes);

    mesh0_->node_get_coordinates(nodes[0], &v0);
    mesh0_->node_get_coordinates(nodes[1], &v1);
    double f0 = poly.Value(v0);
    double f1 = poly.Value(v1);

    WhetStone::VectorPolynomial vf(d_ - 1, 1);
    vf[0].Reshape(d_ - 1, order);
    if (order == 1) {
      vf[0](0, 0) = (f0 + f1) / 2; 
      vf[0](1, 0) = f1 - f0; 
    } else if (order == 2) {
      double f2 = poly.Value((v0 + v1) / 2);
      vf[0](0, 0) = f2;
      vf[0](1, 0) = f1 - f0;
      vf[0](2, 0) = -4 * f2 + 2 * f0 + 2 * f1;
    } else {
      AMANZI_ASSERT(0);
    }

    mesh1_->node_get_coordinates(nodes[0], &v0);
    mesh1_->node_get_coordinates(nodes[1], &v1);

    std::vector<AmanziGeometry::Point> tau;
    tau.push_back(v1 - v0);
    vf[0].InverseChangeCoordinates(xf, tau);
    vvf.push_back(vf);
  }

  VectorPolynomial moments(d_, 1);

  if (order == 2) {
    NumericalIntegration numi(mesh1_);
    double mass = numi.IntegratePolynomialCell(c, poly);

    moments[0](0, 0) = mass / mesh1_->cell_volume(c);
  }

  VectorPolynomial vc(d_, 0);
  MFD3D_LagrangeSerendipity mfd(mesh1_);
  mfd.set_order(order);
  mfd.L2Cell(c, vvf, moments, vc);
  poly = vc[0];
}

}  // namespace WhetStone
}  // namespace Amanzi

