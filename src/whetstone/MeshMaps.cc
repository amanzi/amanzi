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
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity on face f
* NOTE: 2D algorithm for P1 velocity
****************************************************************** */
void MeshMaps::VelocityFace(int f, VectorPolynomial& v) const
{
  AmanziMesh::Entity_ID_List nodes;
  AmanziGeometry::Point x0, x1;

  const AmanziGeometry::Point& xf0 = mesh0_->face_centroid(f);
  const AmanziGeometry::Point& xf1 = mesh1_->face_centroid(f);

  // velocity order 0
  v.resize(d_);
  for (int i = 0; i < d_; ++i) {
    v[i].Reshape(d_, 1);
    v[i](0, 0) = xf1[i] - xf0[i];
  }

  // velocity order 1
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
  v[0](0, 0) -= b * xf0;
  v[1](0, 0) -= (b^xf0)[0];
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

    for (int i = 0; i < nk; ++i) {
      b(i) = 0.0;
      for (int n = 0; n < nx; ++n) {
        b(i) += x2[n][k] * psi(i, n);
      }
    }

    A.Multiply(b, u, false);

    for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
      int n = it.MonomialOrder();
      int m = it.MonomialPosition();
      int i = it.PolynomialPosition();
      v[k].monomials(n).coefs()[m] = u(i);
    }
  }

  return 0;
}


/* ******************************************************************
* Elliptic projector on linear polynomials at time 0.
****************************************************************** */
void MeshMaps::EllipticProjectorP1(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& u) const
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh0_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double vol = mesh0_->cell_volume(c);

  // create zero vector polynomial
  u.resize(d_);
  for (int i = 0; i < d_; ++i) { 
    u[i].Reshape(d_, 1, true);
  }

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh0_->face_normal(f);

    for (int i = 0; i < d_; ++i) {
      double tmp = vf[n][i].Value(xf) * dirs[n] / vol;

      for (int j = 0; j < d_; ++j) {
        u[i].monomials(1).coefs()[j] += tmp * normal[j];
      }
    }
  }

  // fix the constant value
  const AmanziGeometry::Point& xc0 = mesh0_->cell_centroid(c);
  const AmanziGeometry::Point& xc1 = mesh1_->cell_centroid(c);
  AmanziGeometry::Point zero(d_);

  for (int i = 0; i < d_; ++i) {
    u[i](0, 0) = xc1[i] - xc0[i];
    u[i].set_origin(xc0);
    u[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Error calculation requires geometric center.
****************************************************************** */
AmanziGeometry::Point MeshMaps::cell_geometric_center(int id, int c) const
{
  Entity_ID_List nodes;
  AmanziGeometry::Point v(d_), xg(d_);

  mesh0_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  for (int i = 0; i < nnodes; ++i) {
    if (id == 0) {
      mesh0_->node_get_coordinates(nodes[i], &v);
    } else {
      mesh1_->node_get_coordinates(nodes[i], &v);
    }
    xg += v;
  } 
  xg /= nnodes;
  
  return xg;
}

}  // namespace WhetStone
}  // namespace Amanzi

