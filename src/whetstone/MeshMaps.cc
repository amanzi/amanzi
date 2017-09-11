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
* Calculate mesh velocity on face f.
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

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi

