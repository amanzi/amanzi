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

  Derham complex: mimetic inner products on edges.
*/

#include "Mesh.hh"

#include "DeRham_Edge.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
* Non-symmetric tensor is not yet used.
****************************************************************** */
int
DeRham_Edge::L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  int ok;
  if (d_ == 2) {
    ok = L2consistency2D_(c, T, N, Mc);
  } else {
    ok = L2consistency3D_(c, T, N, Mc);
  }

  return ok;
}


/* ******************************************************************
* 2D consistency condition for the mass matrix. Only the upper
* triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
****************************************************************** */
int
DeRham_Edge::L2consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  Tensor Tinv(T);
  Tinv.Inverse();
  // multiply by two 90 degree rotation matrices
  if (Tinv.rank() == 2) {
    double tmp = Tinv(0, 1);
    Tinv(0, 1) = -Tinv(1, 0);
    Tinv(1, 0) = -tmp;

    tmp = Tinv(0, 0);
    Tinv(0, 0) = Tinv(1, 1);
    Tinv(1, 1) = tmp;
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    double a1 = mesh_->getFaceArea(f);
    v2 = Tinv * (xf - xc);

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& yf = mesh_->getFaceCentroid(f);
      double a2 = mesh_->getFaceArea(f);

      v1 = yf - xc;
      Mc(i, j) = (v1 * v2) * (a1 * a2) / volume;
    }
  }

  // Rows of matrix N are oriented tangent vectors. Since N goes to
  // the Gramm-Schmidt orthogonalizetion procedure, we can skip scaling
  // tensorial factor T.
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double a = mesh_->getFaceArea(f);

    for (int k = 0; k < d_; k++) {
      a = -a;
      N(i, k) = normal[1 - k] * dirs[i] / a;
    }
  }

  return 0;
}


/* ******************************************************************
* 3D consistency condition for the mass matrix. Only the upper
* triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
****************************************************************** */
int
DeRham_Edge::L2consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc)
{
  const auto& [faces, fdirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();
  
  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  N.Reshape(nedges, d_);
  Mc.Reshape(nedges, nedges);

  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);
  AmanziGeometry::Point vv[3];

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);

    v1 = xc - xf;

    double a1 = normal * v1;
    for (int k = 0; k < d_; ++k) {
      vv[k].set(normal[k] * v1);
      vv[k][k] -= a1;
    }

    auto [fedges, edirs] = mesh_->getFaceEdgesAndDirections(f);
    int nfedges = fedges.size();

    auto map = mesh_->getFaceCellEdgeMap(f, c);

    for (int m = 0; m < nfedges; ++m) {
      int e = fedges[m];
      const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);

      v3 = xe - xf;

      double len = mesh_->getEdgeLength(e);
      len /= 2.0 * area * area * fdirs[i] * edirs[m];

      for (int k = 0; k < d_; ++k) { N(map[m], k) -= len * ((vv[k] ^ v3) * normal); }
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T
  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = Tinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

  // Rows of matrix N are simply tangents. Since N goes to the
  // Gramm-Schmidt orthogonalizetion procedure, we can skip scaling
  // tensorial factor T.
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->getEdgeVector(e);
    double len = mesh_->getEdgeLength(e);
    for (int k = 0; k < d_; ++k) N(i, k) = tau[k] / len;
  }

  return 0;
}


/* ******************************************************************
* Mass matrix for edge-based discretization.
****************************************************************** */
int
DeRham_Edge::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
****************************************************************** */
int
DeRham_Edge::L2consistencyInverse(int c,
                                  const Tensor& T,
                                  DenseMatrix& R,
                                  DenseMatrix& Wc,
                                  bool symmetry)
{
  int ok;
  if (d_ == 2) {
    ok = L2consistencyInverse2D_(c, T, R, Wc);
  } else {
    ok = L2consistencyInverse3D_(c, T, R, Wc);
  }

  return ok;
}


/* ******************************************************************
* 2D consistency condition for inverse of the mass matrix. Only upper
* triangular part of matrix Wc = N (N^T R)^{-1} N^T is calculated.
****************************************************************** */
int
DeRham_Edge::L2consistencyInverse2D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_);
  Wc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_);
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  // Since N is scaled by T, N = N0 * T, we use tensor T in the
  // inverse L2 consistency term.
  Tensor Trot(T);
  if (Trot.rank() == 2) {
    Trot(0, 1) = -Trot(0, 1);
    Trot(1, 0) = -Trot(1, 0);

    double tmp = Trot(0, 0);
    Trot(0, 0) = Trot(1, 1);
    Trot(1, 1) = tmp;
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double a1 = mesh_->getFaceArea(f);

    v1 = Trot * normal;

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->getFaceNormal(f);
      double a2 = mesh_->getFaceArea(f);
      Wc(i, j) = (v1 * v2) / (a1 * a2 * dirs[i] * dirs[j] * volume);
    }
  }

  // matrix R
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    double len = mesh_->getFaceArea(f);

    v1 = xf - xc;
    for (int k = 0; k < d_; k++) {
      len = -len;
      R(i, k) = v1[1 - k] * len;
    }
  }

  return 0;
}


/* ******************************************************************
* 3D consistency condition for inverse of the mass matrix. Only upper
* triangular part of matrix Wc = N (N^T R)^{-1} N^T is calculated.
* Recall that N^T R = |c| T.
****************************************************************** */
int
DeRham_Edge::L2consistencyInverse3D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc)
{
  const auto& [faces,fdirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  R.Reshape(nedges, d_);
  Wc.Reshape(nedges, nedges);

  AmanziGeometry::Point v3(d_), v4(d_), tau(d_), p1(d_), p2(d_);
  AmanziGeometry::Point vv[3];

  // Since the matrix N is the matrix of scaled tangent vextors,
  // N = N0 T, we use the tensor T.
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nedges; i++) {
    int e1 = edges[i];
    const AmanziGeometry::Point& v1 = mesh_->getEdgeVector(e1);
    double a1 = mesh_->getEdgeLength(e1);
    v3 = T * v1;

    for (int j = i; j < nedges; j++) {
      int e2 = edges[j];
      const AmanziGeometry::Point& v2 = mesh_->getEdgeVector(e2);
      double a2 = mesh_->getEdgeLength(e2);
      Wc(i, j) = (v2 * v3) / (a1 * a2 * volume);
    }
  }

  // Calculate matrix R
  R.PutScalar(0.0);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);

    v4 = xc - xf;

    double a1 = normal * v4;
    for (int k = 0; k < d_; ++k) {
      vv[k].set(normal[k] * v4);
      vv[k][k] -= a1;
    }

    auto [fedges, edirs] = mesh_->getFaceEdgesAndDirections(f);
    int nfedges = fedges.size();

    auto map = mesh_->getFaceCellEdgeMap(f, c);

    for (int m = 0; m < nfedges; ++m) {
      int e = fedges[m];
      const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);

      v3 = xe - xf;

      double len = mesh_->getEdgeLength(e);
      len /= 2.0 * area * area * fdirs[i] * edirs[m];

      for (int k = 0; k < d_; ++k) { R(map[m], k) -= len * ((vv[k] ^ v3) * normal); }
    }
  }

  return 0;
}


/* ******************************************************************
* Inverse of the mass matrix for edge-based discretization.
****************************************************************** */
int
DeRham_Edge::MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, T, R, W, true);
  if (ok) return ok;

  StabilityScalar_(R, W);
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
