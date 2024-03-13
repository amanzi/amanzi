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

  The mimetic finite difference method.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "MFD3D_Electromagnetics.hh"
#include "MFD3D_Diffusion.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
* This is the experimental algorithm.
****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  int ok;
  if (d_ == 2) {
    ok = H1consistency2D_(c, T, N, Ac);
  } else {
    ok = H1consistency3D_(c, T, N, Ac);
  }
  return ok;
}


/* ******************************************************************
* Consistency condition for a stiffness matrix.
****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  const auto& [faces, fdirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nd = 3;
  N.Reshape(nfaces, nd);
  Ac.Reshape(nfaces, nfaces);

  // calculate Ac = R (R^T N)^{+} R^T
  double T00 = T(0, 0);
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    double a1 = mesh_->getFaceArea(f);

    for (int j = i; j < nfaces; ++j) {
      f = faces[j];
      double a2 = mesh_->getFaceArea(f);
      Ac(i, j) = a1 * a2 / (T00 * fdirs[i] * fdirs[j] * volume);
    }
  }

  // Matrix N(:, 1:2) are simply tangents
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    double len = mesh_->getFaceArea(f);

    for (int k = 0; k < d_; ++k) {
      len = -len;
      N(i, k) = normal[1 - k] / len;
    }

    N(i, d_) = (xf - xc) * normal / len;
  }

  return 0;
}


/* ******************************************************************
* Consistency condition for a stiffness matrix.
****************************************************************** */
int
MFD3D_Electromagnetics::H1consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  std::vector<int> map;

  const auto& [faces, fdirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  int nd = 6; // order_ * (order_ + 2) * (order_ + 3) / 2;
  N.Reshape(nedges, nd);
  Ac.Reshape(nedges, nedges);

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

    auto [fedges, edirs] = mesh_->getFaceEdgesAndDirections(f);
    int nfedges = fedges.size();

    map = mesh_->getFaceCellEdgeMap(f, c);

    for (int m = 0; m < nfedges; ++m) {
      int e = fedges[m];
      const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);
      double len = mesh_->getEdgeLength(e);

      len *= 2 * fdirs[i] * edirs[m];
      v2 = xe - xf;

      for (int k = 0; k < d_; ++k) N(map[m], k) += len * v2[k];
    }
  }

  // calculate Ac = R (R^T N)^{+} R^T
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Ac(i, j) = (v2 * v3) / (4 * volume);
    }
  }

  // Matrix N(:, 1:3) are simply tangents
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);
    AmanziGeometry::Point tau = mesh_->getEdgeVector(e);
    double len = mesh_->getEdgeLength(e);

    tau /= len;
    v1 = xe - xc;
    v3 = tau ^ v1;

    for (int k = 0; k < d_; ++k) N(i, k) = tau[k];
    for (int k = 0; k < d_; ++k) N(i, d_ + k) = v3[k];
  }

  return 0;
}


/* ******************************************************************
* Mass matrix optimized for sparsity.
****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M);
  if (ok) return ok;

  ok = StabilityOptimized_(T, N, M);
  return ok;
}


/* ******************************************************************
* Inverse of matrix matrix for edge-based discretization optimized
* for starcity.
****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixInverseOptimized(int c, const Tensor& T, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, T, R, W);
  if (ok) return ok;

  ok = StabilityOptimized_(T, R, W);
  return ok;
}


/* ******************************************************************
* A simple mass matrix for testing.
****************************************************************** */
int
MFD3D_Electromagnetics::MassMatrixDiagonal(int c, const Tensor& T, DenseMatrix& M)
{
  double volume = mesh_->getCellVolume(c);

  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  M.PutScalar(0.0);
  for (int n = 0; n < nedges; n++) { M(n, n) = d_ * volume / (nedges * T(0, 0)); }
  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix M, C;

  StiffnessMatrix(c, T, A, M, C);
  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. Curls in 2D and 3D are
* defined using exterior face normals.
****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrix(int c,
                                        const Tensor& T,
                                        DenseMatrix& A,
                                        DenseMatrix& M,
                                        DenseMatrix& C)
{
  MFD3D_Diffusion diffusion(mesh_);
  diffusion.MassMatrixScaledArea(c, T, M);

  // populate curl matrix
  CurlMatrix(c, C);
  int nfaces = C.NumRows();
  int nedges = C.NumCols();

  A.Reshape(nedges, nedges);
  DenseMatrix MC(nfaces, nedges);

  MC.Multiply(M, C, false);
  A.Multiply(C, MC, true);

  return 0;
}


/* ******************************************************************
* Stiffness matrix: the algorithm based on new consistency condition
****************************************************************** */
int
MFD3D_Electromagnetics::StiffnessMatrix_GradCorrection(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  WhetStone::AddGradient(mesh_, c, N);

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Curl matrix acts onto the space of total fluxes; hence, there is
* no face area scaling below.
****************************************************************** */
void
MFD3D_Electromagnetics::CurlMatrix(int c, DenseMatrix& C)
{
  std::vector<int> map;

  const auto& edges = mesh_->getCellEdges(c);
  int nedges = edges.size();

  const auto& [faces, fdirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  C.Reshape(nfaces, nedges);
  C.PutScalar(0.0);

  if (d_ == 2) {
    auto nodes = mesh_->getCellNodes(c);

    for (int i = 0; i < nfaces; ++i) {
      int j = (i + 1) % nfaces;
      C(i, i) = -1.0;
      C(i, j) = 1.0;
    }
  } else {
    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];

      map = mesh_->getFaceCellEdgeMap(f, c);
      auto [fedges, edirs] = mesh_->getFaceEdgesAndDirections(f);
      int nfedges = fedges.size();

      for (int j = 0; j < nfedges; ++j) {
        int e = fedges[j];
        double len = mesh_->getEdgeLength(e);
        C(i, map[j]) = len * edirs[j] * fdirs[i];
      }
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
