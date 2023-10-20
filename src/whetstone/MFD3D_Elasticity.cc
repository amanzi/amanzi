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

  The mimetic finite difference method for elasticity.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_Elasticity.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for mass matrix in mechanics.
* Only the upper triangular part of Ac is calculated.
* Requires mesh_get_edges to complete the implementation.
****************************************************************** */
int
MFD3D_Elasticity::L2consistency(int c,
                                const Tensor& T,
                                DenseMatrix& N,
                                DenseMatrix& Mc,
                                bool symmetry)
{
  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    v1 = fm - cm;
    double a = normal * v1;
    for (int k = 0; k < d_; k++) v1[k] = a - normal[k] * v1[k];
  }

  return 0;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in mechanics.
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int
MFD3D_Elasticity::H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  double volume = mesh_->getCellVolume(c);

  int nrows = d_ * nnodes;
  N.Reshape(nrows, d_ * (d_ + 1));
  Ac.Reshape(nrows, nrows);

  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  // convolution of tensors
  std::vector<Tensor> vE, vTE;

  for (int k = 0; k < d_; k++) {
    Tensor E(d_, 2);
    E(k, k) = 1.0;
    vE.push_back(E);
    vTE.push_back(T * E);
  }

  for (int k = 0; k < d_; k++) {
    for (int l = k + 1; l < d_; l++) {
      Tensor E(d_, 2);
      E(k, l) = E(l, k) = 1.0;
      vE.push_back(E);
      vTE.push_back(T * E);
    }
  }

  // calculate exact integration matrix
  int modes = d_ * (d_ + 1) / 2;
  DenseMatrix coefM(modes, modes);

  for (int i = 0; i < modes; ++i) {
    for (int j = i; j < modes; ++j) {
      coefM(i, j) = DotTensor(vE[i], vTE[j]) * volume;
      coefM(j, i) = coefM(i, j);
    }
  }

  // to calculate matrix R, we use temporary matrix N
  N.PutScalar(0.0);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    double area = mesh_->getFaceArea(f);

    auto face_nodes = mesh_->getFaceNodes(f);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs[i];
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        p = mesh_->getNodeCoordinate(v);
        pnext = mesh_->getNodeCoordinate(vnext);
        pprev = mesh_->getNodeCoordinate(vprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      for (int k = 0; k < modes; k++) {
        v1 = vTE[k] * normal;
        for (int l = 0; l < d_; l++) N(d_ * pos + l, k) += v1[l] * u;
      }
    }
  }

  // calculate R coefM^{-1} R^T
  DenseVector a1(modes), a2(modes), a3(modes);
  coefM.Inverse();

  for (int i = 0; i < nrows; i++) {
    for (int k = 0; k < modes; ++k) a1(k) = N(i, k);
    coefM.Multiply(a1, a3, false);

    for (int j = i; j < nrows; j++) {
      for (int k = 0; k < modes; ++k) a2(k) = N(j, k);
      Ac(i, j) = a2 * a3;
    }
  }

  // calculate matrix N
  N.PutScalar(0.0);
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh_->getNodeCoordinate(v);
    v1 = p - cm;

    int md = 0;
    for (int k = 0; k < d_; k++) {
      N(d_ * i + k, md) = v1[k];
      md++;
    }
    for (int k = 0; k < d_; k++) {
      for (int l = k + 1; l < d_; l++) {
        N(d_ * i + k, md) = v1[l];
        N(d_ * i + l, md) = v1[k];
        md++;
      }
    }
    for (int k = 0; k < d_; k++) { // additional columns correspod to kernel
      N(d_ * i + k, md) = 1.0;
      md++;
    }
    for (int k = 0; k < d_; k++) {
      for (int l = k + 1; l < d_; l++) {
        N(d_ * i + k, md) = v1[l];
        N(d_ * i + l, md) = -v1[k];
        md++;
      }
    }
  }
  return 0;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrixOptimized(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  StabilityOptimized_(T, N, A);
  return 0;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
* For education purpose only: ther are no M-matrices in elasticity.
****************************************************************** */
int
MFD3D_Elasticity::StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE;
  ok = StabilityMMatrix_(c, N, A, objective);
  return ok;
}

} // namespace WhetStone
} // namespace Amanzi
