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

  The mimetic finite difference method for elasticity with weak symmetry.
  Stress: 1 normal component per face (d dofs)
  Displacement: 1 value per cell (d dofs)
  Rotations:    1 value per node (d(d-1)/2 dofs)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_ElasticityWeakSym_CurvedFace.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for the mass matrix is block diagonal.
****************************************************************** */
int
MFD3D_ElasticityWeakSym_CurvedFace::L2consistency(
  int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nrows = d_ * nfaces;
  N.Reshape(nrows, d_ * d_);
  Mc.Reshape(nrows, nrows);

  const auto& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  Tensor Tinv(T);
  Tinv.Inverse();

  // convolution of tensors
  std::vector<Tensor> vE, vTE;
  std::vector<int> basis_row, basis_index;

  for (int k = 0; k < d_; k++) {
    for (int l = 0; l < d_; l++) {
      Tensor E(d_, 2);
      E(k, l) = 1.0;
      vE.push_back(E);
      vTE.push_back(Tinv * E);

      basis_row.push_back(k);
      basis_index.push_back(l);
    }
  }

  // calculate exact integration matrix
  int modes = d_ * d_;
  DenseMatrix coefM(modes, modes);

  for (int i = 0; i < modes; ++i) {
    for (int j = i; j < modes; ++j) {
      coefM(i, j) = DotTensor(vE[i], vTE[j]) * volume;
      coefM(j, i) = coefM(i, j);
    }
  }
  coefM.Inverse();

  // compute matrix R, we reuse memory for N
  N.PutScalar(0.0);
  int row0(0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->getFaceArea(f);
    const AmanziGeometry::Point xf = (*bf_)[f];

    for (int m = 0; m < modes; ++m) {
      int k0 = basis_row[m];
      int l0 = basis_index[m];
      N(row0 + k0, m) = (xf[l0] - xc[l0]) * area;
    }
    row0 += d_;
  }

  // compute Mc
  DenseVector a1(modes), a2(modes), a3(modes);
  for (int i = 0; i < nrows; ++i) {
    for (int m = 0; m < modes; ++m) a1(m) = N(i, m);
    coefM.Multiply(a1, a3, false);

    for (int j = 0; j < nrows; ++j) {
      for (int m = 0; m < modes; ++m) a2(m) = N(j, m);
      Mc(i, j) = a3 * a2;
    }
  }

  // compute matrix N
  AmanziGeometry::Point v2(d_);
  N.PutScalar(0.0);

  int row(0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);
    for (int m = 0; m < modes; ++m) {
      v2 = vTE[m] * normal;
      for (int k = 0; k < d_; ++k) {
        N(row + k, m) = v2[k] * dirs[i] / area;
      }
    }
    row += d_;
  }

  return 0;
}


/* ******************************************************************
* Rotations operator
****************************************************************** */
void
MFD3D_ElasticityWeakSym_CurvedFace::RotationMatrix(int c, DenseMatrix& G)
{
  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

  int nd = d_ * (d_ - 1) / 2;
  G.Reshape(d_ * nfaces, nd);
  G.PutScalar(0.0);

  int row0(0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->getFaceArea(f);
    AmanziGeometry::Point xf = (*bf_)[f];

    for (int i0 = 0; i0 < nd; ++i0) {
      int i1 = (i0 + 1) % d_;
      G(row0 + i0, i0) = area * (xf[i1] - xc[i1]);
      G(row0 + i1, i0) = area * (xc[i0] - xf[i0]);
    }
    row0 += d_;
  }
}

} // namespace WhetStone
} // namespace Amanzi
