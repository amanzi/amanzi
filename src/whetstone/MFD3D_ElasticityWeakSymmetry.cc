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
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_ElasticityWeakSymmetry.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for the mass matrix is block diagonal.
****************************************************************** */
int
MFD3D_ElasticityWeakSymmetry::L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nrows = d_ * d_ * nfaces;
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
  int index[3] = { 0, 0, 0 };
  NumericalIntegration numi(mesh_);
  std::vector<const PolynomialBase*> polys(2);

  N.PutScalar(0.0);
  int row0(0), row1(d_ * nfaces), offset(d_ * nfaces);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->getFaceArea(f);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    AmanziGeometry::Point normal = mesh_->getFaceNormal(f);

    normal *= dirs[n];
    auto coordsys = std::make_shared<AmanziGeometry::SurfaceCoordinateSystem>(xf, normal);

    for (int m = 0; m < modes; ++m) {
      // zero moments
      int k0 = basis_row[m];
      int l0 = basis_index[m];
      N(row0 + k0, m) = (xf[l0] - xc[l0]) * area;

      // first moments
      set_index_(d_, l0, index);
      Polynomial cmono(d_, index, 1.0);
      cmono.set_origin(xc);
      polys[0] = &cmono;

      for (int l = 0; l < d_ - 1; ++l) {
        set_index_(d_ - 1, l, index);
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, *coordsys->tau());
        polys[1] = &fmono;

        N(l * offset + row1 + k0, m) = numi.IntegratePolynomialsFace(f, polys);
      }
    }
    row0 += d_;
    row1 += d_;
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
      for (int k = 0; k < d_; ++k) { N(row + k, m) = v2[k] * dirs[i] / area; }
    }
    row += d_;
  }

  return 0;
}


/* ******************************************************************
* Mass matrix in space of fluxes.
****************************************************************** */
int
MFD3D_ElasticityWeakSymmetry::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int
MFD3D_ElasticityWeakSymmetry::DivergenceMatrix(int c, DenseMatrix& B)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  B.Reshape(d_ * d_ * nfaces, d_);
  B.PutScalar(0.0);

  int row(0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    double area = mesh_->getFaceArea(f);
    for (int k = 0; k < d_; ++k) B(row++, k) = area;
  }

  return 0;
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
void
MFD3D_ElasticityWeakSymmetry::RotationMatrix(int c, DenseMatrix& G)
{
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

  int nd = d_ * (d_ - 1) / 2;
  G.Reshape(d_ * d_ * nfaces, nd);
  G.PutScalar(0.0);

  int index[3] = { 0, 0, 0 };
  NumericalIntegration numi(mesh_);
  std::vector<const PolynomialBase*> polys(2);

  int row0(0), row1(d_ * nfaces), offset(d_ * nfaces);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->getFaceArea(f);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    AmanziGeometry::Point normal = mesh_->getFaceNormal(f);

    normal *= dirs[n];
    auto coordsys = std::make_shared<AmanziGeometry::SurfaceCoordinateSystem>(xf, normal);

    for (int i0 = 0; i0 < nd; ++i0) {
      int i1 = (i0 + 1) % d_;

      // zero moments
      G(row0 + i0, i0) = area * (xf[i1] - xc[i1]);
      G(row0 + i1, i0) = area * (xc[i0] - xf[i0]);

      // first moments
      for (int l = 0; l < d_ - 1; ++l) {
        set_index_(d_ - 1, l, index);
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, *coordsys->tau());
        polys[1] = &fmono;

        set_index_(d_, i0, index);
        Polynomial cmono1(d_, index, 1.0);
        cmono1.set_origin(xc);
        polys[0] = &cmono1;

        G(l * offset + row1 + i0, i0) = numi.IntegratePolynomialsFace(f, polys);

        set_index_(d_, i1, index);
        Polynomial cmono2(d_, index, -1.0);
        cmono2.set_origin(xc);
        polys[0] = &cmono2;

        G(l * offset + row1 + i1, i0) = numi.IntegratePolynomialsFace(f, polys);
      }
    }
    row0 += d_;
    row1 += d_;
  }
}


/* ******************************************************************
* Lame stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int
MFD3D_ElasticityWeakSymmetry::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix M, G, B;

  MassMatrix(c, T, M);
  RotationMatrix(c, G);
  DivergenceMatrix(c, B);

  int nm = M.NumCols();
  int ng = G.NumCols();
  int nb = B.NumCols();

  int nrows = nm + nb + ng;
  A.Reshape(nrows, nrows);

  // stress continuity contraint
  DenseVector D(nm);
  int nfaces = nm / (d_ * d_);

  int row(0);
  for (int n = 0; n < nfaces; ++n) {
    for (int k = 0; k < d_; ++k) {
      D(row) = -B(row, k);
      row++;
    }
  }

  int offset(d_ * nfaces);
  for (int n = 0; n < offset; ++n) {
    for (int k = 1; k < d_; ++k) D(k * offset + n) = D(n);
  }

  // elliminate stresses
  DenseMatrix Minv(M);
  Minv.Inverse();

  // -- equation for lambdas
  TripleMatrixProduct_(D, Minv, D, A, 0, 0);
  TripleMatrixProduct_(D, Minv, B, A, 0, nm);
  TripleMatrixProduct_(D, Minv, G, A, 0, nm + nb);

  // -- equation for cell-centered displacements
  TripleMatrixProduct_(B, Minv, B, A, nm, nm);
  TripleMatrixProduct_(B, Minv, G, A, nm, nm + nb);

  // -- equation for rotations
  TripleMatrixProduct_(G, Minv, G, A, nm + nb, nm + nb);

  // symmetrize matrix
  for (int i = 0; i < nrows; ++i) {
    for (int j = i + 1; j < nrows; ++j) A(j, i) = A(i, j);
  }

  return 0;
}


/* ******************************************************************
* Helper function: seting index space for monomial
****************************************************************** */
void
MFD3D_ElasticityWeakSymmetry::set_index_(int d, int k, int* index)
{
  for (int i = 0; i < d; ++i) index[i] = 0;
  index[k] = 1;
}


/* ******************************************************************
* Helper function: triple product of matrices, including diagonals
****************************************************************** */
void
MFD3D_ElasticityWeakSymmetry::TripleMatrixProduct_(const DenseVector& ML,
                                                   const DenseMatrix& Minv,
                                                   const DenseVector& MR,
                                                   DenseMatrix& A,
                                                   int i0,
                                                   int j0)
{
  int nm = Minv.NumCols();

  for (int i = 0; i < nm; ++i) {
    for (int j = 0; j < nm; ++j) { A(i + i0, j + j0) = Minv(i, j) * ML(i) * MR(j); }
  }
}


void
MFD3D_ElasticityWeakSymmetry::TripleMatrixProduct_(const DenseVector& ML,
                                                   const DenseMatrix& Minv,
                                                   const DenseMatrix& MR,
                                                   DenseMatrix& A,
                                                   int i0,
                                                   int j0)
{
  int nm = Minv.NumCols();
  int nr = MR.NumCols();

  for (int i = 0; i < nm; ++i) {
    for (int j = 0; j < nr; ++j) {
      double sum(0.0);
      for (int k = 0; k < nm; ++k) { sum += Minv(i, k) * ML(i) * MR(k, j); }
      A(i + i0, j + j0) = sum;
    }
  }
}


void
MFD3D_ElasticityWeakSymmetry::TripleMatrixProduct_(const DenseMatrix& ML,
                                                   const DenseMatrix& Minv,
                                                   const DenseMatrix& MR,
                                                   DenseMatrix& A,
                                                   int i0,
                                                   int j0)
{
  int nm = Minv.NumCols();
  int nl = ML.NumCols();
  int nr = MR.NumCols();

  for (int i = 0; i < nl; ++i) {
    for (int j = 0; j < nr; ++j) {
      double sum(0.0);
      for (int k = 0; k < nm; ++k) {
        for (int l = 0; l < nm; ++l) { sum += Minv(k, l) * ML(k, i) * MR(l, j); }
      }
      A(i + i0, j + j0) = sum;
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
