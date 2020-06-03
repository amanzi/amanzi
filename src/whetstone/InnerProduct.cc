/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for mimetic inner products.
*/

#include "DenseMatrix.hh"
#include "InnerProduct.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Simplest stability term is added to the consistency term. 
****************************************************************** */
void InnerProduct::StabilityScalar_(DenseMatrix& N, DenseMatrix& M)
{
  GrammSchmidt_(N);
  CalculateStabilityScalar_(M);

  int nrows = M.NumRows();
  int ncols = N.NumCols();

  for (int i = 0; i < nrows; i++) {  // add projector ss * (I - N^T N) to matrix M
    M(i, i) += scalar_stability_;

    for (int j = i; j < nrows; j++) {
      double s = 0.0;
      for (int k = 0; k < ncols; k++)  s += N(i, k) * N(j, k);
      M(i, j) -= s * scalar_stability_;
    }
  }

  for (int i = 0; i < nrows; i++) {  // symmetrization
    for (int j = i+1; j < nrows; j++) M(j, i) = M(i, j);
  }
}


/* ******************************************************************
* A simple optimization procedure that returns a diagonal mass
* matrix for a 2D and 3D orthogonal cells and diagonal tensors. 
* The algorithm minimizes off-diagonal entries in the mass matrix.
****************************************************************** */
int InnerProduct::StabilityOptimized_(const Tensor& T, DenseMatrix& N, DenseMatrix& M)
{
  int nrows = N.NumRows();
  int ncols = N.NumCols();

  // find correct scaling of a stability term
  double eigmin = M(0, 0);
  // T.spectral_bounds(&lower, &upper);
  for (int k = 1; k < nrows; k++) eigmin = std::min(eigmin, M(k, k));

  // find null space of N^T
  DenseMatrix U(nrows, nrows);
  int info, ldv = 1, size = 5 * ncols + 3 * nrows;
  double V, S[nrows], work[size];

  DGESVD_F77("A", "N", &nrows, &ncols, N.Values(), &nrows,  // N = u s v
             S, U.Values(), &nrows, &V, &ldv, work, &size, &info);

  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // calculate vectors C and C0
  int mrows = nrows * (nrows - 1) / 2;
  int mcols = nrows - ncols;
  int nparam = (mcols + 1) * mcols / 2;
  DenseMatrix C(mrows, nparam);
  DenseVector F(mrows);

  int m, n = 0;
  for (int k = ncols; k < nrows; k++) {
    m = 0;  // calculate diagonal entries of M_kk = U_k * U_k^T
    for (int i = 0; i < nrows; i++) 
      for (int j = i+1; j < nrows; j++) C(m++, n) = U(i, k) * U(j, k);
    n++; 
  }

  for (int k = ncols; k < nrows; k++) {
    for (int l = k+1; l < nrows; l++) {
      m = 0;  // calculate off-diagonal entries of M_kk + M_ll - M_kl - M_lk 
      for (int i = 0; i < nrows; i++) { 
        for (int j = i+1; j < nrows; j++) {
          C(m, n) = C(m, k-ncols) + C(m, l-ncols) - U(i, k) * U(j, l) - U(i, l) * U(j, k);
          m++;
        }
      }
      n++;
    }
  }

  m = 0;
  for (int i = 0; i < nrows; i++) { 
    for (int j = i+1; j < nrows; j++) F(m++) = -M(i, j);
  }

  // Form a linear system for parameters
  DenseMatrix A(nparam, nparam);
  DenseVector G(nparam);

  A.Multiply(C, C, true);  // A = C^T C
  C.Multiply(F, G, true);

  // Find parameters
  int nrhs = 1;
  DPOSV_F77("U", &nparam, &nrhs, A.Values(), &nparam, G.Values(), &nparam, &info);
  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // project solution on the positive quadrant and convert to matrix
  DenseMatrix P(mcols, mcols);
  P.PutScalar(0.0);

  for (int loop = 0; loop < 3; loop++) {
    if (loop == 1) {   
      for (int i = 0; i < mcols; i++) G(i) = std::max(G(i), 0.0);
    } else if (loop == 2) {
      for (int i = mcols; i < nparam; i++) G(i) = std::max(G(i), 0.0);
    }

    for (int k = 0; k < mcols; k++) P(k, k) = G(k);

    n = mcols;
    for (int k = 0; k < mcols; k++) {
      for (int l = k+1; l < mcols; l++) {
        P(k, k) += G(n);
        P(l, l) += G(n);
        P(l, k) = P(k, l) = -G(n);
        n++;  
      }
    }

    // check SPD property (we use allocated memory)
    DenseMatrix Ptmp(P);
    DSYEV_F77("N", "U", &mcols, Ptmp.Values(), &mcols, S, work, &size, &info); 
    if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    if (S[0] > eigmin) {
      break;
    } else if (loop == 2) {
      for (int k = 0; k < mcols; k++) if (P(k, k) == 0.0) P(k, k) = eigmin;
    }
  }

  // add stability term U G U^T
  DenseMatrix UP(nrows, mcols);
  UP.PutScalar(0.0);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < mcols; j++) {
      double& entry = UP(i, j);
      for (int k = 0; k < mcols; k++) entry += U(i, k+ncols) * P(k, j);
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) {
      double& entry = M(i, j);
      for (int k = 0; k < mcols; k++) entry += UP(i, k) * U(j, k+ncols);
      M(j, i) = M(i, j); 
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}



/* ******************************************************************
* Simple stability term for nonsymmetric tensors.
****************************************************************** */
void InnerProduct::StabilityScalarNonSymmetric_(DenseMatrix& N, DenseMatrix& M)
{
  GrammSchmidt_(N);
  CalculateStabilityScalar_(M);

  int nrows = M.NumRows();
  int ncols = N.NumCols();

  // add projector ss * (I - N^T N) to matrix M
  for (int i = 0; i < nrows; i++) {  
    M(i, i) += scalar_stability_;

    for (int j = i; j < nrows; j++) {
      double s = 0.0;
      for (int k = 0; k < ncols; k++)  s += N(i, k) * N(j, k);

      s *= scalar_stability_;
      M(i, j) -= s;
      if (i - j) M(j, i) -= s;
    }
  }
}


/* ******************************************************************
* Calculate stability factor using matrix and optional scaling.
****************************************************************** */
double InnerProduct::CalculateStabilityScalar_(DenseMatrix& Mc)
{
  int nrows = Mc.NumRows();

  scalar_stability_ = 0.0;
  for (int i = 0; i < nrows; i++) scalar_stability_ += Mc(i, i);
  scalar_stability_ /= double(nrows) / 2.0;

  if (stability_method_ == WHETSTONE_STABILITY_GENERIC_SCALED) {
    scalar_stability_ *= scaling_factor_;
  }

  return scalar_stability_;
}


/* ******************************************************************
* Conventional Gramm-Schmidt orthogonalization of colums of matrix N. 
****************************************************************** */
void InnerProduct::GrammSchmidt_(DenseMatrix& N)
{
  int nrows = N.NumRows();
  int ncols = N.NumCols();

  int i, j, k;
  for (i = 0; i < ncols; i++) {
    double l22 = 0.0;
    for (k = 0; k < nrows; k++) l22 += N(k, i) * N(k, i);

    l22 = 1.0 / sqrt(l22);
    for (k = 0; k < nrows; k++) N(k, i) *= l22;

    for (j = i+1; j < ncols; j++) {
      double s = 0.0;
      for (k = 0; k < nrows; k++) s += N(k, i) * N(k, j);
      for (k = 0; k < nrows; k++) N(k, j) -= s * N(k, i);  // orthogonolize i and j
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi


