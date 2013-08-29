/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: ara-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructors
****************************************************************** */
MFD3D::MFD3D(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
  mesh_(mesh) {
  stability_method_ = WHETSTONE_STABILITY_GENERIC;
}


/* ******************************************************************
* Calculate stability factor using matrix and optional scaling.
****************************************************************** */
double MFD3D::CalculateStabilityScalar(DenseMatrix& Mc)
{
  int nrows = Mc.NumRows();

  scalar_stability_ = 0.0;
  for (int i = 0; i < nrows; i++) scalar_stability_ += Mc(i, i);
  scalar_stability_ /= nrows;

  if (stability_method_ == WHETSTONE_STABILITY_GENERIC_SCALED) {
    scalar_stability_ *= scaling_factor_;
  }

  return scalar_stability_;
}


/* ******************************************************************
* Set up positive scaling factor for a scalar stability term.
* Warning: Ignores silently negative factors.
****************************************************************** */
void MFD3D::ModifyStabilityScalingFactor(double factor)
{
  if (factor > 0.0) {
    stability_method_ = WHETSTONE_STABILITY_GENERIC_SCALED;
    scaling_factor_ = factor;
  }
}


/* ******************************************************************
* Simplest stability term is added to the consistency term. 
****************************************************************** */
void MFD3D::StabilityScalar(int cell, DenseMatrix& N,
                            DenseMatrix& Mc, DenseMatrix& M)
{
  GrammSchmidt(N);
  CalculateStabilityScalar(Mc);

  int nrows = Mc.NumRows();
  int ncols = N.NumCols();

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(i, j) = Mc(i, j);
  }

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
* WARNING: the routine is used for inverse of mass matrix only.
****************************************************************** */
int MFD3D::StabilityOptimized(const Tensor& T, DenseMatrix& N,
                              DenseMatrix& Mc, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = N.NumRows();
  int ncols = N.NumCols();

  // find correct scaling of a stability term
  double lower, upper, eigmin = Mc(0, 0);
  // T.spectral_bounds(&lower, &upper);
  for (int k = 1; k < nrows; k++) eigmin = std::min(eigmin, Mc(k, k));

  // find null space of N^T
  DenseMatrix U(nrows, nrows);
  int info, ldv = 1, size = 5 * d + 3 * nrows;
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
    for (int j = i+1; j < nrows; j++) F(m++) = -Mc(i, j);
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

  int status = WHETSTONE_ELEMENTAL_MATRIX_OK;
  for (int loop = 0; loop < 3; loop++) {
    if (loop == 1) {   
      for (int i = 0; i < mcols; i++) G(i) = std::max(G(i), 0.0);
      status = WHETSTONE_ELEMENTAL_MATRIX_PASSED;
    } else if (loop == 2) {
      for (int i = mcols; i < nparam; i++) G(i) = std::max(G(i), 0.0);
      status = WHETSTONE_ELEMENTAL_MATRIX_PASSED;
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
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(i, j) = Mc(i, j);
  }

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

  return status;
}


/* ******************************************************************
* Conventional Gramm-Schimdt orthogonalization of colums of matrix N. 
****************************************************************** */
void MFD3D::GrammSchmidt(DenseMatrix& N)
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


/* ******************************************************************
* Extension of Mesh API. 
****************************************************************** */
int MFD3D::cell_get_face_adj_cell(const int cell, const int face)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(face, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  if (ncells == 2) {
    int c2 = cells[0];
    if (cell == c2) c2 = cells[1];
    return c2;
  }
  return -1;
}


/* ******************************************************************
* Returns position of the number v in the list of nodes.  
****************************************************************** */
int MFD3D::FindPosition_(int v, AmanziMesh::Entity_ID_List nodes)
{
  for (int i = 0; i < nodes.size(); i++) {
    if (nodes[i] == v) return i;
  }
  return -1;
}

}  // namespace WhetStone
}  // namespace Amanzi



