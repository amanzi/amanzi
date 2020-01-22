/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Dense matrices and operations with them.
*/

#include <vector>
#include <cmath>

#include "dbc.hh"

#include "lapack.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * No memory check is performed: invalid read is possible.
 ****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol, double* data)
{
  m_ = mrow;
  n_ = ncol;
  //access_ = data_access;

  //if (access_ == WHETSTONE_DATA_ACCESS_COPY) {
  mem_ = m_ * n_;
  Kokkos::resize(data_,mem_); 
  //data_ = new double[mem_];
  for (int i = 0; i < mem_; i++) data_[i] = data[i];
  //} else {
  //  mem_ = 0;
  //  data_ = data;
  //}
}
DenseMatrix::DenseMatrix(int mrow, int ncol,
  const Kokkos::View<double*>& data)
{
  m_ = mrow;
  n_ = ncol;
  //access_ = data_access;

  //if (access_ == WHETSTONE_DATA_ACCESS_COPY) {
  mem_ = m_ * n_;
  Kokkos::resize(data_,mem_); 
  Kokkos::deep_copy(data_,data); 
  //data_ = new double[mem_];
  //for (int i = 0; i < mem_; i++) data_[i] = data[i];
  //} else {
  //  mem_ = 0;
  //  data_ = data;
  //}
}

#if 0 

/* ******************************************************************
 * Copy constructor creates an new matrix.
 ****************************************************************** */
DenseMatrix::DenseMatrix(const DenseMatrix& B)
{
  m_ = B.NumRows();
  n_ = B.NumCols();
  //access_ = WHETSTONE_DATA_ACCESS_COPY;

  mem_ = m_ * n_;
  Kokkos::resize(data_,mem_);
  Kokkos::deep_copy(data_,B.Values());  
  //data_ = new double[mem_];
  //const double* dataB = B.Values();
  //for (int i = 0; i < mem_; i++) data_[i] = dataB[i];
}

#endif 
/* ******************************************************************
 * Copy constructor creates a new matrix from the submatrix of B in
 * rows m1 to m2-1 and columns n1 to n2-1.
 ****************************************************************** */
DenseMatrix::DenseMatrix(const DenseMatrix& B, int m1, int m2, int n1, int n2)
{
  m_ = m2 - m1;
  n_ = n2 - n1;
  //access_ = WHETSTONE_DATA_ACCESS_COPY;

  mem_ = m_ * n_;
  Kokkos::resize(data_,mem_); 
  //data_ = new double[mem_];

  int mB = B.NumRows();
  int nB = B.NumCols();
  //const double* dataB = B.Values();

  for (int j = n1; j < n2; ++j) {
    //const double* tmpB = dataB + j * mB + m1;
    //double* tmpA = data_ + (j - n1) * m_;
    for (int i = 0; i < m_; ++i) {
      data_((j-n1)*m_+i) = j*mB+m1+i;
      //*tmpA = *tmpB;
      //tmpA++;
      //tmpB++;
    }
  }
}
/* ******************************************************************
 * Smart memory management. Data destroyed in general.
 ****************************************************************** */
void
DenseMatrix::reshape(int mrow, int ncol)
{
  //AMANZI_ASSERT(access_ == WHETSTONE_DATA_ACCESS_COPY);

  m_ = mrow;
  n_ = ncol;

  //if (mem_ != m_ * n_) {
  //  if (data_ != NULL) { delete[] data_; }
  mem_ = m_ * n_;
  Kokkos::resize(data_,mem_); 
  //data_ = new double[mem_];
  //}
}

/* ******************************************************************
 * Second level routine: submatrix in rows [ib,ie) and columns [jb,je)
 ****************************************************************** */
DenseMatrix
DenseMatrix::SubMatrix(int ib, int ie, int jb, int je)
{
  int mrows(ie - ib), ncols(je - jb);
  DenseMatrix tmp(mrows, ncols);
  Kokkos::View<double*> dataB = tmp.Values(); 
  //double* dataB = tmp.Values();

  for (int j = jb; j < je; ++j) {
    //const double* dataA = data_ + j * m_ + ib;
    for (int i = ib; i < ie; ++i) {
      dataB(j*ie+i) = data_(j*m_+ib+i); 
      //*dataB = *dataA;
      //dataA++;
      //dataB++;
    }
  }

  return tmp;
}

/* ******************************************************************
 * Second level routine: transpose
 ****************************************************************** */
void
DenseMatrix::Transpose(const DenseMatrix& A)
{
  const Kokkos::View<double*> dataA = A.Values(); 
  //const double* dataA = A.Values();
  int mrowsA = A.NumRows(), ncolsA = A.NumCols();

  reshape(ncolsA, mrowsA);

  for (int j = 0; j < ncolsA; ++j) {
    for (int i = 0; i < mrowsA; ++i) {
      data_(i*ncolsA+j) = dataA(j*mrowsA+i); 
      //*(data_ + i * ncolsA + j) = *dataA;
      //dataA++;
    }
  }
}


int
DenseMatrix::Transpose()
{
  if (m_ != n_) return 1;

  for (int i = 0; i < m_; ++i) {
    for (int j = i + 1; j < n_; ++j) {
      double& aij = data_(i * n_ + j);
      double& aji = data_(j * n_ + i);

      double tmp = aij;
      aij = aji;
      aji = tmp;
    }
  }
  return 0;
}

/* ******************************************************************
 * Second level routine: inversion
 ****************************************************************** */
int
DenseMatrix::Inverse()
{
  if (n_ != m_) return 911;

  int ierr;
  int iwork[n_];
  DGETRF_F77(&n_, &n_, &data_(0), &n_, iwork, &ierr);
  if (ierr) return ierr;

  int lwork = n_ * n_;
  double dwork[lwork];
  DGETRI_F77(&n_, &data_(0), &n_, iwork, dwork, &lwork, &ierr);
  return ierr;
}


/* ******************************************************************
 * Second level routine: inversion of symmetric poitive definite
 ****************************************************************** */
int
DenseMatrix::InverseSPD()
{
  if (n_ != m_) return 911;

  int ierr;
  DPOTRF_F77("U", &n_, &data_(0), &n_, &ierr);
  if (ierr) return ierr;

  DPOTRI_F77("U", &n_, &data_(0), &n_, &ierr);
  for (int i = 0; i < n_; i++)
    for (int j = i + 1; j < n_; j++) data_[i * n_ + j] = data_[j * n_ + i];
  return ierr;
}


/* ******************************************************************
 * Second level routine: calculates matrix D such that (*this)^T D = 0.
 * The matrix (*this) must have a full rank and have more rows than
 * columns.
 ****************************************************************** */
int
DenseMatrix::NullSpace(DenseMatrix& D)
{
  // We can treat only one type of rectangular matrix.
  if (m_ <= n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // D must have proper size.
  int m = D.NumRows(), n = D.NumCols();
  if (m != m_ || n != m_ - n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // Allocate memory for Lapack routine.
  int ldv(1), lwork, info;
  lwork = std::max(m_ + 3 * n_, 5 * n_);
  double U[m_ * m_], V, S[n_], work[lwork];

  DGESVD_F77(
    "A", "N", &m_, &n_, &data_(0), &m_, S, U, &m, &V, &ldv, work, &lwork, &info);

  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  Kokkos::View<double*> data = D.Values();
  int offset = m_ * n_;
  for (int i = 0; i < m * n; i++) data[i] = U[offset + i];

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
