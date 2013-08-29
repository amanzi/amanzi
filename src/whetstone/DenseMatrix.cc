/*
 This is the discretization component of the Amanzi code. 

 Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
 Amanzi is released under the three-clause BSD License. 
 The terms of use and "as is" disclaimer for this license are 
 provided in the top-level COPYRIGHT file.

Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)

Usage: 
*/

#include <vector>

#include "lapack.hh"
#include "DenseMatrix.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor.
****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol) 
{ 
  m_ = mrow;
  n_ = ncol;
  data_ = new double[m_ * n_]; 
  access_ = WHETSTONE_DATA_ACCESS_COPY;
}


/* ******************************************************************
* No memory check is performed: invalid read is possible.
****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol, double* data, int data_access)
{
  m_ = mrow;
  n_ = ncol;
  access_ = data_access;

  if (access_ == WHETSTONE_DATA_ACCESS_COPY) {
    data_ = new double[m_ * n_]; 
    for (int i = 0; i < mrow * ncol; i++) data_[i] = data[i];
  } else {
    data_ = data;
  }
}


/* ******************************************************************
* Copy constructor creates an new matrix.
****************************************************************** */
DenseMatrix::DenseMatrix(const DenseMatrix& B)
{
  m_ = B.NumRows();
  n_ = B.NumCols();
  access_ = WHETSTONE_DATA_ACCESS_COPY;

  data_ = new double[m_ * n_];
  const double* dataB = B.Values();
  for (int i = 0; i < m_ * n_; i++) data_[i] = dataB[i];
}


/* ******************************************************************
* Second level routine: inversion
****************************************************************** */
int DenseMatrix::Inverse() 
{
  if (n_ != m_) return 911;

  int ierr;
  int iwork[n_];
  DGETRF_F77(&n_, &n_, data_, &n_, iwork, &ierr);
  if (ierr) return ierr;
 
  int lwork = n_ * n_;
  double dwork[lwork];
  DGETRI_F77(&n_, data_, &n_, iwork, dwork, &lwork, &ierr);
  return ierr;
}

}  // namespace WhetStone
}  // namespace Amanzi
