/*
This is the mimetic discretization component of the Amanzi code. 

Tensors of rank 1 are numbers in all dimensions.
Tensors of rank 2 are square matrices in all dimensions.
Only symmetric tensors of rank 4 are are considered here.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cmath>

#include "Teuchos_LAPACK.hpp"

#include "Point.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {


/* ******************************************************************
* Constructors
****************************************************************** */
Tensor::Tensor(const Tensor& T)
{
  int d = T.dimension();
  int rank = T.rank();
  double* data = T.data();

  if (d && rank) {
    data_ = NULL;
    int mem = init(d, rank);
    for (int i = 0; i < mem; i++) data_[i] = data[i];
  } else {
    d_ = rank_ = size_ = 0;
    data_ = NULL;
  }
}


/* ******************************************************************
* Initialization of a tensor of rank 1, 2 or 4. 
****************************************************************** */
int Tensor::init(const int d, const int rank)
{
  size_ = WHETSTONE_TENSOR_SIZE[d-1][rank-1];
  int mem = size_ * size_;

  if (data_) delete[] data_;
  data_ = new double[mem];

  d_ = d;
  rank_ = rank;
  for (int i = 0; i < mem; i++) data_[i] = 0.0;

  return mem;
}


/* ******************************************************************
* Trace operation with tensors of rank 1 and 2
****************************************************************** */
double Tensor::trace() const
{
  double s = 0.0;
  if (rank_ <= 2) {
    for (int i = 0; i < size_; i++) s += (*this)(i, i);
  }
  return s;
}


/* ******************************************************************
* Inverse operation with tensors of rank 1 and 2
****************************************************************** */
void Tensor::inverse()
{
  if (size_ == 1) {
    data_[0] = 1.0 / data_[0];

  } else if (size_ == 2) {  // We use inverse formula based on minors
    double det = data_[0] * data_[3] - data_[1] * data_[2];

    double a = data_[0];
    data_[0] = data_[3] / det;
    data_[3] = a / det;

    data_[1] /= -det;
    data_[2] /= -det;

  } else if (rank_ <= 2) {
    Teuchos::LAPACK<int, double> lapack;
    int info;
    int ipiv[size_];
    double work[size_];
    lapack.GETRF(size_, size_, data_, size_, ipiv, &info);
    lapack.GETRI(size_, data_, size_, ipiv, work, size_, &info);
  }
}


/* ******************************************************************
* Transpose operator for non-symmetric tensors.
****************************************************************** */
void Tensor::transpose()
{
  if (rank_ == 2 && d_ == 2) {
    double tmp = data_[1];
    data_[1] = data_[2];
    data_[2] = tmp;
  } else if (rank_ == 2 && d_ == 3) {
    double tmp = data_[1];
    data_[1] = data_[3];
    data_[3] = tmp;

    tmp = data_[2];
    data_[2] = data_[6];
    data_[6] = tmp;

    tmp = data_[5];
    data_[5] = data_[7];
    data_[7] = tmp;   
  }
}


/* ******************************************************************
* Determinant of second-order tensors.
****************************************************************** */
double Tensor::determinant()
{
  double det = 0.0;
  if (rank_ == 2 && d_ == 2) {
    det = data_[0] * data_[3] - data_[1] * data_[2];
  } else if (rank_ == 2 && d_ == 3) {
    det = data_[0] * data_[4] * data_[8] 
        + data_[2] * data_[3] * data_[7] 
        + data_[1] * data_[5] * data_[6] 
        - data_[2] * data_[4] * data_[6] 
        - data_[1] * data_[3] * data_[8] 
        - data_[0] * data_[5] * data_[7]; 
  }
  return det;
}


/* ******************************************************************
* Spectral bounds of symmetric tensors of rank 1 and 2
****************************************************************** */
void Tensor::spectral_bounds(double* lower, double* upper) const
{
  if (size_ == 1) {
    *lower = data_[0];
    *upper = data_[0];

  } else if (size_ == 2) {
    double a = data_[0] - data_[3];
    double c = data_[1];
    double D = sqrt(a * a + 4 * c * c);
    double trace = data_[0] + data_[3];

    *lower = (trace - D) / 2;
    *upper = (trace + D) / 2;
  } else if (rank_ <= 2) {
    Teuchos::LAPACK<int, double> lapack;
    int ipiv[size_], lwork(3*size_), info;
    double S[size_], work[lwork];
    
    Tensor T(*this);
    lapack.SYEV('N', 'U', size_, T.data(), size_, S, work, lwork, &info);
    *lower = S[0];
    *upper = S[size_ - 1];
  }
}


/* ******************************************************************
* Multiplication by a constant.
****************************************************************** */
Tensor& Tensor::operator*=(const double& c)
{
  if (rank_ <= 2) {
    for (int i = 0; i < size_*size_; i++) data_[i] *= c;
  }
  return *this;
}


/* ******************************************************************
* Copy operator.
****************************************************************** */
Tensor& Tensor::operator=(const Tensor& T)
{
  int d = T.dimension();
  int rank = T.rank();
  double* data = T.data();

  int mem = init(d, rank);
  for (int i = 0; i < mem; i++) data_[i] = data[i];
  return *this;
}


/* ******************************************************************
* First convolution operation for tensors of rank 1 and 2. 
****************************************************************** */
AmanziGeometry::Point operator*(const Tensor& T, const AmanziGeometry::Point& p)
{
  int rank = T.rank();
  int d = T.dimension();
  double* data = T.data();

  AmanziGeometry::Point p2(p.dim());
  if (rank == 1) {
    p2 = data[0] * p;
    return p2;

  } else if (rank == 2) {
    for (int i = 0; i < d; i++) {
      p2[i] = 0.0;
      for (int j = 0; j < d; j++) {
        p2[i] += (*data) * p[j];
        data++;
      }
    }
    return p2;

  } else if (rank == 4) {
    return p;  // undefined operation (lipnikov@lanl.gov)
  }
  return p;
}


/* ******************************************************************
* Second convolution operation for tensors of rank 1, 2, and 4
****************************************************************** */
Tensor operator*(const Tensor& T1, const Tensor& T2)
{
  int d = T1.dimension();  // the dimensions should be equals
  int rank1 = T1.rank(), rank2 = T2.rank();
  double *data1 = T1.data(), *data2 = T2.data();

  Tensor T3;

  if (d == 2 && rank1 == 4 && rank2 == 2) {
    double a0, b0, c0;
    a0 = T2(0, 0);
    b0 = T2(1, 1);
    c0 = T2(0, 1);

    T3.init(d, rank2);
    T3(0, 0) = T1(0, 0) * a0 + T1(0, 1) * b0 + T1(0, 2) * c0;
    T3(1, 1) = T1(1, 0) * a0 + T1(1, 1) * b0 + T1(1, 2) * c0;
    T3(1, 0) = T3(0, 1) = T1(2, 0) * a0 + T1(2, 1) * b0 + T1(2, 2) * c0;

  } else if (rank1 == 1) {
    int mem = T3.init(d, rank2);
    double *data3 = T3.data();
    for (int i = 0; i < mem; i++) data3[i] = data2[i] * data1[0];

  } else if (rank2 == 1) {
    int mem = T3.init(d, rank1);
    double *data3 = T3.data();
    for (int i = 0; i < mem; i++) data3[i] = data1[i] * data2[0];

  } else if (rank2 == 2) {
    T3.init(d, 2);
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        double& entry = T3(i, j);
        for (int k = 0; k < d; k++) entry += T1(i, k) * T2(k, j);
      }
    }
  }

  return T3;
}


/* ******************************************************************
* Miscaleneous routines: populate tensors of rank 2
****************************************************************** */
int Tensor::add_column(const int column, const AmanziGeometry::Point& p)
{
  if (rank_ == 2) {
    for (int i = 0; i < d_; i++) (*this)(i, column) = p[i];
    return 1;
  }
  return -1;
}


int Tensor::add_row(const int row, const AmanziGeometry::Point& p)
{
  if (rank_ == 2) {
    for (int i = 0; i < d_; i++) (*this)(row, i) = p[i];
    return 1;
  }
  return 0;
}


/* ******************************************************************
* Miscaleneous routines: print
****************************************************************** */
std::ostream& operator<<(std::ostream& os, const Tensor& T)
{
  int d = T.dimension();
  int rank = T.rank();
  int size = T.size();

  os << "Tensor dimension=" << d << "  rank=" << rank << std::endl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) os << T(i, j) << " ";
    os << std::endl;
  }
  return os;
}

}  // namespace WhetStone
}  // namespace Amanzi

