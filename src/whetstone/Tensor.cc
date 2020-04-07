/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tensors of rank 1 are numbers in all dimensions.
  Tensors of rank 2 are square matrices in all dimensions.
  Only symmetric tensors of rank 4 are are considered here.
*/

#include <iostream>
#include <cmath>

#include "Point.hh"
#include "lapack.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor
****************************************************************** */
Tensor::Tensor(const Tensor& T)
{
  int d = T.dimension();
  int rank = T.rank();
  double* data = T.data();

  if (d && rank) {
    data_ = NULL;
    int mem = Init(d, rank);
    for (int i = 0; i < mem; i++) data_[i] = data[i];
  } else {
    d_ = rank_ = size_ = 0;
    data_ = NULL;
  }
}


/* ******************************************************************
* Constructor.
* Warining: no check of data validity is performed. 
****************************************************************** */
Tensor::Tensor(int d, int rank, const double* data)
{
  size_ = WHETSTONE_TENSOR_SIZE[d - 1][rank - 1];
  int mem = size_ * size_;

  data_ = new double[mem];

  d_ = d;
  rank_ = rank;
  for (int i = 0; i < mem; i++) data_[i] = data[i];
}


/* ******************************************************************
* Initialization of a tensor of rank 1, 2 or 4. 
****************************************************************** */
int Tensor::Init(int d, int rank)
{
  size_ = WHETSTONE_TENSOR_SIZE[d - 1][rank - 1];
  int mem = size_ * size_;

  if (data_) delete[] data_;
  data_ = new double[mem];

  d_ = d;
  rank_ = rank;
  for (int i = 0; i < mem; i++) data_[i] = 0.0;

  return mem;
}


/* ******************************************************************
* Assign constant value to the tensor entries 
****************************************************************** */
void Tensor::PutScalar(double val)
{
  if (! data_) return;

  size_ = WHETSTONE_TENSOR_SIZE[d_ - 1][rank_ - 1];
  int mem = size_ * size_;
  for (int i = 0; i < mem; i++) data_[i] = val;
}


/* ******************************************************************
* Trace operation with tensors of rank 1 and 2
****************************************************************** */
double Tensor::Trace() const
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
void Tensor::Inverse()
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

  } else {
    int info, ipiv[size_];
    double work[size_];
    DGETRF_F77(&size_, &size_, data_, &size_, ipiv, &info);
    DGETRI_F77(&size_, data_, &size_, ipiv, work, &size_, &info);
  }
}


/* ******************************************************************
* Pseudo-inverse operation with tensors of rank 1 and 2
* The algorithm is based on eigenvector decomposition. All eigenvalues
* below the tolerance times the largest eigenvale value are neglected.
****************************************************************** */
void Tensor::PseudoInverse()
{
  if (size_ == 1) {
    if (data_[0] != 0.0) data_[0] = 1.0 / data_[0];

  } else {
    int n = size_; 
    int lwork(3 * n), info;
    double S[n], work[lwork];

    Tensor T(*this);
    DSYEV_F77("V", "U", &n, T.data(), &n, S, work, &lwork, &info);

    // pseudo-invert diagonal matrix S
    double norm_inf(fabs(S[0]));
    for (int i = 1; i < n; i++) {
      norm_inf = std::max(norm_inf, fabs(S[i]));
    } 

    double eps = norm_inf * 1e-15;
    for (int i = 0; i < n; i++) {
      double tmp(fabs(S[i]));
      if (tmp > eps) { 
        S[i] = 1.0 / S[i];
      } else {
        S[i] = 0.0;
      }
    }

    // calculate pseudo inverse pinv(A) = V * pinv(S) * V^t
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        double tmp(0.0);
        for (int k = 0; k < n; k++) {
          tmp += T(i, k) * S[k] * T(j, k);
        }
        (*this)(i, j) = tmp;
        (*this)(j, i) = tmp;
      }
    }
  }
}


/* ******************************************************************
* Matrix of co-factors
****************************************************************** */
Tensor Tensor::Cofactors() const
{
  Tensor C(d_, rank_);
  double* dataC = C.data();

  if (size_ == 1) {
    dataC[0] = 1.0;

  } else if (size_ == 2) {
    dataC[0] = data_[3];
    dataC[3] = data_[0];

    dataC[1] = -data_[2];
    dataC[2] = -data_[1];

  } else if (size_ == 3) { 
    dataC[0] = data_[4] * data_[8] - data_[5] * data_[7];
    dataC[1] = data_[5] * data_[6] - data_[3] * data_[8];
    dataC[2] = data_[3] * data_[7] - data_[4] * data_[6];

    dataC[3] = data_[2] * data_[7] - data_[1] * data_[8];
    dataC[4] = data_[0] * data_[8] - data_[2] * data_[6];
    dataC[5] = data_[1] * data_[6] - data_[0] * data_[7];

    dataC[6] = data_[1] * data_[5] - data_[2] * data_[4];
    dataC[7] = data_[2] * data_[3] - data_[0] * data_[5];
    dataC[8] = data_[0] * data_[4] - data_[1] * data_[3];
  }

  return C;
}


/* ******************************************************************
* Transpose operator for non-symmetric tensors.
****************************************************************** */
void Tensor::Transpose()
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
double Tensor::Det() const
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
* Symmetrizing the tensors of rank 2. 
****************************************************************** */
void Tensor::SymmetricPart()
{
  if (rank_ == 2 && d_ == 2) {
    double tmp = (data_[1] + data_[2]) / 2;
    data_[1] = tmp;
    data_[2] = tmp;
  } else if (rank_ == 2 && d_ == 3) {
    double tmp = (data_[1] + data_[3]) / 2;
    data_[1] = tmp;
    data_[3] = tmp;
    tmp = (data_[2] + data_[6]) / 2;
    data_[2] = tmp;
    data_[6] = tmp;
    tmp = (data_[5] + data_[7]) / 2;
    data_[5] = tmp;
    data_[7] = tmp;
  }
}


/* ******************************************************************
* Check that matrix is zero.
****************************************************************** */
bool Tensor::isZero()
{
  for (int i = 0; i < size_ * size_; i++) {
    if (data_[i] != 0.0) return false;
  }
  return true;
}


/* ******************************************************************
* Spectral bounds of symmetric tensors of rank 1 and 2
****************************************************************** */
void Tensor::SpectralBounds(double* lower, double* upper) const
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
    int n = size_; 
    int lwork(3 * n), info;
    double S[n], work[lwork];
    
    Tensor T(*this);
    DSYEV_F77("N", "U", &n, T.data(), &n, S, work, &lwork, &info);
    *lower = S[0];
    *upper = S[n - 1];
  }
}


/* ******************************************************************
* Elementary operations with a constant. Since we use Voigt notation, 
* the identity tensor equals the identity matrix.
****************************************************************** */
Tensor& Tensor::operator*=(double c)
{
  for (int i = 0; i < size_*size_; i++) data_[i] *= c;
  return *this;
}


Tensor& Tensor::operator+=(double c)
{
  for (int i = 0; i < size_*size_; i += size_ + 1) data_[i] += c;
  return *this;
}


/* ******************************************************************
* Elementary operations with a tensor.
****************************************************************** */
Tensor& Tensor::operator+=(const Tensor& T)
{
  double* data = T.data();
  for (int i = 0; i < size_*size_; ++i) data_[i] += data[i];
  return *this;
}

Tensor& Tensor::operator-=(const Tensor& T)
{
  double* data = T.data();
  for (int i = 0; i < size_*size_; ++i) data_[i] -= data[i];
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

  int mem = Init(d, rank);
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
    p2.set(0.0);
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        p2[j] += (*data) * p[i];
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

  if (rank1 == 4 && rank2 == 2) {
    int n = T1.size();
    double *tmp1 = new double[n];
    double *tmp2 = new double[n];

    tmp1[0] = data2[0];
    tmp1[1] = data2[d + 1];
    if (d == 2) {
      tmp1[2] = data2[1];
      tmp1[3] = data2[2];
    } else if (d == 3) {
      tmp1[2] = data2[8];
      tmp1[3] = data2[1];
      tmp1[4] = data2[3];
      tmp1[5] = data2[5];
      tmp1[6] = data2[7];
      tmp1[7] = data2[2];
      tmp1[8] = data2[6];
    }

    T3.Init(d, rank2);

    for (int i = 0; i < n; i++) {
      double s(0.0);
      for (int j = 0; j < n; j++) s += T1(i, j) * tmp1[j];
      tmp2[i] = s;
    }

    T3(0, 0) = tmp2[0];
    T3(1, 1) = tmp2[1];
    if (d == 2) {
      T3(1, 0) = tmp2[2];
      T3(0, 1) = tmp2[3];
    } else if (d == 3) {
      T3(2, 3) = tmp2[2];
      T3(1, 0) = tmp2[3];
      T3(0, 1) = tmp2[4];
      T3(2, 1) = tmp2[5];
      T3(1, 2) = tmp2[6];
      T3(2, 0) = tmp2[7];
      T3(0, 2) = tmp2[8];
    }

    delete [] tmp1;
    delete [] tmp2;

  } else if (rank1 == 1) {
    int mem = T3.Init(d, rank2);
    double *data3 = T3.data();
    for (int i = 0; i < mem; i++) data3[i] = data2[i] * data1[0];

  } else if (rank2 == 1) {
    int mem = T3.Init(d, rank1);
    double *data3 = T3.data();
    for (int i = 0; i < mem; i++) data3[i] = data1[i] * data2[0];

  } else if (rank2 == 2) {
    T3.Init(d, 2);
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
* Dot product of tensors of equal rank.
****************************************************************** */
double DotTensor(const Tensor& T1, const Tensor& T2)
{
  double *data1 = T1.data(), *data2 = T2.data();
  int mem = T1.size() * T1.size();

  double s(0.0); 
  for (int i = 0; i < mem; i++ ) s += data1[i] * data2[i];
  return s;
}


/* ******************************************************************
* Miscaleneous routines: diagonal tensor
****************************************************************** */
void Tensor::MakeDiagonal(double s)
{
  if (! data_) return;

  int mem = size_ * size_;
  for (int i = 1; i < mem; i++) data_[i] = 0.0;
  for (int i = 0; i < mem; i += d_ + 1) data_[i] = s; 
}


/* ******************************************************************
* Miscaleneous routines: populate tensors of rank 2
****************************************************************** */
int Tensor::SetColumn(int column, const AmanziGeometry::Point& p)
{
  if (rank_ == 2) {
    for (int i = 0; i < d_; i++) (*this)(i, column) = p[i];
    return 1;
  }
  return -1;
}


int Tensor::SetRow(int row, const AmanziGeometry::Point& p)
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


/* ******************************************************************
* Convert tensor to a vector and reverse. Used for parallel 
* distribution of tensors. We assume that size of v sufficient to 
* contain tensor of rank 2.
****************************************************************** */
void TensorToVector(const Tensor& T, DenseVector& v) {
  const double* data1 = T.data(); 
  double* data2 = v.Values(); 
  
  if (T.rank() == 2) {
    int mem = T.size() * T.size();
    for (int i = 0; i < mem; ++i) data2[i] = data1[i]; 
  } else if (T.rank() == 1) {
    int d = T.dimension();
    int mem = WHETSTONE_TENSOR_SIZE[d - 1][1];  // rank 2
    v.PutScalar(0.0);
    for (int i = 0; i < mem * mem; i += d + 1) data2[i] = data1[0]; 
  }
}

void VectorToTensor(const DenseVector& v, Tensor& T) {
  AMANZI_ASSERT(v.NumRows() == T.size() * T.size());

  const double* data1 = v.Values(); 
  double* data2 = T.data(); 
  for (int i = 0; i < v.NumRows(); ++i) {
    data2[i] = data1[i]; 
  }
}


/* ******************************************************************
* Rotations
****************************************************************** */
Tensor RotationMatrix90(const AmanziGeometry::Point& u, bool ccw)
{
  AMANZI_ASSERT(u.dim() == 3);
  int sign = (ccw) ? -1 : 1;

  Tensor T(3, 2);
  T(0, 0) = u[0] * u[0];
  T(0, 1) = u[0] * u[1] + u[2] * sign;
  T(0, 2) = u[0] * u[2] - u[1] * sign;

  T(1, 0) = u[1] * u[0] - u[2] * sign;
  T(1, 1) = u[1] * u[1];
  T(1, 2) = u[1] * u[2] + u[0] * sign;

  T(2, 0) = u[2] * u[0] + u[1] * sign;
  T(2, 1) = u[2] * u[1] - u[0] * sign;
  T(2, 2) = u[2] * u[2];

  return T;
}

}  // namespace WhetStone
}  // namespace Amanzi

