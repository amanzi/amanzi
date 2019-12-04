/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tensors of rank 1 are numbers in all dimensions.
  General tensors of rank 2 are square matrices in all dimensions.
  Only symmetric tensors of rank 4 are considered here.
*/

#ifndef AMANZI_TENSOR_HH_
#define AMANZI_TENSOR_HH_

#include <iostream>
#include <cmath>

#include "Point.hh"

#include "DenseVector.hh"

namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_TENSOR_SIZE[3][4] = { { 1, 1, 0, 1 },
                                          { 1, 2, 0, 3 },
                                          { 1, 3, 0, 6 } };

class Tensor {
 public:
  Tensor()
  {
    d_ = rank_ = size_ = 0;
    data_ = NULL;
  }
  Tensor(const Tensor& T);
  Tensor(int d, int rank)
  {
    data_ = NULL;
    Init(d, rank);
  }
  Tensor(int d, int rank, Kokkos::View<double*> data);
  ~Tensor()
  {
    if (data_) delete[] data_;
  }

  // primary members
  int Init(int d, int rank);
  void PutScalar(double val);
  double Trace() const;
  double Det() const;
  void Inverse();
  void PseudoInverse();
  void Transpose();
  Tensor Cofactors() const;
  void SymmetricPart();
  bool isZero();
  void SpectralBounds(double* lower, double* upper) const;

  // elementary operators
  Tensor& operator*=(double c);
  Tensor& operator+=(double c);
  Tensor& operator-=(const Tensor& T);
  Tensor& operator=(const Tensor& T);
  friend AmanziGeometry::Point
  operator*(const Tensor& T, const AmanziGeometry::Point& p);
  friend Tensor operator*(const Tensor& T1, const Tensor& T2);
  friend double DotTensor(const Tensor& T1, const Tensor& T2);

  // initialization
  void MakeDiagonal(double s);
  int SetColumn(int column, const AmanziGeometry::Point& p);
  int SetRow(int row, const AmanziGeometry::Point& p);

  // access members
  double& operator()(int i, int j) { return data_[j * size_ + i]; }
  double& operator()(int i, int j) const { return data_[j * size_ + i]; }

  int dimension() const { return d_; }
  int rank() const { return rank_; }
  int size() const { return size_; }
  double* data() { return data_; }
  double* data() const { return data_; }

  // miscaleneous
  friend std::ostream& operator<<(std::ostream& os, const Tensor& T);

 private:
  int d_, rank_, size_;
  Kokkos::View<double*> data_;
};


struct TensorArray {
  Kokkos::View<double**> data;
  int dim;
  int rank;

  Kokkos::View<double*> operator()(int i) {
    return Kokkos::subview(data, i, Kokkos::ALL);
  }

  WhetStone::Tensor getTensor(int i) {
    return WhetStone::Tensor(dim, rank, this->operator()(i));
  }
};






// non-member functions
// -- comparison operators
inline bool
operator==(const Tensor& T1, const Tensor& T2)
{
  if (T1.rank() != T1.rank()) return false;
  if (T1.size() != T2.size()) return false;

  double* data1 = T1.data();
  double* data2 = T2.data();
  for (int i = 0; i != T1.size(); ++i)
    if (data1[i] != data2[i]) return false;
  return true;
}

inline bool
operator!=(const Tensor& T1, const Tensor& T2)
{
  return !(T1 == T2);
}

// -- expanding tensor to a constant size vector and reverse.
void
TensorToVector(const Tensor& T, DenseVector& v);
void
VectorToTensor(const DenseVector& v, Tensor& T);

} // namespace WhetStone
} // namespace Amanzi

#endif
