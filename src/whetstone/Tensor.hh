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

#include <Kokkos_Core.hpp>

const int WHETSTONE_TENSOR_SIZE[3][4] = { {1, 1, 0, 1},
                                          {1, 2, 0, 3},
                                          {1, 3, 0, 6 }};

namespace Amanzi {
namespace WhetStone {

class Tensor {
 public:

  KOKKOS_INLINE_FUNCTION Tensor(): d_(0), rank_(0), size_(0){}
  //Tensor(const Tensor& T);
  Tensor(const int& d, const int& rank){ Init(d, rank); }

#if 0 
  KOKKOS_INLINE_FUNCTION Tensor(int d, int rank, const Kokkos::View<double*> data)
  {
    size_ = WHETSTONE_TENSOR_SIZE[d - 1][rank - 1];
    int mem = size_ * size_;
    data_ = data; 
    d_ = d;
    rank_ = rank;
  }
#endif 
  
  /*******************************************************************
  * Initialization of a tensor of rank 1, 2 or 4.
  ****************************************************************** */
  int Init(int d, int rank)
  {
    size_ = WHETSTONE_TENSOR_SIZE[d - 1][rank - 1];
    int mem = size_ * size_;
    Kokkos::resize(data_,mem); 
    d_ = d;
    rank_ = rank;
    for (int i = 0; i < mem; i++) data_[i] = 0.0;
    return mem;
  }

  /* ******************************************************************
  * Assign constant value to the tensor entries
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    if (data_.extent(0) != size_*size_) return;
    size_ = WHETSTONE_TENSOR_SIZE[d_ - 1][rank_ - 1];
    int mem = size_ * size_;
    for (int i = 0; i < mem; i++) data_[i] = val;
  }

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
  //Tensor& operator=(const Tensor& T);
  //friend AmanziGeometry::Point
  //operator*(const Tensor& T, const AmanziGeometry::Point& p);
  friend Tensor operator*(const Tensor& T1, const Tensor& T2);
  friend double DotTensor(const Tensor& T1, const Tensor& T2);

  // initialization
  void MakeDiagonal(double s);
  int SetColumn(int column, const AmanziGeometry::Point& p);
  int SetRow(int row, const AmanziGeometry::Point& p);

  // access members
  KOKKOS_INLINE_FUNCTION double& operator()(int i, int j) { return data_[j * size_ + i]; }
  KOKKOS_INLINE_FUNCTION double& operator()(int i, int j) const { return data_[j * size_ + i]; }

  KOKKOS_INLINE_FUNCTION int dimension() const { return d_; }
  KOKKOS_INLINE_FUNCTION int rank() const { return rank_; }
  KOKKOS_INLINE_FUNCTION int size() const { return size_; }
  KOKKOS_INLINE_FUNCTION Kokkos::View<double*> data() { return data_; }
  KOKKOS_INLINE_FUNCTION Kokkos::View<double*> data() const { return data_; }
  KOKKOS_INLINE_FUNCTION double* data_ptr() { return &data_[0];}
  KOKKOS_INLINE_FUNCTION double* data_ptr() const { return &data_[0];}

  // miscaleneous
  friend std::ostream& operator<<(std::ostream& os, const Tensor& T);

 private:
  int d_, rank_, size_;
  Kokkos::View<double*> data_;
};

// non-member functions
// -- comparison operators
inline bool
operator==(const Tensor& T1, const Tensor& T2)
{
  if (T1.rank() != T1.rank()) return false;
  if (T1.size() != T2.size()) return false;

  Kokkos::View<double*> data1 = T1.data();
  Kokkos::View<double*> data2 = T2.data();
  for (int i = 0; i != T1.size(); ++i)
    if (data1[i] != data2[i]) return false;
  return true;
}

inline bool
operator!=(const Tensor& T1, const Tensor& T2)
{
  return !(T1 == T2);
}

KOKKOS_INLINE_FUNCTION AmanziGeometry::Point 
operator*(const Tensor& T, const AmanziGeometry::Point& p)
{
  int rank = T.rank();
  int d = T.dimension();
  Kokkos::View<double*> data = T.data(); 

  AmanziGeometry::Point p2(p.dim());
  if (rank == 1) {
    p2 = data[0] * p;
    return p2;

  } else if (rank == 2) {
    p2.set(0.0);
    int idx = 0; 
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        p2[j] += data[idx++] * p[i];
      }
    }
    return p2;

  } else if (rank == 4) {
    return p; // undefined operation (lipnikov@lanl.gov)
  }
  return p;
}


// -- expanding tensor to a constant size vector and reverse.
void
TensorToVector(const Tensor& T, DenseVector& v);
void
VectorToTensor(const DenseVector& v, Tensor& T);

} // namespace WhetStone
} // namespace Amanzi

#endif
