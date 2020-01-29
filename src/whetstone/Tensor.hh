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
#include "exceptions.hh"
#include "errors.hh"

namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_TENSOR_SIZE[3][4] = {{1, 1, 0, 1},
                                         {1, 2, 0, 3},
                                         {1, 3, 0, 6 }};

class Tensor {
 public:


  KOKKOS_INLINE_FUNCTION Tensor(): d_(0), rank_(0), size_(0) {}

  Tensor(const int& d, const int& rank) {
    Init(d, rank);
  }

  Tensor(const Tensor& other) : Tensor(other.d_, other.rank_) {
    Kokkos::deep_copy(data_, other.data_);
  }

  // Tensor(Tensor&& other) :
  //     d_(other.d_),
  //     rank_(other.rank_),
  //     size_(other.size_),
  //     data_(other.data_) {}
  
  // deep copy
  KOKKOS_INLINE_FUNCTION
  void assign(const Tensor& other) {
    if (this != &other) {
      assert(size_ == other.size_);
      for(int i = 0 ; i < size_; ++i){
        data_[i] = other.data_[i]; 
      }
    }
  }

  // this may change to a shallow copy, but would require careful refactor of
  // existing users.
  Tensor& operator=(const Tensor& other) = default;
  
  /*******************************************************************
  * Initialization of a tensor of rank 1 (scalar), 2 (matrix) or 4.
  ****************************************************************** */
  int Init(int d, int rank)
  {
    int mem;
    if (d <= 0 || rank <= 0) {
      d_ = rank_ = size_ = mem = 0;
      return 0;
    } else if (d > 3) {
      Errors::Message msg("Tensor dimension exceeds limit (3).");
      Exceptions::amanzi_throw(msg);
    } else if (rank > 4) {
      Errors::Message msg("Tensor rank exceeds limit (4).");
      Exceptions::amanzi_throw(msg);
    }
    d_ = d;
    rank_ = rank;
    size_ = WHETSTONE_TENSOR_SIZE[d_ - 1][rank_ - 1];
    mem = size_ * size_;
    Kokkos::resize(data_, mem); 
    Kokkos::deep_copy(data_, 0.);
    return mem;
  }

  /* ******************************************************************
  * Assign constant value to the tensor entries
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    // Need to have a copy for the device 
    // \TODO find another way 
    if (data_.extent(0) != size_*size_) return;
    int mem = size_ * size_;
    for(int i = 0 ; i < mem; ++i){
      data_[i] = val;
    }
  }

  /* ******************************************************************
  * Trace operation with tensors of rank 1 and 2
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION double Trace() const 
  {
    double s = 0.0;
    if (rank_ <= 2) {
      for (int i = 0; i < size_; i++) s += (*this)(i, i);
    }
    return s;
  }


  /* ******************************************************************
  * Transpose operator for non-symmetric tensors.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void Transpose()
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
  KOKKOS_INLINE_FUNCTION double Det() const
  {
    double det = 0.0;
    if (rank_ == 2 && d_ == 2) {
      det = data_[0] * data_[3] - data_[1] * data_[2];
    } else if (rank_ == 2 && d_ == 3) {
      det = data_[0] * data_[4] * data_[8] + data_[2] * data_[3] * data_[7] +
            data_[1] * data_[5] * data_[6] - data_[2] * data_[4] * data_[6] -
            data_[1] * data_[3] * data_[8] - data_[0] * data_[5] * data_[7];
    }
    return det;
  }


  /* ******************************************************************
  * Check that matrix is zero.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION bool isZero()
  {
    for (int i = 0; i < size_ * size_; i++) {
      if (data_[i] != 0.0) return false;
    }
    return true;
  }

  void Inverse();
  void PseudoInverse();
  Tensor Cofactors() const;
  void SymmetricPart();
  void SpectralBounds(double* lower, double* upper) const;

  /* ******************************************************************
  * Elementary operations with a constant. Since we use Voigt notation,
  * the identity tensor equals the identity matrix.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION Tensor& operator*=(double c)
  {
    for (int i = 0; i < size_ * size_; i++) data_[i] *= c;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION Tensor& operator+=(double c)
  {
    for (int i = 0; i < size_ * size_; i += size_ + 1) data_[i] += c;
    return *this;
  }

  /* ******************************************************************
  * Elementary operations with a tensor.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION Tensor& operator-=(const Tensor& T)
  {
    Kokkos::View<double*> data = T.data(); 
    for (int i = 0; i < size_ * size_; ++i) data_[i] -= data[i];
    return *this;
  }

  //Tensor& operator=(const Tensor& T);
  //friend AmanziGeometry::Point
  //operator*(const Tensor& T, const AmanziGeometry::Point& p);
  friend Tensor operator*(const Tensor& T1, const Tensor& T2);

  /* ******************************************************************
  * Dot product of tensors of equal rank.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION double DotTensor(
    const Tensor& T1, const Tensor& T2)
  {
    Kokkos::View<double*> data1 = T1.data(); 
    Kokkos::View<double*> data2 = T2.data(); 
    int mem = T1.size() * T1.size();

    double s(0.0);
    for (int i = 0; i < mem; i++) s += data1[i] * data2[i];
    return s;
  }


  /* ******************************************************************
  * Miscaleneous routines: diagonal tensor
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void MakeDiagonal(double s)
  {
    if (data_.extent(0) != size_*size_) return;

    int mem = size_ * size_;
    for (int i = 1; i < mem; i++) data_[i] = 0.0;
    for (int i = 0; i < mem; i += d_ + 1) data_[i] = s;
  }

  /* ******************************************************************
  * Miscaleneous routines: populate tensors of rank 2
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION int SetColumn(
    int column, const AmanziGeometry::Point& p)
  {
    if (rank_ == 2) {
      for (int i = 0; i < d_; i++) (*this)(i, column) = p[i];
      return 1;
    }
    return -1;
  }


  KOKKOS_INLINE_FUNCTION int SetRow(
    int row, const AmanziGeometry::Point& p)
  {
    if (rank_ == 2) {
      for (int i = 0; i < d_; i++) (*this)(row, i) = p[i];
      return 1;
    }
    return 0;
  }
  
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
