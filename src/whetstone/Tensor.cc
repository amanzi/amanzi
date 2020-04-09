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
  Only symmetric tensors of rank 4 are considered here.
*/

#include <iostream>
#include <cmath>

#include "Point.hh"
#include "lapack.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {


/* ******************************************************************
 * Second convolution operation for tensors of rank 1, 2, and 4
 ****************************************************************** */
template<class MEMSPACE>
Tensor<MEMSPACE> operator*(const Tensor<MEMSPACE>& T1, const Tensor<MEMSPACE>& T2)
{
  int d = T1.dimension(); // the dimensions should be equals
  int rank1 = T1.rank(), rank2 = T2.rank();
  Kokkos::View<double*,MEMSPACE> data1 = T1.data();
  Kokkos::View<double*,MEMSPACE> data2 = T2.data();

  Tensor<MEMSPACE> T3;

  if (rank1 == 4 && rank2 == 2) {
    int n = d * (d + 1) / 2;
    Kokkos::View<double*,MEMSPACE> tmp1("tmp1",n);
    Kokkos::View<double*,MEMSPACE> tmp2("tmp2",n);

    tmp1[0] = data2[0];
    tmp1[1] = data2[d + 1];
    if (d == 2) {
      tmp1[2] = data2[1] + data2[2];
    } else if (d == 3) {
      tmp1[2] = data2[8];
      tmp1[3] = data2[1] + data2[3];
      tmp1[4] = data2[5] + data2[7];
      tmp1[5] = data2[2] + data2[6];
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
      T3(1, 0) = T3(0, 1) = tmp2[2];
    } else if (d == 3) {
      T3(2, 3) = tmp2[2];
      T3(0, 1) = T3(1, 0) = tmp2[3];
      T3(1, 2) = T3(2, 1) = tmp2[4];
      T3(0, 2) = T3(2, 0) = tmp2[5];
    }

  } else if (rank1 == 1) {
    int mem = T3.Init(d, rank2);
    Kokkos::View<double*,MEMSPACE> data3 = T3.data();
    for (int i = 0; i < mem; i++) data3[i] = data2[i] * data1[0];

  } else if (rank2 == 1) {
    int mem = T3.Init(d, rank1);
    Kokkos::View<double*,MEMSPACE> data3 = T3.data();
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
 * Miscaleneous routines: print
 ****************************************************************** */
template<class MEMSPACE>
std::ostream&
operator<<(std::ostream& os, const Tensor<MEMSPACE>& T)
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
template<class MEMSPACE>
void
TensorToVector(const Tensor<MEMSPACE>& T, DenseVector& v)
{
  const Kokkos::View<double*,MEMSPACE> data1 = T.data();
  Kokkos::View<double*,MEMSPACE> data2 = v.Values();

  if (T.rank() == 2) {
    int mem = T.size() * T.size();
    for (int i = 0; i < mem; ++i) data2[i] = data1[i];
  } else if (T.rank() == 1) {
    int d = T.dimension();
    int mem = 2;//T.WHETSTONE_TENSOR_SIZE[d - 1][1]; // rank 2
    v.putScalar(0.0);
    for (int i = 0; i < mem * mem; i += d + 1) data2[i] = data1[0];
  }
}

template<class MEMSPACE>
void
VectorToTensor(const DenseVector& v, Tensor<MEMSPACE>& T)
{
  AMANZI_ASSERT(v.NumRows() == T.size() * T.size());

  const Kokkos::View<double*,MEMSPACE> data1 = v.Values();
  Kokkos::View<double*,MEMSPACE> data2 = T.data();
  for (int i = 0; i < v.NumRows(); ++i) { data2[i] = data1[i]; }
}

template<class MEMSPACE>
Tensor<MEMSPACE> Tensor_ONE()
{
  Kokkos::View<double*,MEMSPACE> t_identity("identity", 1);
  t_identity(0) = 1.0;
  return Tensor<MEMSPACE>(t_identity, 1,1,1);
}


} // namespace WhetStone
} // namespace Amanzi
