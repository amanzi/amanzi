/*
Tensors of rank 1 are numbers in all dimensions.
Tensors of rank 2 are square matrices in all dimensions.
Only symmetric tensors of rank 4 are are considered here.
*/

#include  <iostream>
#include  <math.h>

#include "Teuchos_LAPACK.hpp"

#include "Point.hh"
#include "tensor.hpp"

using namespace std;

namespace Amanzi {
namespace WhetStone {


/* ******************************************************************
* Constructors
****************************************************************** */
Tensor::Tensor(const Tensor& T)
{ 
  int d = T.get_dimension();
  int rank = T.get_rank();
  double* data = T.get_data();

  if (d && rank) {
    data_ = NULL;
    int mem = init(d, rank); 
    for (int i=0; i<mem; i++) data_[i] = data[i];
  } else {
    d_ = rank_ = 0;
    data_ = NULL;
  }
}


/* ******************************************************************
* Initialization of a tensor of rank 1, 2 or 4. 
****************************************************************** */
int Tensor::init(const int d, const int rank) 
{
  int mem = WHETSTONE_TENSOR_SIZE[d-1][rank-1] * WHETSTONE_TENSOR_SIZE[d-1][rank-1]; 
  if (data_) delete[] data_;
  data_ = new double[mem]; 
 
  d_ = d;
  rank_ = rank;
  for (int i=0; i<mem; i++) data_[i] = 0.0;

  return mem;
}


/* ******************************************************************
* Trace operation with tensors of rank 1 and 2
****************************************************************** */
double Tensor::trace()
{
  double s = 0.0;
  if (rank_ <= 2) {
    int size = WHETSTONE_TENSOR_SIZE[d_-1][rank_-1];
    double s = 0.0;
    for (int i=0; i<size; i++) s += (*this)(i, i);
  }
  return s; 
}


/* ******************************************************************
* Inverse operation with tensors of rank 1 and 2
****************************************************************** */
void Tensor::inverse()
{
  int size = WHETSTONE_TENSOR_SIZE[d_-1][rank_-1];
  if (size == 1) { 
    data_[0] = 1.0 / data_[0];
  }
  else if (size == 2) {  // We use inverse formula based on minors
    double det = data_[0] * data_[3] - data_[1] * data_[2];

    double a = data_[0];
    data_[0] = data_[3] / det;
    data_[3] = a / det;

    a = data_[1];
    data_[1] = -data_[2] / det;
    data_[2] = -a / det;
  }
  else if (rank_ <= 2) {
    Teuchos::LAPACK<int, double> lapack;
    int info;
    int ipiv[size];
    double work[size];
    lapack.GETRF(size, size, data_, size, ipiv, &info);
    lapack.GETRI(size, data_, size, ipiv, work, size, &info); 
  }
}


/* ******************************************************************
* Trace operation with tensors of rank 1 and 2
****************************************************************** */
Tensor& Tensor::operator*=(const double& c)
{
  if (rank_ <= 2) {
    int size = WHETSTONE_TENSOR_SIZE[d_-1][rank_-1];
    for (int i=0; i<size; i++) (*this)(i, i) *= c;
  }
  return *this; 
}


/* ******************************************************************
* Rescale the tensor 
****************************************************************** */
Tensor operator*(Tensor& T, const double& c)
{
  int rank = T.get_rank();
  int d = T.get_dimension();
  double* data = T.get_data();

  Tensor T1(d, rank);
  double* data1 = T1.get_data();
  
  int size = WHETSTONE_TENSOR_SIZE[d-1][rank-1];
  for (int i=0; i<size; i++) data1[i] = data[i] * c;

  return T1;
}


/* ******************************************************************
* First convolution operation for tensors of rank 1 and 2. 
****************************************************************** */
AmanziGeometry::Point operator*(Tensor& T, const AmanziGeometry::Point& p)
{
  int rank = T.get_rank();
  int d = T.get_dimension();
  double* data = T.get_data();

  AmanziGeometry::Point p2(p.dim());
  if (rank == 1) { 
    p2 = data[0] * p;
    return p2;
  } 
  else if (rank == 2) {
    for (int i=0; i<d; i++) {
      p2[i] = 0.0;
      for (int j=0; j<d; j++) {
        p2[i] += (*data) * p[j];
        data++;
      }
    }
    return p2;
  }
  else if (rank == 4) {
    return p;  // undefined operation (lipnikov@lanl.gov) 
  }
}


/* ******************************************************************
* Second convolution operation for tensors of rank 1, 2, and 4
****************************************************************** */
Tensor operator*(Tensor& T1, Tensor& T2)
{
  int d = T1.get_dimension();  // the dimensions should be equals
  int rank1 = T1.get_rank(), rank2 = T2.get_rank();
  double *data1 = T1.get_data(), *data2 = T2.get_data();

  if (d == 2 && rank1 == 4 && rank2 == 2) {
    double a0, b0, c0;
    a0 = T2(0, 0);  
    b0 = T2(1, 1);  
    c0 = T2(0, 1);

    Tensor T3(d, rank2);
    T3(0, 0) = T1(0, 0) * a0 + T1(0, 1) * b0 + T1(0, 2) * c0;
    T3(1, 1) = T1(1, 0) * a0 + T1(1, 1) * b0 + T1(1, 2) * c0;
    T3(1, 0) = T3(0, 1) = T1(2, 0) * a0 + T1(2, 1) * b0 + T1(2, 2) * c0;
    return T3;
  }
}


/* ******************************************************************
* Miscaleneous routines: populate tensors of rank 2
****************************************************************** */
int Tensor::add_column(const int column, const AmanziGeometry::Point& p)
{
  if (rank_ == 2) {
    for (int i=0; i<d_; i++) (*this)(i, column) = p[i];
    return 1;
  }
  return -1;
}


int Tensor::add_row(const int row, const AmanziGeometry::Point& p)
{
  if (rank_ == 2) {
    for (int i=0; i<d_; i++) (*this)(row, i) = p[i];
    return 1;
  }
  return 0;
}


/* ******************************************************************
* Miscaleneous routines: print
****************************************************************** */
std::ostream& operator<<(std::ostream& os, const Tensor& T)
{
  int d = T.get_dimension();
  int rank = T.get_rank();
  int size = WHETSTONE_TENSOR_SIZE[d-1][rank-1];

  os << "Tensor dimension=" << d << "  rank=" << rank << endl; 
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) os << T(i, j) << " ";
    os << endl;
  }
  return os;
}

}  // namespace WhetStone
}  // namespace Amanzi

 
