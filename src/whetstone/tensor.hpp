#ifndef   __TENSOR_HXX
#define   __TENSOR_HXX

/*
Tensors of rank 1 are numbers in all dimensions.
Tensors of rank 2 are square matrices in all dimensions.
Only symmetric tensors of rank 4 are considered here.
*/

#include  <iostream>
#include  <math.h>

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_TENSOR_SIZE[3][3] = {1, 1, 1, 
                                         1, 2, 3,
                                         1, 3, 6};

class Tensor {
 public:
  Tensor() { d_ = rank_ = 0; data_ = NULL; }
  Tensor(const Tensor& T);
  Tensor(const int d, const int rank) { data_ = NULL; init(d, rank); }
  ~Tensor() { if (data_) delete[] data_; }

  // primary members
  int init(const int d, const int rank);
  double trace();
  double det();
  void inverse();
  void transpose();

  // elementary operators
  friend Tensor operator*(Tensor& T, const double& c);
  friend AmanziGeometry::Point operator*(Tensor& T, const AmanziGeometry::Point& p);
  friend Tensor operator*(Tensor& T1, Tensor& T2);
  Tensor& operator*=(const double& c);

  // access members
  double& operator()(int i, int j) { return data_[i * WHETSTONE_TENSOR_SIZE[d_-1][rank_-1] + j]; }
  double& operator()(int i, int j) const { return data_[i * WHETSTONE_TENSOR_SIZE[d_-1][rank_-1] + j]; }
  int add_column(const int column, const AmanziGeometry::Point& p); 
  int add_row(const int row, const AmanziGeometry::Point& p); 

  int get_dimension() const { return d_; }
  int get_rank() const { return rank_; }
  double* get_data() { return data_; }
  double* get_data() const { return data_; }

  // miscaleneous
  friend std::ostream& operator<<(std::ostream& os, const Tensor& T);
 
 private:
  int d_, rank_;
  double* data_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
