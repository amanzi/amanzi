#ifndef   __TENSOR_HH__
#define   __TENSOR_HH__

/*
Tensors of rank 1 are numbers in all dimensions.
Tensors of rank 2 are square matrices in all dimensions.
Only symmetric tensors of rank 4 are considered here.
*/

#include <iostream>
#include <cmath>

#include "Point.hh"


namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_TENSOR_SIZE[3][3] = {1, 1, 1, 
                                         1, 2, 3,
                                         1, 3, 6};

class Tensor {
 public:
  Tensor() { 
    d_ = rank_ = size_ = 0; 
    data_ = NULL;
  }
  Tensor(const Tensor& T);
  Tensor(int d, int rank) {
    data_ = NULL;
    init(d, rank);
  }
  Tensor(int d, int rank, const double* data);
  ~Tensor() { if (data_) delete[] data_; }

  // primary members
  int init(int d, int rank);
  double trace() const;
  double det();
  void inverse();
  void transpose();
  double determinant();
  void spectral_bounds(double* lower, double* upper) const;

  // elementary operators
  Tensor& operator*=(const double& c);
  Tensor& operator=(const Tensor& T);
  friend AmanziGeometry::Point operator*(const Tensor& T, const AmanziGeometry::Point& p);
  friend Tensor operator*(const Tensor& T1, const Tensor& T2);

  // access members
  double& operator()(int i, int j) { return data_[i * size_ + j]; }
  double& operator()(int i, int j) const { return data_[i * size_ + j]; }
  int add_column(const int column, const AmanziGeometry::Point& p); 
  int add_row(const int row, const AmanziGeometry::Point& p); 

  int dimension() const { return d_; }
  int rank() const { return rank_; }
  int size() const { return size_; }
  double* data() { return data_; }
  double* data() const { return data_; }

  // miscaleneous
  friend std::ostream& operator<<(std::ostream& os, const Tensor& T);
 
 private:
  int d_, rank_, size_;
  double* data_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
