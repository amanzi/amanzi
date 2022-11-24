/*
  Geometry

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Special implementation of tensors of rank 2.
*/

#ifndef AMANZI_TENSOR_SIMPLE_HH_
#define AMANZI_TENSOR_SIMPLE_HH_

#include <cmath>

#include "Point.hh"

namespace Amanzi {
namespace AmanziGeometry {

class TensorSimple {
 public:
  TensorSimple()
  {
    d_ = 0;
    for (int i = 0; i < 9; ++i) data_[i] = 0.0;
  }
  ~TensorSimple(){};

  // primary members
  void set(int d) { d_ = d; }

  void Inverse()
  {
    if (d_ == 2) {
      double det = data_[0] * data_[3] - data_[1] * data_[2];

      double a = data_[0];
      data_[0] = data_[3] / det;
      data_[3] = a / det;

      data_[1] /= -det;
      data_[2] /= -det;
    } else if (d_ == 3) {
      double det, copy[9];
      for (int i = 0; i < 9; ++i) copy[i] = data_[i];

      det = data_[0] * data_[4] * data_[8] + data_[2] * data_[3] * data_[7] +
            data_[1] * data_[5] * data_[6] - data_[2] * data_[4] * data_[6] -
            data_[1] * data_[3] * data_[8] - data_[0] * data_[5] * data_[7];

      data_[0] = (copy[4] * copy[8] - copy[5] * copy[7]) / det;
      data_[1] = (copy[2] * copy[7] - copy[1] * copy[8]) / det;
      data_[2] = (copy[1] * copy[5] - copy[2] * copy[4]) / det;
      data_[3] = (copy[5] * copy[6] - copy[3] * copy[8]) / det;
      data_[4] = (copy[0] * copy[8] - copy[2] * copy[6]) / det;
      data_[5] = (copy[2] * copy[3] - copy[0] * copy[5]) / det;
      data_[6] = (copy[3] * copy[7] - copy[4] * copy[6]) / det;
      data_[7] = (copy[1] * copy[6] - copy[0] * copy[7]) / det;
      data_[8] = (copy[0] * copy[4] - copy[1] * copy[4]) / det;
    }
  }

  // elementary operators
  friend Point operator*(const TensorSimple& T, const Point& p)
  {
    int d = p.dim();
    const double* data = &T(0, 0);
    Point p2(d);

    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        p2[j] += (*data) * p[i];
        data++;
      }
    }
    return p2;
  }

  // access members
  int dim() const { return d_; }
  double& operator()(int i, int j) { return data_[j * d_ + i]; }
  const double& operator()(int i, int j) const { return data_[j * d_ + i]; }

  // I/O
  friend std::ostream& operator<<(std::ostream& os, const TensorSimple& T)
  {
    int d = T.dim();
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) os << T(i, j) << " ";
      os << std::endl;
    }
    return os;
  }

 private:
  int d_;
  double data_[9];
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
