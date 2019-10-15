/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "FunctionBilinear.hh"
#include "errors.hh"

namespace Amanzi {

FunctionBilinear::FunctionBilinear(const Kokkos::View<double*>& x,
                                   const Kokkos::View<double*>& y,
                                   const Kokkos::View<double**>& v,
                                   const int xi, const int yi)
  : x_(x), y_(y), v_(v), xi_(xi), yi_(yi)
{
  check_args(x, y, v);
}


void
FunctionBilinear::check_args(const Kokkos::View<double*>& x,
                             const Kokkos::View<double*>& y,
                             const Kokkos::View<double**>& v) const
{
  if (x.extent(0) != v.extent(0)) {
    Errors::Message m;
    m << "the number of x values and row in v differ";
    Exceptions::amanzi_throw(m);
  }
  if (y.extent(0) != v.extent(1)) {
    Errors::Message m;
    m << "the number of y values and columns in v differ";
    Exceptions::amanzi_throw(m);
  }
  if (x.extent(0) < 2) {
    Errors::Message m;
    m << "at least two rows values must be given";
    Exceptions::amanzi_throw(m);
  }
  if (y.extent(0) < 2) {
    Errors::Message m;
    m << "at least two column values must be given";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 1; j < x.extent(0); ++j) {
    if (x[j] <= x[j - 1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
  for (int j = 1; j < y.extent(0); ++j) {
    if (y[j] <= y[j - 1]) {
      Errors::Message m;
      m << "y values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
}

double
FunctionBilinear::operator()(const Kokkos::View<double*>& x) const
{
  double v;
  int nx = x_.extent(0);
  int ny = y_.extent(0);
  double xv = x[xi_];
  double yv = x[yi_];
  // if xv and yv are out of bounds
  if (xv <= x_[0] && yv <= y_[0]) {
    v = v_(0, 0);
  } else if (xv >= x_[nx - 1] && yv <= y_[0]) {
    v = v_(nx - 1, 0);
  } else if (xv >= x_[nx - 1] && yv >= y_[ny - 1]) {
    v = v_(nx - 1, ny - 1);
  } else if (xv <= x_[0] && yv >= y_[ny - 1]) {
    v = v_(0, ny - 1);
  } else {
    // binary search to find interval containing xv
    int j1 = 0, j2 = nx - 1;
    while (j2 - j1 > 1) {
      int j = (j1 + j2) / 2;
      if (xv >= x_[j]) { // right continuous
                         // if (xv > x_[j]) { // left continuous
        j1 = j;
      } else {
        j2 = j;
      }
    }
    // binary search to find interval containing yv
    int k1 = 0, k2 = ny - 1;
    while (k2 - k1 > 1) {
      int k = (k1 + k2) / 2;
      if (yv >= y_[k]) { // right continuous
                         // if (yv > y_[k]) { // left continuous
        k1 = k;
      } else {
        k2 = k;
      }
    }
    // if only xv is out of bounds, linear interpolation
    if (xv <= x_[0] && yv > y_[0] && yv < y_[ny - 1]) {
      v = v_(0, k1) +
          ((v_(0, k2) - v_(0, k1)) / (y_[k2] - y_[k1])) * (yv - y_[k1]);
    } else if (xv > x_[nx - 1] && yv > y_[0] && yv < y_[ny - 1]) {
      v =
        v_(nx - 1, k1) +
        ((v_(nx - 1, k2) - v_(nx - 1, k1)) / (y_[k2] - y_[k1])) * (yv - y_[k1]);
      // if only yv is out of bounds, linear interpolation
    } else if (yv <= y_[0] && xv > x_[0] && xv < x_[nx - 1]) {
      v = v_(j1, 0) +
          ((v_(j2, 0) - v_(j1, 0)) / (x_[j2] - x_[j1])) * (xv - x_[j1]);
    } else if (yv > y_[ny - 1] && xv > x_[0] && xv < x_[nx - 1]) {
      v =
        v_(j1, ny - 1) +
        ((v_(j2, ny - 1) - v_(j1, ny - 1)) / (x_[j2] - x_[j1])) * (xv - x_[j1]);
    } else {
      // bilinear interpolation
      v = v_(j1, k1) * (x_[j2] - xv) * (y_[k2] - yv) +
          v_(j2, k1) * (xv - x_[j1]) * (y_[k2] - yv) +
          v_(j1, k2) * (x_[j2] - xv) * (yv - y_[k1]) +
          v_(j2, k2) * (xv - x_[j1]) * (yv - y_[k1]);
      v = v / ((x_[j2] - x_[j1]) * (y_[k2] - y_[k1]));
    }
  }
  return v;
}


} // namespace Amanzi
