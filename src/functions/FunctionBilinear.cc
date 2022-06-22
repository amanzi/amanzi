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

FunctionBilinear::FunctionBilinear(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                                   const Kokkos::View<double*,Kokkos::HostSpace>& y,
                                   const Kokkos::View<double**,Kokkos::HostSpace>& v,
                                   const int xi, const int yi)
  : xi_(xi), yi_(yi)
{
  Kokkos::resize(x_,x.extent(0)); 
  Kokkos::resize(y_,y.extent(0)); 
  Kokkos::resize(v_,v.extent(0),v.extent(1)); 

  Kokkos::deep_copy(x_.view_host(),x);
  Kokkos::deep_copy(y_.view_host(),y);
  Kokkos::deep_copy(v_.view_host(),v);

  Kokkos::deep_copy(x_.view_device(),x_.view_host());
  Kokkos::deep_copy(y_.view_device(),y_.view_host());
  Kokkos::deep_copy(v_.view_device(),v_.view_host());

  check_args(x, y, v);
}


void
FunctionBilinear::check_args(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                             const Kokkos::View<double*,Kokkos::HostSpace>& y,
                             const Kokkos::View<double**,Kokkos::HostSpace>& v) const
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
FunctionBilinear::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& x) const
{
  double v;
  int nx = x_.extent(0);
  int ny = y_.extent(0);
  double xv = x[xi_];
  double yv = x[yi_];
  auto vv = v_.view_host(); 
  auto vx = x_.view_host(); 
  auto vy = y_.view_host();
  // if xv and yv are out of bounds
  if (xv <= vx[0] && yv <= vy[0]) {
    v = vv(0, 0);
  } else if (xv >= vx[nx - 1] && yv <= vy[0]) {
    v = vv(nx - 1, 0);
  } else if (xv >= vx[nx - 1] && yv >= vy[ny - 1]) {
    v = vv(nx - 1, ny - 1);
  } else if (xv <= vx[0] && yv >= vy[ny - 1]) {
    v = vv(0, ny - 1);
  } else {
    // binary search to find interval containing xv
    int j1 = 0, j2 = nx - 1;
    while (j2 - j1 > 1) {
      int j = (j1 + j2) / 2;
      if (xv >= vx[j]) { // right continuous
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
      if (yv >= vy[k]) { // right continuous
                         // if (yv > y_[k]) { // left continuous
        k1 = k;
      } else {
        k2 = k;
      }
    }
    // if only xv is out of bounds, linear interpolation
    if (xv <= vx[0] && yv > vy[0] && yv < vy[ny - 1]) {
      v = vv(0, k1) +
          ((vv(0, k2) - vv(0, k1)) / (vy[k2] - vy[k1])) * (yv - vy[k1]);
    } else if (xv > vx[nx - 1] && yv > vy[0] && yv < vy[ny - 1]) {
      v =
        vv(nx - 1, k1) +
        ((vv(nx - 1, k2) - vv(nx - 1, k1)) / (vy[k2] - vy[k1])) * (yv - vy[k1]);
      // if only yv is out of bounds, linear interpolation
    } else if (yv <= vy[0] && xv > vx[0] && xv < vx[nx - 1]) {
      v = vv(j1, 0) +
          ((vv(j2, 0) - vv(j1, 0)) / (vx[j2] - vx[j1])) * (xv - vx[j1]);
    } else if (yv > vy[ny - 1] && xv > vx[0] && xv < vx[nx - 1]) {
      v =
        vv(j1, ny - 1) +
        ((vv(j2, ny - 1) - vv(j1, ny - 1)) / (vx[j2] - vx[j1])) * (xv - vx[j1]);
    } else {
      // bilinear interpolation
      v = vv(j1, k1) * (vx[j2] - xv) * (vy[k2] - yv) +
          vv(j2, k1) * (xv - vx[j1]) * (vy[k2] - yv) +
          vv(j1, k2) * (vx[j2] - xv) * (yv - vy[k1]) +
          vv(j2, k2) * (xv - vx[j1]) * (yv - vy[k1]);
      v = v / ((vx[j2] - vx[j1]) * (vy[k2] - vy[k1]));
    }
  }
  return v;
}


} // namespace Amanzi
