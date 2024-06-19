/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon ecoon@lanl.gov
*/

/* -------------------------------------------------------------------------
  Spline

  Simple spline functor which is used for smoothing.

------------------------------------------------------------------------- */


#ifndef AMANZI_UTILS_SPLINE_HH_
#define AMANZI_UTILS_SPLINE_HH_

#include "Kokkos_Core.hpp"
#include "dbc.hh"

namespace Amanzi {
namespace Utils {

class Spline {
 public:
  Spline() {};

  Spline(double x1, double y1, double dy1, double x2, double y2, double dy2)
  {
    Setup(x1, y1, dy1, x2, y2, dy2);
  }

  void Setup(double x1, double y1, double dy1, double x2, double y2, double dy2)
  {
    x1_ = x1;
    y1_ = y1;
    dy1_ = dy1;
    x2_ = x2;
    y2_ = y2;
    dy2_ = dy2;

    AMANZI_ASSERT(x1_ < x2_);
    if (y2_ >= y1_) {
      AMANZI_ASSERT(dy1_ >= 0.);
      AMANZI_ASSERT(dy2_ >= 0.);
    }
    if (y1_ >= y2_) {
      AMANZI_ASSERT(dy1_ <= 0.);
      AMANZI_ASSERT(dy2_ <= 0.);
    }
  }


  KOKKOS_INLINE_FUNCTION
  double operator()(double x) const { return Value(x); }

  KOKKOS_INLINE_FUNCTION
  double Value(double x) const {
    assert(x1_ <= x && x <= x2_);
    double t = T(x);
    return Kokkos::pow(1 - t, 2) * ((1 + 2 * t) * y1_ + t * (x2_ - x1_) * dy1_) +
      Kokkos::pow(t, 2) * ((3 - 2 * t) * y2_ + (t - 1) * (x2_ - x1_) * dy2_);
  }

  KOKKOS_INLINE_FUNCTION
  double Derivative(double x) const {
    assert(x1_ <= x && x <= x2_);
    double t = T(x);
    double dtdx = 1. / (x2_ - x1_);
    double dydt =
      (6 * Kokkos::pow(t, 2) - 6 * t) * y1_ + (3 * Kokkos::pow(t, 2) - 4 * t + 1) * (x2_ - x1_) * dy1_ +
      (-6 * Kokkos::pow(t, 2) + 6 * t) * y2_ + (3 * Kokkos::pow(t, 2) - 2 * t) * (x2_ - x1_) * dy2_;
    return dydt * dtdx;
  }

  KOKKOS_INLINE_FUNCTION
  double T(double x) const {
    return (x - x1_) / (x2_ - x1_);
  }

 private:
  double x1_, x2_, y1_, y2_, dy1_, dy2_;
};


} // namespace Utils
} // namespace Amanzi

#endif
