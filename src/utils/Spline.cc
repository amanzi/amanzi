/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/* -------------------------------------------------------------------------
  Spline

  Author: Ethan Coon coonet@ornl.gov

 Cubic Hermite spline functor which is used for smoothing.  NOTE: this is NOT
 guaranteed to be monotonic!  In practice I hope it is usually!

------------------------------------------------------------------------- */

#include "Teuchos_SerialDenseMatrix.hpp"
#include "dbc.hh"

#include "Spline.hh"


namespace Amanzi {
namespace Utils {

void
Spline::Setup(double x1, double y1, double dy1, double x2, double y2,
              double dy2)
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


double
Spline::Value(double x)
{
  AMANZI_ASSERT(x1_ <= x <= x2_);
  double t = T(x);
  return std::pow(1 - t, 2) * ((1 + 2 * t) * y1_ + t * (x2_ - x1_) * dy1_) +
         std::pow(t, 2) * ((3 - 2 * t) * y2_ + (t - 1) * (x2_ - x1_) * dy2_);

  // above is a bit cleaner
  // return (2*std::pow(t,3) - 3*std::pow(t,2) + 1 )* y1_
  //     + (std::pow(t,3) - 2*std::pow(t,2) + t) * (x2_ - x1_) * dy1_
  //     + (-2*std::pow(t,3) + 3*std::pow(t,2)) * y2_
  //     + (std::pow(t,3) - std::pow(t,2)) * (x2_ - x1_) * dy2_;
}

double
Spline::Derivative(double x)
{
  AMANZI_ASSERT(x1_ <= x <= x2_);
  double t = T(x);
  double dtdx = 1. / (x2_ - x1_);
  double dydt = (6 * std::pow(t, 2) - 6 * t) * y1_ +
                (3 * std::pow(t, 2) - 4 * t + 1) * (x2_ - x1_) * dy1_ +
                (-6 * std::pow(t, 2) + 6 * t) * y2_ +
                (3 * std::pow(t, 2) - 2 * t) * (x2_ - x1_) * dy2_;
  return dydt * dtdx;
}


} // namespace Utils
} // namespace Amanzi
