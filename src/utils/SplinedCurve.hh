/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/* -------------------------------------------------------------------------
  Spline

  Spline fit of a curve, given points, values, and derivatives of those
  values.  If needed, this could easily provide a constructor without
  derivatives, and use an algorithm to construct these derivatives.

  Provides options for monotonicity preservation.

------------------------------------------------------------------------- */


#ifndef AMANZI_UTILS_SPLINED_CURVE_HH_
#define AMANZI_UTILS_SPLINED_CURVE_HH_

#include <vector>

#include "errors.hh"

namespace Amanzi {
namespace Utils {

class SplinedCurve {
 public:
  enum struct SplineExtrapolation_t {
    CONSTANT, // outside the interval, use the endpoint value
    LINEAR,   // outside the interval, linearly extrapolate
    SPLINE,   // outside the interval continue the spline
    THROW     // outside the interval, throw and error
  };
  typedef std::pair<SplineExtrapolation_t, SplineExtrapolation_t> SplineEndpoints_t;

  class Amanzi_exception_SplineOutOfRange : public Exceptions::Amanzi_exception {};

  SplinedCurve() {}

  SplinedCurve(std::vector<double> x,
               std::vector<double> y,
               std::vector<double> dydx,
               SplineEndpoints_t endpoints = SplineEndpoints_t(SplineExtrapolation_t::THROW,
                                                               SplineExtrapolation_t::THROW),
               bool enforce_monotonicity = false);

  void Setup(std::vector<double> x,
             std::vector<double> y,
             std::vector<double> dydx,
             SplineEndpoints_t endpoints = SplineEndpoints_t(SplineExtrapolation_t::THROW,
                                                             SplineExtrapolation_t::THROW),
             bool enforce_monotonicity = false);

  double operator()(double x) { return Value(x); }
  double Value(double x);
  double Derivative(double x);

 private:
  double T_(double x, int i);
  int interval_(double x);
  void Setup_();

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> dy_;
  std::vector<double> dx_;

  SplineExtrapolation_t left_, right_;

  bool mono_;
};


inline int
SplinedCurve::interval_(double x)
{
  return std::lower_bound(x_.begin(), x_.end(), x) - x_.begin() - 1;
}

inline double
SplinedCurve::T_(double x, int i)
{
  return (x - x_[i]) / dx_[i];
}


inline double
SplinedCurve::Value(double x)
{
  int i = interval_(x);
  if (i == -1) {
    // first make sure it isn't the left end, which returns -1 thanks to lower_bound
    if (x == x_[0]) {
      i = 0;
    } else {
      // off the left end
      switch (left_) {
      case SplineExtrapolation_t::CONSTANT:
        return y_[0];
      case SplineExtrapolation_t::LINEAR:
        return y_[0] - (x - x_[0]) * dy_[0];
      case SplineExtrapolation_t::SPLINE:
        ++i;
        break;
      case SplineExtrapolation_t::THROW:
        Exceptions::amanzi_throw(Amanzi_exception_SplineOutOfRange());
      }
    }
  } else if (i == dx_.size()) {
    // off the right end
    switch (right_) {
    case SplineExtrapolation_t::CONSTANT:
      return y_[i];
    case SplineExtrapolation_t::LINEAR:
      return y_[i] + (x - x_[i]) * dy_[i];
    case SplineExtrapolation_t::SPLINE:
      --i;
      break;
    case SplineExtrapolation_t::THROW:
      Exceptions::amanzi_throw(Amanzi_exception_SplineOutOfRange());
    }
  }

  double t = T_(x, i);
  return (1 - t) * (1 - t) * ((1 + 2 * t) * y_[i] + t * dx_[i] * dy_[i]) +
         t * t * ((3 - 2 * t) * y_[i + 1] + (t - 1) * dx_[i] * dy_[i + 1]);
}

inline double
SplinedCurve::Derivative(double x)
{
  int i = interval_(x);
  if (i == -1) {
    // first make sure it isn't the left end, which returns -1 thanks to lower_bound
    if (x == x_[0]) {
      i = 0;
    } else {
      // off the left end
      switch (left_) {
      case SplineExtrapolation_t::CONSTANT:
        return 0.;
      case SplineExtrapolation_t::LINEAR:
        return dy_[0];
      case SplineExtrapolation_t::SPLINE:
        ++i;
        break;
      case SplineExtrapolation_t::THROW:
        Exceptions::amanzi_throw(Amanzi_exception_SplineOutOfRange());
      }
    }
  } else if (i == dx_.size()) {
    // off the right end
    switch (right_) {
    case SplineExtrapolation_t::CONSTANT:
      return 0.;
    case SplineExtrapolation_t::LINEAR:
      return dy_[i];
    case SplineExtrapolation_t::SPLINE:
      --i;
      break;
    case SplineExtrapolation_t::THROW:
      Exceptions::amanzi_throw(Amanzi_exception_SplineOutOfRange());
    }
  }

  double t = T_(x, i);
  double t_sq = t * t;
  double dydt = (6 * t_sq - 6 * t) * y_[i] + (3 * t_sq - 4 * t + 1) * dx_[i] * dy_[i] +
                (-6 * t_sq + 6 * t) * y_[i + 1] + (3 * t_sq - 2 * t) * dx_[i] * dy_[i + 1];
  return dydt / dx_[i];
}


} // namespace Utils
} // namespace Amanzi


#endif
