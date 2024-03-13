/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Cubic spline

*/

#include <algorithm>
#include <iostream>

#include "SplinedCurve.hh"

namespace Amanzi {
namespace Utils {

SplinedCurve::SplinedCurve(std::vector<double> x,
                           std::vector<double> y,
                           std::vector<double> dydx,
                           SplineEndpoints_t endpoints,
                           bool enforce_monotonicity)
  : x_(x),
    y_(y),
    dy_(dydx),
    left_(endpoints.first),
    right_(endpoints.second),
    mono_(enforce_monotonicity)
{
  Setup_();
};


SplinedCurve::SplinedCurve(std::vector<double> x,
                           std::vector<double> y,
                           SplineEndpoints_t endpoints,
                           bool enforce_monotonicity)
  : x_(x), y_(y), left_(endpoints.first), right_(endpoints.second), mono_(enforce_monotonicity)
{
  // approximate derivatives
  int n = x.size();
  dy_.resize(n, 0.0);
  for (int i = 1; i < n - 1; ++i) {
    double dxl = x[i] - x[i - 1];
    double dxr = x[i + 1] - x[i];

    double gl = (y[i] - y[i - 1]) / dxl;
    double gr = (y[i + 1] - y[i]) / dxr;

    // weighting gives one additional order of gradient approximation
    if (gl * gr > 0.0) dy_[i] = (gl * dxr + gr * dxl) / (dxl + dxr);
  }

  Setup_();
};


void
SplinedCurve::Setup(std::vector<double> x,
                    std::vector<double> y,
                    std::vector<double> dydx,
                    SplineEndpoints_t endpoints,
                    bool enforce_monotonicity)
{
  x_ = x;
  y_ = y;
  dy_ = dydx;
  left_ = endpoints.first;
  right_ = endpoints.second;
  mono_ = enforce_monotonicity;
  Setup_();
}


double
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


double
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


void
SplinedCurve::Setup_()
{
  if (x_.size() < 2) {
    Errors::Message msg("SplinedCurve must get arrays of length 2 or greater.");
    Exceptions::amanzi_throw(msg);
  }

  if (x_.size() != y_.size() || x_.size() != dy_.size()) {
    Errors::Message msg("SplinedCurve must get arrays of equal length.");
    Exceptions::amanzi_throw(msg);
  }

  dx_.resize(x_.size() - 1, 0.);
  for (int i = 0; i != dx_.size(); ++i) {
    dx_[i] = x_[i + 1] - x_[i];
    if (dx_[i] <= 0.) {
      Errors::Message msg("SplinedCurve independent variable x must be monotonically increasing.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // enforce monotonicity using the approach of Hyman '83 Accurate
  // Mono. Pres. Cubic Interp.
  if (mono_) {
    // determine whether monotonically increasing, decreasing, or zero
    double mean_slope = (y_[y_.size() - 1] - y_[0]);

    // calculate the finite difference slopes
    std::vector<double> dy_fd(dx_.size(), 0.);
    for (int i = 0; i != dx_.size(); ++i) { dy_fd[i] = (y_[i + 1] - y_[i]) / dx_[i]; }

    // enforce sufficient monotonicity constraint on slopes
    if (mean_slope > 0.) {
      double min_slope = *std::min_element(dy_.begin(), dy_.end());
      if (min_slope < 0.) {
        Errors::Message msg("SplinedCurve requested monotonicity_preserving with postive mean "
                            "slope but slopes provided are not uniformly non-negative.");
        Exceptions::amanzi_throw(msg);
      }

      // dy = min(3*min(S-1/2, S+1/2), dy)
      for (int i = 1; i != dy_.size() - 1; ++i) {
        dy_[i] = std::min(3 * std::min(dy_fd[i - 1], dy_fd[i]), dy_[i]);
      }

    } else if (mean_slope < 0.) {
      double max_slope = *std::max_element(dy_.begin(), dy_.end());
      if (max_slope > 0.) {
        Errors::Message msg("SplinedCurve requested monotonicity_preserving with negative mean "
                            "slope but slopes provided are not uniformly non-positive.");
        Exceptions::amanzi_throw(msg);
      }

      // dy = max(3*max(S-1/2, S+1/2), dy)
      for (int i = 1; i != dy_.size() - 1; ++i) {
        dy_[i] = std::max(3 * std::max(dy_fd[i - 1], dy_fd[i]), dy_[i]);
      }

    } else {
      double min_slope = *std::min_element(dy_.begin(), dy_.end());
      double max_slope = *std::max_element(dy_.begin(), dy_.end());
      if (min_slope != 0.0 || max_slope != 0.0) {
        Errors::Message msg("SplinedCurve requested monotonicity_preserving with zero mean slope "
                            "but slopes provided are not uniformly zero.");
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}

} // namespace Utils
} // namespace Amanzi
