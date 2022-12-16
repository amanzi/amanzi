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

#include <algorithm>

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

      // dy = min( 3*min(S-1/2, S+1/2), dy)
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

      // dy = max( 3*max(S-1/2, S+1/2), dy)
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
