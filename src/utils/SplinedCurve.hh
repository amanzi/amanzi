/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Spline
/*!

Spline fit of a curve, given points, values, and derivatives (optional) of those
values. Missing derivatives are constructed using a simple algorithm.
The class orovides options for monotonicity preservation.

*/


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

  SplinedCurve(std::vector<double> x,
               std::vector<double> y,
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

} // namespace Utils
} // namespace Amanzi


#endif
