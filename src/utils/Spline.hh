/* -------------------------------------------------------------------------
  Spline

  Author: Ethan Coon ecoon@lanl.gov

  Simple spline functor which is used for smoothing.

------------------------------------------------------------------------- */


#ifndef AMANZI_UTILS_SPLINE_HH_
#define AMANZI_UTILS_SPLINE_HH_


#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Utils {

class Spline {

 public:
  Spline(double x1, double y1, double dy1,
         double x2, double y2, double dy2);

  double operator()(double x) { return Value(x); }

  double Value(double x);
  double Derivative(double x);

 private:
  double T(double x);

  double x1_, x2_, y1_, y2_, dy1_, dy2_;

};


inline
double
Spline::T(double x) {
  return (x - x1_) / (x2_ - x1_);
}

} // namespace Utils
} // namespace Amanzi


#endif
