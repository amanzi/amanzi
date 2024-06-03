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


#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Utils {

class Spline {
 public:
  Spline(){};

  Spline(double x1, double y1, double dy1, double x2, double y2, double dy2)
  {
    Setup(x1, y1, dy1, x2, y2, dy2);
  }

  void Setup(double x1, double y1, double dy1, double x2, double y2, double dy2);

  double operator()(double x) const { return Value(x); }

  double Value(double x) const;
  double Derivative(double x) const;

 private:
  double T(double x) const;

  double x1_, x2_, y1_, y2_, dy1_, dy2_;
};


inline double
Spline::T(double x) const
{
  return (x - x1_) / (x2_ - x1_);
}

} // namespace Utils
} // namespace Amanzi

#endif
