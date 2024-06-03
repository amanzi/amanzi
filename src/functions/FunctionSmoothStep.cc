/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include "FunctionSmoothStep.hh"
#include "errors.hh"

namespace Amanzi {

FunctionSmoothStep::FunctionSmoothStep(double x0, double y0, double x1, double y1)
{
  x0_ = x0;
  y0_ = y0;
  x1_ = x1;
  y1_ = y1;
  if (x0 >= x1) {
    Errors::Message m;
    m << "require x0 < x1";
    Exceptions::amanzi_throw(m);
  }
}

double
FunctionSmoothStep::operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const
{
  double y;
  if (x[0] <= x0_) {
    y = y0_;
  } else if (x[0] >= x1_) {
    y = y1_;
  } else {
    double s = (x[0] - x0_) / (x1_ - x0_);
    y = y0_ + (y1_ - y0_) * s * s * (3 - 2 * s);
  }
  return y;
}

} // namespace Amanzi
