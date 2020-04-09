/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "FunctionSquareDistance.hh"
#include "errors.hh"
#include <cmath>

namespace Amanzi {

FunctionSquareDistance::FunctionSquareDistance(
  const Kokkos::View<double*>& x0, const Kokkos::View<double*>& metric)
{
  if (x0.extent(0) != metric.extent(0)) {
    Errors::Message m;
    m << "Mismatch of metric and point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  x0_ = x0;
  metric_ = metric;
}

double
FunctionSquareDistance::operator()(const Kokkos::View<double*>& x) const
{
  double tmp(0.), y(0.0);
  if (x.extent(0) < x0_.extent(0)) {
    Errors::Message m;
    m << "FunctionSquareDistance expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < x0_.extent(0); ++j) {
    tmp = x[j] - x0_[j];
    y += metric_[j] * tmp * tmp;
  }
  return y;
}

} // namespace Amanzi
