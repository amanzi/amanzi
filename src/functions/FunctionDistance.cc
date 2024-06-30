/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>
#include "errors.hh"
#include "FunctionDistance.hh"
#include "ViewUtils.hh"

namespace Amanzi {

FunctionDistance::FunctionDistance(const Kokkos::View<const double*, Kokkos::HostSpace>& x0,
        const Kokkos::View<const double*, Kokkos::HostSpace>& metric,
        bool squared)
  : squared_(squared)
{
  if (x0.extent(0) != metric.extent(0)) {
    Errors::Message m;
    m << "Mismatch of metric and point dimensions.";
    Exceptions::amanzi_throw(m);
  }

  x0_ = asDualView(x0);
  metric_ = asDualView(metric);
}


double
FunctionDistance::operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const
{
  double tmp(0.0), y(0.0);
  if (x.extent(0) < x0_.extent(0)) {
    assert(false && "FunctionDistance expects higher-dimension argument.");
    // Errors::Message m;
    // m << "FunctionDistance expects higher-dimensional argument.";
    // Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < x0_.extent(0); ++j) {
    tmp = x[j] - x0_.view_host()[j];
    y += metric_.view_host()[j] * tmp * tmp;
  }
  y = sqrt(y);
  return y;
}

} // namespace Amanzi
