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
FunctionDistance::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionDistanceFunctor(x0_.view_host(), metric_.view_host(), squared_, x);
  return f(0);
}


void
FunctionDistance::apply(const Kokkos::View<const double**>& in,
           Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  AMANZI_ASSERT(in.extent(1) == out.extent(0));
  if (in.extent(0) < x0_.extent(0)) {
    Errors::Message m;
    m << "FunctionDistance expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }

  {
    auto f = Impl::FunctionDistanceFunctor(x0_.view_device(), metric_.view_device(), squared_, in);

    if (ids) {
      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "FunctionBilinear::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
          out(ids_loc(i)) = f(ids_loc(i));
        });

    } else {
      Kokkos::parallel_for(
        "FunctionBilinear::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(i) = f(i);
        });
    }
  }
}

} // namespace Amanzi
