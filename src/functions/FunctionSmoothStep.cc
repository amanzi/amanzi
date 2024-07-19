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
FunctionSmoothStep::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionSmoothStepFunctor(x0_, y0_, x1_, y1_, x);
  return f(0);
}


void
FunctionSmoothStep::apply(const Kokkos::View<const double**>& in,
           Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionSmoothStepFunctor(x0_, y0_, x1_, y1_, in);

  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionSmoothStep::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(i);
      });
  } else {
    Kokkos::parallel_for(
      "FunctionSmoothStep::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}


} // namespace Amanzi
