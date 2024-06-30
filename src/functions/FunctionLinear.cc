/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include "FunctionLinear.hh"
#include "errors.hh"

namespace Amanzi {

FunctionLinear::FunctionLinear(double y0,
                               const Kokkos::View<const double*, Kokkos::HostSpace>& grad,
                               const Kokkos::View<const double*, Kokkos::HostSpace>& x0)
{
  if (grad.extent(0) < 1) {
    Errors::Message m;
    m << "At least one value required for the gradient vector";
    Exceptions::amanzi_throw(m);
  }
  if (x0.extent(0) != grad.extent(0)) {
    Errors::Message m;
    m << "Mismatch of gradient and point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  y0_ = y0;
  x0_ = asDualView(x0);
  grad_ = asDualView(grad);
}

FunctionLinear::FunctionLinear(double y0, const Kokkos::View<const double*, Kokkos::HostSpace>& grad)
{
  if (grad.extent(0) < 1) {
    Errors::Message m;
    m << "at least one value required for the gradient vector";
    Exceptions::amanzi_throw(m);
  }
  y0_ = y0;

  std::vector<double> x(grad.extent(0), 0.0);
  x_ = asDualView(x);

  grad_ = asDualView(grad);
}

double
FunctionLinear::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionLinearFunctor(y0_, grad_.view_host(), x0_.view_host(), x);
  return f(0);
}


void
FunctionLinear::apply(const Kokkos::View<const double**>& in,
           Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionLinearFunctor(y0_, grad_.view_device(), x0_.view_device(), in);
  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionBilinear::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(ids_loc[i]);
      });
  } else {
    Kokkos::parallel_for(
      "FunctionBilinear::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}


} // namespace Amanzi
