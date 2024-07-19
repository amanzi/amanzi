/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>
#include "errors.hh"
#include "ViewUtils.hh"
#include "FunctionMonomial.hh"

namespace Amanzi {

FunctionMonomial::FunctionMonomial(double c,
                                   const Kokkos::View<const double*, Kokkos::HostSpace>& x0,
                                   const Kokkos::View<const int*, Kokkos::HostSpace>& p)
{
  if (x0.extent(0) != p.extent(0)) {
    Errors::Message m;
    m << "Mismatch of multi-index and reference point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  c_ = c;
  x0_ = asDualView(x0);
  p_ = asDualView(p);
}


double
FunctionMonomial::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& in) const {
  auto f = Impl::FunctionMonomialFunctor(c_, x0_.view_host(), p_.view_host(), in);
  return f(0);
}


void
FunctionMonomial::apply(const Kokkos::View<const double**>& in,
           Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  if (in.extent(0) < x0_.extent(0)) {
    Errors::Message msg("FunctionMonomial expects higher-dimensional argument.");
    Exceptions::amanzi_throw(msg);
  }
  auto f = Impl::FunctionMonomialFunctor(c_, x0_.view_device(), p_.view_device(), in);

  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionBilinear::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(i);
      });
  } else {
    Kokkos::parallel_for(
      "FunctionBilinear::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}


} // namespace Amanzi
