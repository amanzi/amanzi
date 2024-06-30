/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include <algorithm>

#include "FunctionPolynomial.hh"
#include "errors.hh"

namespace Amanzi {

FunctionPolynomial::FunctionPolynomial(const Kokkos::View<const double*, Kokkos::HostSpace>& c,
                                       const Kokkos::View<const int*, Kokkos::HostSpace>& p,
                                       double x0)
{
  if (c.extent(0) < 1) {
    Errors::Message m;
    m << "at least one value requred for the coefficient vector";
    Exceptions::amanzi_throw(m);
  }
  if (p.extent(0) != c.extent(0)) {
    Errors::Message m;
    m << "the number of values for the coefficient and exponent vectors differ";
    Exceptions::amanzi_throw(m);
  }

  x0_ = x0;

  // Find min and max
  pmin_ = p(0);
  pmax_ = p(0);
  for (int i = 1; i < p.extent(0); ++i) {
    pmin_ = std::min(pmin_, p(i));
    pmax_ = std::max(pmax_, p(i));
  }
  pmin_ = std::min(0, pmin_);
  pmax_ = std::max(0, pmax_);

  int n = pmax_ - pmin_ + 1;
  std::vector<double> c(n, 0.);
  for (int j = 0; j < c.size(); ++j) c_[p[j] - pmin_] += c[j];
  c_ = asDualView(c);
}

double
FunctionPolynomial::operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionPolynomialFunctor(c_.view_host(), pmin_, pmax_, x0_, x);
  return f(0);
}


void
FunctionPolynomial::apply(const Kokkos::View<const double**>& in,
                          Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionPolynomialFunctor(c_.view_device(), pmin_, pmax_, x0_, in);

  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionPolynomial::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(ids_loc(i));
      });
  } else {
    Kokkos::parallel_for(
      "FunctionPolynomial::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}

} // namespace Amanzi
