/*
  Copyright 2010-201x held jointly by participating institutions.
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

FunctionPolynomial::FunctionPolynomial(const Kokkos::View<double*,Kokkos::HostSpace>& c,
                                       const Kokkos::View<int*,Kokkos::HostSpace>& p, double x0)
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

  // Minimum and maximum powers.
  // pmin_ = std::min(0, *(std::min_element(p.begin(), p.end())));
  // pmax_ = std::max(0, *(std::max_element(p.begin(), p.end())));

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
  Kokkos::resize(c_, n);
  for (int i = 0; i < n; ++i) { c_.view_host()(i) = 0.0; }
  // c_.resize(n);
  // c_.assign(n, 0.0);
  for (int j = 0; j < c.extent(0); ++j) c_.view_host()[p[j] - pmin_] += c[j];
  x0_ = x0;
  Kokkos::deep_copy(c_.view_device(),c_.view_host());
}

double
FunctionPolynomial::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& x) const
{
  // Polynomial terms with non-negative exponents
  double y = c_.view_host()[pmax_ - pmin_];
  if (pmax_ > 0) {
    double z = x[0] - x0_;
    for (int j = pmax_; j > 0; --j) y = c_.view_host()[j - 1 - pmin_] + z * y;
  }
  // Polynomial terms with negative exponents.
  if (pmin_ < 0) {
    double w = c_.view_host()[0];
    double z = 1.0 / (x[0] - x0_);
    for (int j = pmin_; j < -1; ++j) w = c_.view_host()[j + 1 - pmin_] + z * w;
    y += z * w;
  }
  return y;
}
} // namespace Amanzi
