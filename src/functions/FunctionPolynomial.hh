/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! FunctionPolynomial: a polynomial
/*!

A generic polynomial function is given by the following expression:

.. math::
  f(x) = \sum_{{j=0}}^n c_j (x - x_0)^{{p_j}}

where :math:`c_j` are coefficients of monomials,
:math:`p_j` are integer exponents, and :math:`x_0` is the reference point.

.. _function-polynomial-spec:
.. admonition:: function-polynomial-spec

   * `"coefficients`" ``[Array(double)]`` c_j polynomial coefficients
   * `"exponents`" ``[Array(int)]`` p_j polynomail exponents
   * `"reference point`" ``[double]`` x0 to which polynomial argument is normalized.

Example:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{1.0, 1.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 4}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>

*/

#pragma once

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionPolynomial : public Function {
 public:
  FunctionPolynomial(const Kokkos::View<const double*, Kokkos::HostSpace>& c,
                     const Kokkos::View<const int*, Kokkos::HostSpace>& p,
                     double x0 = 0.0);

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionPolynomial>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>&) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  int pmin_;
  int pmax_;
  double x0_;
  Kokkos::DualView<double*> c_;
};


namespace Impl {

template <class DoubleView_type,
          class InView_type>
class FunctionPolynomialFunctor {
 public:
  FunctionPolynomialFunctor(const DoubleView_type& c,
                            int pmin,
                            int pmax,
                            double x0,
                            const InView_type& in)
    : c_(c), pmin_(pmin), pmax_(pmax), x0_(x0), in_(in) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const int i) const
  {
    // Polynomial terms with non-negative exponents
    double y = c_[pmax_ - pmin_];
    if (pmax_ > 0) {
      double z = in_(0, i) - x0_;
      for (int j = pmax_; j > 0; --j) y = c_[j - 1 - pmin_] + z * y;
    }

    // Polynomial terms with negative exponents.
    if (pmin_ < 0) {
      double w = c_[0];
      double z = 1.0 / (in_(0, i) - x0_);
      for (int j = pmin_; j < -1; ++j) w = c_[j + 1 - pmin_] + z * w;
      y += z * w;
    }
    return y;
  }

 private:
  DoubleView_type c_;
  double x0_;
  int pmax_, pmin_;
  InView_type in_;
};


} // namespace Impl
} // namespace Amanzi


