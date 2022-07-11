/*
  Copyright 2010-201x held jointly by participating institutions.
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

* `"coefficients`" ``[Array(double)]`` c_j polynomial coefficients
* `"exponents`" ``[Array(int)]`` p_j polynomail exponents
* `"reference point`" ``[double]`` x0 to which polynomial argument is
normalized.

Example:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{1.0, 1.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 4}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>

*/

#ifndef AMANZI_POLYNOMIAL_FUNCTION_HH_
#define AMANZI_POLYNOMIAL_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionPolynomial : public Function {
 public:
  FunctionPolynomial(const Kokkos::View<double*,Kokkos::HostSpace>& c,
                     const Kokkos::View<int*,Kokkos::HostSpace>& p, double x0 = 0.0);
  ~FunctionPolynomial() {}
  FunctionPolynomial* Clone() const { return new FunctionPolynomial(*this); }
  double operator()(const Kokkos::View<double*,Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    auto vc = c_.view_device(); 
    // Polynomial terms with non-negative exponents
    double y = vc[pmax_ - pmin_];
    if (pmax_ > 0) {
      double z = x(0, i) - x0_;
      for (int j = pmax_; j > 0; --j) y = vc[j - 1 - pmin_] + z * y;
    }
    // Polynomial terms with negative exponents.
    if (pmin_ < 0) {
      double w = vc[0];
      double z = 1.0 / (x(0, i) - x0_);
      for (int j = pmin_; j < -1; ++j) w = vc[j + 1 - pmin_] + z * w;
      y += z * w;
    }
    return y;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::parallel_for(
      "FunctionPolynomial::apply",
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }


 private:
  int pmin_;
  int pmax_;
  double x0_;
  Kokkos::DualView<double*> c_;
};

} // namespace Amanzi

#endif // AMANZI_POLYNOMIAL_FUNCTION_HH_
