/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! FunctionMonomial: a multivariate monomial function.
/*!

A multi-variable monomial function is given by the following expression:

.. math::
  f(x) = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

with the constant factor :math:`c`, the reference point :math:`x_0`, and
integer exponents :math:`p_j`.
Note that the first parameter in :math:`x` can be time.

.. _function-monomial-spec:
.. admonition:: function-monomial-spec

   * `"c`" ``[double]`` c in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`
   * `"x0`" ``[Array(double)]`` x0 in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`
   * `"exponents`" ``[Array(int)]`` p in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`

Conditions:

.. code-block:: python

  len(x0) == len(exponents)

Here is an example of monomial of degree 6 in three variables:

.. code-block:: xml

  <ParameterList name="function-monomial">
    <Parameter name="c" type="double" value="1.0"/>
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 3, 1}"/>
  </ParameterList>

*/

#pragma once

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionMonomial : public Function {
 public:
  FunctionMonomial(double c,
                   const Kokkos::View<const double*, Kokkos::HostSpace>& x0,
                   const Kokkos::View<const int*, Kokkos::HostSpace>& p);
  ~FunctionMonomial() {}
  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionMonomial>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& in) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  double c_;
  Kokkos::DualView<const double*> x0_;
  Kokkos::DualView<const int*> p_;
};


namespace Impl {

template<typename View_type, typename IntView_type, typename InView_type>
class FunctionMonomialFunctor {
 public:
  FunctionMonomialFunctor(double c,
                          const View_type& x0,
                          const IntView_type& p,
                          const InView_type& in)
    : c_(c), x0_(x0), p_(p), in_(in) {}

  KOKKOS_INLINE_FUNCTION
  double
  operator()(const int i) const
  {
    double y = c_;
    for (int j = 0; j < x0_.extent(0); ++j)
      y *= Kokkos::pow(in_(j, i) - x0_[j], p_[j]);
    return y;
  }

 private:
  View_type x0_;
  IntView_type p_;
  double c_;

  InView_type in_;
};

} // namespace Impl
} // namespace Amanzi

