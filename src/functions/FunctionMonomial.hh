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

#ifndef AMANZI_MONOMIAL_FUNCTION_HH_
#define AMANZI_MONOMIAL_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionMonomial : public Function {
 public:
  FunctionMonomial(double c,
                   const Kokkos::View<double*, Kokkos::HostSpace>& x0,
                   const Kokkos::View<int*, Kokkos::HostSpace>& p);
  ~FunctionMonomial() {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionMonomial>(*this); }
  double operator()(const Kokkos::View<double*, Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double y = c_;
    if (x.extent(0) < x0_.extent(0)) {
      assert(false && "FunctionMonomial expects higher-dimensional argument.");
    }
    for (int j = 0; j < x0_.extent(0); ++j)
      y *= pow(x(j, i) - x0_.view_device()[j], p_.view_device()[j]);
    return y;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out, const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
  {
    if (ids) {
      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "FunctionBilinear::apply1", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(ids_loc(i)) = apply_gpu(in, i);
        });
    } else {
      assert(in.extent(1) == out.extent(0));
      Kokkos::parallel_for(
        "FunctionBilinear::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(i) = apply_gpu(in, i);
        });
    }
  }

 private:
  double c_;
  Kokkos::DualView<double*> x0_;
  Kokkos::DualView<int*> p_;
};

} // namespace Amanzi

#endif
