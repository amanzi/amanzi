/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! FunctionLinear: a multivariate linear function.
/*!

A multi-variable linear function is formally defined by

.. math::
  f(x) = y_0 + \sum_{{j=0}}^{{n-1}} g_j (x_j - x_{{0,j}})

with the constant term "math:`y_0` and gradient :math:`g_0,\, g_1\,...,
g_{{n-1}}`. If the reference point :math:`x_0` is specified, it must have the
same number of values as the gradient.  Otherwise, it defaults to zero. Note
that one of the parameters in a multi-valued linear function can be time.

.. _function-linear-spec:
.. admonition:: function-linear-spec

   * `"y0`" ``[double]`` y_0 in f = y0 + g * (x - x0)
   * `"gradient`" ``[Array(double)]`` g in f = y0 + g * (x - x0)
   * `"x0`" ``[Array(double)]`` x0 in f = y0 + g * (x - x0)

Conditions:

.. code-block:: python

  len(x0) == len(gradient)


Example:

.. code-block:: xml

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value="1.0"/>
    <Parameter name="gradient" type="Array(double)" value="{1.0, 2.0, 3.0}"/>
    <Parameter name="x0" type="Array(double)" value="{2.0, 3.0, 1.0}"/>
  </ParameterList>

*/

#ifndef AMANZI_LINEAR_FUNCTION_HH_
#define AMANZI_LINEAR_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionLinear : public Function {
 public:
  FunctionLinear(double y0, const Kokkos::View<double*, Kokkos::HostSpace>& grad);
  FunctionLinear(double y0,
                 const Kokkos::View<double*, Kokkos::HostSpace>& grad,
                 const Kokkos::View<double*, Kokkos::HostSpace>& x0);
  ~FunctionLinear() {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionLinear>(*this); }
  double operator()(const Kokkos::View<double*, Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double y = y0_;
    if (x.extent(0) < grad_.extent(0)) {
      assert(false && "FunctionLinear expects higher-dimensional argument.");
    }
    for (int j = 0; j < grad_.extent(0); ++j)
      y += grad_.view_device()[j] * (x(j, i) - x0_.view_device()[j]);
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
  double y0_;
  Kokkos::DualView<double*> grad_, x0_;
};

} // namespace Amanzi

#endif
