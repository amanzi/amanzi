/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//! FunctionSmoothStep: a smoothed discontinuity.

/*!

A smooth :math:`C^2` function `f(x)` on interval :math:`[x_0,\,x_1]` is
defined such that `f(x) = y_0` for `x < x0`, `f(x) = y_1` for `x > x_1`, and
monotonically increasing for :math:`x \in [x_0, x_1]` through cubic
interpolation.

* `"x0`" ``[double]`` First fitting point
* `"y0`" ``[double]`` First fitting value
* `"x1`" ``[double]`` Second fitting point
* `"y1`" ``[double]`` Second fitting value

Example:

.. code-block:: xml

  <ParameterList name="function-smooth-step">
    <Parameter name="x0" type="double" value="0.0"/>
    <Parameter name="y0" type="double" value="0.0"/>
    <Parameter name="x1" type="double" value="1.0"/>
    <Parameter name="y1" type="double" value="2.0"/>
  </ParameterList>

*/

#ifndef AMANZI_SMOOTH_STEP_FUNCTION_HH_
#define AMANZI_SMOOTH_STEP_FUNCTION_HH_

#include "Function.hh"

namespace Amanzi {

class FunctionSmoothStep : public Function {
 public:
  FunctionSmoothStep(double x0, double y0, double x1, double y1);
  ~FunctionSmoothStep(){};
  FunctionSmoothStep* Clone() const { return new FunctionSmoothStep(*this); }
  double operator()(const Kokkos::View<double*>&) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double y;
    if (x(0, i) <= x0_) {
      y = y0_;
    } else if (x(0, i) >= x1_) {
      y = y1_;
    } else {
      double s = (x(0, i) - x0_) / (x1_ - x0_);
      y = y0_ + (y1_ - y0_) * s * s * (3 - 2 * s);
    }
    return y;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::parallel_for(
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }

 private:
  double x0_, y0_, x1_, y1_;
};

} // namespace Amanzi

#endif // AMANZI_SMOOTH_STEP_FUNCTION_HH_
