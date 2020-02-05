/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! FunctionDistance: distance from a reference point.

/*!

A distance function calculates distance from reference point :math:`x_0`
using by the following expression:

.. math::
  f(x) = \sqrt( \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2 )

Note that the first parameter in :math:`x` can be time.

* `"x0`" ``[Array(double)]`` Point from which distance is measured.
* `"metric`" ``[Array(double)]`` Linear scaling metric, typically all 1s.

Here is an example of a distance function using isotropic metric:

Example:
.. code-block:: xml

  <ParameterList name="function-distance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
  </ParameterList>

*/

#ifndef AMANZI_DISTANCE_FUNCTION_HH_
#define AMANZI_DISTANCE_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionDistance : public Function {
 public:
  FunctionDistance(const Kokkos::View<double*>& x0,
                   const Kokkos::View<double*>& metric);
  ~FunctionDistance() {}
  FunctionDistance* Clone() const { return new FunctionDistance(*this); }
  double operator()(const Kokkos::View<double*>&) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double tmp(0.0), y(0.0);
    if (x.extent(0) < x0_.extent(0)) {
      assert(false && "FunctionDistance expects higher-dimension argument.");
      // Errors::Message m;
      // m << "FunctionDistance expects higher-dimensional argument.";
      // Exceptions::amanzi_throw(m);
    }
    for (int j = 0; j < x0_.extent(0); ++j) {
      tmp = x(j, i) - x0_[j];
      y += metric_[j] * tmp * tmp;
    }
    y = sqrt(y);
    return y;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    assert(in.extent(1) == out.extent(0));
    Kokkos::parallel_for(
      "FunctionDistance::apply",
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }

 private:
  Kokkos::View<double*> x0_, metric_;
};

} // namespace Amanzi

#endif
