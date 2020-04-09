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
  f(x) = \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2

Note that the first parameter in :math:`x` can be time.
Here is an example of a distance function using isotropic metric:

Example:
.. code-block:: xml

  <ParameterList name="function-squaredistance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
  </ParameterList>

*/

#ifndef AMANZI_SQUAREDISTANCE_FUNCTION_HH_
#define AMANZI_SQUAREDISTANCE_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionSquareDistance : public Function {
 public:
  FunctionSquareDistance(const Kokkos::View<double*>& x0,
                         const Kokkos::View<double*>& metric);
  ~FunctionSquareDistance() {}
  FunctionSquareDistance* Clone() const
  {
    return new FunctionSquareDistance(*this);
  }
  double operator()(const Kokkos::View<double*>& x) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double tmp(0.), y(0.0);
    if (x.extent(0) < x0_.extent(0)) {
      assert(false &&
             "FunctionSquareDistance expects higher-dimensional argument.");
    }
    for (int j = 0; j < x0_.extent(0); ++j) {
      tmp = x(j, i) - x0_[j];
      y += metric_[j] * tmp * tmp;
    }
    return y;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::parallel_for(
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }


 private:
  Kokkos::View<double*> x0_, metric_;
};

} // namespace Amanzi

#endif
