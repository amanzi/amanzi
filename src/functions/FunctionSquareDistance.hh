/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
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
  FunctionSquareDistance(const Kokkos::View<double*, Kokkos::HostSpace>& x0,
                         const Kokkos::View<double*, Kokkos::HostSpace>& metric);
  ~FunctionSquareDistance() {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionSquareDistance>(*this); }
  double operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const;

  KOKKOS_INLINE_FUNCTION double apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double tmp(0.), y(0.0);
    if (x.extent(0) < x0_.extent(0)) {
      assert(false && "FunctionSquareDistance expects higher-dimensional argument.");
    }
    for (int j = 0; j < x0_.extent(0); ++j) {
      tmp = x(j, i) - x0_.view_device()[j];
      y += metric_.view_device()[j] * tmp * tmp;
    }
    return y;
  }


  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out, const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
  {
    if (ids) {
      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "FunctionSquareDistance::apply1", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(ids_loc(i)) = apply_gpu(in, i);
        });
    } else {
      assert(in.extent(1) == out.extent(0));
      Kokkos::parallel_for(
        "FunctionSquareDistance::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(i) = apply_gpu(in, i);
        });
    }
  }


 private:
  Kokkos::DualView<double*> x0_, metric_;
};

} // namespace Amanzi

#endif
