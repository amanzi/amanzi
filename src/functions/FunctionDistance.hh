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

#pragma once

#include "dbc.hh"
#include "Function.hh"

namespace Amanzi {

class FunctionDistance : public Function {
 public:
  FunctionDistance(const Kokkos::View<const double*, Kokkos::HostSpace>& x0,
                   const Kokkos::View<const double*, Kokkos::HostSpace>& metric,
                   bool squared);

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionDistance>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>&) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  Kokkos::DualView<const double*> x0_, metric_;
  bool squared_;
};


namespace Impl {

template <class ParView_type,
          class InView_type>
class FunctionDistanceFunctor {
 public:
  FunctionDistanceFunctor(const ParView_type& x0,
                          const ParView_type& metric,
                          bool squared,
                          const InView_type& in)
    : x0_(x0), metric_(metric), in_(in), squared_(squared_)  {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const int i) const {
    double y(0.);
    for (int j = 0; j < x0_.extent(0); ++j) {
      double tmp = in_(j, i) - x0_[j];
      y += metric_[j] * tmp * tmp;
    }
    return squared_ ? y : Kokkos::sqrt(y);
  }

 private:
  ParView_type x0_, metric_;
  bool squared_;
  InView_type in_;

};

} // namespace Impl
} // namespace Amanzi


