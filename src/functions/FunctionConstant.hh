/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! FunctionConstant: Implements the Function interface using a constant value.
/*!

Constant function is defined as :math:`f(x) = a`, for all :math:`x`.

* `"value`" ``[double]`` The constant to be applied.

Example:

.. code-block:: xml

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value="1.0"/>
  </ParameterList>

*/

#ifndef AMANZI_CONSTANT_FUNCTION_HH_
#define AMANZI_CONSTANT_FUNCTION_HH_

#include "Function.hh"

namespace Amanzi {

class FunctionConstant : public Function {
 public:
  FunctionConstant(double c) : c_(c) {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionConstant>(*this); }
  double operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const { return c_; }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out, const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
  {
    if (ids) {
      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "FunctionConstant::apply", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(ids_loc(i)) = c_; });
    } else {
      Kokkos::deep_copy(out, c_);
    }
  }

 private:
  double c_;
};

} // namespace Amanzi

#endif // AMANZI_CONSTANT_FUNCTION_HH_
