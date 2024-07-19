/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Functions

*/

/*!

This function parses a string expression. The function has min(N, D + 1)
arguments t, x, y, and z. The argument t is required. D is the space dimension,
and N is the user specified number of arguments which could be less than D + 1.

Example:

.. code-block:: xml

  <ParameterList name="function-exprtk">
    <Parameter name="number of arguments" type="int" value="3"/>
    <Parameter name="formula" type="string" value="t + x + 2 * y"/>
  </ParameterList>
*/

#ifndef AMANZI_EXPRTK_FUNCTION_HH_
#define AMANZI_EXPRTK_FUNCTION_HH_

#include <memory>
#include <vector>

#include "ExprTK.hh"

#include "Function.hh"

namespace Amanzi {

class FunctionExprTK : public Function {
 public:
  FunctionExprTK(int n, const std::string& formula);

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionExprTK>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override
  {
    // NOTE ExprTK cannot be used on device!
    Kokkos::View<double**, Kokkos::HostSpace> in_host(
      "ExpTK work space in", in.extent(0), in.extent(1));
    Kokkos::deep_copy(in_host, in);
    Kokkos::View<double*, Kokkos::HostSpace> out_host("ExpTK work space out", in.extent(1));
    for (int i = 0; i != out_host.extent(0); ++i) {
      out_host(i) = operator()(Kokkos::subview(in_host, Kokkos::ALL, Kokkos::pair{i,i+1}));
    }

    if (ids) {
      Kokkos::View<double*> out_dev("ExpTK work space dev", in.extent(1));
      Kokkos::deep_copy(out_dev, out_host);

      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "ExpTK copy to out", in.extent(1), KOKKOS_CLASS_LAMBDA(const int& i) {
          out(ids_loc(i)) = out_dev(i);
        });
    } else {
      Kokkos::deep_copy(out, out_host);
    }
  }

 private:
  std::shared_ptr<Utils::ExprTK> exprtk_;
};

} // namespace Amanzi

#endif
