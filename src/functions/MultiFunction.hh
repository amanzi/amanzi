/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Function from R^d to R^n.

/*!

  A MultiFunction is simply an array of functions, which allow Functions to
  be used for MultiVectors.

  Factory for vector functions which are composed of multiple scalar functions.
  The expected plist is of the form:

  <ParameterList name="constuctor plist">
    <Parameter name="number of dofs">
  <ParameterList name="dof 1 function">
      <ParameterList name="function-constant">
        ...
      </ParameterList>
    </ParameterList>

  <ParameterList name="dof 2 function">
      <ParameterList name="function-linear">
        ...
      </ParameterList>
    </ParameterList>

    ...
  </ParameterList>

  Where each of the "Function X" lists are valid input to the
  function-factory Create() method (see ./function-factory.hh).
*/


#ifndef AMANZI_MULTIVECTOR_FUNCTION_HH_
#define AMANZI_MULTIVECTOR_FUNCTION_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

class MultiFunction {
 public:
  MultiFunction(const std::vector<Teuchos::RCP<const Function>>& functions);
  MultiFunction(const Teuchos::RCP<const Function>& function);
  MultiFunction(Teuchos::ParameterList& plist);

  ~MultiFunction();

  int size() const;
  Kokkos::View<double*> operator()(const Kokkos::View<double*>& xt) const;

  //
  // NOTE: this requirement of the out to be LayoutLeft is because of the
  // expectation that out_i is NOT LayoutStride.  In an ideal world, out_i
  // COULD be layout stride, in which case the single-function apply methods
  // would have to be templated on view type (i.e. could take either
  // LayoutStride or not (contiguous in memory).  But single-function apply
  // must be virtual, and in C++11 at least, we can't have both worlds.
  //
  // So for now, we require out to be LayoutLeft to enforce that out_i is
  // contiguous.  Likely this is important for performance anyway, so I doubt
  // we're losing much generality, and may even be making performance more
  // robust.
  void apply(const Kokkos::View<double**>& in,
             Kokkos::View<double**, Kokkos::LayoutLeft>& out) const
  {
    for (int i = 0; i < size(); ++i) {
      Kokkos::View<double*> out_i = Kokkos::subview(out, Kokkos::ALL, i);
      functions_[i]->apply(in, out_i);
    }
  }

 private:
  std::vector<Teuchos::RCP<const Function>> functions_;
  Kokkos::View<double*> values_;
};

} // namespace Amanzi

#endif
