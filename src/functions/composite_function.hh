/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Factory for vector functions which are composed of multiple scalar functions.

The expected plist is of the form:

<ParameterList name="constuctor plist">
  <Parameter name="Number of DoFs">

  <ParameterList name="Function 1">
    <ParameterList name="function-constant">
      ...
    </ParameterList>
  </ParameterList>

  <ParameterList name="Function 2">
    <ParameterList name="function-linear">
      ...
    </ParameterList>
  </ParameterList>

  ...
</ParameterList>


Where each of the "Function X" lists are valid input to the
function-factory Create() method (see ./function-factory.hh).

------------------------------------------------------------------------- */


#ifndef AMANZI_COMPOSITE_FUNCTION_HH_
#define AMANZI_COMPOSITE_FUNCTION_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "function.hh"
#include "function-factory.hh"
#include "vector_function.hh"
#include "vector_function_factory.hh"

namespace Amanzi {

class CompositeFunction : public VectorFunction {

public:
  CompositeFunction(const std::vector<Teuchos::RCP<const Function> >& functions);

  CompositeFunction(const Teuchos::RCP<const Function>& function);

  ~CompositeFunction();

  VectorFunction* Clone() const;

  int size() const;

  double* operator() (const double* xt) const;

private:
  std::vector<Teuchos::RCP<const Function> > functions_;
  double* values_;

  static RegisteredVectorFunctionFactory reg_;

};

// Non-member constructor/factory
Teuchos::RCP<VectorFunction> CreateCompositeFunction(Teuchos::ParameterList& plist);

} // namespace

#endif
