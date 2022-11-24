/* -------------------------------------------------------------------------
ATS & Amanzi

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Function from R^d to R^n.
------------------------------------------------------------------------- */

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
  double* operator()(const std::vector<double>& xt) const;

 private:
  std::vector<Teuchos::RCP<const Function>> functions_;
  double* values_;
};

} // namespace Amanzi

#endif
