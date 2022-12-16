/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/* -------------------------------------------------------------------------
ATS & Amanzi

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
