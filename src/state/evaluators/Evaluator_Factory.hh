/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Self-registering factory for FIELD_EVALUATOR implementations.
*/

#ifndef AMANZI_STATE_EVALUATOR_FACTORY_HH_
#define AMANZI_STATE_EVALUATOR_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Evaluator.hh"
#include "Factory.hh"

namespace Amanzi {

class Evaluator_Factory : public Utils::Factory<Evaluator> {
 public:
  Teuchos::RCP<Evaluator> createEvaluator(Teuchos::ParameterList& plist);
};

} // namespace Amanzi

#endif
