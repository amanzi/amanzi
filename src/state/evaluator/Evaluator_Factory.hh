/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for FIELD_EVALUATOR implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_EVALUATOR_FACTORY_HH_
#define AMANZI_STATE_EVALUATOR_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Evaluator.hh"
#include "factory.hh"

namespace Amanzi {

class Evaluator_Factory : public Utils::Factory<Evaluator> {

public:
  Teuchos::RCP<Evaluator> createEvaluator(Teuchos::ParameterList &plist);
};

} // namespace Amanzi

#endif
