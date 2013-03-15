/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for FIELD_EVALUATOR implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_FIELD_EVALUATOR_FACTORY_HH_
#define _PK_FLOW_FIELD_EVALUATOR_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "FieldEvaluator.hh"
#include "factory.hh"

namespace Amanzi {

class FieldEvaluator_Factory : public Utils::Factory<FieldEvaluator> {

public:
  Teuchos::RCP<FieldEvaluator> createFieldEvaluator(Teuchos::ParameterList& plist);
};

} // namespace

#endif
