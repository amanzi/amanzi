/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for FIELD_EVALUATOR implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

// explicity instantitate the static data of Factory<FieldEvaluator>
template<> Utils::Factory<FieldEvaluator>::map_type* Utils::Factory<FieldEvaluator>::map_;

// method for instantiating FieldEvaluator implementations
Teuchos::RCP<FieldEvaluator> FieldEvaluator_Factory::createFieldEvaluator(Teuchos::ParameterList& plist) {
  std::string field_evaluator_typename = plist.get<std::string>("field evaluator type");
  return Teuchos::rcp(CreateInstance(field_evaluator_typename, plist));
};

} // namespace

