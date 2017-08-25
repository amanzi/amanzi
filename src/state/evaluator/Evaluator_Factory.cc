/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for FIELD_EVALUATOR implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "Evaluator_Factory.hh"

namespace Amanzi {

namespace Utils {
// explicity instantitate the static data of Factory<Evaluator>
template<> Factory<Evaluator>::map_type* Factory<Evaluator>::map_;
} // namespace Utils

// method for instantiating Evaluator implementations
Teuchos::RCP<Evaluator> Evaluator_Factory::createEvaluator(Teuchos::ParameterList& plist) {
  std::string field_evaluator_typename = plist.get<std::string>("field evaluator type");
  return Teuchos::rcp(CreateInstance(field_evaluator_typename, plist));
};

} // namespace

