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

#include <string>

#include "Evaluator_Factory.hh"

namespace Amanzi {

namespace Utils {
template <>
Factory<Evaluator>::map_type* Factory<Evaluator>::map_;
} // namespace Utils

// method for instantiating Evaluator implementations
Teuchos::RCP<Evaluator>
Evaluator_Factory::createEvaluator(Teuchos::ParameterList& plist)
{
  std::string type = plist.get<std::string>("evaluator type");
  return Teuchos::rcp(CreateInstance(type, plist));
};

} // namespace Amanzi
