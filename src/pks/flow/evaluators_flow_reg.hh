/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Flow PK

  Self-registering factory for flow evaluators.
*/

#include "ModelEvaluator.hh"
#include "SpecificStorage.hh"
#include "WaterStorage.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<Evaluator, WaterStorage> WaterStorage::reg_("water storage");
template<>
Utils::RegisteredFactory<Evaluator, ModelEvaluator<SpecificStorage>> ModelEvaluator<SpecificStorage>::reg_("specific storage");

} // namespace Flow
} // namespace Amanzi
