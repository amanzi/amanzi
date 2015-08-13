/*
  Litter drainage rate.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "litter_drainage_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LitterDrainageEvaluator> LitterDrainageEvaluator::factory_("litter drainage");

} // namespace
} // namespace
} // namespace
