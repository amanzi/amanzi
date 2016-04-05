/*
  Drainage rate.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "drainage_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,DrainageEvaluator> DrainageEvaluator::factory_("interception drainage");

} // namespace
} // namespace
} // namespace
