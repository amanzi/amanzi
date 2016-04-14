/*
  Litter drainage rate.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "interception_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,InterceptionEvaluator> InterceptionEvaluator::factory_("interception/throughfall");

} // namespace
} // namespace
} // namespace
