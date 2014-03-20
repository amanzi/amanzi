/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  SnowDistributionEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "explicit_snow_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ExplicitSnowDistributionEvaluator> ExplicitSnowDistributionEvaluator::factory_("explicit snow distribution");

} // namespace
} // namespace
} // namespace
