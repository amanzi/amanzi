/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  SnowDistributionEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "snow_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SnowDistributionEvaluator> SnowDistributionEvaluator::factory_("snow distribution");

} // namespace
} // namespace
} // namespace
