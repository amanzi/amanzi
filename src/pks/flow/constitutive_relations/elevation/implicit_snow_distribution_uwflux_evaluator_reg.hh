/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  SnowDistributionUWFluxEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "implicit_snow_distribution_uwflux_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ImplicitSnowDistributionUWFluxEvaluator> ImplicitSnowDistributionUWFluxEvaluator::factory_("implicit snow distribution, upwind flux");

} // namespace
} // namespace
} // namespace
