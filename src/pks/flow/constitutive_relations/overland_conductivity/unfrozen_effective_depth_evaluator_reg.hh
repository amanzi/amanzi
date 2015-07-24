/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen effective depth.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,UnfrozenEffectiveDepthEvaluator> UnfrozenEffectiveDepthEvaluator::fac_("unfrozen effective depth");

} //namespace
} //namespace
} //namespace

