/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_fraction_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,UnfrozenFractionEvaluator>
UnfrozenFractionEvaluator::fac_("unfrozen fraction");

} //namespace
} //namespace
} //namespace

