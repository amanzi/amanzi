/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#include "wrm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> WRMEvaluator::factory_("WRM");

} //namespace
} //namespace
} //namespace
