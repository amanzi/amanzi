/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_relation_factory.hh"
#include "viscosity_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");

} // namespace
} // namespace
