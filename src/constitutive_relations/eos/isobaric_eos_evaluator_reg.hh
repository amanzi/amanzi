/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "isobaric_eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,IsobaricEOSEvaluator> IsobaricEOSEvaluator::factory_("isobaric eos");

} // namespace
} // namespace
