/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator evaluates saturation of gas, liquid, and ice from
  capillary pressures for the ice-liquid and liquid-gas pairs.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_permafrost_evaluator.hh"
#include "wrm_partition.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMPermafrostEvaluator> WRMPermafrostEvaluator::factory_("permafrost WRM");

} // namespace
} // namespace
} // namespace



