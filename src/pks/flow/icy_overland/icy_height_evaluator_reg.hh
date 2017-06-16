#include "icy_height_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,IcyHeightEvaluator> IcyHeightEvaluator::factory_("icy ponded depth");

} //namespace
} //namespace
