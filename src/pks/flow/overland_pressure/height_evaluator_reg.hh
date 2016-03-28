#include "height_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,HeightEvaluator> HeightEvaluator::factory_("ponded depth");

} //namespace
} //namespace
} //namespace
