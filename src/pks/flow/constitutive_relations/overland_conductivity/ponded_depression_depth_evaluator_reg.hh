#include "ponded_depression_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,PondedDepressionDepthEvaluator> PondedDepressionDepthEvaluator::factory_("ponded depression depth");

}
}
}
