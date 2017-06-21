#include "overland_conductivity_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,OverlandConductivityEvaluator> OverlandConductivityEvaluator::factory_("overland conductivity");

}
}
