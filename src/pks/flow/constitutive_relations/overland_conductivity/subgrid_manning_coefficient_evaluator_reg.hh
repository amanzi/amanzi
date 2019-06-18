#include "subgrid_manning_coefficient_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SubgridManningCoefficientEvaluator> SubgridManningCoefficientEvaluator::factory_("subgrid manning coefficient");

}
}
}
