#include "WaterContentEvaluator.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<FieldEvaluator,WaterContentEvaluator> WaterContentEvaluator::reg_("water content");

}  // namespace Flow
}  // namespace Amanzi
