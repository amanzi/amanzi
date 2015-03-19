#include "VWContentEvaluator.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<FieldEvaluator,VWContentEvaluator> VWContentEvaluator::reg_("water content");

}  // namespace Flow
}  // namespace Amanzi
