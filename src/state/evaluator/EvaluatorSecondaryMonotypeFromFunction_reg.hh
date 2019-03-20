#include "evaluator/EvaluatorSecondaryMonotypeFromFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction>
    EvaluatorSecondaryMonotypeFromFunction::fac_("secondary variable from function");

} // namespace Amanzi
