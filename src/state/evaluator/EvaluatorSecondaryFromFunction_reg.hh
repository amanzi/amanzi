#include "evaluator/EvaluatorSecondaryFromFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryFromFunction>
    EvaluatorSecondaryFromFunction::fac_("secondary variable from function");

} // namespace Amanzi
