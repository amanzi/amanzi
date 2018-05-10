#include "evaluator/EvaluatorAlgebraicFromFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorAlgebraicFromFunction>
    EvaluatorAlgebraicFromFunction::fac_("secondary variable from function");

} // namespace Amanzi
