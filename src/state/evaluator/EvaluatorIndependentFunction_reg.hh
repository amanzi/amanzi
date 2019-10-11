#include "evaluator/EvaluatorIndependentFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction>
  EvaluatorIndependentFunction::fac_("independent variable");

} // namespace Amanzi
