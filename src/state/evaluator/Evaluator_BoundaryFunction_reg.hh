#include "evaluator/Evaluator_BoundaryFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_BoundaryFunction>
    Evaluator_BoundaryFunction::fac_("boundary function");

} // namespace Amanzi
