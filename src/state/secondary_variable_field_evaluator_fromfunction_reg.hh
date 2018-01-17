#include "secondary_variable_field_evaluator_fromfunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,
                         SecondaryVariableFieldEvaluatorFromFunction>
    SecondaryVariableFieldEvaluatorFromFunction::fac_(
        "secondary variable from function");

} // namespace Amanzi
