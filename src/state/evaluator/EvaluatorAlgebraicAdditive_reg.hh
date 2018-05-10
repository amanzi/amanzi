#include "evaluator/EvaluatorAlgebraicAdditive.hh"
namespace Amanzi {

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorAlgebraicAdditive<CompositeVector,CompositeVectorSpace>>
EvaluatorAlgebraicAdditive<CompositeVector,CompositeVectorSpace>::fac_("additive evaluator");

} // namespace Amanzi
