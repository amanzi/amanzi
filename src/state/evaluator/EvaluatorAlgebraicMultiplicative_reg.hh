#include "evaluator/EvaluatorAlgebraicMultiplicative.hh"
namespace Amanzi {

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorAlgebraicMultiplicative<CompositeVector,CompositeVectorSpace>>
 EvaluatorAlgebraicMultiplicative<CompositeVector,CompositeVectorSpace>::fac_("multiplicative");

} // namespace Amanzi
