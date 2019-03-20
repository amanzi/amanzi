#include "evaluator/EvaluatorSecondaryMonotypeAdditive.hh"
namespace Amanzi {

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeAdditive<CompositeVector,CompositeVectorSpace>>
EvaluatorSecondaryMonotypeAdditive<CompositeVector,CompositeVectorSpace>::fac_("additive");

} // namespace Amanzi
