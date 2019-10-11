#include "evaluator/EvaluatorSecondaryMonotypeMultiplicative.hh"
namespace Amanzi {

template <>
Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeMultiplicative<
                                      CompositeVector, CompositeVectorSpace>>
  EvaluatorSecondaryMonotypeMultiplicative<
    CompositeVector, CompositeVectorSpace>::fac_("multiplicative");

} // namespace Amanzi
