#include "bioturbation_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

// registry of method 
  Utils::RegisteredFactory<FieldEvaluator,BioturbationEvaluator> BioturbationEvaluator::fac_("bioturbation evaluator"); 

} //namespace
} //namespace
} //namespace
