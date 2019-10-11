#include "evaluator/Evaluator_PDE_Accumulation.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Accumulation>
  Evaluator_PDE_Accumulation::fac_("accumulation operator");

} // namespace Amanzi
