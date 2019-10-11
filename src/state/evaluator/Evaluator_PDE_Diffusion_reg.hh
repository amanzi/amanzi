#include "evaluator/Evaluator_PDE_Diffusion.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Diffusion>
  Evaluator_PDE_Diffusion::fac_("diffusion operator");

} // namespace Amanzi
