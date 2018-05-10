#include "evaluator/Evaluator_BCs.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_BCs>
    Evaluator_BCs::fac_("general boundary conditions");

} // namespace Amanzi
