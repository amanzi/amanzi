#include "richards_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,RichardsEnergyEvaluator> RichardsEnergyEvaluator::reg_("richards energy");

} //namespace
} //namespace
} //namespace
