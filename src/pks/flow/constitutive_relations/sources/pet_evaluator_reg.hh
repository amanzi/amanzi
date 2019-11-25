#include "pet_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,PETEvaluator> PETEvaluator::reg_("potential evapotranspiration");

} //namespace
} //namespace
} //namespace
