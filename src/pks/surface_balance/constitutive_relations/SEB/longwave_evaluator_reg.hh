#include "longwave_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LongwaveEvaluator> LongwaveEvaluator::reg_("incoming longwave radiation");

} //namespace
} //namespace
