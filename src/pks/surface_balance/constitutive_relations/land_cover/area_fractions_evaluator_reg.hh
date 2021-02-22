#include "area_fractions_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AreaFractionsEvaluator> AreaFractionsEvaluator::reg_("snow area fraction");

} //namespace
} //namespace
