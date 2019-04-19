#include "area_fractions_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AreaFractionsEvaluator> AreaFractionsEvaluator::reg_("surface balance area fractions");

} //namespace
} //namespace
