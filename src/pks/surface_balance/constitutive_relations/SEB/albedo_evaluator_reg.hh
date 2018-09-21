#include "albedo_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AlbedoEvaluator> AlbedoEvaluator::reg_("albedo");

} //namespace
} //namespace
