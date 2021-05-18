#include "albedo_twocomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AlbedoTwoComponentEvaluator>
  AlbedoTwoComponentEvaluator::reg_("two-component subgrid albedos");

} //namespace
} //namespace
