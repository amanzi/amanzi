#include "albedo_threeecomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AlbedoThreeComponentEvaluator>
AlbedoThreeComponentEvaluator::reg_("three-component subgrid albedos");

} //namespace
} //namespace
