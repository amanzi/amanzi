#include "albedo_subgrid_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AlbedoSubgridEvaluator> AlbedoSubgridEvaluator::reg_("albedo subgrid");

} //namespace
} //namespace
