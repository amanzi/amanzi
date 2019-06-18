#include "area_fractions_subgrid_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,AreaFractionsSubgridEvaluator> AreaFractionsSubgridEvaluator::reg_("surface balance subgrid area fractions");

} //namespace
} //namespace
