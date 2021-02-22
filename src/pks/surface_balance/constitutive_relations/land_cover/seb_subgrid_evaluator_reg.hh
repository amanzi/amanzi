#include "seb_subgrid_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SubgridEvaluator> SubgridEvaluator::reg_("surface balance subgrid");

} //namespace
} //namespace
