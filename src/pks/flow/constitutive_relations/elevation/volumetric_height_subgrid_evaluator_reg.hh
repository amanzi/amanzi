#include "volumetric_height_subgrid_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,VolumetricHeightSubgridEvaluator> VolumetricHeightSubgridEvaluator::reg_("volumetric ponded and snow depths");

} //namespace
} //namespace
