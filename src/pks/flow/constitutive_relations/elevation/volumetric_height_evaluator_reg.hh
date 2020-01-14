#include "volumetric_height_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,VolumetricHeightEvaluator> VolumetricHeightEvaluator::reg_("volumetric ponded depth");

} //namespace
} //namespace
