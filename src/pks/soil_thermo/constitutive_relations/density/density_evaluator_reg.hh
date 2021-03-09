#include "density_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,DensityEvaluator> DensityEvaluator::factory_("density");

} //namespace
} //namespace
