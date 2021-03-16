#include "soil_density_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilDensityEvaluator> SoilDensityEvaluator::factory_("soil density");

} //namespace
} //namespace
