#include "soil_heat_capacity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilHeatCapacityEvaluator> SoilHeatCapacityEvaluator::factory_("soil heat capacity");

} //namespace
} //namespace
