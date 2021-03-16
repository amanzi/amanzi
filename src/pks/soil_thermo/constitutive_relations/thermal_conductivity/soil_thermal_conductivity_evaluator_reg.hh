#include "soil_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilThermalConductivityEvaluator> SoilThermalConductivityEvaluator::factory_("soil thermal conductivity");

} //namespace
} //namespace
