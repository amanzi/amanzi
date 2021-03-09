#include "thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ThermalConductivityEvaluator> ThermalConductivityEvaluator::factory_("soil thermal conductivity");

} //namespace
} //namespace
