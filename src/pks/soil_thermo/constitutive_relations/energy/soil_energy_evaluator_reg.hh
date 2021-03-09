#include "soil_energy_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilEnergyEvaluator> SoilEnergyEvaluator::factory_("soil energy");

} //namespace
} //namespace
