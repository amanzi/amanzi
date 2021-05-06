#include "soil_enthalpy_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilEnthalpyEvaluator> SoilEnthalpyEvaluator::factory_("enthalpy");

} //namespace
} //namespace
