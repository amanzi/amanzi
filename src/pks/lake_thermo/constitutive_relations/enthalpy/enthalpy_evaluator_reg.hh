#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EnthalpyEvaluator> EnthalpyEvaluator::factory_("enthalpy");

} //namespace
} //namespace
